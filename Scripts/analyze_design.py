#!/usr/bin/env python3
"""Analyze sequence diversity for a materialized CombPhage assembly design.

This module replaces the old analyze_fragments.py and analyze_assemblies.py
scripts. It is designed around the outputs of design_assemblies.py and can be
used either directly from the command line or imported by design_assemblies.py
and the Streamlit GUI.

Main analyses
-------------
1. Pairwise nucleotide identity among parental versions of each fragment.
2. Homologous CDS nucleotide and amino-acid identity within each fragment.
3. Identification of the most variable CDS clusters per fragment.
4. Assembly-vs-parent identity for every generated chimera.
5. Optional all-vs-all assembly identity when the number of assemblies is
   small enough to remain useful and computationally reasonable.
6. Publication/debug-friendly TSV tables and PNG visualizations.

Assembly identities are calculated from aligned *core fragment contributions*
(the exact non-duplicated sequence contributions used when chimeras are
constructed) rather than by repeatedly aligning entire chimeric genomes. This
preserves the fragment structure of the design and scales much better for
combinatorial libraries.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Align import PairwiseAligner


DEFAULT_MAX_PAIRWISE_ASSEMBLIES = 250
DEFAULT_TOP_VARIABLE_CDS = 5
DEFAULT_MAX_HEATMAP_ASSEMBLIES = 80


def make_global_aligner() -> PairwiseAligner:
    """Return a deterministic global aligner for nucleotide/protein identity."""
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -5.0
    aligner.extend_gap_score = -0.5
    return aligner


def pairwise_identity_stats(
    seq1: str,
    seq2: str,
    aligner: Optional[PairwiseAligner] = None,
) -> dict:
    """Globally align two sequences and return identity/alignment statistics."""
    seq1 = str(seq1).upper()
    seq2 = str(seq2).upper()

    if not seq1 and not seq2:
        return {
            "identity": 100.0,
            "matches": 0,
            "alignment_length": 0,
            "mismatches": 0,
            "gap_columns": 0,
            "aligned_seq1": "",
            "aligned_seq2": "",
        }

    if not seq1 or not seq2:
        alignment_length = max(len(seq1), len(seq2))
        return {
            "identity": 0.0,
            "matches": 0,
            "alignment_length": alignment_length,
            "mismatches": 0,
            "gap_columns": alignment_length,
            "aligned_seq1": seq1 or ("-" * alignment_length),
            "aligned_seq2": seq2 or ("-" * alignment_length),
        }

    if seq1 == seq2:
        return {
            "identity": 100.0,
            "matches": len(seq1),
            "alignment_length": len(seq1),
            "mismatches": 0,
            "gap_columns": 0,
            "aligned_seq1": seq1,
            "aligned_seq2": seq2,
        }

    if aligner is None:
        aligner = make_global_aligner()

    alignment = aligner.align(seq1, seq2)[0]
    aligned1 = str(alignment[0])
    aligned2 = str(alignment[1])

    matches = 0
    mismatches = 0
    gap_columns = 0

    for a, b in zip(aligned1, aligned2):
        if a == "-" or b == "-":
            gap_columns += 1
        elif a == b:
            matches += 1
        else:
            mismatches += 1

    alignment_length = len(aligned1)
    identity = (
        (matches / alignment_length) * 100.0
        if alignment_length
        else 0.0
    )

    return {
        "identity": float(identity),
        "matches": int(matches),
        "alignment_length": int(alignment_length),
        "mismatches": int(mismatches),
        "gap_columns": int(gap_columns),
        "aligned_seq1": aligned1,
        "aligned_seq2": aligned2,
    }


def write_pairwise_alignment(
    output_path: Path,
    name1: str,
    name2: str,
    stats: dict,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as handle:
        handle.write(f">{name1}\n{stats['aligned_seq1']}\n")
        handle.write(f">{name2}\n{stats['aligned_seq2']}\n")
        handle.write(f"; percent_identity={stats['identity']:.4f}\n")
        handle.write(f"; matches={stats['matches']}\n")
        handle.write(f"; alignment_length={stats['alignment_length']}\n")


def read_genomes_fasta(fasta_path: Path) -> Dict[str, str]:
    sequences = {}
    for record in SeqIO.parse(str(fasta_path), "fasta"):
        sequences[record.id] = str(record.seq).upper()
    if not sequences:
        raise ValueError(f"No sequences found in {fasta_path}")
    return sequences


def find_prepped_genbanks(prepped_dir: Path) -> List[Path]:
    genbank_dir = prepped_dir / "genbanks"
    if not genbank_dir.exists():
        raise FileNotFoundError(
            f"Expected normalized GenBanks at {genbank_dir}"
        )

    extensions = {".gb", ".gbk", ".gbff", ".genbank"}
    files = sorted(
        path for path in genbank_dir.iterdir()
        if path.is_file() and path.suffix.lower() in extensions
    )
    if not files:
        raise FileNotFoundError(
            f"No normalized GenBank files found in {genbank_dir}"
        )
    return files


def parse_normalized_cds(prepped_dir: Path) -> pd.DataFrame:
    """Extract normalized CDS sequences and annotations from prepared GenBanks."""
    rows = []

    for path in find_prepped_genbanks(prepped_dir):
        for record in SeqIO.parse(str(path), "genbank"):
            genome = record.id
            genome_length = len(record.seq)

            for feature in record.features:
                if feature.type != "CDS":
                    continue

                qualifiers = feature.qualifiers
                gene = None
                for key in ("locus_tag", "protein_id", "gene"):
                    values = qualifiers.get(key)
                    if values:
                        gene = str(values[0])
                        break

                if gene is None:
                    continue

                cds_seq = str(feature.extract(record.seq)).upper()
                translation = (
                    str(qualifiers.get("translation", [""])[0]).upper()
                    if qualifiers.get("translation")
                    else ""
                )

                product = str(
                    qualifiers.get("product", ["hypothetical protein"])[0]
                )

                low = int(feature.location.start) + 1
                high = int(feature.location.end)
                midpoint = (low + high) / 2.0
                strand = feature.location.strand

                rows.append({
                    "gene": gene,
                    "genome": genome,
                    "genome_length": genome_length,
                    "start": low,
                    "stop": high,
                    "midpoint": midpoint,
                    "strand": strand,
                    "product": product,
                    "cds_nt": cds_seq,
                    "protein": translation,
                    "cds_nt_length": len(cds_seq),
                    "protein_length": len(translation),
                })

    if not rows:
        raise ValueError(
            f"No CDS features could be parsed from {prepped_dir / 'genbanks'}"
        )

    return pd.DataFrame(rows)


def position_in_interval(
    position: float,
    start: int,
    stop: int,
    seq_len: int,
) -> bool:
    """Test a position against a 1-based inclusive circular interval."""
    position = float(position)
    start = int(start)
    stop = int(stop)

    if start <= stop:
        return start <= position <= stop
    return position >= start or position <= stop


def _interval_span(start: int, stop: int, seq_len: int) -> int:
    if start <= stop:
        return stop - start + 1
    return (seq_len - start + 1) + stop


def load_fragment_regions(
    materialized_dir: Path,
    genome_sequences: Dict[str, str],
) -> pd.DataFrame:
    """Load unique per-genome fragment contribution regions.

    When chimeras were generated, the chimera assembly log contains the exact
    non-overlapping contribution of each source genome to each fragment. Those
    coordinates are preferred for CDS-to-fragment assignment. Otherwise the
    overlapping fragment-manifest coordinates are used as a fallback.
    """
    chimera_log = materialized_dir / "chimeras" / "chimera_assembly_log.tsv"

    if chimera_log.exists() and chimera_log.stat().st_size:
        log = pd.read_csv(chimera_log, sep="\t")
        required = {
            "fragment_number",
            "source_genome",
            "actual_parent_start",
            "actual_parent_stop",
        }
        missing = required - set(log.columns)
        if not missing:
            regions = (
                log[list(required)]
                .drop_duplicates()
                .rename(columns={
                    "actual_parent_start": "region_start",
                    "actual_parent_stop": "region_stop",
                })
            )
            regions["region_source"] = "chimera_core"
            return regions.sort_values(
                ["fragment_number", "source_genome"]
            ).reset_index(drop=True)

    manifest_path = materialized_dir / "fragment_manifest.tsv"
    if not manifest_path.exists():
        raise FileNotFoundError(
            f"Missing fragment manifest: {manifest_path}"
        )

    manifest = pd.read_csv(manifest_path, sep="\t")
    required = {
        "fragment_number",
        "source_genome",
        "parent_start",
        "parent_stop",
    }
    missing = required - set(manifest.columns)
    if missing:
        raise ValueError(
            "fragment_manifest.tsv is missing: "
            + ", ".join(sorted(missing))
        )

    regions = manifest[list(required)].copy()
    regions = regions.rename(columns={
        "parent_start": "region_start",
        "parent_stop": "region_stop",
    })
    regions["region_source"] = "fragment_manifest"
    return regions.sort_values(
        ["fragment_number", "source_genome"]
    ).reset_index(drop=True)


def assign_cds_to_fragments(
    cds_df: pd.DataFrame,
    fragment_regions: pd.DataFrame,
    genome_sequences: Dict[str, str],
) -> pd.DataFrame:
    """Assign each CDS to one fragment based on its genomic midpoint."""
    region_lookup = {}
    for genome, group in fragment_regions.groupby("source_genome"):
        region_lookup[str(genome)] = group.to_dict("records")

    fragment_numbers = []
    assignment_sources = []

    for _, cds in cds_df.iterrows():
        genome = str(cds["genome"])
        seq_len = len(genome_sequences[genome])
        midpoint = float(cds["midpoint"])
        candidates = []

        for region in region_lookup.get(genome, []):
            if position_in_interval(
                midpoint,
                int(region["region_start"]),
                int(region["region_stop"]),
                seq_len,
            ):
                candidates.append(region)

        if not candidates:
            fragment_numbers.append(np.nan)
            assignment_sources.append("unassigned")
            continue

        # Manifest regions can overlap. Prefer the shortest interval containing
        # the CDS midpoint, which minimizes overlap-driven double assignment.
        candidates.sort(
            key=lambda region: _interval_span(
                int(region["region_start"]),
                int(region["region_stop"]),
                seq_len,
            )
        )
        chosen = candidates[0]
        fragment_numbers.append(int(chosen["fragment_number"]))
        assignment_sources.append(str(chosen["region_source"]))

    result = cds_df.copy()
    result["fragment_number"] = fragment_numbers
    result["fragment_assignment_source"] = assignment_sources
    return result


def analyze_fragment_fastas(
    materialized_dir: Path,
    outdir: Path,
    save_alignments: bool = False,
) -> Tuple[pd.DataFrame, pd.DataFrame, Dict[int, dict]]:
    fragments_dir = materialized_dir / "fragments"
    manifest_path = materialized_dir / "fragment_manifest.tsv"

    if not fragments_dir.exists():
        raise FileNotFoundError(f"Missing fragments directory: {fragments_dir}")
    if not manifest_path.exists():
        raise FileNotFoundError(f"Missing fragment manifest: {manifest_path}")

    manifest = pd.read_csv(manifest_path, sep="\t")
    record_to_genome = dict(
        zip(manifest["fasta_record"].astype(str), manifest["source_genome"].astype(str))
    )

    matrix_dir = outdir / "fragment_identity_matrices"
    alignment_dir = outdir / "alignments" / "fragments"
    matrix_dir.mkdir(parents=True, exist_ok=True)
    if save_alignments:
        alignment_dir.mkdir(parents=True, exist_ok=True)

    aligner = make_global_aligner()
    pair_rows = []
    summary_rows = []
    lookup = {}

    fasta_files = sorted(fragments_dir.glob("Fragment_*.fasta"))

    for fasta_path in fasta_files:
        try:
            fragment_number = int(fasta_path.stem.split("_")[-1])
        except ValueError:
            continue

        records = list(SeqIO.parse(str(fasta_path), "fasta"))
        if not records:
            continue

        by_genome = {}
        for record in records:
            genome = record_to_genome.get(record.id)
            if genome is None:
                # Fallback for manually generated files: choose the longest
                # known genome ID that prefixes the record name.
                known = sorted(set(record_to_genome.values()), key=len, reverse=True)
                genome = next(
                    (g for g in known if record.id.startswith(g + "_")),
                    record.id,
                )
            by_genome[genome] = str(record.seq).upper()

        genomes = sorted(by_genome)
        matrix = pd.DataFrame(100.0, index=genomes, columns=genomes)
        identities = []
        lookup[fragment_number] = {}

        for genome in genomes:
            lookup[fragment_number].setdefault(genome, {})
            seq = by_genome[genome]
            lookup[fragment_number][genome][genome] = {
                "identity": 100.0,
                "matches": len(seq),
                "alignment_length": len(seq),
            }

        for genome_a, genome_b in itertools.combinations(genomes, 2):
            stats = pairwise_identity_stats(
                by_genome[genome_a],
                by_genome[genome_b],
                aligner,
            )
            identities.append(stats["identity"])
            matrix.loc[genome_a, genome_b] = stats["identity"]
            matrix.loc[genome_b, genome_a] = stats["identity"]

            compact = {
                "identity": stats["identity"],
                "matches": stats["matches"],
                "alignment_length": stats["alignment_length"],
            }
            lookup[fragment_number].setdefault(genome_a, {})[genome_b] = compact
            lookup[fragment_number].setdefault(genome_b, {})[genome_a] = compact

            pair_rows.append({
                "fragment_number": fragment_number,
                "genome_a": genome_a,
                "genome_b": genome_b,
                "nucleotide_identity": stats["identity"],
                "matches": stats["matches"],
                "alignment_length": stats["alignment_length"],
                "mismatches": stats["mismatches"],
                "gap_columns": stats["gap_columns"],
                "length_a": len(by_genome[genome_a]),
                "length_b": len(by_genome[genome_b]),
            })

            if save_alignments:
                write_pairwise_alignment(
                    alignment_dir
                    / f"Fragment_{fragment_number}"
                    / f"{genome_a}_vs_{genome_b}.fasta",
                    genome_a,
                    genome_b,
                    stats,
                )

        matrix.to_csv(
            matrix_dir / f"Fragment_{fragment_number}_identity.tsv",
            sep="\t",
        )

        lengths = [len(sequence) for sequence in by_genome.values()]
        summary_rows.append({
            "fragment_number": fragment_number,
            "n_genomes": len(genomes),
            "mean_length": float(np.mean(lengths)),
            "min_length": int(min(lengths)),
            "max_length": int(max(lengths)),
            "length_range": int(max(lengths) - min(lengths)),
            "mean_pairwise_nucleotide_identity": (
                float(np.mean(identities)) if identities else 100.0
            ),
            "min_pairwise_nucleotide_identity": (
                float(np.min(identities)) if identities else 100.0
            ),
            "max_pairwise_nucleotide_identity": (
                float(np.max(identities)) if identities else 100.0
            ),
        })

    pair_df = pd.DataFrame(pair_rows)
    summary_df = pd.DataFrame(summary_rows).sort_values("fragment_number")

    pair_df.to_csv(outdir / "fragment_pairwise_identity.tsv", sep="\t", index=False)
    summary_df.to_csv(outdir / "fragment_sequence_summary.tsv", sep="\t", index=False)

    return pair_df, summary_df, lookup


def analyze_cds_diversity(
    prepped_dir: Path,
    cluster_membership_path: Path,
    materialized_dir: Path,
    genome_sequences: Dict[str, str],
    outdir: Path,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    cds = parse_normalized_cds(prepped_dir)
    clusters = pd.read_csv(cluster_membership_path, sep="\t", dtype=str)

    required = {"gene", "cluster"}
    missing = required - set(clusters.columns)
    if missing:
        raise ValueError(
            "cluster_membership.tsv is missing: "
            + ", ".join(sorted(missing))
        )

    merge_cols = ["gene", "cluster"]
    for optional in ("resolution", "source_mmseqs_representative"):
        if optional in clusters.columns:
            merge_cols.append(optional)

    cds = cds.merge(
        clusters[merge_cols].drop_duplicates("gene"),
        on="gene",
        how="left",
    )

    fragment_regions = load_fragment_regions(materialized_dir, genome_sequences)
    cds = assign_cds_to_fragments(cds, fragment_regions, genome_sequences)
    cds.to_csv(outdir / "cds_fragment_assignments.tsv", sep="\t", index=False)

    aligner = make_global_aligner()
    pair_rows = []
    cluster_rows = []

    clustered = cds.dropna(subset=["cluster", "fragment_number"]).copy()
    clustered["fragment_number"] = clustered["fragment_number"].astype(int)

    for cluster, group in clustered.groupby("cluster", sort=False):
        # A resolved CombPhage cluster should normally contain one CDS per
        # genome. If unexpected duplicates remain, keep all pairwise comparisons
        # and expose that fact via n_genes/n_genomes.
        records = group.to_dict("records")
        if len(records) < 2:
            continue

        nt_values = []
        protein_values = []
        fragment_values = []

        for a, b in itertools.combinations(records, 2):
            if a["genome"] == b["genome"]:
                continue

            nt = pairwise_identity_stats(a["cds_nt"], b["cds_nt"], aligner)
            protein_identity = np.nan
            protein_matches = np.nan
            protein_alignment_length = np.nan

            if a["protein"] and b["protein"]:
                prot = pairwise_identity_stats(a["protein"], b["protein"], aligner)
                protein_identity = prot["identity"]
                protein_matches = prot["matches"]
                protein_alignment_length = prot["alignment_length"]
                protein_values.append(prot["identity"])

            nt_values.append(nt["identity"])
            fragment_values.extend([
                int(a["fragment_number"]),
                int(b["fragment_number"]),
            ])

            pair_rows.append({
                "cluster": cluster,
                "fragment_a": int(a["fragment_number"]),
                "fragment_b": int(b["fragment_number"]),
                "gene_a": a["gene"],
                "genome_a": a["genome"],
                "gene_b": b["gene"],
                "genome_b": b["genome"],
                "product_a": a["product"],
                "product_b": b["product"],
                "cds_nucleotide_identity": nt["identity"],
                "cds_nt_matches": nt["matches"],
                "cds_nt_alignment_length": nt["alignment_length"],
                "protein_identity": protein_identity,
                "protein_matches": protein_matches,
                "protein_alignment_length": protein_alignment_length,
            })

        if not nt_values:
            continue

        fragment_counts = pd.Series(fragment_values).value_counts()
        assigned_fragment = int(fragment_counts.index[0])
        fragment_consistent = len(fragment_counts) == 1
        products = sorted(set(str(x) for x in group["product"].dropna()))

        row = {
            "fragment_number": assigned_fragment,
            "fragment_consistent_across_genomes": bool(fragment_consistent),
            "cluster": cluster,
            "products": " | ".join(products),
            "n_genes": int(len(group)),
            "n_genomes": int(group["genome"].nunique()),
            "genes": ";".join(group["gene"].astype(str)),
            "genomes": ";".join(group["genome"].astype(str)),
            "mean_pairwise_cds_nucleotide_identity": float(np.mean(nt_values)),
            "min_pairwise_cds_nucleotide_identity": float(np.min(nt_values)),
            "max_pairwise_cds_nucleotide_identity": float(np.max(nt_values)),
            "cds_nucleotide_divergence": float(100.0 - np.mean(nt_values)),
            "mean_pairwise_protein_identity": (
                float(np.mean(protein_values)) if protein_values else np.nan
            ),
            "min_pairwise_protein_identity": (
                float(np.min(protein_values)) if protein_values else np.nan
            ),
            "max_pairwise_protein_identity": (
                float(np.max(protein_values)) if protein_values else np.nan
            ),
            "protein_divergence": (
                float(100.0 - np.mean(protein_values))
                if protein_values
                else np.nan
            ),
        }

        if "resolution" in group.columns:
            values = sorted(set(group["resolution"].dropna().astype(str)))
            row["resolution"] = ";".join(values)

        cluster_rows.append(row)

    pair_df = pd.DataFrame(pair_rows)
    cluster_df = pd.DataFrame(cluster_rows)

    if not cluster_df.empty:
        cluster_df = cluster_df.sort_values(
            ["fragment_number", "mean_pairwise_cds_nucleotide_identity"]
        )

    fragment_rows = []
    if not cluster_df.empty:
        for fragment_number, group in cluster_df.groupby("fragment_number"):
            most_variable = group.sort_values(
                "mean_pairwise_cds_nucleotide_identity"
            ).iloc[0]
            fragment_rows.append({
                "fragment_number": int(fragment_number),
                "n_compared_cds_clusters": int(len(group)),
                "mean_cds_nucleotide_identity": float(
                    group["mean_pairwise_cds_nucleotide_identity"].mean()
                ),
                "median_cds_nucleotide_identity": float(
                    group["mean_pairwise_cds_nucleotide_identity"].median()
                ),
                "minimum_cds_nucleotide_identity": float(
                    group["min_pairwise_cds_nucleotide_identity"].min()
                ),
                "mean_protein_identity": float(
                    group["mean_pairwise_protein_identity"].mean()
                ) if group["mean_pairwise_protein_identity"].notna().any() else np.nan,
                "most_variable_cluster": most_variable["cluster"],
                "most_variable_product": most_variable["products"],
                "most_variable_cluster_mean_nt_identity": float(
                    most_variable["mean_pairwise_cds_nucleotide_identity"]
                ),
                "most_variable_cluster_mean_protein_identity": (
                    float(most_variable["mean_pairwise_protein_identity"])
                    if pd.notna(most_variable["mean_pairwise_protein_identity"])
                    else np.nan
                ),
            })

    fragment_cds_df = pd.DataFrame(fragment_rows)

    pair_df.to_csv(outdir / "cds_pairwise_identity.tsv", sep="\t", index=False)
    cluster_df.to_csv(outdir / "cds_cluster_summary.tsv", sep="\t", index=False)
    fragment_cds_df.to_csv(outdir / "fragment_cds_summary.tsv", sep="\t", index=False)

    return pair_df, cluster_df, fragment_cds_df


def extract_inclusive_interval(
    sequence: str,
    start: int,
    stop: int,
) -> str:
    """Extract 1-based inclusive coordinates with circular wrap support."""
    start_idx = int(start) - 1
    stop = int(stop)
    if start <= stop:
        return sequence[start_idx:stop]
    return sequence[start_idx:] + sequence[:stop]


def build_core_fragment_sequences(
    materialized_dir: Path,
    genome_sequences: Dict[str, str],
) -> Dict[int, Dict[str, str]]:
    """Return exact source sequences contributed by each fragment to chimeras."""
    regions = load_fragment_regions(materialized_dir, genome_sequences)
    core = {}

    for _, row in regions.iterrows():
        fragment = int(row["fragment_number"])
        genome = str(row["source_genome"])
        if genome not in genome_sequences:
            continue
        sequence = extract_inclusive_interval(
            genome_sequences[genome],
            int(row["region_start"]),
            int(row["region_stop"]),
        )
        core.setdefault(fragment, {})[genome] = sequence

    return core


def build_core_identity_lookup(
    core_sequences: Dict[int, Dict[str, str]],
) -> Dict[int, Dict[str, Dict[str, dict]]]:
    aligner = make_global_aligner()
    lookup = {}

    for fragment, sequence_map in core_sequences.items():
        lookup[fragment] = {}
        genomes = sorted(sequence_map)

        for genome in genomes:
            seq = sequence_map[genome]
            lookup[fragment].setdefault(genome, {})[genome] = {
                "identity": 100.0,
                "matches": len(seq),
                "alignment_length": len(seq),
            }

        for a, b in itertools.combinations(genomes, 2):
            stats = pairwise_identity_stats(
                sequence_map[a], sequence_map[b], aligner
            )
            compact = {
                "identity": stats["identity"],
                "matches": stats["matches"],
                "alignment_length": stats["alignment_length"],
            }
            lookup[fragment].setdefault(a, {})[b] = compact
            lookup[fragment].setdefault(b, {})[a] = compact

    return lookup


def load_assembly_compositions(
    materialized_dir: Path,
    genomes: List[str],
    fragment_numbers: List[int],
) -> pd.DataFrame:
    """Return one source-genome assignment per assembly and fragment."""
    rows = []

    # Parental assemblies are useful references even though they are not written
    # into the chimera directory.
    for genome in genomes:
        for fragment in fragment_numbers:
            rows.append({
                "assembly": genome,
                "assembly_type": "parental",
                "fragment_number": fragment,
                "source_genome": genome,
            })

    chimera_log = materialized_dir / "chimeras" / "chimera_assembly_log.tsv"
    if chimera_log.exists() and chimera_log.stat().st_size:
        log = pd.read_csv(chimera_log, sep="\t")
        required = {"chimera_id", "fragment_number", "source_genome"}
        missing = required - set(log.columns)
        if not missing:
            unique = log[list(required)].drop_duplicates()
            for _, row in unique.iterrows():
                rows.append({
                    "assembly": str(row["chimera_id"]),
                    "assembly_type": "chimera",
                    "fragment_number": int(row["fragment_number"]),
                    "source_genome": str(row["source_genome"]),
                })

    result = pd.DataFrame(rows)
    if result.empty:
        raise ValueError("No assembly compositions could be determined.")

    return result.sort_values(
        ["assembly_type", "assembly", "fragment_number"]
    ).reset_index(drop=True)


def compositional_identity(
    composition_a: Dict[int, str],
    composition_b: Dict[int, str],
    identity_lookup: Dict[int, Dict[str, Dict[str, dict]]],
) -> dict:
    matches = 0
    alignment_length = 0
    compared_fragments = 0

    for fragment in sorted(set(composition_a) & set(composition_b)):
        source_a = composition_a[fragment]
        source_b = composition_b[fragment]
        try:
            stats = identity_lookup[fragment][source_a][source_b]
        except KeyError:
            continue
        matches += int(stats["matches"])
        alignment_length += int(stats["alignment_length"])
        compared_fragments += 1

    identity = (
        (matches / alignment_length) * 100.0
        if alignment_length
        else np.nan
    )

    return {
        "identity": float(identity) if pd.notna(identity) else np.nan,
        "matches": int(matches),
        "alignment_length": int(alignment_length),
        "compared_fragments": int(compared_fragments),
    }


def analyze_assemblies(
    materialized_dir: Path,
    genome_sequences: Dict[str, str],
    outdir: Path,
    max_pairwise_assemblies: int = DEFAULT_MAX_PAIRWISE_ASSEMBLIES,
) -> Tuple[pd.DataFrame, pd.DataFrame, Optional[pd.DataFrame], pd.DataFrame]:
    core_sequences = build_core_fragment_sequences(
        materialized_dir, genome_sequences
    )
    fragment_numbers = sorted(core_sequences)
    genomes = sorted(genome_sequences)
    identity_lookup = build_core_identity_lookup(core_sequences)
    composition_df = load_assembly_compositions(
        materialized_dir,
        genomes,
        fragment_numbers,
    )
    composition_df.to_csv(
        outdir / "assembly_composition.tsv", sep="\t", index=False
    )

    composition_maps = {}
    type_map = {}
    for assembly, group in composition_df.groupby("assembly"):
        composition_maps[assembly] = {
            int(row["fragment_number"]): str(row["source_genome"])
            for _, row in group.iterrows()
        }
        type_map[assembly] = str(group["assembly_type"].iloc[0])

    assemblies = sorted(
        composition_maps,
        key=lambda name: (0 if type_map[name] == "parental" else 1, name),
    )

    parent_rows = []
    for assembly in assemblies:
        for parent in genomes:
            stats = compositional_identity(
                composition_maps[assembly],
                {fragment: parent for fragment in fragment_numbers},
                identity_lookup,
            )
            parent_rows.append({
                "assembly": assembly,
                "assembly_type": type_map[assembly],
                "parent": parent,
                "nucleotide_identity": stats["identity"],
                "matches": stats["matches"],
                "alignment_length": stats["alignment_length"],
                "compared_fragments": stats["compared_fragments"],
            })

    assembly_parent_df = pd.DataFrame(parent_rows)
    assembly_parent_df.to_csv(
        outdir / "assembly_vs_parent_identity.tsv", sep="\t", index=False
    )

    assembly_summary_rows = []
    for assembly, group in assembly_parent_df.groupby("assembly"):
        best = group.sort_values("nucleotide_identity", ascending=False).iloc[0]
        source_pattern = composition_maps[assembly]
        assembly_summary_rows.append({
            "assembly": assembly,
            "assembly_type": type_map[assembly],
            "source_pattern": "-".join(
                source_pattern[fragment] for fragment in fragment_numbers
            ),
            "nearest_parent": best["parent"],
            "nearest_parent_identity": float(best["nucleotide_identity"]),
            "mean_parent_identity": float(group["nucleotide_identity"].mean()),
            "min_parent_identity": float(group["nucleotide_identity"].min()),
            "max_parent_identity": float(group["nucleotide_identity"].max()),
        })

    assembly_summary_df = pd.DataFrame(assembly_summary_rows)
    assembly_summary_df.to_csv(
        outdir / "assembly_summary.tsv", sep="\t", index=False
    )

    pairwise_df = None
    if len(assemblies) <= max_pairwise_assemblies:
        pair_rows = []
        matrix = pd.DataFrame(
            100.0,
            index=assemblies,
            columns=assemblies,
        )

        for a, b in itertools.combinations(assemblies, 2):
            stats = compositional_identity(
                composition_maps[a], composition_maps[b], identity_lookup
            )
            matrix.loc[a, b] = stats["identity"]
            matrix.loc[b, a] = stats["identity"]
            pair_rows.append({
                "assembly_a": a,
                "assembly_b": b,
                "identity": stats["identity"],
                "matches": stats["matches"],
                "alignment_length": stats["alignment_length"],
                "compared_fragments": stats["compared_fragments"],
            })

        pairwise_df = pd.DataFrame(pair_rows)
        pairwise_df.to_csv(
            outdir / "assembly_pairwise_identity.tsv", sep="\t", index=False
        )
        matrix.to_csv(
            outdir / "assembly_pairwise_identity_matrix.tsv", sep="\t"
        )

    return (
        composition_df,
        assembly_parent_df,
        pairwise_df,
        assembly_summary_df,
    )


def _plot_heatmap(
    matrix: pd.DataFrame,
    output_path: Path,
    title: str,
    annotate_limit: int = 12,
) -> None:
    if matrix.empty:
        return

    fig_width = max(6.0, min(16.0, 0.45 * len(matrix.columns) + 3.0))
    fig_height = max(5.0, min(16.0, 0.40 * len(matrix.index) + 2.5))
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    values = matrix.astype(float).values
    finite = values[np.isfinite(values)]
    vmin = max(0.0, float(np.min(finite)) - 0.5) if finite.size else 0.0
    image = ax.imshow(values, aspect="auto", vmin=vmin, vmax=100.0)
    fig.colorbar(image, ax=ax, label="Nucleotide identity (%)")

    ax.set_xticks(range(len(matrix.columns)))
    ax.set_xticklabels(matrix.columns, rotation=90, fontsize=8)
    ax.set_yticks(range(len(matrix.index)))
    ax.set_yticklabels(matrix.index, fontsize=8)
    ax.set_title(title)

    if len(matrix.index) <= annotate_limit and len(matrix.columns) <= annotate_limit:
        for i in range(len(matrix.index)):
            for j in range(len(matrix.columns)):
                value = values[i, j]
                if np.isfinite(value):
                    ax.text(j, i, f"{value:.1f}", ha="center", va="center", fontsize=7)

    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def make_visualizations(
    outdir: Path,
    fragment_summary: pd.DataFrame,
    cds_cluster_summary: pd.DataFrame,
    assembly_parent_df: pd.DataFrame,
    assembly_pairwise_df: Optional[pd.DataFrame],
    assembly_summary_df: pd.DataFrame,
    top_variable_cds: int = DEFAULT_TOP_VARIABLE_CDS,
    max_heatmap_assemblies: int = DEFAULT_MAX_HEATMAP_ASSEMBLIES,
) -> List[Path]:
    plots_dir = outdir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    outputs = []

    # Per-fragment nucleotide identity matrices.
    matrix_dir = outdir / "fragment_identity_matrices"
    for matrix_path in sorted(matrix_dir.glob("Fragment_*_identity.tsv")):
        matrix = pd.read_csv(matrix_path, sep="\t", index_col=0)
        output = plots_dir / f"{matrix_path.stem}_heatmap.png"
        _plot_heatmap(
            matrix,
            output,
            title=matrix_path.stem.replace("_identity", " pairwise identity"),
        )
        outputs.append(output)

    if not fragment_summary.empty:
        fig, ax = plt.subplots(figsize=(8, 5))
        data = fragment_summary.sort_values("fragment_number")
        x = np.arange(len(data))
        ax.bar(x, data["mean_pairwise_nucleotide_identity"])
        ax.set_xticks(x)
        ax.set_xticklabels([f"Fragment {int(v)}" for v in data["fragment_number"]])
        ax.set_ylabel("Mean pairwise nucleotide identity (%)")
        ax.set_ylim(
            max(0.0, float(data["min_pairwise_nucleotide_identity"].min()) - 2.0),
            100.0,
        )
        ax.set_title("Fragment sequence conservation")
        fig.tight_layout()
        output = plots_dir / "fragment_sequence_conservation.png"
        fig.savefig(output, dpi=200)
        plt.close(fig)
        outputs.append(output)

    if not cds_cluster_summary.empty:
        for fragment, group in cds_cluster_summary.groupby("fragment_number"):
            variable = group.nsmallest(
                top_variable_cds,
                "mean_pairwise_cds_nucleotide_identity",
            ).copy()
            variable = variable.sort_values(
                "mean_pairwise_cds_nucleotide_identity",
                ascending=False,
            )
            labels = [
                f"{row.cluster}: {str(row.products)[:55]}"
                for row in variable.itertuples()
            ]
            divergence = 100.0 - variable["mean_pairwise_cds_nucleotide_identity"]

            fig_height = max(4.0, 0.55 * len(variable) + 1.8)
            fig, ax = plt.subplots(figsize=(10, fig_height))
            y = np.arange(len(variable))
            ax.barh(y, divergence)
            ax.set_yticks(y)
            ax.set_yticklabels(labels, fontsize=8)
            ax.set_xlabel("Mean CDS nucleotide divergence (%)")
            ax.set_title(f"Most variable CDS clusters in Fragment {int(fragment)}")
            fig.tight_layout()
            output = plots_dir / f"Fragment_{int(fragment)}_variable_cds.png"
            fig.savefig(output, dpi=200)
            plt.close(fig)
            outputs.append(output)

    if not assembly_parent_df.empty:
        matrix = assembly_parent_df.pivot(
            index="assembly",
            columns="parent",
            values="nucleotide_identity",
        )
        if len(matrix.index) <= max_heatmap_assemblies:
            output = plots_dir / "assembly_vs_parent_identity_heatmap.png"
            _plot_heatmap(
                matrix,
                output,
                title="Assembly vs parental genome identity",
                annotate_limit=15,
            )
            outputs.append(output)

    chimera_summary = assembly_summary_df[
        assembly_summary_df["assembly_type"] == "chimera"
    ] if not assembly_summary_df.empty else pd.DataFrame()

    if not chimera_summary.empty:
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.hist(chimera_summary["nearest_parent_identity"], bins="auto")
        ax.set_xlabel("Identity to nearest parental genome (%)")
        ax.set_ylabel("Number of chimeras")
        ax.set_title("Nearest-parent identity of chimeric assemblies")
        fig.tight_layout()
        output = plots_dir / "chimera_nearest_parent_identity_distribution.png"
        fig.savefig(output, dpi=200)
        plt.close(fig)
        outputs.append(output)

    if assembly_pairwise_df is not None and not assembly_pairwise_df.empty:
        assemblies = sorted(
            set(assembly_pairwise_df["assembly_a"])
            | set(assembly_pairwise_df["assembly_b"])
        )
        if len(assemblies) <= max_heatmap_assemblies:
            matrix = pd.DataFrame(100.0, index=assemblies, columns=assemblies)
            for row in assembly_pairwise_df.itertuples():
                matrix.loc[row.assembly_a, row.assembly_b] = row.identity
                matrix.loc[row.assembly_b, row.assembly_a] = row.identity
            output = plots_dir / "assembly_pairwise_identity_heatmap.png"
            _plot_heatmap(
                matrix,
                output,
                title="Pairwise assembly nucleotide identity",
                annotate_limit=12,
            )
            outputs.append(output)

    return outputs


def merge_fragment_summaries(
    sequence_summary: pd.DataFrame,
    cds_summary: pd.DataFrame,
) -> pd.DataFrame:
    if sequence_summary.empty:
        return cds_summary.copy()
    if cds_summary.empty:
        return sequence_summary.copy()
    return sequence_summary.merge(
        cds_summary,
        on="fragment_number",
        how="outer",
    ).sort_values("fragment_number")



def analyze_fragment_design(
    materialized_dir: str | Path,
    prepped_dir: str | Path,
    cluster_membership: str | Path,
    outdir: str | Path | None = None,
    save_alignments: bool = False,
    top_variable_cds: int = DEFAULT_TOP_VARIABLE_CDS,
) -> dict:
    """Analyze only the parental fragment set for a working assembly design.

    This is intended for rapid design iteration in the GUI. It calculates
    fragment-level nucleotide identity, homologous CDS nucleotide/protein
    identity, most-variable CDSs, and fragment-focused visualizations without
    generating or analyzing the combinatorial chimera library.
    """
    materialized_dir = Path(materialized_dir).resolve()
    prepped_dir = Path(prepped_dir).resolve()
    cluster_membership = Path(cluster_membership).resolve()
    outdir = (
        Path(outdir).resolve()
        if outdir is not None
        else materialized_dir / "analysis"
    )
    outdir.mkdir(parents=True, exist_ok=True)

    genomes_fasta = prepped_dir / "genomes.fasta"
    if not genomes_fasta.exists():
        raise FileNotFoundError(
            f"Missing prepared genome FASTA: {genomes_fasta}"
        )
    if not cluster_membership.exists():
        raise FileNotFoundError(
            f"Missing cluster membership table: {cluster_membership}"
        )

    genome_sequences = read_genomes_fasta(genomes_fasta)

    fragment_pairwise_df, fragment_sequence_summary, _ = analyze_fragment_fastas(
        materialized_dir,
        outdir,
        save_alignments=save_alignments,
    )

    cds_pairwise_df, cds_cluster_summary, fragment_cds_summary = analyze_cds_diversity(
        prepped_dir,
        cluster_membership,
        materialized_dir,
        genome_sequences,
        outdir,
    )

    fragment_summary = merge_fragment_summaries(
        fragment_sequence_summary,
        fragment_cds_summary,
    )
    fragment_summary.to_csv(
        outdir / "fragment_summary.tsv",
        sep="\t",
        index=False,
    )

    # Reuse the common plotting code while intentionally passing empty
    # assembly tables so only fragment/CDS plots are generated.
    plot_files = make_visualizations(
        outdir,
        fragment_summary,
        cds_cluster_summary,
        pd.DataFrame(),
        None,
        pd.DataFrame(),
        top_variable_cds=top_variable_cds,
    )

    most_variable_fragment = None
    if (
        not fragment_summary.empty
        and "mean_pairwise_nucleotide_identity" in fragment_summary.columns
    ):
        variable_row = fragment_summary.sort_values(
            "mean_pairwise_nucleotide_identity"
        ).iloc[0]
        most_variable_fragment = {
            "fragment_number": int(variable_row["fragment_number"]),
            "mean_pairwise_nucleotide_identity": float(
                variable_row["mean_pairwise_nucleotide_identity"]
            ),
        }

    overview = {
        "analysis_level": "fragments",
        "n_genomes": len(genome_sequences),
        "n_fragments": int(fragment_summary["fragment_number"].nunique())
        if not fragment_summary.empty
        else 0,
        "n_compared_cds_clusters": int(len(cds_cluster_summary)),
        "mean_fragment_pairwise_nucleotide_identity": (
            float(fragment_pairwise_df["nucleotide_identity"].mean())
            if not fragment_pairwise_df.empty
            else 100.0
        ),
        "mean_cds_cluster_nucleotide_identity": (
            float(
                cds_cluster_summary[
                    "mean_pairwise_cds_nucleotide_identity"
                ].mean()
            )
            if not cds_cluster_summary.empty
            else None
        ),
        "mean_cds_cluster_protein_identity": (
            float(
                cds_cluster_summary[
                    "mean_pairwise_protein_identity"
                ].mean()
            )
            if (
                not cds_cluster_summary.empty
                and cds_cluster_summary[
                    "mean_pairwise_protein_identity"
                ].notna().any()
            )
            else None
        ),
        "most_variable_fragment": most_variable_fragment,
    }

    overview_file = outdir / "analysis_summary.json"
    overview_file.write_text(json.dumps(overview, indent=2))

    return {
        "outdir": outdir,
        "overview": overview,
        "overview_file": overview_file,
        "fragment_summary": outdir / "fragment_summary.tsv",
        "fragment_pairwise_identity": outdir / "fragment_pairwise_identity.tsv",
        "cds_pairwise_identity": outdir / "cds_pairwise_identity.tsv",
        "cds_cluster_summary": outdir / "cds_cluster_summary.tsv",
        "plot_files": plot_files,
    }


def analyze_materialized_design(
    materialized_dir: str | Path,
    prepped_dir: str | Path,
    cluster_membership: str | Path,
    outdir: str | Path | None = None,
    save_alignments: bool = False,
    top_variable_cds: int = DEFAULT_TOP_VARIABLE_CDS,
    max_pairwise_assemblies: int = DEFAULT_MAX_PAIRWISE_ASSEMBLIES,
    max_heatmap_assemblies: int = DEFAULT_MAX_HEATMAP_ASSEMBLIES,
) -> dict:
    """GUI-friendly public entry point for assembly-design analysis."""
    materialized_dir = Path(materialized_dir).resolve()
    prepped_dir = Path(prepped_dir).resolve()
    cluster_membership = Path(cluster_membership).resolve()
    outdir = (
        Path(outdir).resolve()
        if outdir is not None
        else materialized_dir / "analysis"
    )
    outdir.mkdir(parents=True, exist_ok=True)

    genomes_fasta = prepped_dir / "genomes.fasta"
    if not genomes_fasta.exists():
        raise FileNotFoundError(
            f"Missing prepared genome FASTA: {genomes_fasta}"
        )
    if not cluster_membership.exists():
        raise FileNotFoundError(
            f"Missing cluster membership table: {cluster_membership}"
        )

    genome_sequences = read_genomes_fasta(genomes_fasta)

    fragment_pairwise_df, fragment_sequence_summary, _ = analyze_fragment_fastas(
        materialized_dir,
        outdir,
        save_alignments=save_alignments,
    )

    cds_pairwise_df, cds_cluster_summary, fragment_cds_summary = analyze_cds_diversity(
        prepped_dir,
        cluster_membership,
        materialized_dir,
        genome_sequences,
        outdir,
    )

    (
        assembly_composition_df,
        assembly_parent_df,
        assembly_pairwise_df,
        assembly_summary_df,
    ) = analyze_assemblies(
        materialized_dir,
        genome_sequences,
        outdir,
        max_pairwise_assemblies=max_pairwise_assemblies,
    )

    fragment_summary = merge_fragment_summaries(
        fragment_sequence_summary,
        fragment_cds_summary,
    )
    fragment_summary.to_csv(
        outdir / "fragment_summary.tsv", sep="\t", index=False
    )

    plot_files = make_visualizations(
        outdir,
        fragment_summary,
        cds_cluster_summary,
        assembly_parent_df,
        assembly_pairwise_df,
        assembly_summary_df,
        top_variable_cds=top_variable_cds,
        max_heatmap_assemblies=max_heatmap_assemblies,
    )

    chimera_summary = assembly_summary_df[
        assembly_summary_df["assembly_type"] == "chimera"
    ] if not assembly_summary_df.empty else pd.DataFrame()

    overview = {
        "n_genomes": len(genome_sequences),
        "n_fragments": int(fragment_summary["fragment_number"].nunique())
        if not fragment_summary.empty
        else 0,
        "n_chimeras": int(len(chimera_summary)),
        "n_compared_cds_clusters": int(len(cds_cluster_summary)),
        "mean_fragment_pairwise_nucleotide_identity": (
            float(fragment_pairwise_df["nucleotide_identity"].mean())
            if not fragment_pairwise_df.empty
            else 100.0
        ),
        "mean_cds_cluster_nucleotide_identity": (
            float(cds_cluster_summary["mean_pairwise_cds_nucleotide_identity"].mean())
            if not cds_cluster_summary.empty
            else None
        ),
        "mean_cds_cluster_protein_identity": (
            float(cds_cluster_summary["mean_pairwise_protein_identity"].mean())
            if (
                not cds_cluster_summary.empty
                and cds_cluster_summary["mean_pairwise_protein_identity"].notna().any()
            )
            else None
        ),
        "mean_chimera_nearest_parent_identity": (
            float(chimera_summary["nearest_parent_identity"].mean())
            if not chimera_summary.empty
            else None
        ),
        "all_vs_all_assembly_identity_calculated": assembly_pairwise_df is not None,
    }

    overview_file = outdir / "analysis_summary.json"
    overview_file.write_text(json.dumps(overview, indent=2))

    return {
        "outdir": outdir,
        "overview": overview,
        "overview_file": overview_file,
        "fragment_summary": outdir / "fragment_summary.tsv",
        "fragment_pairwise_identity": outdir / "fragment_pairwise_identity.tsv",
        "cds_pairwise_identity": outdir / "cds_pairwise_identity.tsv",
        "cds_cluster_summary": outdir / "cds_cluster_summary.tsv",
        "assembly_composition": outdir / "assembly_composition.tsv",
        "assembly_vs_parent_identity": outdir / "assembly_vs_parent_identity.tsv",
        "assembly_summary": outdir / "assembly_summary.tsv",
        "assembly_pairwise_identity": (
            outdir / "assembly_pairwise_identity.tsv"
            if assembly_pairwise_df is not None
            else None
        ),
        "plot_files": plot_files,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Analyze nucleotide, CDS, protein, and assembly diversity for a "
            "materialized CombPhage assembly design."
        )
    )
    parser.add_argument(
        "--assembly-dir",
        required=True,
        help=(
            "Materialized design directory containing fragments/, "
            "fragment_manifest.tsv, and optionally chimeras/."
        ),
    )
    parser.add_argument(
        "--prepped-dir",
        required=True,
        help="Prepared CombPhage genome directory containing genomes.fasta and genbanks/.",
    )
    parser.add_argument(
        "--cluster-membership",
        required=True,
        help="Path to feature_overlaps/cluster_membership.tsv",
    )
    parser.add_argument(
        "--outdir",
        default=None,
        help="Analysis output directory (default: <assembly-dir>/analysis)",
    )
    parser.add_argument(
        "--save-alignments",
        action="store_true",
        help="Save individual pairwise fragment alignments in FASTA format.",
    )
    parser.add_argument(
        "--fragments-only",
        action="store_true",
        help=(
            "Analyze only parental fragment/CDS diversity and skip assembly "
            "comparisons. Useful while iterating on an assembly design."
        ),
    )
    parser.add_argument(
        "--top-variable-cds",
        type=int,
        default=DEFAULT_TOP_VARIABLE_CDS,
        help=(
            "Number of most variable CDS clusters to plot per fragment "
            f"(default: {DEFAULT_TOP_VARIABLE_CDS})."
        ),
    )
    parser.add_argument(
        "--max-pairwise-assemblies",
        type=int,
        default=DEFAULT_MAX_PAIRWISE_ASSEMBLIES,
        help=(
            "Maximum total assemblies for full all-vs-all identity calculation "
            f"(default: {DEFAULT_MAX_PAIRWISE_ASSEMBLIES})."
        ),
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.fragments_only:
        result = analyze_fragment_design(
            materialized_dir=args.assembly_dir,
            prepped_dir=args.prepped_dir,
            cluster_membership=args.cluster_membership,
            outdir=args.outdir,
            save_alignments=args.save_alignments,
            top_variable_cds=args.top_variable_cds,
        )
    else:
        result = analyze_materialized_design(
            materialized_dir=args.assembly_dir,
            prepped_dir=args.prepped_dir,
            cluster_membership=args.cluster_membership,
            outdir=args.outdir,
            save_alignments=args.save_alignments,
            top_variable_cds=args.top_variable_cds,
            max_pairwise_assemblies=args.max_pairwise_assemblies,
        )

    print(f"Analysis output: {result['outdir']}")
    print(json.dumps(result["overview"], indent=2))
    for plot in result["plot_files"]:
        print(f"Plot: {plot}")


if __name__ == "__main__":
    main()
