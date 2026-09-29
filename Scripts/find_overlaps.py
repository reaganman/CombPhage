#!/usr/bin/env python3

import argparse
import shutil
import subprocess
from collections import defaultdict
from itertools import combinations, product
from pathlib import Path

import pandas as pd
from Bio import SeqIO


COMPLEMENT = str.maketrans(
    "ACGTRYMKWSBDHVNacgtrymkwsbdhvn",
    "TGCAYRKMWSVHDBNtgcayrkmwsvhdbn",
)


def run_command(command):
    """Run an external command and stop immediately if it fails."""
    print("[cmd]", " ".join(map(str, command)))
    subprocess.run([str(x) for x in command], check=True)


def require_mmseqs(mmseqs_exe):
    """Ensure that the MMseqs2 executable is available."""
    if shutil.which(mmseqs_exe) is None:
        raise RuntimeError(
            f"Could not find MMseqs2 executable '{mmseqs_exe}' in PATH. "
            "Install/activate MMseqs2 or provide --mmseqs-exe."
        )


def mmseqs_db_exists(db_prefix):
    """Return True when an MMseqs2 database prefix appears to exist."""
    db_prefix = Path(db_prefix)
    return Path(str(db_prefix) + ".dbtype").exists()


def ensure_mmseqs_database(proteins_faa, mmseqs_dir, mmseqs_exe="mmseqs"):
    """Create the MMseqs2 protein database if it does not already exist."""
    mmseqs_dir = Path(mmseqs_dir)
    db_dir = mmseqs_dir / "DB"
    db_dir.mkdir(parents=True, exist_ok=True)

    protein_db = db_dir / "proteomesDB"

    if not mmseqs_db_exists(protein_db):
        print("[+] Creating MMseqs2 protein database")
        run_command([
            mmseqs_exe,
            "createdb",
            proteins_faa,
            protein_db,
        ])

    return protein_db


def run_mmseqs_clustering(
    proteins_faa,
    mmseqs_dir,
    mmseqs_exe="mmseqs",
    min_seq_id=0.25,
    coverage=0.8,
    sensitivity=7.5,
    cluster_mode=0,
    cov_mode=1,
):
    """Cluster proteins with MMseqs2 and return the path to clusters.tsv."""
    proteins_faa = Path(proteins_faa)
    mmseqs_dir = Path(mmseqs_dir)
    clusters_tsv = mmseqs_dir / "clusters.tsv"

    require_mmseqs(mmseqs_exe)

    if clusters_tsv.exists():
        print(f"[+] Reusing existing MMseqs2 clusters: {clusters_tsv}")
        return clusters_tsv

    results_dir = mmseqs_dir / "Results"
    tmp_dir = mmseqs_dir / "tmp"
    results_dir.mkdir(parents=True, exist_ok=True)
    tmp_dir.mkdir(parents=True, exist_ok=True)

    protein_db = ensure_mmseqs_database(
        proteins_faa,
        mmseqs_dir,
        mmseqs_exe=mmseqs_exe,
    )
    cluster_db = results_dir / "clusters"

    print("[+] Clustering proteins with MMseqs2")

    run_command([
        mmseqs_exe,
        "cluster",
        protein_db,
        cluster_db,
        tmp_dir,
        "--cluster-mode", str(cluster_mode),
        "--cov-mode", str(cov_mode),
        "-c", str(coverage),
        "-s", str(sensitivity),
        "--min-seq-id", str(min_seq_id),
    ])

    run_command([
        mmseqs_exe,
        "createtsv",
        protein_db,
        protein_db,
        cluster_db,
        clusters_tsv,
    ])

    return clusters_tsv


def run_mmseqs_all_vs_all(
    proteins_faa,
    mmseqs_dir,
    mmseqs_exe="mmseqs",
    sensitivity=7.5,
    max_seqs=10000,
):
    """
    Run an all-vs-all protein search and write pairwise similarities.

    Output columns:
        query, target, pident, alnlen, qcov, tcov, evalue
    """
    proteins_faa = Path(proteins_faa)
    mmseqs_dir = Path(mmseqs_dir)
    allvsall_tsv = mmseqs_dir / "allvsall_pident.tsv"

    require_mmseqs(mmseqs_exe)

    if allvsall_tsv.exists():
        print(f"[+] Reusing existing MMseqs2 all-vs-all output: {allvsall_tsv}")
        return allvsall_tsv

    results_dir = mmseqs_dir / "Results"
    tmp_dir = mmseqs_dir / "tmp"
    results_dir.mkdir(parents=True, exist_ok=True)
    tmp_dir.mkdir(parents=True, exist_ok=True)

    protein_db = ensure_mmseqs_database(
        proteins_faa,
        mmseqs_dir,
        mmseqs_exe=mmseqs_exe,
    )
    allvsall_db = results_dir / "allvsall"

    print("[+] Running MMseqs2 all-vs-all protein search")

    run_command([
        mmseqs_exe,
        "search",
        protein_db,
        protein_db,
        allvsall_db,
        tmp_dir,
        "--search-type", "1",
        "-s", str(sensitivity),
        "-a",
        "--max-seqs", str(max_seqs),
    ])

    run_command([
        mmseqs_exe,
        "convertalis",
        protein_db,
        protein_db,
        allvsall_db,
        allvsall_tsv,
        "--format-output",
        "query,target,pident,alnlen,qcov,tcov,evalue",
    ])

    return allvsall_tsv


def load_prepared_inputs(input_dir):
    """Load the standardized files produced by prepare_genbanks.py."""
    input_dir = Path(input_dir)
    genomes_path = input_dir / "genomes.fasta"
    proteins_path = input_dir / "proteins.faa"
    cds_path = input_dir / "cds.tsv"

    missing = [p for p in (genomes_path, proteins_path, cds_path) if not p.exists()]
    if missing:
        missing_text = "\n".join(f"  {p}" for p in missing)
        raise FileNotFoundError(
            "The prepared input directory is missing required file(s):\n"
            f"{missing_text}\n"
            "Expected output from prepare_genbanks.py."
        )

    genomes = list(SeqIO.parse(genomes_path, "fasta"))
    if not genomes:
        raise ValueError(f"No genome sequences found in {genomes_path}")

    genome_ids = [rec.id for rec in genomes]
    if len(genome_ids) != len(set(genome_ids)):
        raise ValueError("genomes.fasta contains duplicate genome IDs.")

    cds_df = pd.read_csv(cds_path, sep="\t", dtype=str)
    required_columns = {
        "gene",
        "genome",
        "start",
        "stop",
        "frame",
        "cds_index",
    }
    missing_columns = required_columns - set(cds_df.columns)
    if missing_columns:
        raise ValueError(
            "cds.tsv is missing required columns: "
            + ", ".join(sorted(missing_columns))
        )

    if cds_df["gene"].duplicated().any():
        duplicates = cds_df.loc[cds_df["gene"].duplicated(), "gene"].tolist()
        raise ValueError(
            "cds.tsv contains duplicate gene IDs, including: "
            + ", ".join(duplicates[:5])
        )

    unknown_genomes = sorted(set(cds_df["genome"]) - set(genome_ids))
    if unknown_genomes:
        raise ValueError(
            "cds.tsv refers to genome IDs not present in genomes.fasta: "
            + ", ".join(unknown_genomes)
        )

    bad_frames = sorted(set(cds_df["frame"]) - {"+", "-"})
    if bad_frames:
        raise ValueError(
            "cds.tsv contains invalid frame/strand values: "
            + ", ".join(bad_frames)
        )

    return genomes, proteins_path, cds_df


def load_mmseqs_clusters(clusters_path):
    """Read MMseqs2 createtsv output into representative -> member sets."""
    clusters = defaultdict(set)

    with open(clusters_path) as handle:
        for line_no, line in enumerate(handle, start=1):
            line = line.rstrip("\n")
            if not line:
                continue

            fields = line.split("\t")
            if len(fields) < 2:
                raise ValueError(
                    f"Malformed MMseqs2 cluster file at line {line_no}: {line}"
                )

            rep, mem = fields[0], fields[1]
            clusters[rep].add(rep)
            clusters[rep].add(mem)

    return clusters


def diagnostic_rows(cluster_df, representative, classification, genome_ids, min_genomes):
    """Create one diagnostic row per CDS in an ambiguous/incomplete cluster."""
    genome_counts = cluster_df["genome"].value_counts()
    represented = set(genome_counts.index)
    missing_genomes = [g for g in genome_ids if g not in represented]
    duplicated_genomes = sorted(
        genome for genome, count in genome_counts.items() if count > 1
    )

    metadata = {
        "mmseqs_representative": representative,
        "classification": classification,
        "member_count": len(cluster_df),
        "distinct_genomes": int(genome_counts.size),
        "required_genomes": min_genomes,
        "missing_genomes": ";".join(missing_genomes),
        "duplicated_genomes": ";".join(duplicated_genomes),
        "has_within_genome_multiplicity": bool(duplicated_genomes),
    }

    rows = []
    for _, mem in cluster_df.iterrows():
        row = metadata.copy()
        row.update(mem.to_dict())
        rows.append(row)

    return rows


def classify_clusters(clusters, cds_df, genome_ids, min_genomes, reference_genome):
    """
    Classify every MMseqs2 protein family.

    Direct:
      * at least min_genomes distinct genomes
      * no genome contributes more than one CDS

    Ambiguous:
      * at least min_genomes distinct genomes
      * one or more genomes contribute multiple CDSs

    Incomplete:
      * fewer than min_genomes distinct genomes

    Unclustered CDSs are CDSs that never appear in clusters.tsv at all.

    Returns the actual direct/ambiguous cluster DataFrames in addition to the
    diagnostic tables so ambiguous families can subsequently be resolved with
    the MMseqs2 all-vs-all alignments.
    """
    cds_by_gene = cds_df.set_index("gene", drop=False)
    direct_clusters = []
    ambiguous_clusters = []
    ambiguous_rows = []
    incomplete_rows = []
    clustered_genes = set()

    for representative, members in clusters.items():
        missing_members = sorted(
            gene for gene in members if gene not in cds_by_gene.index
        )
        if missing_members:
            raise ValueError(
                f"MMseqs2 cluster '{representative}' contains gene IDs absent "
                "from cds.tsv, including: " + ", ".join(missing_members[:5])
            )

        cluster_df = cds_by_gene.loc[sorted(members)].copy()
        clustered_genes.update(cluster_df["gene"].tolist())

        genome_counts = cluster_df["genome"].value_counts()
        distinct_genomes = int(genome_counts.size)
        has_multiplicity = bool((genome_counts > 1).any())

        if distinct_genomes < min_genomes:
            incomplete_rows.extend(
                diagnostic_rows(
                    cluster_df,
                    representative,
                    "incomplete",
                    genome_ids,
                    min_genomes,
                )
            )
        elif has_multiplicity:
            ambiguous_clusters.append((representative, cluster_df))
            ambiguous_rows.extend(
                diagnostic_rows(
                    cluster_df,
                    representative,
                    "ambiguous",
                    genome_ids,
                    min_genomes,
                )
            )
        else:
            direct_clusters.append((representative, cluster_df))

    def cluster_sort_key(item):
        _, cluster_df = item
        ref_rows = cluster_df[cluster_df["genome"] == reference_genome]
        if not ref_rows.empty:
            return int(ref_rows["cds_index"].iloc[0])
        return int(cluster_df["cds_index"].astype(int).min())

    direct_clusters.sort(key=cluster_sort_key)

    ambiguous_df = pd.DataFrame(ambiguous_rows)
    incomplete_df = pd.DataFrame(incomplete_rows)

    unclustered_df = cds_df.loc[
        ~cds_df["gene"].isin(clustered_genes)
    ].copy()
    if not unclustered_df.empty:
        unclustered_df.insert(0, "classification", "unclustered")

    print(f"[+] MMseqs2 clusters read: {len(clusters)}")
    print(f"[+] Direct one-CDS-per-genome clusters: {len(direct_clusters)}")
    print(f"[+] Ambiguous clusters: {len(ambiguous_clusters)}")
    print(
        "[+] Incomplete clusters: "
        f"{incomplete_df['mmseqs_representative'].nunique() if not incomplete_df.empty else 0}"
    )
    print(f"[+] Unclustered CDSs: {len(unclustered_df)}")

    return (
        direct_clusters,
        ambiguous_clusters,
        ambiguous_df,
        incomplete_df,
        unclustered_df,
    )


def load_allvsall_similarity(allvsall_path):
    """
    Load MMseqs2 all-vs-all output and retain the best alignment for each
    unordered protein pair.

    The main resolution metric is balanced coverage:
        min(qcov, tcov)

    pident is retained as a secondary tie-breaker. MMseqs2 qcov/tcov are
    normally fractions (0-1), but percentage-like values are normalized if
    encountered.
    """
    columns = ["query", "target", "pident", "alnlen", "qcov", "tcov", "evalue"]
    df = pd.read_csv(allvsall_path, sep="\t", names=columns, dtype={"query": str, "target": str})

    if df.empty:
        raise ValueError(f"MMseqs2 all-vs-all output is empty: {allvsall_path}")

    for col in ("pident", "alnlen", "qcov", "tcov", "evalue"):
        df[col] = pd.to_numeric(df[col], errors="coerce")

    if df[["pident", "qcov", "tcov"]].isna().any().any():
        raise ValueError(
            "MMseqs2 all-vs-all output contains non-numeric pident/qcov/tcov values."
        )

    pair_metrics = {}

    for row in df.itertuples(index=False):
        if row.query == row.target:
            continue

        qcov = float(row.qcov)
        tcov = float(row.tcov)
        if qcov > 1.0:
            qcov /= 100.0
        if tcov > 1.0:
            tcov /= 100.0

        pident = float(row.pident)
        if pident > 1.0:
            pident /= 100.0

        balanced_cov = min(qcov, tcov)
        key = tuple(sorted((row.query, row.target)))
        metrics = {
            "balanced_coverage": balanced_cov,
            "qcov": qcov,
            "tcov": tcov,
            "pident": pident,
            "evalue": float(row.evalue) if pd.notna(row.evalue) else float("inf"),
        }

        previous = pair_metrics.get(key)
        if previous is None or (
            metrics["balanced_coverage"], metrics["pident"]
        ) > (
            previous["balanced_coverage"], previous["pident"]
        ):
            pair_metrics[key] = metrics

    return pair_metrics


def pair_metrics_for_genes(gene_a, gene_b, pair_metrics):
    """Return pairwise similarity metrics for two different CDS proteins."""
    if gene_a == gene_b:
        return {
            "balanced_coverage": 1.0,
            "qcov": 1.0,
            "tcov": 1.0,
            "pident": 1.0,
            "evalue": 0.0,
        }
    return pair_metrics.get(tuple(sorted((gene_a, gene_b))))


def score_resolution_combination(genes, pair_metrics, min_coverage):
    """
    Score a one-CDS-per-genome candidate using all cross-genome pairs.

    Every pair must be present in the all-vs-all table and satisfy the minimum
    balanced coverage. The primary score is mean balanced coverage. Mean
    identity is used only as a secondary tie-breaker.
    """
    pairwise = []

    for gene_a, gene_b in combinations(genes, 2):
        metrics = pair_metrics_for_genes(gene_a, gene_b, pair_metrics)
        if metrics is None:
            return None

        if metrics["balanced_coverage"] < min_coverage:
            return None

        pairwise.append(metrics)

    if not pairwise:
        return None

    coverages = [m["balanced_coverage"] for m in pairwise]
    identities = [m["pident"] for m in pairwise]

    return {
        "mean_balanced_coverage": sum(coverages) / len(coverages),
        "min_balanced_coverage": min(coverages),
        "mean_pident": sum(identities) / len(identities),
        "pair_count": len(pairwise),
    }


def resolve_ambiguous_cluster(
    representative,
    cluster_df,
    genome_ids,
    pair_metrics,
    min_coverage=0.8,
    max_combinations=100000,
):
    """
    Split an ambiguous MMseqs2 family into one or more one-CDS-per-genome
    subclusters using pairwise protein coverage.

    Resolution strategy:
      1. Build every possible combination containing exactly one CDS from each
         genome represented in the original MMseqs family.
      2. Reject combinations where any protein pair has balanced coverage
         min(qcov, tcov) below min_coverage.
      3. Rank valid combinations primarily by mean balanced coverage, then by
         minimum pairwise coverage and mean identity.
      4. Select the best combination, mark all of its CDSs as used, and repeat
         on the remaining CDSs.

    A CDS can therefore participate in at most one resolved subcluster.
    """
    genome_order = [g for g in genome_ids if g in set(cluster_df["genome"])]
    members_by_genome = {
        genome: cluster_df.loc[cluster_df["genome"] == genome, "gene"].tolist()
        for genome in genome_order
    }

    unused = set(cluster_df["gene"])
    resolved = []
    stop_reason = "resolved_all_possible"

    while all(any(g in unused for g in members_by_genome[genome]) for genome in genome_order):
        candidate_lists = [
            [gene for gene in members_by_genome[genome] if gene in unused]
            for genome in genome_order
        ]

        combination_count = 1
        for candidates in candidate_lists:
            combination_count *= len(candidates)

        if combination_count > max_combinations:
            stop_reason = f"too_many_combinations:{combination_count}"
            break

        valid_candidates = []
        for genes in product(*candidate_lists):
            score = score_resolution_combination(
                genes,
                pair_metrics,
                min_coverage=min_coverage,
            )
            if score is None:
                continue

            valid_candidates.append((genes, score))

        if not valid_candidates:
            stop_reason = "no_valid_full_coverage_combination"
            break

        valid_candidates.sort(
            key=lambda item: (
                item[1]["mean_balanced_coverage"],
                item[1]["min_balanced_coverage"],
                item[1]["mean_pident"],
            ),
            reverse=True,
        )

        best_genes, best_score = valid_candidates[0]
        best_df = cluster_df[cluster_df["gene"].isin(best_genes)].copy()

        resolved.append({
            "source_mmseqs_representative": representative,
            "resolution_index": len(resolved) + 1,
            "cluster_df": best_df,
            **best_score,
        })

        unused.difference_update(best_genes)

    unresolved_df = cluster_df[cluster_df["gene"].isin(sorted(unused))].copy()
    return resolved, unresolved_df, stop_reason


def resolve_ambiguous_clusters(
    ambiguous_clusters,
    genome_ids,
    pair_metrics,
    min_coverage=0.8,
    max_combinations=100000,
):
    """Resolve all ambiguous MMseqs2 families and produce audit tables."""
    resolved_groups = []
    resolved_rows = []
    unresolved_rows = []
    summary_rows = []

    for representative, cluster_df in ambiguous_clusters:
        resolved, unresolved_df, stop_reason = resolve_ambiguous_cluster(
            representative=representative,
            cluster_df=cluster_df,
            genome_ids=genome_ids,
            pair_metrics=pair_metrics,
            min_coverage=min_coverage,
            max_combinations=max_combinations,
        )

        for group in resolved:
            resolved_groups.append(group)
            for _, member in group["cluster_df"].iterrows():
                row = {
                    "source_mmseqs_representative": representative,
                    "resolution_index": group["resolution_index"],
                    "mean_balanced_coverage": group["mean_balanced_coverage"],
                    "min_balanced_coverage": group["min_balanced_coverage"],
                    "mean_pident": group["mean_pident"],
                }
                row.update(member.to_dict())
                resolved_rows.append(row)

        for _, member in unresolved_df.iterrows():
            row = {
                "source_mmseqs_representative": representative,
                "classification": "unresolved_ambiguous_member",
                "resolution_stop_reason": stop_reason,
            }
            row.update(member.to_dict())
            unresolved_rows.append(row)

        summary_rows.append({
            "source_mmseqs_representative": representative,
            "original_member_count": len(cluster_df),
            "distinct_genomes": cluster_df["genome"].nunique(),
            "resolved_subclusters": len(resolved),
            "resolved_members": sum(len(group["cluster_df"]) for group in resolved),
            "unresolved_members": len(unresolved_df),
            "resolution_stop_reason": stop_reason,
        })

    resolved_df = pd.DataFrame(resolved_rows)
    unresolved_df = pd.DataFrame(unresolved_rows)
    summary_df = pd.DataFrame(summary_rows)

    print(f"[+] Resolved subclusters recovered: {len(resolved_groups)}")
    print(f"[+] Unresolved ambiguous CDSs remaining: {len(unresolved_df)}")

    return resolved_groups, resolved_df, unresolved_df, summary_df

def position_at_offset(anchor_coord, strand, sense, offset):
    """Return zero-based genomic position at a gene-relative offset."""
    if strand == "+":
        return anchor_coord - offset if sense == "up" else anchor_coord + offset
    return anchor_coord + offset if sense == "up" else anchor_coord - offset


def oriented_base(seq, position, strand):
    """Return the nucleotide in CDS-oriented 5'->3' space."""
    base = str(seq[position])
    return base if strand == "+" else base.translate(COMPLEMENT)


def shared_overlap(cluster_df, genome_dict, ref_point, sense):
    """
    Find the identical nucleotide run shared by every cluster member around a
    CDS start or stop boundary.
    """
    base_mem = cluster_df.iloc[0]
    base_seq = genome_dict[base_mem["genome"]]
    base_anchor = int(base_mem[ref_point]) - 1
    base_strand = base_mem["frame"]

    overlap_len = 0

    while True:
        base_pos = position_at_offset(
            base_anchor,
            base_strand,
            sense,
            overlap_len,
        )

        if base_pos < 0 or base_pos >= len(base_seq):
            break

        base_nt = oriented_base(base_seq, base_pos, base_strand)
        all_match = True

        for _, mem in cluster_df.iterrows():
            seq = genome_dict[mem["genome"]]
            anchor = int(mem[ref_point]) - 1
            strand = mem["frame"]
            pos = position_at_offset(anchor, strand, sense, overlap_len)

            if pos < 0 or pos >= len(seq):
                all_match = False
                break

            if oriented_base(seq, pos, strand) != base_nt:
                all_match = False
                break

        if not all_match:
            break

        overlap_len += 1

    return overlap_len


def overlap_bounds(anchor, strand, up_len, down_len):
    """Return zero-based inclusive genomic bounds of a merged overlap."""
    if up_len == 0 and down_len == 0:
        return None

    if strand == "+":
        low = anchor - max(up_len - 1, 0)
        high = anchor + max(down_len - 1, 0)
    else:
        low = anchor - max(down_len - 1, 0)
        high = anchor + max(up_len - 1, 0)

    return low, high


def find_overlaps(
    cluster_df,
    genomes,
    cluster_id,
    reference_genome,
    resolution="direct",
    source_mmseqs_representative="",
):
    """Find shared start/stop overlap regions for one direct/resolved cluster."""
    genome_dict = {rec.id: rec.seq for rec in genomes}

    ref_rows = cluster_df[cluster_df["genome"] == reference_genome]
    if ref_rows.empty:
        base_mem = cluster_df.iloc[0]
    else:
        base_mem = ref_rows.iloc[0]

    rows = []

    for region in ("start", "stop"):
        up_len = shared_overlap(cluster_df, genome_dict, region, "up")
        down_len = shared_overlap(cluster_df, genome_dict, region, "down")

        if up_len == 0 and down_len == 0:
            continue

        base_seq = genome_dict[base_mem["genome"]]
        base_anchor = int(base_mem[region]) - 1
        base_bounds = overlap_bounds(
            base_anchor,
            base_mem["frame"],
            up_len,
            down_len,
        )

        if base_bounds is None:
            continue

        base_low, base_high = base_bounds
        base_low = max(0, base_low)
        base_high = min(len(base_seq) - 1, base_high)

        ov_seq = str(base_seq[base_low:base_high + 1])

        row_data = {
            "cluster": cluster_id,
            "resolution": resolution,
            "source_mmseqs_representative": source_mmseqs_representative,
            "region": region,
            "reference_genome": base_mem["genome"],
            "reference_gene": base_mem["gene"],
            "reference_frame": base_mem["frame"],
            "ov_seq": ov_seq,
            "overlap_length": len(ov_seq),
            "upstream_match_length": up_len,
            "downstream_match_length": down_len,
        }

        for _, mem in cluster_df.iterrows():
            seq = genome_dict[mem["genome"]]
            anchor = int(mem[region]) - 1
            bounds = overlap_bounds(
                anchor,
                mem["frame"],
                up_len,
                down_len,
            )

            if bounds is None:
                continue

            low, high = bounds
            low = max(0, low)
            high = min(len(seq) - 1, high)

            genome_id = mem["genome"]
            row_data[f"{genome_id}_frame"] = mem["frame"]
            row_data[f"{genome_id}_ov_start"] = low + 1
            row_data[f"{genome_id}_ov_stop"] = high + 1

        rows.append(row_data)

    return pd.DataFrame(rows)


def write_diagnostic_table(df, path, cds_columns):
    """Write a diagnostic TSV with stable headers even when it is empty."""
    path = Path(path)

    if df.empty:
        if path.name == "unclustered_cds.tsv":
            columns = ["classification"] + list(cds_columns)
        else:
            columns = [
                "mmseqs_representative",
                "classification",
                "member_count",
                "distinct_genomes",
                "required_genomes",
                "missing_genomes",
                "duplicated_genomes",
                "has_within_genome_multiplicity",
            ] + list(cds_columns)
        pd.DataFrame(columns=columns).to_csv(path, sep="\t", index=False)
    else:
        df.to_csv(path, sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Cluster proteins from prepare_genbanks.py with MMseqs2 and identify "
            "shared nucleotide overlaps around homologous CDS start/stop regions."
        )
    )

    parser.add_argument(
        "--input-dir",
        required=True,
        help=(
            "Directory produced by prepare_genbanks.py containing genomes.fasta, "
            "proteins.faa, and cds.tsv"
        ),
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory for MMseqs2 and overlap results",
    )
    parser.add_argument(
        "--clusters",
        help=(
            "Existing MMseqs2 clusters.tsv. If omitted, this script runs MMseqs2 "
            "on input-dir/proteins.faa. The all-vs-all search is still generated."
        ),
    )
    parser.add_argument(
        "--cluster-threshold",
        type=int,
        help=(
            "Minimum number of distinct genomes represented in a protein cluster. "
            "Default: all input genomes."
        ),
    )
    parser.add_argument(
        "--reference-genome",
        help=(
            "Genome ID used to order clusters and report the representative overlap "
            "sequence. Default: first record in genomes.fasta."
        ),
    )

    # MMseqs2 options
    parser.add_argument("--mmseqs-exe", default="mmseqs")
    parser.add_argument("--min-seq-id", type=float, default=0.25)
    parser.add_argument("--coverage", type=float, default=0.8)
    parser.add_argument("--sensitivity", type=float, default=7.5)
    parser.add_argument("--cluster-mode", type=int, default=0)
    parser.add_argument("--cov-mode", type=int, default=1)
    parser.add_argument("--max-seqs", type=int, default=10000)
    parser.add_argument(
        "--resolve-min-coverage",
        type=float,
        default=0.8,
        help=(
            "Minimum balanced pairwise coverage min(qcov, tcov) required to "
            "place CDSs into a resolved ambiguous subcluster. Default: 0.8"
        ),
    )
    parser.add_argument(
        "--max-resolution-combinations",
        type=int,
        default=100000,
        help=(
            "Maximum one-per-genome combinations evaluated in one resolution "
            "round for an ambiguous MMseqs family. Default: 100000"
        ),
    )
    parser.add_argument(
        "--rerun-mmseqs",
        action="store_true",
        help=(
            "Delete and regenerate MMseqs2 outputs created under OUTDIR/mmseqs. "
            "If --clusters points inside that directory, use a fresh output "
            "directory instead."
        ),
    )

    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    genomes, proteins_faa, cds_df = load_prepared_inputs(input_dir)
    genome_ids = [rec.id for rec in genomes]

    reference_genome = args.reference_genome or genome_ids[0]
    if reference_genome not in genome_ids:
        parser.error(
            f"--reference-genome '{reference_genome}' is not present in genomes.fasta"
        )

    min_genomes = args.cluster_threshold or len(genome_ids)
    if min_genomes < 1 or min_genomes > len(genome_ids):
        parser.error(
            f"--cluster-threshold must be between 1 and {len(genome_ids)}"
        )

    if not 0.0 <= args.resolve_min_coverage <= 1.0:
        parser.error("--resolve-min-coverage must be between 0 and 1")

    if args.max_resolution_combinations < 1:
        parser.error("--max-resolution-combinations must be at least 1")

    print(f"[+] Loaded {len(genomes)} genomes")
    print(f"[+] Loaded {len(cds_df)} CDS annotations")
    print(f"[+] Reference genome: {reference_genome}")
    print(f"[+] Required distinct genomes per cluster: {min_genomes}")

    mmseqs_dir = outdir / "mmseqs"

    if args.rerun_mmseqs and mmseqs_dir.exists():
        if args.clusters:
            cluster_resolved = Path(args.clusters).resolve()
            mmseqs_resolved = mmseqs_dir.resolve()
            if mmseqs_resolved == cluster_resolved.parent or mmseqs_resolved in cluster_resolved.parents:
                parser.error(
                    "--rerun-mmseqs would delete the --clusters file because it is "
                    "inside OUTDIR/mmseqs. Use a fresh --outdir or omit "
                    "--rerun-mmseqs."
                )
        shutil.rmtree(mmseqs_dir)

    if args.clusters:
        clusters_path = Path(args.clusters)
        if not clusters_path.exists():
            parser.error(f"Cluster file does not exist: {clusters_path}")
    else:
        clusters_path = run_mmseqs_clustering(
            proteins_faa=proteins_faa,
            mmseqs_dir=mmseqs_dir,
            mmseqs_exe=args.mmseqs_exe,
            min_seq_id=args.min_seq_id,
            coverage=args.coverage,
            sensitivity=args.sensitivity,
            cluster_mode=args.cluster_mode,
            cov_mode=args.cov_mode,
        )

    allvsall_path = run_mmseqs_all_vs_all(
        proteins_faa=proteins_faa,
        mmseqs_dir=mmseqs_dir,
        mmseqs_exe=args.mmseqs_exe,
        sensitivity=args.sensitivity,
        max_seqs=args.max_seqs,
    )

    raw_clusters = load_mmseqs_clusters(clusters_path)
    (
        direct_clusters,
        ambiguous_clusters,
        ambiguous_df,
        incomplete_df,
        unclustered_df,
    ) = classify_clusters(
        clusters=raw_clusters,
        cds_df=cds_df,
        genome_ids=genome_ids,
        min_genomes=min_genomes,
        reference_genome=reference_genome,
    )

    pair_metrics = load_allvsall_similarity(allvsall_path)
    (
        resolved_groups,
        resolved_df,
        unresolved_ambiguous_df,
        resolution_summary_df,
    ) = resolve_ambiguous_clusters(
        ambiguous_clusters=ambiguous_clusters,
        genome_ids=genome_ids,
        pair_metrics=pair_metrics,
        min_coverage=args.resolve_min_coverage,
        max_combinations=args.max_resolution_combinations,
    )

    # Combine direct and successfully resolved groups, then order them along
    # the reference genome before assigning final CombPhage cluster names.
    usable_clusters = []

    for representative, cluster_df in direct_clusters:
        usable_clusters.append({
            "cluster_df": cluster_df,
            "resolution": "direct",
            "source_mmseqs_representative": representative,
            "resolution_index": "",
            "mean_balanced_coverage": "",
            "min_balanced_coverage": "",
            "mean_pident": "",
        })

    for group in resolved_groups:
        usable_clusters.append({
            "cluster_df": group["cluster_df"],
            "resolution": "resolved",
            "source_mmseqs_representative": group["source_mmseqs_representative"],
            "resolution_index": group["resolution_index"],
            "mean_balanced_coverage": group["mean_balanced_coverage"],
            "min_balanced_coverage": group["min_balanced_coverage"],
            "mean_pident": group["mean_pident"],
        })

    def usable_sort_key(group):
        cluster_df = group["cluster_df"]
        ref_rows = cluster_df[cluster_df["genome"] == reference_genome]
        if not ref_rows.empty:
            return int(ref_rows["cds_index"].iloc[0])
        return int(cluster_df["cds_index"].astype(int).min())

    usable_clusters.sort(key=usable_sort_key)

    named_clusters = [
        (f"cluster{i}", group)
        for i, group in enumerate(usable_clusters, start=1)
    ]

    # MMseqs families partition proteins, and the resolver additionally removes
    # each CDS after assignment. This assertion protects the one-CDS-one-final-
    # cluster invariant against future code changes.
    final_genes = []
    for _, group in named_clusters:
        final_genes.extend(group["cluster_df"]["gene"].tolist())
    if len(final_genes) != len(set(final_genes)):
        raise RuntimeError(
            "Internal error: at least one CDS was assigned to multiple final clusters."
        )

    membership_rows = []
    overlap_frames = []

    for cluster_id, group in named_clusters:
        cluster_df = group["cluster_df"]
        for _, mem in cluster_df.iterrows():
            membership_rows.append({
                "cluster": cluster_id,
                "resolution": group["resolution"],
                "source_mmseqs_representative": group["source_mmseqs_representative"],
                "resolution_index": group["resolution_index"],
                "mean_balanced_coverage": group["mean_balanced_coverage"],
                "min_balanced_coverage": group["min_balanced_coverage"],
                "mean_pident": group["mean_pident"],
                "gene": mem["gene"],
                "genome": mem["genome"],
                "cds_index": int(mem["cds_index"]),
                "frame": mem["frame"],
            })

        overlap_df = find_overlaps(
            cluster_df,
            genomes,
            cluster_id=cluster_id,
            reference_genome=reference_genome,
            resolution=group["resolution"],
            source_mmseqs_representative=group["source_mmseqs_representative"],
        )
        if not overlap_df.empty:
            overlap_frames.append(overlap_df)

    membership_path = outdir / "cluster_membership.tsv"
    membership_df = pd.DataFrame(
        membership_rows,
        columns=[
            "cluster", "resolution", "source_mmseqs_representative",
            "resolution_index", "mean_balanced_coverage",
            "min_balanced_coverage", "mean_pident", "gene", "genome",
            "cds_index", "frame"
        ],
    )
    membership_df.to_csv(membership_path, sep="\t", index=False)

    overlaps_path = outdir / "cluster_overlap_summary.tsv"
    if overlap_frames:
        final_overlap_df = pd.concat(overlap_frames, ignore_index=True)
    else:
        final_overlap_df = pd.DataFrame(
            columns=[
                "cluster",
                "resolution",
                "source_mmseqs_representative",
                "region",
                "reference_genome",
                "reference_gene",
                "ov_seq",
                "overlap_length",
                "upstream_match_length",
                "downstream_match_length",
            ]
        )

    final_overlap_df.to_csv(overlaps_path, sep="\t", index=False)

    ambiguous_path = outdir / "ambiguous_clusters.tsv"
    resolved_path = outdir / "resolved_clusters.tsv"
    unresolved_path = outdir / "unresolved_ambiguous_members.tsv"
    resolution_summary_path = outdir / "cluster_resolution_summary.tsv"
    incomplete_path = outdir / "incomplete_clusters.tsv"
    unclustered_path = outdir / "unclustered_cds.tsv"

    write_diagnostic_table(ambiguous_df, ambiguous_path, cds_df.columns)
    write_diagnostic_table(incomplete_df, incomplete_path, cds_df.columns)
    write_diagnostic_table(unclustered_df, unclustered_path, cds_df.columns)

    if resolved_df.empty:
        pd.DataFrame(columns=[
            "source_mmseqs_representative", "resolution_index",
            "mean_balanced_coverage", "min_balanced_coverage", "mean_pident"
        ] + list(cds_df.columns)).to_csv(resolved_path, sep="\t", index=False)
    else:
        resolved_df.to_csv(resolved_path, sep="\t", index=False)

    if unresolved_ambiguous_df.empty:
        pd.DataFrame(columns=[
            "source_mmseqs_representative", "classification",
            "resolution_stop_reason"
        ] + list(cds_df.columns)).to_csv(unresolved_path, sep="\t", index=False)
    else:
        unresolved_ambiguous_df.to_csv(unresolved_path, sep="\t", index=False)

    if resolution_summary_df.empty:
        pd.DataFrame(columns=[
            "source_mmseqs_representative", "original_member_count",
            "distinct_genomes", "resolved_subclusters", "resolved_members",
            "unresolved_members", "resolution_stop_reason"
        ]).to_csv(resolution_summary_path, sep="\t", index=False)
    else:
        resolution_summary_df.to_csv(resolution_summary_path, sep="\t", index=False)

    direct_count = sum(1 for _, group in named_clusters if group["resolution"] == "direct")
    resolved_count = sum(1 for _, group in named_clusters if group["resolution"] == "resolved")

    print(f"[✓] Saved cluster membership:        {membership_path}")
    print(f"[✓] Saved overlap summary:           {overlaps_path}")
    print(f"[✓] Saved raw ambiguous clusters:    {ambiguous_path}")
    print(f"[✓] Saved resolved clusters:         {resolved_path}")
    print(f"[✓] Saved unresolved ambiguous CDSs: {unresolved_path}")
    print(f"[✓] Saved resolution summary:        {resolution_summary_path}")
    print(f"[✓] Saved incomplete clusters:       {incomplete_path}")
    print(f"[✓] Saved unclustered CDSs:          {unclustered_path}")
    print(f"[✓] Saved all-vs-all search:         {allvsall_path}")
    print(f"[✓] Final direct clusters:           {direct_count}")
    print(f"[✓] Final resolved clusters:         {resolved_count}")
    print(f"[✓] Overlap candidates found:        {len(final_overlap_df)}")


if __name__ == "__main__":
    main()
