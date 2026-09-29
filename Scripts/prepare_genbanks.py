#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


GENBANK_EXTENSIONS = {".gb", ".gbk", ".gbff", ".genbank"}


def find_genbank_files(inputs):
    """Collect GenBank files from files and/or directories supplied by the user."""
    genbank_files = []

    for input_path in inputs:
        path = Path(input_path)

        if not path.exists():
            raise FileNotFoundError(f"Input does not exist: {path}")

        if path.is_file():
            if path.suffix.lower() not in GENBANK_EXTENSIONS:
                raise ValueError(
                    f"Unsupported GenBank extension: {path}\n"
                    f"Expected one of: {', '.join(sorted(GENBANK_EXTENSIONS))}"
                )
            genbank_files.append(path)

        elif path.is_dir():
            files = sorted(
                p for p in path.iterdir()
                if p.is_file() and p.suffix.lower() in GENBANK_EXTENSIONS
            )

            if not files:
                raise ValueError(f"No GenBank files found in directory: {path}")

            genbank_files.extend(files)

    seen = set()
    unique_files = []

    for path in genbank_files:
        resolved = path.resolve()
        if resolved not in seen:
            seen.add(resolved)
            unique_files.append(path)

    return unique_files


def clean_protein_sequence(sequence):
    """Remove whitespace and a terminal stop-codon marker."""
    sequence = str(sequence).replace(" ", "").replace("\n", "")

    if sequence.endswith("*"):
        sequence = sequence[:-1]

    return sequence


def translate_feature(feature, record, default_table=11):
    """
    Return the protein translation for a CDS.

    Prefer the GenBank /translation qualifier. If absent,
    derive the translation from the nucleotide sequence.
    """
    if "translation" in feature.qualifiers:
        protein = feature.qualifiers["translation"][0]
        return clean_protein_sequence(protein), "genbank"

    cds_sequence = feature.extract(record.seq)

    codon_start = int(feature.qualifiers.get("codon_start", ["1"])[0])
    cds_sequence = cds_sequence[codon_start - 1:]

    transl_table = int(
        feature.qualifiers.get("transl_table", [str(default_table)])[0]
    )

    remainder = len(cds_sequence) % 3
    if remainder:
        cds_sequence = cds_sequence[:-remainder]

    if len(cds_sequence) == 0:
        raise ValueError("CDS sequence is empty after codon adjustment.")

    protein = cds_sequence.translate(
        table=transl_table,
        to_stop=False
    )

    return clean_protein_sequence(protein), "derived"


def remove_duplicate_cds_features(cds_features, genome_id):
    """
    Remove exact duplicate CDS annotations.

    Two CDS features are duplicates only when they have identical:
      - genomic start
      - genomic end
      - strand

    Overlapping CDSs with different boundaries are retained.
    """
    seen = set()
    unique_features = []
    duplicate_features = []

    for feature in cds_features:
        key = (
            int(feature.location.start),
            int(feature.location.end),
            feature.location.strand,
        )

        if key in seen:
            duplicate_features.append(feature)
            continue

        seen.add(key)
        unique_features.append(feature)

    if duplicate_features:
        print(
            f"[!] {genome_id}: skipped "
            f"{len(duplicate_features)} duplicate CDS feature(s)"
        )

        for feature in duplicate_features:
            left = int(feature.location.start) + 1
            right = int(feature.location.end)
            strand = "+" if feature.location.strand == 1 else "-"

            locus_tag = feature.qualifiers.get("locus_tag", [""])[0]
            label = f" ({locus_tag})" if locus_tag else ""

            print(
                f"    duplicate CDS: "
                f"{left}-{right} ({strand}){label}"
            )

    return unique_features, len(duplicate_features)


def prepare_record(record, default_table=11):
    """
    Normalize a GenBank record for CombPhage.

    Returns:
        normalized_record
        protein_records
        cds_rows
        duplicate_count
    """
    genome_id = record.id

    if not genome_id:
        raise ValueError("GenBank record has no record ID.")

    if len(record.seq) == 0:
        raise ValueError(
            f"Genome {genome_id} contains no nucleotide sequence."
        )

    cds_features = [
        feature
        for feature in record.features
        if feature.type == "CDS"
    ]

    if not cds_features:
        raise ValueError(
            f"Genome {genome_id} contains no CDS features."
        )

    original_cds_count = len(cds_features)

    cds_features, duplicate_count = remove_duplicate_cds_features(
        cds_features,
        genome_id
    )

    if not cds_features:
        raise ValueError(
            f"Genome {genome_id} has no unique CDS features after "
            "duplicate removal."
        )

    cds_features.sort(
        key=lambda feature: (
            int(feature.location.start),
            int(feature.location.end),
        )
    )

    protein_records = []
    cds_rows = []

    for cds_index, feature in enumerate(cds_features, start=1):

        strand_value = feature.location.strand

        if strand_value not in (1, -1):
            raise ValueError(
                f"{genome_id}: CDS at {feature.location} "
                "has no valid strand information."
            )

        frame = "+" if strand_value == 1 else "-"

        left = int(feature.location.start) + 1
        right = int(feature.location.end)

        if frame == "+":
            cds_start = left
            cds_stop = right
        else:
            cds_start = right
            cds_stop = left

        gene_id = f"{genome_id}_CDS_{cds_index:04d}"

        original_locus_tag = feature.qualifiers.get("locus_tag", [""])[0]
        original_protein_id = feature.qualifiers.get("protein_id", [""])[0]
        original_gene = feature.qualifiers.get("gene", [""])[0]
        product = feature.qualifiers.get(
            "product", ["hypothetical protein"]
        )[0]

        transl_table = int(
            feature.qualifiers.get(
                "transl_table",
                [str(default_table)]
            )[0]
        )

        protein, translation_source = translate_feature(
            feature,
            record,
            default_table=default_table
        )

        if original_locus_tag:
            old_tags = feature.qualifiers.get("old_locus_tag", [])
            if original_locus_tag not in old_tags:
                feature.qualifiers["old_locus_tag"] = (
                    old_tags + [original_locus_tag]
                )

        feature.qualifiers["locus_tag"] = [gene_id]
        feature.qualifiers["protein_id"] = [gene_id]
        feature.qualifiers["translation"] = [protein]

        protein_record = SeqRecord(
            Seq(protein),
            id=gene_id,
            description=product
        )
        protein_records.append(protein_record)

        cds_rows.append(
            {
                "gene": gene_id,
                "genome": genome_id,
                "start": cds_start,
                "stop": cds_stop,
                "frame": frame,
                "cds_index": cds_index,
                "product": product,
                "translation_source": translation_source,
                "transl_table": transl_table,
                "original_locus_tag": original_locus_tag,
                "original_protein_id": original_protein_id,
                "original_gene": original_gene,
                "location": str(feature.location),
            }
        )

    print(
        f"    {genome_id}: {len(record.seq):,} bp, "
        f"{original_cds_count} CDS input, "
        f"{len(cds_features)} unique CDS"
    )

    return record, protein_records, cds_rows, duplicate_count


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Prepare annotated GenBank files for CombPhage. "
            "Produces standardized genome FASTA, protein FASTA, "
            "CDS table, and normalized GenBank files."
        )
    )

    parser.add_argument(
        "-i",
        "--genbanks",
        nargs="+",
        required=True,
        help=(
            "GenBank files and/or directories containing GenBank files "
            "(.gb, .gbk, .gbff, .genbank)"
        ),
    )

    parser.add_argument(
        "-o",
        "--outdir",
        required=True,
        help="Output directory",
    )

    parser.add_argument(
        "--default-transl-table",
        type=int,
        default=11,
        help=(
            "Translation table to use when a CDS does not specify "
            "/transl_table. Default: 11"
        ),
    )

    args = parser.parse_args()

    output_dir = Path(args.outdir)
    normalized_genbank_dir = output_dir / "genbanks"

    output_dir.mkdir(parents=True, exist_ok=True)
    normalized_genbank_dir.mkdir(parents=True, exist_ok=True)

    genbank_files = find_genbank_files(args.genbanks)

    print(f"[+] Found {len(genbank_files)} GenBank file(s).")

    genome_records = []
    protein_records = []
    cds_rows = []

    genome_ids = set()
    total_duplicates = 0

    for genbank_file in genbank_files:
        print(f"[+] Reading {genbank_file}")

        records = list(SeqIO.parse(genbank_file, "genbank"))

        if not records:
            raise ValueError(
                f"No GenBank records found in {genbank_file}"
            )

        for record in records:

            if record.id in genome_ids:
                raise ValueError(
                    f"Duplicate genome ID encountered: {record.id}\n"
                    "Every genome supplied to CombPhage must have "
                    "a unique GenBank record ID."
                )

            genome_ids.add(record.id)

            normalized_record, proteins, rows, duplicate_count = (
                prepare_record(
                    record,
                    default_table=args.default_transl_table
                )
            )

            total_duplicates += duplicate_count

            genome_records.append(
                SeqRecord(
                    normalized_record.seq,
                    id=normalized_record.id,
                    description=""
                )
            )

            protein_records.extend(proteins)
            cds_rows.extend(rows)

            normalized_path = (
                normalized_genbank_dir
                / f"{normalized_record.id}.gbk"
            )

            SeqIO.write(
                normalized_record,
                normalized_path,
                "genbank"
            )

    genomes_path = output_dir / "genomes.fasta"
    SeqIO.write(genome_records, genomes_path, "fasta")

    proteins_path = output_dir / "proteins.faa"
    SeqIO.write(protein_records, proteins_path, "fasta")

    cds_path = output_dir / "cds.tsv"
    cds_df = pd.DataFrame(cds_rows)
    cds_df.to_csv(cds_path, sep="\t", index=False)

    print()
    print("[✓] CombPhage inputs prepared successfully")
    print(f"    Genomes:                {len(genome_records)}")
    print(f"    Unique CDSs:            {len(cds_rows)}")
    print(f"    Duplicate CDSs skipped: {total_duplicates}")
    print()
    print(f"    {genomes_path}")
    print(f"    {proteins_path}")
    print(f"    {cds_path}")
    print(f"    {normalized_genbank_dir}/")


if __name__ == "__main__":
    main()