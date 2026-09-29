#!/usr/bin/env python3

import argparse
from pathlib import Path

from Bio import Entrez, SeqIO


def read_accessions(accessions, accession_file=None):
    """
    Combine accessions provided directly on the command line
    with accessions provided in a text file.
    """

    combined = []

    if accessions:
        combined.extend(accessions)

    if accession_file:
        with open(accession_file) as f:
            for line in f:
                line = line.strip()

                # Ignore blank lines and comments
                if not line or line.startswith("#"):
                    continue

                # Allow either one accession per line or whitespace-separated
                combined.extend(line.split())

    # Remove duplicates while preserving order
    seen = set()
    unique_accessions = []

    for accession in combined:
        if accession not in seen:
            seen.add(accession)
            unique_accessions.append(accession)

    return unique_accessions


def fetch_genbank(accession, output_dir):
    """
    Download one GenBank record from NCBI and save it to disk.
    """

    print(f"[+] Fetching {accession}...")

    with Entrez.efetch(
        db="nuccore",
        id=accession,
        rettype="gbwithparts",
        retmode="text"
    ) as handle:

        records = list(SeqIO.parse(handle, "genbank"))

    if not records:
        raise ValueError(
            f"No GenBank record returned for accession {accession}"
        )

    if len(records) > 1:
        raise ValueError(
            f"Accession {accession} returned {len(records)} records; "
            "expected exactly one."
        )

    record = records[0]

    # Use the accession/version actually returned by NCBI
    filename = f"{record.id}.gbk"
    output_path = output_dir / filename

    SeqIO.write(
        record,
        output_path,
        "genbank"
    )

    print(
        f"    Saved {record.id} "
        f"({len(record.seq):,} bp, "
        f"{sum(f.type == 'CDS' for f in record.features)} CDS)"
    )

    return output_path


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Download GenBank records from NCBI for a list of "
            "nucleotide accessions."
        )
    )

    parser.add_argument(
        "-a",
        "--accessions",
        nargs="+",
        help=(
            "NCBI nucleotide accessions, e.g. "
            "MZ501081.1 MZ501078.1 PQ850602.1"
        )
    )

    parser.add_argument(
        "-f",
        "--accession-file",
        help=(
            "Text file containing accessions. Accessions may be "
            "one per line or whitespace-separated."
        )
    )

    parser.add_argument(
        "-o",
        "--outdir",
        required=True,
        help="Directory in which downloaded GenBank files will be saved."
    )

    parser.add_argument(
        "--email",
        required=True,
        help=(
            "Email address supplied to NCBI Entrez."
        )
    )

    parser.add_argument(
        "--api-key",
        help=(
            "Optional NCBI API key."
        )
    )

    args = parser.parse_args()

    if not args.accessions and not args.accession_file:
        parser.error(
            "Provide accessions with --accessions and/or --accession-file."
        )

    #
    # Configure NCBI Entrez
    #
    Entrez.email = args.email
    Entrez.tool = "CombPhage"

    if args.api_key:
        Entrez.api_key = args.api_key

    #
    # Create output directory
    #
    output_dir = Path(args.outdir)
    output_dir.mkdir(parents=True, exist_ok=True)

    #
    # Collect accession list
    #
    accessions = read_accessions(
        args.accessions,
        args.accession_file
    )

    if not accessions:
        parser.error("No valid accessions were provided.")

    print(f"[+] Downloading {len(accessions)} GenBank record(s)")
    print(f"[+] Output directory: {output_dir}")
    print()

    successful = []
    failed = []

    #
    # Download records
    #
    for accession in accessions:
        try:
            output_path = fetch_genbank(
                accession,
                output_dir
            )

            successful.append(
                (accession, output_path)
            )

        except Exception as exc:
            print(
                f"[!] Failed to fetch {accession}: {exc}"
            )
            failed.append(
                (accession, str(exc))
            )

    #
    # Summary
    #
    print()
    print("Download summary")
    print("----------------")
    print(f"Successful: {len(successful)}")
    print(f"Failed:     {len(failed)}")

    if successful:
        print()
        print("Downloaded files:")

        for accession, path in successful:
            print(f"  {accession} -> {path}")

    if failed:
        print()
        print("Failed accessions:")

        for accession, error in failed:
            print(f"  {accession}: {error}")


if __name__ == "__main__":
    main()