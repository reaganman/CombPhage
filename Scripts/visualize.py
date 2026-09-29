#!/usr/bin/env python3

import argparse
import shutil
import subprocess
from pathlib import Path

import pandas as pd


DEFAULT_GREY = "#cccccc"
START_COLOUR = "#15ed2b"
STOP_COLOUR = "#ed1515"
BOTH_COLOUR = "#d367eb"


def run_command(cmd):
    print("\n[+] Running:")
    print("    " + " ".join(str(x) for x in cmd))
    subprocess.run(cmd, check=True)


def validate_inputs(genbanks_dir, clusters_file, overlaps_file):
    if not genbanks_dir.exists() or not genbanks_dir.is_dir():
        raise FileNotFoundError(
            f"GenBank directory does not exist or is not a directory: {genbanks_dir}"
        )

    gb_extensions = {".gb", ".gbk", ".gbff", ".genbank"}
    gb_files = [
        p for p in genbanks_dir.iterdir()
        if p.is_file() and p.suffix.lower() in gb_extensions
    ]

    if not gb_files:
        raise FileNotFoundError(f"No GenBank files found in: {genbanks_dir}")

    if not clusters_file.exists():
        raise FileNotFoundError(f"Cluster membership file not found: {clusters_file}")

    if not overlaps_file.exists():
        raise FileNotFoundError(f"Filtered overlap file not found: {overlaps_file}")

    return gb_files


def colour_for_regions(regions):
    has_start = "start" in regions
    has_stop = "stop" in regions

    if has_start and has_stop:
        return BOTH_COLOUR
    if has_start:
        return START_COLOUR
    if has_stop:
        return STOP_COLOUR
    return DEFAULT_GREY


def load_gene_info(clusters_file, overlaps_file):
    """
    Build a mapping:

        gene -> {
            "cluster": cluster_id,
            "regions": {"start", "stop"}
        }

    Only clusters that have at least one filtered overlap receive regions.
    """
    clusters_df = pd.read_csv(clusters_file, sep="\t", dtype=str)
    overlaps_df = pd.read_csv(overlaps_file, sep="\t", dtype=str)

    missing = {"cluster", "gene"} - set(clusters_df.columns)
    if missing:
        raise ValueError(
            f"{clusters_file} is missing required column(s): {', '.join(sorted(missing))}"
        )

    missing = {"cluster", "region"} - set(overlaps_df.columns)
    if missing:
        raise ValueError(
            f"{overlaps_file} is missing required column(s): {', '.join(sorted(missing))}"
        )

    overlaps_df = overlaps_df[
        overlaps_df["region"].isin({"start", "stop"})
    ].copy()

    cluster_to_regions = (
        overlaps_df.groupby("cluster")["region"]
        .apply(lambda values: set(values))
        .to_dict()
    )

    gene_info = {}

    for _, row in clusters_df.iterrows():
        gene = str(row["gene"])
        cluster = str(row["cluster"])

        # Each CDS should belong to at most one final CombPhage cluster.
        # If that invariant is ever violated, fail loudly rather than
        # assigning an arbitrary label.
        if gene in gene_info and gene_info[gene]["cluster"] != cluster:
            raise ValueError(
                f"Gene {gene} appears in more than one final cluster: "
                f"{gene_info[gene]['cluster']} and {cluster}."
            )

        gene_info[gene] = {
            "cluster": cluster,
            "regions": cluster_to_regions.get(cluster, set()),
        }

    return gene_info, clusters_df


def create_feature_override_table(base_annotation_file, clusters_file, overlaps_file, output_file):
    annotations_df = pd.read_csv(base_annotation_file, sep="\t", dtype=str)

    if "feature_id" not in annotations_df.columns:
        raise ValueError(
            f"{base_annotation_file} does not contain a feature_id column."
        )

    gene_info, clusters_df = load_gene_info(clusters_file, overlaps_file)

    rows = []
    counts = {"start": 0, "stop": 0, "both": 0, "grey": 0, "labelled": 0}
    feature_ids = set(annotations_df["feature_id"].dropna().astype(str))

    for feature_id in annotations_df["feature_id"].astype(str):
        info = gene_info.get(feature_id)

        if info is not None and info["regions"]:
            cluster = info["cluster"]
            regions = info["regions"]
            colour = colour_for_regions(regions)
            name = cluster
            show_label = 1
            counts["labelled"] += 1

            if regions == {"start"}:
                counts["start"] += 1
            elif regions == {"stop"}:
                counts["stop"] += 1
            elif "start" in regions and "stop" in regions:
                counts["both"] += 1

        else:
            colour = DEFAULT_GREY
            name = ""
            show_label = 0
            counts["grey"] += 1

        rows.append(
            {
                "feature_id": feature_id,
                "fill_colour": colour,
                "name": name,
                "show_label": show_label,
            }
        )

    pd.DataFrame(rows).to_csv(output_file, sep="\t", index=False)

    cluster_genes = set(clusters_df["gene"].dropna().astype(str))
    missing_feature_ids = sorted(cluster_genes - feature_ids)

    print("\n[+] Feature-colour override table created")
    print(f"    Total features:             {len(rows)}")
    print(f"    Start overlap only:         {counts['start']}")
    print(f"    Stop overlap only:          {counts['stop']}")
    print(f"    Start + stop:               {counts['both']}")
    print(f"    No filtered overlap:        {counts['grey']}")
    print(f"    Labelled overlap CDSs:      {counts['labelled']}")

    if missing_feature_ids:
        print(
            f"\n[!] {len(missing_feature_ids)} cluster-member CDS ID(s) were not found "
            "in LoVis4u's feature table."
        )
        print(
            "    Check that cluster_membership.tsv and the normalized GenBanks "
            "use the same CDS identifiers."
        )
        for feature_id in missing_feature_ids[:10]:
            print(f"      {feature_id}")
        if len(missing_feature_ids) > 10:
            print(f"      ... and {len(missing_feature_ids) - 10} more")


def remove_output_dir(path, force):
    if path.exists():
        if not force:
            raise FileExistsError(
                f"Output directory already exists: {path}\n"
                "Use --force to remove and regenerate it."
            )
        shutil.rmtree(path)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Run LoVis4u, colour CombPhage CDSs based on filtered start/stop "
            "overlaps, and rerun LoVis4u with the updated feature colours."
        )
    )

    parser.add_argument("--genbanks", required=True, help="Directory containing normalized GenBank files.")
    parser.add_argument(
        "--overlaps-dir",
        required=True,
        help=(
            "CombPhage results directory containing cluster_membership.tsv "
            "and filtered_overlaps.tsv."
        ),
    )
    parser.add_argument(
        "--clusters",
        help="Optional path to cluster_membership.tsv. Default: <overlaps-dir>/cluster_membership.tsv",
    )
    parser.add_argument(
        "--overlaps",
        help="Optional path to filtered_overlaps.tsv. Default: <overlaps-dir>/filtered_overlaps.tsv",
    )
    parser.add_argument("--config", default="A4L", help="LoVis4u configuration name/path. Default: A4L")
    parser.add_argument("--lovis4u-bin", default="lovis4u", help="LoVis4u executable name/path. Default: lovis4u")
    parser.add_argument(
        "--reuse-base",
        action="store_true",
        help="Reuse an existing initial LoVis4u run in <overlaps-dir>/lovis4u.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Delete existing LoVis4u output directories before rerunning.",
    )

    args = parser.parse_args()

    if args.force and args.reuse_base:
        parser.error("--force and --reuse-base cannot be used together.")

    genbanks_dir = Path(args.genbanks).resolve()
    overlaps_dir = Path(args.overlaps_dir).resolve()

    clusters_file = (
        Path(args.clusters).resolve()
        if args.clusters
        else overlaps_dir / "cluster_membership.tsv"
    )
    overlaps_file = (
        Path(args.overlaps).resolve()
        if args.overlaps
        else overlaps_dir / "filtered_overlaps.tsv"
    )

    base_dir = overlaps_dir / "lovis4u"
    final_dir = overlaps_dir / "lovis4u_overlaps"
    base_feature_table = base_dir / "feature_annotation_table.tsv"
    locus_table = base_dir / "locus_annotation_table.tsv"
    feature_override = base_dir / "feature_overlap_colours.tsv"

    overlaps_dir.mkdir(parents=True, exist_ok=True)

    gb_files = validate_inputs(genbanks_dir, clusters_file, overlaps_file)

    print(f"[+] GenBank directory:   {genbanks_dir}")
    print(f"[+] GenBank files:       {len(gb_files)}")
    print(f"[+] Cluster membership:  {clusters_file}")
    print(f"[+] Filtered overlaps:    {overlaps_file}")

    lovis_bin = shutil.which(args.lovis4u_bin)
    if lovis_bin is None:
        explicit = Path(args.lovis4u_bin)
        if explicit.exists():
            lovis_bin = str(explicit.resolve())
        else:
            raise FileNotFoundError(
                f"Could not find LoVis4u executable: {args.lovis4u_bin}\n"
                "Activate the CombPhage conda environment and verify that `lovis4u -h` works."
            )

    if args.reuse_base:
        print("\n[+] Reusing existing initial LoVis4u run.")
        if not base_feature_table.exists():
            raise FileNotFoundError(f"Missing: {base_feature_table}")
        if not locus_table.exists():
            raise FileNotFoundError(f"Missing: {locus_table}")
    else:
        remove_output_dir(base_dir, args.force)
        run_command([
            lovis_bin,
            "-gb", str(genbanks_dir),
            "-o", str(base_dir),
            "-hl",
            "-sxa",
            "-align",
            "-cl-off",
            "-c", args.config,
        ])

        if not base_feature_table.exists():
            raise FileNotFoundError(f"LoVis4u did not create: {base_feature_table}")
        if not locus_table.exists():
            raise FileNotFoundError(f"LoVis4u did not create: {locus_table}")

    create_feature_override_table(
        base_feature_table,
        clusters_file,
        overlaps_file,
        feature_override,
    )

    print(f"\n[+] Saved feature override table: {feature_override}")

    remove_output_dir(final_dir, args.force)

    run_command([
        lovis_bin,
        "-gb", str(genbanks_dir),
        "-o", str(final_dir),
        "-hl",
        "-sxa",
        "-cl-off",
        "-c", args.config,
        "--locus-annotation-file", str(locus_table),
        "--feature-annotation-file", str(feature_override),
    ])

    print("\n[✓] LoVis4u overlap visualization complete")
    print(f"    Base output:   {base_dir}")
    print(f"    Colour table:  {feature_override}")
    print(f"    Final output:  {final_dir}")


if __name__ == "__main__":
    main()
