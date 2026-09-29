#!/usr/bin/env python3
import argparse
import pandas as pd
from itertools import combinations
from Bio import SeqIO
import ast
import os



def assign_color(feature_id, overlap):
        if feature_id :
            return "#cccccc"  # grey
        if "start" in regions and "stop" in regions:
            return "#d367eb"  # purple
        elif "start" in regions:
            return "#15ed2b"  # green
        elif "stop" in regions:
            return "#ed1515"  # red
        return "#cccccc"  # fallback grey



def main():
    parser = argparse.ArgumentParser(
        description="Update lovis4u annotation table to color code CDSs where overlaps were found"
    )
    parser.add_argument("--clusters", required=True, help="cluster_membership.tsv")
    parser.add_argument("--overlaps", required=True, help="optimized_overlaps.tsv")
    parser.add_argument("--annotations", required=True, help="feature_annotations_table.tsv from lovis4u output")
    parser.add_argument("--out_dir", required=False, help="output directory (default is same directory as annotations input)")
    args = parser.parse_args()

    if args.out_dir:
        out_file = os.path.splitext(os.path.basename(args.annotations))[0] + "_overlaps.tsv"
        out_path = os.path.join(args.out_dir, out_file)
    else:
        out_path = os.path.splitext(args.annotations)[0] + "_overlaps.tsv"

    # load clusters
    clusters_df = pd.read_csv(args.clusters, sep="\t")

    # load overlaps
    overlaps_df = pd.read_csv(args.overlaps, sep="\t")

    # load lovis annotations table
    annotations_df = pd.read_csv(args.annotations, sep="\t")

    # update annotations and colors
    for i, row in annotations_df.iterrows():
        try:
            gene_id = row["feature_id"]
            cluster_value = clusters_df.query("gene == @gene_id")["cluster"].iloc[0]
            cluster_overlap_regions = overlaps_df.query("cluster == @cluster_value")["region"].tolist()

            print(f"Cluster value: {cluster_value}\tOverlap regions: {cluster_overlap_regions}\n")
            if "start" in cluster_overlap_regions and "stop" in cluster_overlap_regions:
                color = "#d367eb"
            elif "start" in cluster_overlap_regions:
                color = "#15ed2b"
            elif "stop" in cluster_overlap_regions:
                color = "#ed1515"
            else:
                color = "#cccccc"
                #cluster_value = " "
        except Exception:
            color = "#cccccc"
            #cluster_value = " "

        # Update df
        #annotations_df.loc[i, "name"] = cluster_value
        annotations_df.loc[i, "fill_colour"] = color


    # Save updated annotations
    annotations_df.to_csv(out_path, sep="\t", index=False)


if __name__ == "__main__":
    main()
