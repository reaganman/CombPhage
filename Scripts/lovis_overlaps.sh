#!/bin/bash
RESULTS_DIR=$1

ml anaconda
conda activate lovis4u_stable

# Initial lovis4u run to generate base feature annotation table
lovis4u -gff $RESULTS_DIR/pharokka_results/single_gffs/ \
    -o $RESULTS_DIR/lovis4u/ \
    -scc -hl -sxa -c A4L 

# update annotations table to color code overlaps
python Scripts/update_annotations_table.py \
    --clusters $RESULTS_DIR/cluster_membership.tsv \
    --overlaps $RESULTS_DIR/optimized_overlaps.tsv \
    --annotations $RESULTS_DIR/lovis4u/feature_annotation_table.tsv

# Re-run lovis4u
OUT_DIR="$RESULTS_DIR/lovis4u_overlaps"
mkdir -p $OUT_DIR
ANNO_FILE="$RESULTS_DIR/lovis4u/feature_annotation_table_overlaps.tsv"

echo $ANNO_FILE

lovis4u -gff $RESULTS_DIR/pharokka_results/single_gffs/ \
    -o "$OUT_DIR" \
    -scc -sxa -hl -c A4L \
    --feature-annotation-file "$ANNO_FILE"

