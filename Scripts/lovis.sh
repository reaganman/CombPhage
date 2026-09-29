#!/bin/bash
GBNKS_DIR=$1
OVERLAPS_DIR=$2

ml anaconda 
conda activate combphage

lovis4u_linux

# Initial lovis4u run to generate base feature annotation table
lovis4u -gb $GBNKS_DIR/ \
    -o $OVERLAPS_DIR/lovis4u/ \
    -hl -sxa -align -cl-off -c A4L 

# update annotations table to color code overlaps
python Scripts/update_annotations.py \
    --clusters $OVERLAPS_DIR/cluster_membership.tsv \
    --overlaps $OVERLAPS_DIR/filtered_overlaps.tsv \
    --annotations $OVERLAPS_DIR/lovis4u/feature_annotation_table.tsv

# Re-run lovis4u
OUT_DIR="$OVERLAPS_DIR/lovis4u_overlaps"
mkdir -p $OUT_DIR
ANNO_FILE="$OVERLAPS_DIR/lovis4u/feature_annotation_table_overlaps.tsv"
LOC_FILE="$OVERLAPS_DIR/lovis4u/locus_annotation_table.tsv"

echo $ANNO_FILE

lovis4u -gb $GBNKS_DIR -o "$OUT_DIR" \
    -sxa -hl -align -cl-off -c A4L \
    --feature-annotation-file "$ANNO_FILE"


