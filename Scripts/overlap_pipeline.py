conda activate combphage

python Scripts/fetch_genbanks.py -a MZ501081.1 MZ501078.1 V01146.1 -o genbanks --email email@email.com

python Scripts/prepare_genbanks.py -i genbanks/ -o prepped_genbanks

python Scripts/find_overlaps.py --input-dir prepped_genbanks/ --outdir feature_overlaps

python Scripts/filter_overlaps.py --overlaps feature_overlaps/cluster_overlap_summary.tsv --min_size 25

python Scripts/visualize.py --genbanks prepped_genbanks/genbanks/ --overlaps-dir feature_overlaps