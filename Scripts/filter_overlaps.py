#!/usr/bin/env python3

import argparse
import os
import re

import pandas as pd
from Bio.Seq import Seq


def check_overlap_gc(ov_seq, min_size=20):
    """
    Check that the GC content of ov_seq is between 40-60%. If yes, return
    (sequence, GC%). If not, search for subsequences >= min_size within ov_seq
    that meet the GC range and return the longest one and its GC%.
    """
    if not ov_seq:
        return None, None

    seq = ov_seq.upper()

    def gc_content(s):
        return sum(base in {"G", "C"} for base in s) / len(s)

    gc = gc_content(seq)
    if 0.4 <= gc <= 0.6:
        return seq, round(gc * 100, 2)

    best_subseq = None
    best_gc = 0

    for i in range(len(seq)):
        for j in range(i + min_size, len(seq) + 1):
            subseq = seq[i:j]
            gc_sub = gc_content(subseq)
            if 0.4 <= gc_sub <= 0.6:
                if best_subseq is None or len(subseq) > len(best_subseq) or (
                    len(subseq) == len(best_subseq) and gc_sub > best_gc
                ):
                    best_subseq = subseq
                    best_gc = gc_sub

    if best_subseq:
        return best_subseq, round(best_gc * 100, 2)

    return None, None


def find_homopolymers(ov_seq, max_homo_size=4):
    """
    Return homopolymer runs longer than max_homo_size.

    For example, when max_homo_size=4, AAAA is allowed and AAAAA is flagged.
    Coordinates are 0-based and stop is exclusive.
    """
    ov_seq = ov_seq.upper()
    homos = []
    i = 0
    n = len(ov_seq)

    while i < n:
        j = i + 1
        while j < n and ov_seq[j] == ov_seq[i]:
            j += 1

        if j - i > max_homo_size:
            homos.append((i, j, ov_seq[i:j]))

        i = j

    return homos


def recommend_no_homo(ov_seq, homos, max_homo_size=4, min_len=1):
    """
    Produce continuous subsequences that contain no homopolymer longer than
    max_homo_size. Up to max_homo_size bases from a flagged run are retained
    at fragment boundaries.
    """
    ov_seq = ov_seq.upper()
    if not homos:
        return [ov_seq] if len(ov_seq) >= min_len else []

    fragments = []
    last_cut = 0

    for start, stop, _seq in sorted(homos, key=lambda h: h[0]):
        run_len = stop - start
        if run_len > max_homo_size:
            left_frag = ov_seq[last_cut:start + max_homo_size]
            if len(left_frag) >= min_len:
                fragments.append(left_frag)

            last_cut = max(0, stop - max_homo_size)

    if len(ov_seq) - last_cut >= min_len:
        fragments.append(ov_seq[last_cut:])

    return fragments


def find_repeats(seq, motif_len=2, min_repeats=3):
    """
    Find tandem repeats of motifs of length motif_len repeated >= min_repeats.
    Returns [(start, stop, motif, repeat_count, sequence), ...].
    """
    seq = seq.upper()
    pattern = re.compile(rf'((\w{{{motif_len}}})\2{{{min_repeats-1},}})')
    results = []

    for match in pattern.finditer(seq):
        start, stop = match.span()
        full = match.group(1)
        motif = match.group(2)
        repeat_count = len(full) // motif_len
        results.append((start, stop, motif, repeat_count, full))

    return results


def find_palindromes(seq, min_len=4, max_len=10, min_size=20):
    """
    Find reverse-complement palindromes, remove nested hits, and return the
    longest non-palindromic subsequence of at least min_size when available.
    """
    seq = seq.upper()
    palindromes = []

    for i in range(len(seq)):
        for length in range(min_len, max_len + 1):
            frag = seq[i:i + length]
            if len(frag) == length and frag == str(Seq(frag).reverse_complement()):
                palindromes.append((i, i + length, frag))

    if not palindromes:
        return [], seq if len(seq) >= min_size else None

    palindromes.sort(key=lambda x: (x[0], -(x[1] - x[0])))
    filtered = []
    current_end = -1

    for start, stop, frag in palindromes:
        if start >= current_end:
            filtered.append((start, stop, frag))
            current_end = stop

    non_pal_regions = []
    prev_end = 0

    for start, stop, _frag in filtered:
        if start > prev_end:
            non_pal_regions.append(seq[prev_end:start])
        prev_end = stop

    if prev_end < len(seq):
        non_pal_regions.append(seq[prev_end:])

    longest_non_pal = max(
        (frag for frag in non_pal_regions if len(frag) >= min_size),
        key=len,
        default=None,
    )

    return filtered, longest_non_pal


def optimize_overlap(ov_seq, min_size=20, max_homo_size=4):
    """
    Return (optimized_sequence, actions, rejection_reason).

    actions records any trimming that occurred. If optimization fails,
    optimized_sequence is None and rejection_reason describes why.
    """
    best_ov = ov_seq.upper()
    actions = []

    homos = find_homopolymers(best_ov, max_homo_size=max_homo_size)
    if homos:
        print(f"{len(homos)} homopolymer(s) > {max_homo_size} bp found: {homos}")
        no_homo = recommend_no_homo(
            best_ov,
            homos,
            max_homo_size=max_homo_size,
            min_len=min_size,
        )

        if not no_homo:
            return None, actions, f"no >= {min_size} bp region after homopolymer filtering"

        trimmed = max(no_homo, key=len)
        if trimmed != best_ov:
            actions.append("homopolymer_trim")
            best_ov = trimmed
            print(
                f"Recommended overlap without homopolymers > {max_homo_size}:\n"
                f"{best_ov}  Length: {len(best_ov)}"
            )

    di_repeats = find_repeats(best_ov, 2, 3)
    if di_repeats:
        print(f"Dinucleotide repeat region(s) found in {best_ov}\n{di_repeats}")
        return None, actions, "dinucleotide_repeat"

    iv_repeats, non_pal = find_palindromes(best_ov, 6, 10, min_size)
    if iv_repeats:
        print(f"Inverted repeat(s) found in {best_ov}\n{iv_repeats}")
        if non_pal is None:
            return None, actions, f"no >= {min_size} bp region after palindrome filtering"

        if non_pal != best_ov:
            actions.append("palindrome_trim")
            best_ov = non_pal
            print(
                f"Recommended overlap without inverted repeat(s): "
                f"{best_ov}\tLength: {len(best_ov)}"
            )

    gc_seq, gc_percent = check_overlap_gc(best_ov, min_size)
    if gc_seq is None:
        print(f"Could not find optimized overlap (40-60% GC) for {ov_seq}\n")
        return None, actions, "no_region_with_40_60_percent_gc"

    if gc_seq != best_ov:
        actions.append("gc_trim")
        best_ov = gc_seq

    print(f"Optimized overlap for {ov_seq}\n{best_ov}\tGC%: {gc_percent}")

    if len(best_ov) < min_size:
        return None, actions, f"below_min_size_after_filtering:{len(best_ov)}<{min_size}"

    return best_ov, actions, None


def adjust_coordinates_strand_aware(row, optimized_row, left_trim, right_trim):
    """
    Adjust per-genome genomic overlap coordinates after trimming.

    *_ov_start and *_ov_stop are stored as low/high genomic coordinates.
    If a genome's CDS is on the same strand as the reference CDS, left-side
    trimming removes bases from the low-coordinate side. If it is on the
    opposite strand, left-side trimming removes bases from the high-coordinate
    side instead.
    """
    reference_frame = row.get("reference_frame", None)

    if pd.isna(reference_frame) or reference_frame not in {"+", "-"}:
        raise ValueError(
            "Input overlap table must contain a valid 'reference_frame' column "
            "for strand-aware coordinate trimming."
        )

    start_cols = [col for col in row.index if col.endswith("_ov_start")]

    for start_col in start_cols:
        genome_id = start_col[:-len("_ov_start")]
        stop_col = f"{genome_id}_ov_stop"
        frame_col = f"{genome_id}_frame"

        if stop_col not in row.index or frame_col not in row.index:
            raise ValueError(
                f"Missing {stop_col} or {frame_col}; cannot perform "
                f"strand-aware trimming for {genome_id}."
            )

        if pd.isna(row[start_col]) or pd.isna(row[stop_col]):
            continue

        genome_frame = row[frame_col]
        if genome_frame not in {"+", "-"}:
            raise ValueError(
                f"Invalid frame '{genome_frame}' for genome {genome_id}."
            )

        if genome_frame == reference_frame:
            optimized_row[start_col] = int(row[start_col]) + left_trim
            optimized_row[stop_col] = int(row[stop_col]) - right_trim
        else:
            optimized_row[start_col] = int(row[start_col]) + right_trim
            optimized_row[stop_col] = int(row[stop_col]) - left_trim

        if optimized_row[start_col] > optimized_row[stop_col]:
            raise ValueError(
                f"Coordinate trimming inverted the overlap coordinates for "
                f"{genome_id}: {optimized_row[start_col]} > {optimized_row[stop_col]}"
            )


def make_report_row(row, status, reason, original_seq, final_seq=None,
                    left_trim=None, right_trim=None, actions=None):
    """Create one compact audit row for overlap_filter_report.tsv."""
    report = {
        "cluster": row.get("cluster", ""),
        "region": row.get("region", ""),
        "resolution": row.get("resolution", ""),
        "source_mmseqs_representative": row.get(
            "source_mmseqs_representative", ""
        ),
        "reference_genome": row.get("reference_genome", ""),
        "reference_gene": row.get("reference_gene", ""),
        "status": status,
        "reason": reason,
        "filter_actions": ";".join(actions or []),
        "original_ov_seq": original_seq,
        "original_overlap_length": len(original_seq),
        "final_ov_seq": final_seq if final_seq is not None else "",
        "final_overlap_length": len(final_seq) if final_seq is not None else "",
        "left_trim": left_trim if left_trim is not None else "",
        "right_trim": right_trim if right_trim is not None else "",
    }
    return report


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Select optimal overlap sequences from those identified by "
            "find_cluster_overlaps.py"
        )
    )
    parser.add_argument("--overlaps", required=True, help="cluster_overlap_summary.tsv")
    parser.add_argument(
        "--outdir",
        required=False,
        help="Directory to save filtered overlaps. Default: directory of --overlaps.",
    )
    parser.add_argument(
        "--min_size",
        type=int,
        default=20,
        help="Minimum overlap size in bp. Default: 20.",
    )
    parser.add_argument(
        "--max_homo_size",
        type=int,
        default=4,
        help=(
            "Maximum allowed homopolymer run length. Runs longer than this "
            "are trimmed. Default: 4."
        ),
    )

    args = parser.parse_args()

    outdir = (
        args.outdir
        if args.outdir
        else os.path.dirname(os.path.abspath(args.overlaps)) or "."
    )
    os.makedirs(outdir, exist_ok=True)

    outfile = os.path.join(outdir, "filtered_overlaps.tsv")
    report_file = os.path.join(outdir, "overlap_filter_report.tsv")

    ov_df = pd.read_csv(args.overlaps, sep="\t")

    required_cols = {"ov_seq", "reference_frame"}
    missing = required_cols - set(ov_df.columns)
    if missing:
        raise ValueError(
            "Input overlap table is missing required column(s): "
            + ", ".join(sorted(missing))
        )

    optimized_overlaps_rows = []
    report_rows = []

    for _, row in ov_df.iterrows():
        ov = str(row["ov_seq"]).upper()

        if len(ov) < args.min_size:
            report_rows.append(
                make_report_row(
                    row=row,
                    status="rejected",
                    reason=f"below_min_size:{len(ov)}<{args.min_size}",
                    original_seq=ov,
                )
            )
            continue

        print(f"Analyzing: {ov}")
        optimized_ov, actions, rejection_reason = optimize_overlap(
            ov,
            min_size=args.min_size,
            max_homo_size=args.max_homo_size,
        )

        if optimized_ov is None:
            print(f"Could not optimize {ov}\n")
            report_rows.append(
                make_report_row(
                    row=row,
                    status="rejected",
                    reason=rejection_reason,
                    original_seq=ov,
                    actions=actions,
                )
            )
            continue

        delta_start = ov.find(optimized_ov)
        if delta_start < 0:
            raise ValueError(
                "Optimized overlap is not a subsequence of the original overlap: "
                f"{optimized_ov} not found in {ov}"
            )

        left_trim = delta_start
        right_trim = len(ov) - delta_start - len(optimized_ov)

        print(
            f"Optimized overlap left trim:\t{left_trim} "
            f"right trim:\t{right_trim}\tSize:\t{len(optimized_ov)}\n"
        )

        opt_ov_row = row.copy()

        # Preserve original sequence/length before replacing the overlap.
        opt_ov_row["original_ov_seq"] = ov
        opt_ov_row["original_overlap_length"] = len(ov)
        opt_ov_row["ov_seq"] = optimized_ov
        opt_ov_row["overlap_length"] = len(optimized_ov)

        adjust_coordinates_strand_aware(
            row=row,
            optimized_row=opt_ov_row,
            left_trim=left_trim,
            right_trim=right_trim,
        )

        optimized_overlaps_rows.append(opt_ov_row)

        status = "trimmed" if optimized_ov != ov else "accepted"
        reason = ";".join(actions) if actions else "passed_without_trimming"

        report_rows.append(
            make_report_row(
                row=row,
                status=status,
                reason=reason,
                original_seq=ov,
                final_seq=optimized_ov,
                left_trim=left_trim,
                right_trim=right_trim,
                actions=actions,
            )
        )

    # Preserve the expected columns even if no overlaps pass filtering.
    output_columns = list(ov_df.columns)
    for extra_col in ["original_ov_seq", "original_overlap_length"]:
        if extra_col not in output_columns:
            output_columns.append(extra_col)

    optimized_ov_df = pd.DataFrame(
        optimized_overlaps_rows,
        columns=output_columns,
    )
    optimized_ov_df.to_csv(outfile, sep="\t", index=False)

    report_df = pd.DataFrame(report_rows)
    report_df.to_csv(report_file, sep="\t", index=False)

    accepted_count = sum(report_df["status"].isin(["accepted", "trimmed"])) if not report_df.empty else 0
    rejected_count = sum(report_df["status"] == "rejected") if not report_df.empty else 0

    print(f"\n[✓] Saved {accepted_count} filtered overlaps to: {outfile}")
    print(f"[✓] Saved filter report for {len(report_df)} input overlaps to: {report_file}")
    print(f"    Accepted/trimmed: {accepted_count}")
    print(f"    Rejected:         {rejected_count}")


if __name__ == "__main__":
    main()
