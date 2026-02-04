import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
import pandas as pd
import ast
import os
from collections import defaultdict
from BCBio import GFF
import itertools




def which_genes():
    """
    Determine which genes (CDSs) are found in each selected fragment
    and calculate pairwise similarity? (metric?) with homologs on
     sibling frags if present
    """
    pass

def split_gff(gff_path, out_path, frag_start, frag_stop, frag_num):
    """
    Create a new GFF file containing only features completely enclosed
    between frag_start and frag_stop, and write only custom header lines.
    Keeps original source, score, and phase.
    """
    filtered_rows = []
    seq_id = None

    # --- Read the original GFF lines manually ---
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#"):  # skip header/meta lines
                continue
            parts = line.strip().split("\t")
            if len(parts) != 9:
                continue  # skip malformed lines

            seqid, source, ftype, start, end, score, strand, phase, attributes = parts
            start, end = int(start), int(end)

            # keep only features fully inside fragment range
            if start >= frag_start and end <= frag_stop:
                filtered_rows.append((seqid, source, ftype, start, end, score, strand, phase, attributes))
                seq_id = seqid

    # --- Write new filtered GFF ---
    with open(out_path, "w") as out_handle:
        # Custom header lines only
        out_handle.write("##gff-version 3\n")
        out_handle.write(
            f"##sequence-region {seq_id} {frag_start} {frag_stop} fragment={frag_num}\n"
        )

        # Write filtered feature lines exactly as in original
        for row in filtered_rows:
            out_handle.write("\t".join(map(str, row)) + "\n")

    print(f"✅ Filtered GFF written to {out_path}")



def make_frags(chosen_seqs, n_frags, topology, overlaps):
    """
    Split each genome into n_frags fragments based on overlap sequences,
    including the overlap region on both sides for each fragment.
    Handles circular topology by including overlaps at both ends of fragment 1.
    """

    fragments = [[] for _ in range(n_frags)]
    boundaries = {}

    for seq in chosen_seqs:
        full_seq = str(seq.seq)
        seq_id = seq.id

        overlap_positions = []
        overlap_lengths = []

        # find positions and lengths of each overlap
        for i in range(1, n_frags):
            ov_seq = overlaps.get(f"{i}-{i+1}")
            if ov_seq:
                pos = full_seq.find(ov_seq)
                if pos == -1:
                    print(f"⚠️ Overlap {ov_seq[:15]}... not found in {seq_id}")
                    continue
                overlap_positions.append(pos)
                overlap_lengths.append(len(ov_seq))

        # handle circular topology: add last overlap at start for fragment 1
        if topology == "circular":
            if f"{n_frags}-1" in overlaps:
                ov_seq = overlaps[f"{n_frags}-1"]
                pos = full_seq.find(ov_seq)
                if pos == -1:
                    print(f"⚠️ Circular overlap {ov_seq[:15]}... not found in {seq_id}")
                else:
                    # Insert at the beginning for fragment 1
                    overlap_positions = [pos] + overlap_positions
                    overlap_lengths = [len(ov_seq)] + overlap_lengths

        # construct fragment boundaries
        cut_points = []
        prev_end = 0
        for idx, pos in enumerate(overlap_positions):
            end = pos + overlap_lengths[idx]  # include overlap in current frag
            cut_points.append((prev_end, end))
            prev_end = pos  # next frag starts at overlap start
        cut_points.append((prev_end, len(full_seq)))  # final fragment

        # adjust fragment 1 for circular topology
        if topology == "circular":
            start1, end1 = cut_points[0]
            # fragment 1 should start at last overlap start
            new_start1 = overlap_positions[0]  # last fragment overlap
            cut_points[0] = (new_start1, end1)

        # make SeqRecords
        genome_boundaries = []
        for i, (start, end) in enumerate(cut_points):
            frag_seq = full_seq[start:end]
            frag_id = f"{seq_id}_frag{i+1}"
            frag_record = SeqRecord(
                Seq(frag_seq),
                id=frag_id,
                description=f"Fragment {i+1} from {seq_id} ({start+1}:{end})"
            )
            fragments[i].append(frag_record)
            genome_boundaries.append((start+1, end))

        boundaries[seq_id] = genome_boundaries

    return fragments, boundaries




def write_concatenated_assemblies_no_overlap(fragments, overlap_lengths, outdir):
    """
    Write one FASTA per fragment combination, concatenating fragments
    while removing duplicated overlap regions.
    """
    os.makedirs(outdir, exist_ok=True)

    for idx, combo in enumerate(itertools.product(*fragments), start=1):
        assembled_seq = []
        frag_ids = []

        for i, frag in enumerate(combo):
            seq = str(frag.seq)

            # Trim left overlap for all but first fragment
            if i > 0:
                trim = overlap_lengths[i - 1]
                seq = seq[trim:]

            assembled_seq.append(seq)
            frag_ids.append(frag.id)

        full_seq = "".join(assembled_seq)

        record = SeqRecord(
            Seq(full_seq),
            id=f"assembly_{idx}",
            description=" | ".join(frag_ids)
        )

        out_path = os.path.join(outdir, f"assembly_{idx}.fasta")
        SeqIO.write(record, out_path, "fasta")

        print(f"Wrote {out_path}")




def write_assembly_manifest_tsv(
    fragments,
    boundaries,
    overlap_lengths,
    assembly_manifest_path,
    module_table_path,
):
    """
    Writes two TSVs:

    1) Assembly manifest (assembly coordinates respect overlap trimming)
    2) Module table (untrimmed source + assembly coordinates, long format)
    """

    n_frags = len(fragments)

    # -----------------------
    # Assembly manifest header
    # -----------------------
    manifest_header = ["assembly"]
    for i in range(1, n_frags + 1):
        manifest_header.extend([
            f"frag_{i}_source",
            f"frag_{i}_source_start",
            f"frag_{i}_source_stop",
            f"frag_{i}_assembly_start",
            f"frag_{i}_assembly_stop",
        ])

    with open(assembly_manifest_path, "w") as mf, open(module_table_path, "w") as mt:
        # Write headers
        mf.write("\t".join(manifest_header) + "\n")
        mt.write(
            "\t".join([
                "assembly",
                "fragment",
                "source_genome",
                "source_start",
                "source_stop",
                "assembly_start",
                "assembly_stop",
                "trim_left",
                "trim_right",
            ]) + "\n"
        )

        # Iterate over all assemblies
        for asm_idx, combo in enumerate(itertools.product(*fragments), start=1):
            assembly_id = f"assembly_{asm_idx}"
            manifest_row = [assembly_id]

            assembly_pos = 1  # 1-based assembly coordinates

            for i, frag in enumerate(combo):
                frag_id = frag.id
                genome_id, frag_num = frag_id.rsplit("_frag", 1)
                frag_num = int(frag_num)

                # Source genome coordinates (1-based, untrimmed)
                src_start, src_stop = boundaries[genome_id][frag_num - 1]

                frag_len = len(frag.seq)

                # Overlap trimming logic
                trim_left = overlap_lengths[i - 1] if i > 0 else 0
                trim_right = overlap_lengths[i] if i < n_frags - 1 else 0

                effective_len = frag_len - trim_left - trim_right

                asm_start = assembly_pos
                asm_stop = assembly_pos + effective_len - 1

                # ---- Assembly manifest (unchanged semantics) ----
                manifest_row.extend([
                    genome_id,
                    src_start,
                    src_stop,
                    asm_start,
                    asm_stop,
                ])

                # ---- Module table (untrimmed coordinates) ----
                mt.write(
                    "\t".join(map(str, [
                        assembly_id,
                        f"frag_{i+1}",
                        genome_id,
                        src_start,
                        src_stop,
                        asm_start,
                        asm_stop,
                        trim_left,
                        trim_right,
                    ])) + "\n"
                )

                assembly_pos = asm_stop + 1

            mf.write("\t".join(map(str, manifest_row)) + "\n")

    print(f"Assembly manifest written to {assembly_manifest_path}")
    print(f"Module table written to {module_table_path}")


def is_primer_ok(primer_seq, max_homopolymer=4, max_self_comp=4):
    """
    Check primer quality:
    - No runs of a single nucleotide longer than max_homopolymer
    - Self-complementarity below max_self_comp (approximate)
    """
    seq_str = str(primer_seq).upper()
    
    # check homopolymers
    for base in "ATGC":
        if base * max_homopolymer in seq_str:
            return False
    
    # simple self-complementarity check: max contiguous reverse complement matches
    rev_comp = str(Seq(seq_str).reverse_complement())
    max_match = 0
    for i in range(len(seq_str)):
        match_len = 0
        for j in range(len(seq_str)-i):
            if seq_str[i+j] == rev_comp[j]:
                match_len += 1
                max_match = max(max_match, match_len)
            else:
                match_len = 0
    if max_match >= max_self_comp:
        return False
    
    return True



def design_primers_for_fragment(frag_record, target_tm=60, tm_tol=1, max_primer_len=30, min_primer_len=18):
    """
    Returns forward and reverse primers adjusted to reach target Tm.
    Extends primers up to max_primer_len if necessary.
    """
    seq = frag_record.seq
    frag_len = len(seq)

    # --- Forward primer ---
    f_len = min(min_primer_len, frag_len)
    forward = seq[:f_len]
    tm_fwd = mt.Tm_NN(forward)

    while tm_fwd < target_tm - tm_tol and f_len < min(max_primer_len, frag_len):
        f_len += 1
        forward = seq[:f_len]
        tm_fwd = mt.Tm_NN(forward)

    # --- Reverse primer ---
    r_len = min(min_primer_len, frag_len)
    reverse = seq[-r_len:].reverse_complement()
    tm_rev = mt.Tm_NN(reverse)

    while tm_rev < target_tm - tm_tol and r_len < min(max_primer_len, frag_len):
        r_len += 1
        reverse = seq[-r_len:].reverse_complement()
        tm_rev = mt.Tm_NN(reverse)

    return str(forward), float(tm_fwd), str(reverse), float(tm_rev)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Generate assembly fragments from user indicated genomes "
            "with user indicated overlaps\n"
            "Example: (from results dir) \npython ../Scripts/make_assembly_frags.py "
            "--genomes ../../../fastas/bas14-15-16-18.fasta "
            "--gff_dir ./pharokka_results/single_gffs/ "
            "--kmers optimized_overlaps.tsv "
            "--fasta_outdir fragment_fastas "
            "--gff_outdir fragment_gffs"
                    ), formatter_class=argparse.RawDescriptionHelpFormatter)   
    parser.add_argument("--genomes", required=True, help="Combined genome FASTA")
    parser.add_argument("--gff_dir", required=True, help="Path to directory with GFFs from pharroka output")
    parser.add_argument("--kmers", required=True, help="optimized_overlaps.csv")
    parser.add_argument("--fasta_outdir", required=True, help="Output directory to save fragment FASTAS")
    parser.add_argument("--gff_outdir", required=True, help="Output directory to save fragment GFFs")
    args = parser.parse_args()

    def score_overlap(seq, other_overlaps=[]):
        """Lower score = better: self + cross complementarity"""
        self_comp = max_self_complementarity(seq)
        cross_comp = max([max_self_complementarity(seq, ov) for ov in other_overlaps] + [0])
        return self_comp + cross_comp

    def shared_by_all(row_accs):
        return all(acc in row_accs for acc in chosen_accessions)


    # --- load genome sequences and kmers df ---
    genome_seqs = list(SeqIO.parse(args.genomes, "fasta"))
    accessions = [seq.id for seq in genome_seqs]
    kmers_df = pd.read_csv(args.kmers)


    # --- prompt user for which accessions to make fragments from and how many frags ---
    print("Which genomes do you want to make fragments with?")
    for i, acc in enumerate(accessions, start=1):
        print(f"{i}. {acc}")

    selected = input("\nEnter the numbers of the genomes you want (comma-separated, e.g. 1,3,5): ")

    # --- parse user input ---
    try:
        selected_indices = [int(x.strip()) - 1 for x in selected.split(",")]
        chosen_accessions = [accessions[i] for i in selected_indices if 0 <= i < len(accessions)]
    except ValueError:
        print("Invalid input. Please enter numbers separated by commas.")
        return

    chosen_seqs = [seq for seq in genome_seqs if seq.id in chosen_accessions]

    n_frags = input("\nHow many assembly fragments? (integer, e.g. 3,4,5): ")
    n_frags = int(n_frags)
    print(f"\nMaking {n_frags} fragments with: {', '.join(chosen_accessions)}")


    # --- prompt user for assembly topology ---
    while True:
        topology = input("Select topology (1 for linear, 2 for circular assembly): ").strip()
        
        if topology == "1":
            topology = "linear"
            break
        elif topology == "2":
            topology = "circular"
            break
        else:
            print("Invalid choice. Please enter 1 or 2.\n")

    # --- prompt user for overlaps ---
    overlaps_user = {}
    for i in range(1, n_frags):
        ov = input(f"\nEnter overlap between fragment {i} and {i+1}: ").strip()
        overlaps_user[f"{i}-{i+1}"] = ov

    if topology == "circular":
        key = f"{n_frags}-1"
        ov = input(f"\nEnter overlap between fragment {n_frags} and 1: ").strip()
        overlaps_user[key] = ov

    overlap_lengths = []
    for i in range(1, n_frags):
        ov = overlaps_user.get(f"{i}-{i+1}")
        if ov is None:
            raise ValueError(f"Missing overlap for fragments {i}-{i+1}")
        overlap_lengths.append(len(ov))

    # --- use chosen overlaps for fragment creation ---
    fragments, boundaries = make_frags(chosen_seqs, n_frags, topology, overlaps_user)


    # make GFF for each fragment
    # make sure output directory exists
    os.makedirs(args.gff_outdir, exist_ok=True)
    for seq_id, frags in boundaries.items():
        gff_path = os.path.join(args.gff_dir, f"{seq_id}.gff")
        for frag_num, (frag_start, frag_stop) in enumerate(frags, start=1):
            out_path = os.path.join(args.gff_outdir, f"{seq_id}_frag{frag_num}.gff")
            split_gff(gff_path, out_path, frag_start, frag_stop, frag_num)



    # write each fragment group to its own FASTA
    # make sure output directory exists
    os.makedirs(args.fasta_outdir, exist_ok=True)
    for i, frag_group in enumerate(fragments, start=1):
        out_path = os.path.join(args.fasta_outdir, f"fragment_{i}.fasta")
        SeqIO.write(frag_group, out_path, "fasta")
        print(f"Wrote {out_path}")


    # --- make concatenated genome assemblies without duplicated overlaps ---
    assembly_outdir = os.path.join(args.fasta_outdir, "concatenated_assemblies")
    write_concatenated_assemblies_no_overlap(fragments, overlap_lengths, assembly_outdir)


    manifest_path = os.path.join(args.fasta_outdir, "assembly_manifest.tsv")
    module_table_path = os.path.join(args.fasta_outdir, "module_table.tsv")
    write_assembly_manifest_tsv(
        fragments=fragments,
        boundaries=boundaries,
        overlap_lengths=overlap_lengths,
        assembly_manifest_path=manifest_path, 
        module_table_path=module_table_path
    )


    # --- design primers ---
    design_primers = input("Design PCR primers? (y/n): ").strip().lower()
    if design_primers == "y":
        target_tm = float(input("Enter target Tm (°C, e.g., 60): "))

        primers_out = os.path.join(args.outdir, "primers.txt")
        unique_primers_out = os.path.join(args.outdir, "unique_primers.txt")
        unique_primers = defaultdict(lambda: {"fragments": [], "Tm": None})

        with open(primers_out, "w") as f:
            f.write("Fragment\tForward_Primer\tTm_Fwd\tReverse_Primer\tTm_Rev\n")
            for frag_group in fragments:
                for frag in frag_group:
                    fwd, tm_fwd, rev, tm_rev = design_primers_for_fragment(frag, target_tm=target_tm)
                    f.write(f"{frag.id}\t{fwd}\t{tm_fwd:.1f}\t{rev}\t{tm_rev:.1f}\n")

                    # add to unique primers
                    if fwd not in unique_primers:
                        unique_primers[fwd]["Tm"] = tm_fwd
                    unique_primers[fwd]["fragments"].append(frag.id)

                    if rev not in unique_primers:
                        unique_primers[rev]["Tm"] = tm_rev
                    unique_primers[rev]["fragments"].append(frag.id)

        print(f"✅ Primers saved to {primers_out}")

        # --- write unique primer table ---
        with open(unique_primers_out, "w") as f:
            f.write("Primer_Sequence\tTm\tFragments\n")
            for primer, info in unique_primers.items():
                frag_list = ",".join(info["fragments"])
                f.write(f"{primer}\t{info['Tm']:.1f}\t{frag_list}\n")
        print(f"✅ Unique primers saved to {unique_primers_out}")


    
if __name__ == "__main__":
    main()