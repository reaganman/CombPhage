#!/usr/bin/env python3

import argparse
import pandas as pd
import re
import os

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# --------------------------------------------------
# Helper Functions
# --------------------------------------------------

def calculate_identity(seq1, seq2):
    """
    Calculates identity ignoring positions where both sequences have gaps.
    Gaps vs Nucleotides are counted as mismatches.
    """
    if len(seq1) != len(seq2) or len(seq1) == 0:
        return 0.0
    
    seq1 = seq1.upper()
    seq2 = seq2.upper()

    matches = 0
    informative_positions = 0

    for a, b in zip(seq1, seq2):
        # If both are gaps, skip this position entirely (artificial inflation)
        if a == "-" and b == "-":
            continue
        
        informative_positions += 1
        
        # If they match and are not both gaps, count as match
        if a == b:
            matches += 1
            
    if informative_positions == 0:
        return 0.0

    return (matches / informative_positions) * 100

def parse_fragment_info(record):
    desc = record.description
    m = re.search(r"from (\S+) \((\d+):(\d+)\)", desc)
    if not m:
        raise ValueError(f"Cannot parse fragment header: {desc}")
    return m.group(1), int(m.group(2)), int(m.group(3))

def build_pos_map(aligned_seq):
    pos_map = {}
    genome_pos = 0
    for aln_pos, base in enumerate(aligned_seq):
        if base != "-":
            genome_pos += 1
            pos_map[genome_pos] = aln_pos
    return pos_map

def place_fragment_in_alignment(
    frag_seq, parent_id, start, stop, trim_left, trim_right, 
    parent_pos_map, parent_aln_seq, aln_length
):
    actual_start = start + trim_left
    actual_stop = stop - trim_right
    
    try:
        aln_start = parent_pos_map[actual_start]
        aln_stop = parent_pos_map[actual_stop]
    except KeyError:
        raise ValueError(f"{parent_id}: coord {actual_start}-{actual_stop} not in alignment")

    aligned = ["-"] * aln_length
    frag_ptr = trim_left 
    
    for pos in range(aln_start, aln_stop + 1):
        parent_base = parent_aln_seq[pos]
        if parent_base == "-":
            aligned[pos] = "-"
        else:
            if frag_ptr < len(frag_seq):
                aligned[pos] = frag_seq[frag_ptr]
                frag_ptr += 1
    return "".join(aligned)

def build_assemblies(module_df, frag_info, parent_pos_maps, aligned_parental_records):
    parent_seqs = {rec.id: str(rec.seq) for rec in aligned_parental_records}
    aln_length = len(aligned_parental_records[0].seq)
    assemblies = {}
    assembly_sources = {} 

    for asm_id, group in module_df.groupby("assembly"):
        assembly = ["-"] * aln_length
        sources = set()

        for _, row in group.iterrows():
            frag_id = row["fragment"]
            parent = row["source_genome"]
            sources.add(parent)
            
            info = frag_info[(parent, frag_id)]
            aligned_frag = place_fragment_in_alignment(
                info["seq"], parent, row["source_start"], row["source_stop"],
                row["trim_left"], row["trim_right"], parent_pos_maps[parent],
                parent_seqs[parent], aln_length
            )

            for i, base in enumerate(aligned_frag):
                if base == "-": continue
                if assembly[i] == "-":
                    assembly[i] = base
                elif assembly[i] != base:
                    raise ValueError(f"Conflict in {asm_id} at {i}: {assembly[i]} vs {base}")

        assemblies[asm_id] = "".join(assembly)
        assembly_sources[asm_id] = sources

    return assemblies, assembly_sources

# --------------------------------------------------
# Main
# --------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--frag_fasta", required=True)
    parser.add_argument("--parental_msa", required=True)
    parser.add_argument("--module_table", required=True)
    parser.add_argument("--out_fasta", required=True)
    parser.add_argument("--identity_out", required=True)
    parser.add_argument("--frag_identity_dir", required=True, help="Directory for fragment-specific matrices")
    parser.add_argument("--aligned_frag_fasta", required=False)
    args = parser.parse_args()

    # Ensure output directory exists
    if not os.path.exists(args.frag_identity_dir):
        os.makedirs(args.frag_identity_dir)

    frag_records = list(SeqIO.parse(args.frag_fasta, "fasta"))
    aligned_parental_records = list(SeqIO.parse(args.parental_msa, "fasta"))
    parent_seqs = {rec.id: str(rec.seq) for rec in aligned_parental_records}
    module_df = pd.read_csv(args.module_table, sep="\t")

    # 1. Map Fragments and Track by ID
    frag_info = {}
    frags_by_num = {} # { 'frag_1': { 'ParentA': 'aligned_seq', ... }, ... }

    parent_pos_maps = {rec.id: build_pos_map(str(rec.seq)) for rec in aligned_parental_records}

    for rec in frag_records:
        source, start, stop = parse_fragment_info(rec)
        m = re.search(r"frag(\d+)", rec.id)
        frag_num_id = f"frag_{m.group(1)}"
        
        # Store basic info
        frag_info[(source, frag_num_id)] = {"seq": str(rec.seq), "start": start, "stop": stop}
        
        # Calculate aligned version for fragment-to-fragment comparison
        aln_frag = place_fragment_in_alignment(
            str(rec.seq), source, start, stop, 0, 0,
            parent_pos_maps[source], parent_seqs[source], len(parent_seqs[source])
        )
        
        if frag_num_id not in frags_by_num:
            frags_by_num[frag_num_id] = {}
        frags_by_num[frag_num_id][source] = aln_frag

    # 2. Build Assemblies
    assemblies, assembly_sources = build_assemblies(
        module_df, frag_info, parent_pos_maps, aligned_parental_records
    )

    # 3. Pairwise Global Identity (Matrices)
    identity_entries = []
    assembly_records = []

    for asm_id in sorted(assemblies):
        sources = assembly_sources[asm_id]
        is_recomb = len(sources) > 1
        desc = "aligned recombinant assembly" if is_recomb else "aligned parental assembly"
        assembly_records.append(SeqRecord(Seq(assemblies[asm_id]), id=asm_id, description=desc))

        if is_recomb:
            for p_rec in aligned_parental_records:
                score = calculate_identity(assemblies[asm_id], str(p_rec.seq))
                identity_entries.append({"Query": asm_id, "Parent": p_rec.id, "Identity": round(score, 3)})

    for i, p1 in enumerate(aligned_parental_records):
        for p2 in aligned_parental_records:
            score = calculate_identity(str(p1.seq), str(p2.seq))
            identity_entries.append({"Query": p1.id, "Parent": p2.id, "Identity": round(score, 3)})

    SeqIO.write(assembly_records, args.out_fasta, "fasta")
    pd.DataFrame(identity_entries).pivot(index="Query", columns="Parent", values="Identity").to_csv(args.identity_out, sep="\t")

    # 4. Fragment-Specific Pairwise Identities
    for frag_id, parents_dict in frags_by_num.items():
        frag_matrix_data = []
        parent_list = sorted(parents_dict.keys())
        
        for p1 in parent_list:
            for p2 in parent_list:
                score = calculate_identity(parents_dict[p1], parents_dict[p2])
                frag_matrix_data.append({"Parent_A": p1, "Parent_B": p2, "Identity": round(score, 3)})
        
        if frag_matrix_data:
            frag_df = pd.DataFrame(frag_matrix_data).pivot(index="Parent_A", columns="Parent_B", values="Identity")
            out_path = os.path.join(args.frag_identity_dir, f"{frag_id}_identity.tsv")
            frag_df.to_csv(out_path, sep="\t")

    # 5. Optional Aligned Frag Output
    if args.aligned_frag_fasta:
        aligned_frag_records = []
        for frag_id, parents_dict in frags_by_num.items():
            for parent, seq in parents_dict.items():
                rec = SeqRecord(Seq(seq), id=f"{frag_id}_{parent}", description=f"aligned {frag_id} from {parent}")
                aligned_frag_records.append(rec)
        SeqIO.write(aligned_frag_records, args.aligned_frag_fasta, "fasta")

if __name__ == "__main__":
    main()



"""
python Scripts/analyze_assemblies.py --frag_fasta results_bas14-15-16-18/assmeblies/all_frags.fasta --parental_msa ~/phage_genomics/fastas/bas14-15-16-18_msa.fsata --module_table results_bas14-15-16-18/assmeblies/module_table.tsv --out_fasta results_bas14-15-16-18/assmeblies/aligned_assemblies.fasta --aligned_frag_fasta test_frags_msa.fasta --identity_out test_id.tsv --frag_identity_dir test_frag_ids
"""