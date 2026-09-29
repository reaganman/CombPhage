import argparse
import pandas as pd
import numpy as np
import itertools
import math
import sys

# --- CONFIGURATION ---
IDEAL_OVERLAP_LENGTH = 40
IDEAL_GC_CONTENT = 0.50

# Constraints
MAX_FRAGMENT_SIZE = 10000

# Weights for scoring
WEIGHT_GROUP_STD_DEV = 1.0
WEIGHT_INDIVIDUAL_PENALTY = 20.0 

SHORT_OVERLAP_MULTIPLIER = 50.0 
GC_MULTIPLIER = 200.0
LARGE_FRAGMENT_MULTIPLIER = 2.0 # penalty for exceeding MAX_FRAGMENT_SIZE

# Number of top combinations to save
TOP_N_RESULTS = 3

def calculate_gc(seq):
    if not seq:
        return 0
    seq = seq.upper()
    return (seq.count('G') + seq.count('C')) / len(seq)

def calculate_individual_penalty(seq):
    """Penalize overlaps that deviate from ideal length and ideal GC content."""
    length = len(seq)
    gc = calculate_gc(seq)
    
    if length < IDEAL_OVERLAP_LENGTH:
        len_penalty = (IDEAL_OVERLAP_LENGTH - length) * SHORT_OVERLAP_MULTIPLIER
    else:
        len_penalty = length - IDEAL_OVERLAP_LENGTH
        
    gc_penalty = abs(gc - IDEAL_GC_CONTENT) * GC_MULTIPLIER 
    
    return len_penalty + gc_penalty

def optimize_overlap(ov, genomes):
    """
    Scans overlaps longer than IDEAL_OVERLAP_LENGTH with a sliding window.
    Updates the coordinates for all genomes based on the best window.
    """
    seq = ov['seq']
    orig_length = len(seq)
    
    if orig_length <= IDEAL_OVERLAP_LENGTH:
        ov['penalty'] = calculate_individual_penalty(seq)
        ov['orig_coords'] = ov['coords']
        ov['new_coords'] = ov['coords']
        return ov
        
    best_penalty = float('inf')
    best_subseq_data = None
    
    window_size = IDEAL_OVERLAP_LENGTH
    for i in range(orig_length - window_size + 1):
        subseq = seq[i : i + window_size]
        gc = calculate_gc(subseq)
        
        gc_penalty = abs(gc - IDEAL_GC_CONTENT) * 100 
        penalty = gc_penalty 
        
        if penalty < best_penalty:
            best_penalty = penalty
            
            new_coords = {}
            for g in genomes:
                new_start = ov['coords'][g]['start'] + i
                new_stop = new_start + window_size
                new_coords[g] = {
                    'start': new_start,
                    'stop': new_stop,
                    'center': new_start + (window_size / 2.0)
                }
            
            best_subseq_data = {
                'id': ov['id'],
                'orig_coords': ov['coords'],
                'new_coords': new_coords,
                'seq': subseq,
                'length': window_size,
                'gc': gc,
                'penalty': penalty
            }
            
    return best_subseq_data

def evaluate_combination(combo, num_fragments, is_circular, genomes, genome_lengths):
    all_std_devs = []
    all_frag_sizes = {}
    total_large_frag_penalty = 0
    
    # 1. Calculate Fragment Sizes for EACH genome
    for g in genomes:
        positions = [ov['new_coords'][g]['center'] for ov in combo]
        positions.sort() 
        
        fragment_sizes = []
        if is_circular:
            for i in range(len(positions)):
                if i == len(positions) - 1:
                    frag_len = (genome_lengths[g] - positions[i]) + positions[0]
                else:
                    frag_len = positions[i+1] - positions[i]
                fragment_sizes.append(frag_len)
        else:
            fragment_sizes.append(positions[0]) 
            for i in range(len(positions) - 1):
                fragment_sizes.append(positions[i+1] - positions[i])
            fragment_sizes.append(genome_lengths[g] - positions[-1])
            
        # Check for fragments exceeding the maximum allowed size
        for frag_len in fragment_sizes:
            if frag_len > MAX_FRAGMENT_SIZE:
                # Add a compounding penalty based on how many bp over the limit it is
                total_large_frag_penalty += (frag_len - MAX_FRAGMENT_SIZE) * LARGE_FRAGMENT_MULTIPLIER
            
        all_std_devs.append(np.std(fragment_sizes))
        all_frag_sizes[g] = fragment_sizes
        
    group_penalty = np.mean(all_std_devs)
    
    # 2. Get Individual Penalties (Pre-calculated)
    avg_indiv_penalty = np.mean([ov['penalty'] for ov in combo])
    
    # 3. Total Score (Includes the new large fragment penalty)
    total_score = (WEIGHT_GROUP_STD_DEV * group_penalty) + (WEIGHT_INDIVIDUAL_PENALTY * avg_indiv_penalty) + total_large_frag_penalty
    
    return total_score, all_std_devs, all_frag_sizes


def parse_fasta_lengths(fasta_path):
    """
    Parse a multi-FASTA file and return a dict mapping each sequence header
    (the full description line, minus the leading '>') to its sequence length.
    """
    lengths = {}
    current_header = None
    current_length = 0

    with open(fasta_path, 'r') as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_header is not None:
                    lengths[current_header] = current_length
                current_header = line[1:]  # strip '>'
                current_length = 0
            else:
                current_length += len(line)

    if current_header is not None:
        lengths[current_header] = current_length

    return lengths


def match_genome_to_fasta(genome_name, fasta_lengths):
    """
    Find the best FASTA header for a genome name derived from the TSV columns.
    Tries, in order:
      1. Exact match against the full header.
      2. Exact match against the first whitespace-delimited token of the header.
      3. Case-insensitive substring match against the full header.
    Returns the matched length or None if no match is found.
    """
    genome_lower = genome_name.lower()

    for header, length in fasta_lengths.items():
        if header == genome_name:
            return length

    for header, length in fasta_lengths.items():
        first_token = header.split()[0]
        if first_token == genome_name:
            return length

    for header, length in fasta_lengths.items():
        if genome_lower in header.lower():
            return length

    return None


def main():
    parser = argparse.ArgumentParser(
        description="Optimize Gibson assembly overlap positions across multiple genomes."
    )
    parser.add_argument("tsv", help="Path to the overlap candidate TSV file.")
    parser.add_argument("fasta", help="Path to the FASTA file containing genome sequences.")
    args = parser.parse_args()

    try:
        df = pd.read_csv(args.tsv, sep='\t')
    except Exception as e:
        print(f"Error reading TSV file: {e}")
        sys.exit(1)
        
    start_cols = [c for c in df.columns if 'start' in c.lower()]
    genomes = [c.rsplit('_', 2)[0] for c in start_cols] 
    print(f"\nIdentified {len(genomes)} genomes in the TSV: {', '.join(genomes)}")
    
    raw_overlaps = []
    for index, row in df.iterrows():
        coords = {}
        for g in genomes:
            start_col = f"{g}_ov_start"
            stop_col = f"{g}_ov_stop"
            coords[g] = {
                'start': row[start_col],
                'stop': row[stop_col],
                'center': (row[start_col] + row[stop_col]) / 2.0
            }
            
        seq = str(row['ov_seq'])
        raw_overlaps.append({
            'id': f"{row['cluster']}_{row['region']}",
            'coords': coords,
            'seq': seq,
            'length': len(seq),
            'gc': calculate_gc(seq)
        })
        
    print("Optimizing individual candidate sequences...")
    overlaps = [optimize_overlap(ov, genomes) for ov in raw_overlaps]
        
    primary_genome = genomes[0]
    overlaps = sorted(overlaps, key=lambda x: x['new_coords'][primary_genome]['center'])

    num_fragments = int(input("\nHow many fragments do you want to assemble? "))
    
    topo_input = input("Is the assembly linear or circular? (l/c): ").strip().lower()
    is_circular = True if topo_input.startswith('c') else False

    # --- FASTA-based genome length resolution ---
    try:
        fasta_lengths = parse_fasta_lengths(args.fasta)
    except Exception as e:
        print(f"Error reading FASTA file: {e}")
        sys.exit(1)

    print(f"  Parsed {len(fasta_lengths)} sequence(s) from FASTA.")

    genome_lengths = {}
    print("\nResolving genome lengths from FASTA:")
    for g in genomes:
        matched_length = match_genome_to_fasta(g, fasta_lengths)
        max_coord = max(ov['orig_coords'][g]['stop'] for ov in overlaps)

        if matched_length is not None:
            genome_lengths[g] = matched_length
            print(f"  {g}: {matched_length:,} bp (matched from FASTA)")
        else:
            genome_lengths[g] = max_coord
            print(f"  WARNING: '{g}' not found in FASTA. Falling back to max overlap coordinate: {max_coord:,} bp")
    # --- end FASTA block ---

    num_overlaps_needed = num_fragments if is_circular else num_fragments - 1
    
    if num_overlaps_needed > len(overlaps):
        print(f"\nError: You want {num_fragments} fragments, requiring {num_overlaps_needed} overlaps, but only have {len(overlaps)} candidates.")
        sys.exit(1)
        
    total_combinations = math.comb(len(overlaps), num_overlaps_needed)
    print(f"\nAnalyzing {total_combinations} possible combinations...")
    
    top_results = []
    
    for combo in itertools.combinations(overlaps, num_overlaps_needed):
        score, std_devs, frag_sizes = evaluate_combination(combo, num_fragments, is_circular, genomes, genome_lengths)
        
        top_results.append({
            'score': score,
            'combo': combo,
            'stats': (std_devs, frag_sizes)
        })
        
        top_results = sorted(top_results, key=lambda x: x['score'])[:TOP_N_RESULTS]
            
    print("\n" + "="*70)
    print(f"            TOP {TOP_N_RESULTS} OPTIMAL OVERLAP GROUPS")
    print("="*70)
    
    for rank, result in enumerate(top_results, 1):
        print(f"\n--- RANK {rank} (Score: {result['score']:.2f}) ---")
        std_devs, frag_sizes = result['stats']
        for i, g in enumerate(genomes):
            print(f"  {g} - Fragment Size Std Dev: {std_devs[i]:.1f}")
            
    print("\nGenerating TSV files...")
    
    def save_tsv(rank, combo):
        filename = f"optimized_overlaps_rank{rank}.tsv"
        header = ["cluster_region", "optimized_seq", "length", "gc_content"]
        for g in genomes:
            header.extend([f"{g}_start", f"{g}_stop"])
            
        rows = []
        for ov in combo:
            row_data = [
                ov['id'],
                ov['seq'],
                str(ov['length']),
                f"{ov['gc']:.3f}"
            ]
            for g in genomes:
                row_data.extend([
                    str(ov['new_coords'][g]['start']),
                    str(ov['new_coords'][g]['stop'])
                ])
            rows.append("\t".join(row_data))
            
        with open(filename, 'w') as f:
            f.write("\t".join(header) + "\n")
            f.write("\n".join(rows) + "\n")
            
        print(f"  Saved: {filename}")

    # Initial TSV Generation
    for rank, result in enumerate(top_results, 1):
        save_tsv(rank, result['combo'])

    # --- INTERACTIVE OVERLAP EDITING BLOCK ---
    while True:
        edit_choice = input("\nWould you like to manually change any of the overlap sets? (y/n): ").strip().lower()
        if edit_choice != 'y':
            break
            
        rank_to_edit = input(f"Which set of overlaps would you like to edit? (1-{len(top_results)}): ").strip()
        try:
            rank_idx = int(rank_to_edit) - 1
            if rank_idx < 0 or rank_idx >= len(top_results):
                print("Invalid selection. Please choose a valid rank number.")
                continue
        except ValueError:
            print("Invalid input. Please enter a number.")
            continue

        target_combo = list(top_results[rank_idx]['combo'])
        print(f"\nCurrent overlap regions in Rank {rank_idx + 1}:")
        for ov in target_combo:
            print(f" - {ov['id']}")

        old_cluster = input("\nEnter the cluster_region you want to change: ").strip()
        
        # Verify the old cluster region exists in the selected combo
        old_idx = -1
        old_ov = None
        for i, ov in enumerate(target_combo):
            if ov['id'] == old_cluster:
                old_idx = i
                old_ov = ov
                break
                
        if old_idx == -1:
            print(f"Error: '{old_cluster}' was not found in the Rank {rank_idx + 1} set.")
            continue

        # Find nearby candidate cluster regions based on distance to the old overlap's center
        old_center = old_ov['orig_coords'][primary_genome]['center']
        current_ids = set([o['id'] for o in target_combo])
        
        available_raw = [r for r in raw_overlaps if r['id'] not in current_ids]
        
        def get_distance(raw_ov):
            return abs(raw_ov['coords'][primary_genome]['center'] - old_center)
            
        available_raw.sort(key=get_distance)
        
        print("\nNearby candidate cluster_regions (Top 10 closest):")
        for r in available_raw[:10]:
            dist = get_distance(r)
            print(f" -> {r['id']} (Distance: {dist:.1f} bp)")
            print(f"    Raw Seq: {r['seq']}")
            
        new_cluster = input("\nEnter the ID of the candidate you want to use: ").strip()
        
        chosen_raw = next((r for r in available_raw if r['id'] == new_cluster), None)
        if not chosen_raw:
            print(f"Error: '{new_cluster}' is not a valid candidate ID from the list.")
            continue
            
        new_seq = input("Enter the specific subsequence you want to use from this candidate's raw seq: ").strip().upper()

        # Validate that the requested sequence is actually a subsequence of the raw overlap
        raw_seq = chosen_raw['seq'].upper()
        subseq_idx = raw_seq.find(new_seq)
        
        if subseq_idx == -1:
            print(f"\nError: The sequence provided is not a valid subsequence of '{new_cluster}'.")
            continue
            
        # Dynamically build the new overlap object using the specific substring index
        window_size = len(new_seq)
        new_coords = {}
        for g in genomes:
            new_start = chosen_raw['coords'][g]['start'] + subseq_idx
            new_stop = new_start + window_size
            new_coords[g] = {
                'start': new_start,
                'stop': new_stop,
                'center': new_start + (window_size / 2.0)
            }
            
        custom_gc = calculate_gc(new_seq)
        custom_penalty = calculate_individual_penalty(new_seq)
        
        valid_replacement = {
            'id': chosen_raw['id'],
            'orig_coords': chosen_raw['coords'],
            'new_coords': new_coords,
            'seq': new_seq,
            'length': window_size,
            'gc': custom_gc,
            'penalty': custom_penalty
        }

        # Execute Swap
        target_combo[old_idx] = valid_replacement
        top_results[rank_idx]['combo'] = tuple(target_combo)
        
        # Re-evaluate the score internally for accurate data
        new_score, new_std_devs, new_frag_sizes = evaluate_combination(
            top_results[rank_idx]['combo'], num_fragments, is_circular, genomes, genome_lengths
        )
        top_results[rank_idx]['score'] = new_score
        top_results[rank_idx]['stats'] = (new_std_devs, new_frag_sizes)
        
        print(f"\nSuccessfully swapped '{old_cluster}' for '{new_cluster}'.")
        print(f"New combo score: {new_score:.2f}. Updating output file...")
        
        save_tsv(rank_idx + 1, top_results[rank_idx]['combo'])

if __name__ == "__main__":
    main() 