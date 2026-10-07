#!/usr/bin/env python3
from __future__ import annotations

import argparse
import heapq
import itertools
import json
import math
import shutil
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from Bio import SeqIO

DEFAULT_IDEAL_OVERLAP_LENGTH = 40
DEFAULT_IDEAL_GC = 0.50
DEFAULT_MAX_FRAGMENT_SIZE = 10000
DEFAULT_TOP_N = 3
DEFAULT_WEIGHT_GROUP_STD_DEV = 1.0
DEFAULT_WEIGHT_INDIVIDUAL_PENALTY = 20.0
DEFAULT_SHORT_OVERLAP_MULTIPLIER = 50.0
DEFAULT_LONG_OVERLAP_MULTIPLIER = 1.0
DEFAULT_GC_MULTIPLIER = 200.0
DEFAULT_LARGE_FRAGMENT_MULTIPLIER = 2.0
DEFAULT_MAX_COMBINATIONS = 2000000
DEFAULT_MAX_CHIMERAS = 10000


def calculate_gc(seq: str) -> float:
    seq = str(seq).upper()
    if not seq:
        return 0.0
    return (seq.count('G') + seq.count('C')) / len(seq)


def normalize_frame(value):
    if value is None or pd.isna(value):
        return None
    text = str(value).strip().lower()
    if text in {'+', '+1', '1', 'forward', 'fwd'}:
        return '+'
    if text in {'-', '-1', 'reverse', 'rev'}:
        return '-'
    return None


def calculate_individual_penalty(
    seq: str,
    ideal_length: int = DEFAULT_IDEAL_OVERLAP_LENGTH,
    ideal_gc: float = DEFAULT_IDEAL_GC,
    short_overlap_multiplier: float = DEFAULT_SHORT_OVERLAP_MULTIPLIER,
    long_overlap_multiplier: float = DEFAULT_LONG_OVERLAP_MULTIPLIER,
    gc_multiplier: float = DEFAULT_GC_MULTIPLIER,
) -> float:
    length = len(seq)
    gc = calculate_gc(seq)

    if length < ideal_length:
        length_penalty = (ideal_length - length) * short_overlap_multiplier
    else:
        length_penalty = (length - ideal_length) * long_overlap_multiplier

    gc_penalty = abs(gc - ideal_gc) * gc_multiplier
    return length_penalty + gc_penalty


def read_genomes_fasta(fasta_path: str | Path) -> Dict[str, str]:
    sequences = {}
    for record in SeqIO.parse(str(fasta_path), 'fasta'):
        if record.id in sequences:
            raise ValueError(f'Duplicate FASTA record ID: {record.id}')
        sequences[record.id] = str(record.seq).upper()

    if not sequences:
        raise ValueError(f'No FASTA sequences found in {fasta_path}')
    return sequences


def detect_genomes(overlap_df: pd.DataFrame) -> List[str]:
    suffix = '_ov_start'
    genomes = [
        col[:-len(suffix)]
        for col in overlap_df.columns
        if col.endswith(suffix)
    ]
    genomes = list(dict.fromkeys(genomes))

    if not genomes:
        raise ValueError("No columns ending in '_ov_start' were found.")

    for genome in genomes:
        if f'{genome}_ov_stop' not in overlap_df.columns:
            raise ValueError(f'Missing column: {genome}_ov_stop')

    return genomes


def interval_length(start: int, stop: int, seq_len: int) -> int:
    start = int(start)
    stop = int(stop)

    if start < stop:
        return stop - start
    if start == stop:
        return seq_len
    return (seq_len - start + 1) + (stop - 1)


def interval_center(start: int, stop: int, seq_len: int) -> float:
    length = interval_length(start, stop, seq_len)
    return ((start - 1 + length / 2.0) % seq_len) + 1.0


def advance_start(start: int, amount: int, seq_len: int) -> int:
    return ((int(start) - 1 + int(amount)) % seq_len) + 1


def retreat_stop(stop: int, amount: int, seq_len: int) -> int:
    stop = int(stop)
    amount = int(amount)
    candidate = stop - amount
    if 1 <= candidate <= seq_len + 1:
        return candidate
    return ((stop - 1 - amount) % seq_len) + 1


def optimize_overlap_candidate(
    candidate: dict,
    genome_lengths: Dict[str, int],
    ideal_length: int,
    ideal_gc: float,
    short_overlap_multiplier: float = DEFAULT_SHORT_OVERLAP_MULTIPLIER,
    long_overlap_multiplier: float = DEFAULT_LONG_OVERLAP_MULTIPLIER,
    gc_multiplier: float = DEFAULT_GC_MULTIPLIER,
) -> dict:
    sequence = candidate['original_seq']
    original_length = len(sequence)

    if original_length <= ideal_length:
        best_offset = 0
        best_seq = sequence
        best_penalty = calculate_individual_penalty(
            sequence,
            ideal_length=ideal_length,
            ideal_gc=ideal_gc,
            short_overlap_multiplier=short_overlap_multiplier,
            long_overlap_multiplier=long_overlap_multiplier,
            gc_multiplier=gc_multiplier,
        )
    else:
        best_offset = 0
        best_seq = sequence[:ideal_length]
        best_penalty = float('inf')

        for offset in range(original_length - ideal_length + 1):
            subseq = sequence[offset:offset + ideal_length]
            penalty = calculate_individual_penalty(
                subseq,
                ideal_length=ideal_length,
                ideal_gc=ideal_gc,
                short_overlap_multiplier=short_overlap_multiplier,
                long_overlap_multiplier=long_overlap_multiplier,
                gc_multiplier=gc_multiplier,
            )
            if penalty < best_penalty:
                best_offset = offset
                best_seq = subseq
                best_penalty = penalty

    left_trim = best_offset
    right_trim = original_length - left_trim - len(best_seq)
    reference_frame = candidate.get('reference_frame')
    optimized_coords = {}

    for genome, coord in candidate['original_coords'].items():
        seq_len = genome_lengths[genome]
        genome_frame = coord.get('frame')

        same_orientation = True
        if reference_frame is not None and genome_frame is not None:
            same_orientation = reference_frame == genome_frame

        if same_orientation:
            new_start = advance_start(coord['start'], left_trim, seq_len)
            new_stop = retreat_stop(coord['stop'], right_trim, seq_len)
        else:
            new_start = advance_start(coord['start'], right_trim, seq_len)
            new_stop = retreat_stop(coord['stop'], left_trim, seq_len)

        optimized_coords[genome] = {
            'start': int(new_start),
            'stop': int(new_stop),
            'center': interval_center(int(new_start), int(new_stop), seq_len),
            'frame': genome_frame,
        }

    result = dict(candidate)
    result.update({
        'seq': best_seq,
        'length': len(best_seq),
        'gc': calculate_gc(best_seq),
        'penalty': float(best_penalty),
        'left_trim': left_trim,
        'right_trim': right_trim,
        'coords': optimized_coords,
    })
    return result


def load_overlap_candidates(
    overlap_path: str | Path,
    genome_sequences: Dict[str, str],
    ideal_length: int,
    ideal_gc: float,
    short_overlap_multiplier: float = DEFAULT_SHORT_OVERLAP_MULTIPLIER,
    long_overlap_multiplier: float = DEFAULT_LONG_OVERLAP_MULTIPLIER,
    gc_multiplier: float = DEFAULT_GC_MULTIPLIER,
) -> Tuple[List[dict], List[str]]:
    df = pd.read_csv(overlap_path, sep='\t')
    required = {'cluster', 'region', 'ov_seq'}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            'Filtered overlap table is missing: ' + ', '.join(sorted(missing))
        )

    genomes = detect_genomes(df)
    missing_genomes = [g for g in genomes if g not in genome_sequences]
    if missing_genomes:
        raise ValueError(
            'Genome IDs not found in genomes.fasta: ' + ', '.join(missing_genomes)
        )

    genome_lengths = {g: len(genome_sequences[g]) for g in genomes}
    candidates = []

    for row_index, row in df.iterrows():
        seq = str(row['ov_seq']).upper()
        if not seq or seq == 'NAN':
            continue

        cluster = str(row['cluster'])
        region = str(row['region'])
        reference_frame = (
            normalize_frame(row['reference_frame'])
            if 'reference_frame' in df.columns
            else None
        )

        coords = {}
        for genome in genomes:
            start = int(row[f'{genome}_ov_start'])
            stop = int(row[f'{genome}_ov_stop'])
            frame_col = f'{genome}_frame'
            coords[genome] = {
                'start': start,
                'stop': stop,
                'center': interval_center(start, stop, genome_lengths[genome]),
                'frame': (
                    normalize_frame(row[frame_col])
                    if frame_col in df.columns
                    else None
                ),
            }

        candidate = {
            'id': f'{cluster}_{region}',
            'cluster': cluster,
            'region': region,
            'source_row': int(row_index),
            'original_seq': seq,
            'original_length': len(seq),
            'original_gc': calculate_gc(seq),
            'original_coords': coords,
            'reference_frame': reference_frame,
        }
        candidates.append(
            optimize_overlap_candidate(
                candidate,
                genome_lengths,
                ideal_length,
                ideal_gc,
                short_overlap_multiplier,
                long_overlap_multiplier,
                gc_multiplier,
            )
        )

    if not candidates:
        raise ValueError('No usable overlap candidates were found.')

    return candidates, genomes


def prepare_design_context(
    overlaps: str | Path,
    fasta: str | Path,
    ideal_overlap_length: int = DEFAULT_IDEAL_OVERLAP_LENGTH,
    ideal_gc: float = DEFAULT_IDEAL_GC,
    short_overlap_multiplier: float = DEFAULT_SHORT_OVERLAP_MULTIPLIER,
    long_overlap_multiplier: float = DEFAULT_LONG_OVERLAP_MULTIPLIER,
    gc_multiplier: float = DEFAULT_GC_MULTIPLIER,
):
    """
    Load and optimize candidate overlaps once for reuse by a GUI or other
    interactive client.

    The returned context can be passed to rank_designs(), design_from_ids(),
    and materialize_design() without rereading the input files.
    """
    genome_sequences = read_genomes_fasta(fasta)
    candidates, genomes = load_overlap_candidates(
        overlaps,
        genome_sequences,
        ideal_overlap_length,
        ideal_gc,
        short_overlap_multiplier,
        long_overlap_multiplier,
        gc_multiplier,
    )
    genome_lengths = {
        genome: len(genome_sequences[genome])
        for genome in genomes
    }
    candidate_map = {candidate['id']: candidate for candidate in candidates}

    if len(candidate_map) != len(candidates):
        duplicates = pd.Series(
            [candidate['id'] for candidate in candidates]
        ).value_counts()
        duplicates = duplicates[duplicates > 1].index.tolist()
        raise ValueError(
            'Candidate overlap IDs are not unique: ' + ', '.join(duplicates)
        )

    return {
        'genome_sequences': genome_sequences,
        'genomes': genomes,
        'genome_lengths': genome_lengths,
        'candidates': candidates,
        'candidate_map': candidate_map,
        'primary_genome': genomes[0],
        'ideal_overlap_length': ideal_overlap_length,
        'ideal_gc': ideal_gc,
        'short_overlap_multiplier': short_overlap_multiplier,
        'long_overlap_multiplier': long_overlap_multiplier,
        'gc_multiplier': gc_multiplier,
    }


def rank_designs(
    context: dict,
    num_fragments: int,
    topology: str,
    top_n: int = DEFAULT_TOP_N,
    max_fragment_size: int = DEFAULT_MAX_FRAGMENT_SIZE,
    max_combinations: int = DEFAULT_MAX_COMBINATIONS,
    weight_group_std_dev: float = DEFAULT_WEIGHT_GROUP_STD_DEV,
    weight_individual_penalty: float = DEFAULT_WEIGHT_INDIVIDUAL_PENALTY,
    large_fragment_multiplier: float = DEFAULT_LARGE_FRAGMENT_MULTIPLIER,
):
    """Rank overlap combinations from a prepared design context."""
    return find_top_designs(
        context['candidates'],
        context['genomes'],
        context['genome_lengths'],
        num_fragments,
        topology,
        top_n=top_n,
        max_fragment_size=max_fragment_size,
        max_combinations=max_combinations,
        weight_group_std_dev=weight_group_std_dev,
        weight_individual_penalty=weight_individual_penalty,
        large_fragment_multiplier=large_fragment_multiplier,
    )



def trim_candidate(
    candidate: dict,
    genome_lengths: Dict[str, int],
    left_trim: int = 0,
    right_trim: int = 0,
    ideal_overlap_length: int = DEFAULT_IDEAL_OVERLAP_LENGTH,
    ideal_gc: float = DEFAULT_IDEAL_GC,
) -> dict:
    """
    Return a copy of an optimized overlap after additional manual trimming.

    ``left_trim`` and ``right_trim`` are defined relative to the displayed
    overlap sequence. Coordinate updates are strand-aware in the same way as
    the initial automatic overlap optimization.
    """
    left_trim = int(left_trim)
    right_trim = int(right_trim)

    if left_trim < 0 or right_trim < 0:
        raise ValueError('Manual trim amounts cannot be negative.')

    sequence = str(candidate['seq']).upper()
    if left_trim + right_trim >= len(sequence):
        raise ValueError(
            'Manual trimming must leave at least one nucleotide in the overlap.'
        )

    end = len(sequence) - right_trim if right_trim else len(sequence)
    trimmed_seq = sequence[left_trim:end]
    reference_frame = candidate.get('reference_frame')
    trimmed_coords = {}

    for genome, coord in candidate['coords'].items():
        seq_len = int(genome_lengths[genome])
        genome_frame = coord.get('frame')

        same_orientation = True
        if reference_frame is not None and genome_frame is not None:
            same_orientation = reference_frame == genome_frame

        if same_orientation:
            new_start = advance_start(coord['start'], left_trim, seq_len)
            new_stop = retreat_stop(coord['stop'], right_trim, seq_len)
        else:
            new_start = advance_start(coord['start'], right_trim, seq_len)
            new_stop = retreat_stop(coord['stop'], left_trim, seq_len)

        trimmed_coords[genome] = {
            'start': int(new_start),
            'stop': int(new_stop),
            'center': interval_center(int(new_start), int(new_stop), seq_len),
            'frame': genome_frame,
        }

    result = dict(candidate)
    result.update({
        'seq': trimmed_seq,
        'length': len(trimmed_seq),
        'gc': calculate_gc(trimmed_seq),
        'penalty': calculate_individual_penalty(
            trimmed_seq,
            ideal_length=ideal_overlap_length,
            ideal_gc=ideal_gc,
        ),
        'coords': trimmed_coords,
        'manual_left_trim': left_trim,
        'manual_right_trim': right_trim,
        'left_trim': int(candidate.get('left_trim', 0)) + left_trim,
        'right_trim': int(candidate.get('right_trim', 0)) + right_trim,
    })
    return result



def candidate_from_subsequence(
    candidate: dict,
    genome_lengths: Dict[str, int],
    subsequence: str,
    ideal_overlap_length: int = DEFAULT_IDEAL_OVERLAP_LENGTH,
    ideal_gc: float = DEFAULT_IDEAL_GC,
    min_length: int = 1,
    short_overlap_multiplier: float = DEFAULT_SHORT_OVERLAP_MULTIPLIER,
    long_overlap_multiplier: float = DEFAULT_LONG_OVERLAP_MULTIPLIER,
    gc_multiplier: float = DEFAULT_GC_MULTIPLIER,
) -> dict:
    """Return a candidate using an exact subsequence of its full overlap.

    The user-facing sequence is selected from ``candidate['original_seq']``
    rather than from the automatically optimized sequence.  Coordinates are
    recalculated from ``original_coords`` so any valid contiguous subsequence
    of the full filtered overlap can be used without accumulating coordinate
    shifts from previous edits.

    If the requested sequence occurs more than once in the full overlap, the
    location is ambiguous and the edit is rejected rather than silently using
    the first occurrence.
    """
    full_sequence = str(candidate['original_seq']).upper()
    subsequence = ''.join(str(subsequence).split()).upper()
    min_length = int(min_length)

    if not subsequence:
        raise ValueError('Enter a non-empty overlap subsequence.')

    if len(subsequence) < min_length:
        raise ValueError(
            f'Overlap subsequence must be at least {min_length} bp long.'
        )

    if len(subsequence) > len(full_sequence):
        raise ValueError(
            'Overlap subsequence cannot be longer than the full overlap.'
        )

    offsets = [
        offset
        for offset in range(len(full_sequence) - len(subsequence) + 1)
        if full_sequence.startswith(subsequence, offset)
    ]

    if not offsets:
        raise ValueError(
            'The working sequence must be an exact contiguous subsequence of '
            'the full overlap shown above.'
        )

    if len(offsets) > 1:
        raise ValueError(
            f'This sequence occurs {len(offsets)} times in the full overlap, '
            'so its genomic position is ambiguous. Choose a longer or unique '
            'subsequence.'
        )

    left_trim = int(offsets[0])
    right_trim = len(full_sequence) - left_trim - len(subsequence)
    reference_frame = candidate.get('reference_frame')
    edited_coords = {}

    for genome, coord in candidate['original_coords'].items():
        seq_len = int(genome_lengths[genome])
        genome_frame = coord.get('frame')

        same_orientation = True
        if reference_frame is not None and genome_frame is not None:
            same_orientation = reference_frame == genome_frame

        if same_orientation:
            new_start = advance_start(coord['start'], left_trim, seq_len)
            new_stop = retreat_stop(coord['stop'], right_trim, seq_len)
        else:
            # The displayed sequence is oriented relative to the reference.
            # Therefore the sequence's left and right offsets are reversed in
            # genomic coordinates for homologous features on the other strand.
            new_start = advance_start(coord['start'], right_trim, seq_len)
            new_stop = retreat_stop(coord['stop'], left_trim, seq_len)

        edited_coords[genome] = {
            'start': int(new_start),
            'stop': int(new_stop),
            'center': interval_center(int(new_start), int(new_stop), seq_len),
            'frame': genome_frame,
        }

    result = dict(candidate)
    result.update({
        'seq': subsequence,
        'length': len(subsequence),
        'gc': calculate_gc(subsequence),
        'penalty': calculate_individual_penalty(
            subsequence,
            ideal_length=ideal_overlap_length,
            ideal_gc=ideal_gc,
            short_overlap_multiplier=short_overlap_multiplier,
            long_overlap_multiplier=long_overlap_multiplier,
            gc_multiplier=gc_multiplier,
        ),
        'coords': edited_coords,
        'manual_subsequence': subsequence,
        'manual_left_trim': left_trim,
        'manual_right_trim': right_trim,
        'left_trim': left_trim,
        'right_trim': right_trim,
    })
    return result

def design_from_candidates(
    context: dict,
    candidates,
    topology: str,
    max_fragment_size: int = DEFAULT_MAX_FRAGMENT_SIZE,
    rank='edited',
    weight_group_std_dev: float = DEFAULT_WEIGHT_GROUP_STD_DEV,
    weight_individual_penalty: float = DEFAULT_WEIGHT_INDIVIDUAL_PENALTY,
    large_fragment_multiplier: float = DEFAULT_LARGE_FRAGMENT_MULTIPLIER,
):
    """Build and score a custom design from candidate dictionaries.

    Unlike :func:`design_from_ids`, this accepts manually trimmed candidate
    copies, making it the preferred entry point for GUI-edited designs.
    """
    candidates = list(candidates)

    if not candidates:
        raise ValueError('At least one overlap must be selected.')

    candidate_ids = [candidate['id'] for candidate in candidates]
    if len(candidate_ids) != len(set(candidate_ids)):
        raise ValueError(
            'Each junction must use a different overlap candidate.'
        )

    stats = evaluate_combination(
        tuple(candidates),
        context['genomes'],
        context['genome_lengths'],
        topology,
        max_fragment_size=max_fragment_size,
        weight_group_std_dev=weight_group_std_dev,
        weight_individual_penalty=weight_individual_penalty,
        large_fragment_multiplier=large_fragment_multiplier,
    )

    return {
        'rank': rank,
        'combo': tuple(candidates),
        **stats,
    }

def design_from_ids(
    context: dict,
    candidate_ids,
    topology: str,
    max_fragment_size: int = DEFAULT_MAX_FRAGMENT_SIZE,
    rank='edited',
    weight_group_std_dev: float = DEFAULT_WEIGHT_GROUP_STD_DEV,
    weight_individual_penalty: float = DEFAULT_WEIGHT_INDIVIDUAL_PENALTY,
    large_fragment_multiplier: float = DEFAULT_LARGE_FRAGMENT_MULTIPLIER,
):
    """
    Build and score a custom assembly design from selected candidate IDs.

    This is the main entry point for GUI editing. Candidate order does not
    affect scoring or extraction; materialization sorts junctions by genomic
    position in the primary genome.
    """
    candidate_ids = list(candidate_ids)

    if not candidate_ids:
        raise ValueError('At least one overlap must be selected.')

    if len(candidate_ids) != len(set(candidate_ids)):
        raise ValueError(
            'Each junction must use a different overlap candidate.'
        )

    missing = [
        candidate_id for candidate_id in candidate_ids
        if candidate_id not in context['candidate_map']
    ]
    if missing:
        raise ValueError(
            'Unknown overlap candidate(s): ' + ', '.join(missing)
        )

    combo = [
        context['candidate_map'][candidate_id]
        for candidate_id in candidate_ids
    ]

    return design_from_candidates(
        context,
        combo,
        topology=topology,
        max_fragment_size=max_fragment_size,
        rank=rank,
        weight_group_std_dev=weight_group_std_dev,
        weight_individual_penalty=weight_individual_penalty,
        large_fragment_multiplier=large_fragment_multiplier,
    )


def candidate_summary_table(context: dict) -> pd.DataFrame:
    """Return a GUI-friendly summary of every optimized candidate."""
    primary = context['primary_genome']
    rows = []

    for candidate in sorted(
        context['candidates'],
        key=lambda item: item['coords'][primary]['center'],
    ):
        rows.append({
            'candidate_id': candidate['id'],
            'cluster': candidate['cluster'],
            'region': candidate['region'],
            'primary_position': candidate['coords'][primary]['center'],
            'optimized_seq': candidate['seq'],
            'optimized_length': candidate['length'],
            'gc_content': candidate['gc'],
            'overlap_penalty': candidate['penalty'],
            'original_length': candidate['original_length'],
        })

    return pd.DataFrame(rows)


def calculate_fragment_sizes(combo, genomes, genome_lengths, topology):
    """
    Calculate the exact physical fragment lengths that will be written to FASTA.

    Fragment boundaries deliberately retain the shared assembly overlap on BOTH
    adjacent fragments. For example, an internal fragment runs from the start
    of its left overlap through the exclusive stop of its right overlap. Thus
    the displayed/scored fragment size includes both terminal overlap sequences
    exactly as they occur on the physical fragment molecule.

    The overlap order is defined once from the primary genome and then reused
    for every homologous genome, matching materialize_design().
    """
    if not genomes:
        return {}

    primary_genome = genomes[0]
    ordered_overlaps = sorted(
        combo,
        key=lambda ov: ov['coords'][primary_genome]['center'],
    )

    all_sizes = {}

    for genome in genomes:
        seq_len = genome_lengths[genome]
        bounds = build_fragment_bounds(
            ordered_overlaps,
            genome,
            seq_len,
            topology,
        )

        all_sizes[genome] = [
            float(interval_length(start, stop, seq_len))
            for start, stop in bounds
        ]

    return all_sizes


def evaluate_combination(
    combo,
    genomes,
    genome_lengths,
    topology,
    max_fragment_size=DEFAULT_MAX_FRAGMENT_SIZE,
    weight_group_std_dev=DEFAULT_WEIGHT_GROUP_STD_DEV,
    weight_individual_penalty=DEFAULT_WEIGHT_INDIVIDUAL_PENALTY,
    large_fragment_multiplier=DEFAULT_LARGE_FRAGMENT_MULTIPLIER,
):
    fragment_sizes = calculate_fragment_sizes(
        combo, genomes, genome_lengths, topology
    )
    std_devs = {
        g: float(np.std(fragment_sizes[g]))
        for g in genomes
    }
    group_std_dev = float(np.mean(list(std_devs.values())))
    average_overlap_penalty = float(
        np.mean([ov['penalty'] for ov in combo])
    )

    large_fragment_penalty = 0.0
    for genome in genomes:
        for size in fragment_sizes[genome]:
            if size > max_fragment_size:
                large_fragment_penalty += (
                    size - max_fragment_size
                ) * large_fragment_multiplier

    score = (
        weight_group_std_dev * group_std_dev
        + weight_individual_penalty * average_overlap_penalty
        + large_fragment_penalty
    )

    return {
        'score': float(score),
        'group_std_dev': group_std_dev,
        'average_overlap_penalty': average_overlap_penalty,
        'large_fragment_penalty': float(large_fragment_penalty),
        'std_devs': std_devs,
        'fragment_sizes': fragment_sizes,
    }


def find_top_designs(
    candidates,
    genomes,
    genome_lengths,
    num_fragments,
    topology,
    top_n=DEFAULT_TOP_N,
    max_fragment_size=DEFAULT_MAX_FRAGMENT_SIZE,
    max_combinations=DEFAULT_MAX_COMBINATIONS,
    weight_group_std_dev=DEFAULT_WEIGHT_GROUP_STD_DEV,
    weight_individual_penalty=DEFAULT_WEIGHT_INDIVIDUAL_PENALTY,
    large_fragment_multiplier=DEFAULT_LARGE_FRAGMENT_MULTIPLIER,
):
    if topology not in {'linear', 'circular'}:
        raise ValueError("topology must be 'linear' or 'circular'")
    if num_fragments < 2:
        raise ValueError('num_fragments must be at least 2')

    overlaps_needed = num_fragments if topology == 'circular' else num_fragments - 1
    if overlaps_needed > len(candidates):
        raise ValueError(
            f'{num_fragments} {topology} fragments require {overlaps_needed} '
            f'overlaps, but only {len(candidates)} candidates are available.'
        )

    total_combinations = math.comb(len(candidates), overlaps_needed)
    if total_combinations > max_combinations:
        raise ValueError(
            f'Design requires evaluating {total_combinations:,} combinations, '
            f'exceeding the configured limit of {max_combinations:,}.'
        )

    heap = []
    counter = 0

    for combo in itertools.combinations(candidates, overlaps_needed):
        stats = evaluate_combination(
            combo,
            genomes,
            genome_lengths,
            topology,
            max_fragment_size=max_fragment_size,
            weight_group_std_dev=weight_group_std_dev,
            weight_individual_penalty=weight_individual_penalty,
            large_fragment_multiplier=large_fragment_multiplier,
        )
        result = {'combo': combo, **stats}
        entry = (-stats['score'], counter, result)
        counter += 1

        if len(heap) < top_n:
            heapq.heappush(heap, entry)
        elif stats['score'] < -heap[0][0]:
            heapq.heapreplace(heap, entry)

    results = [entry[2] for entry in heap]
    results.sort(key=lambda x: x['score'])
    for rank, result in enumerate(results, start=1):
        result['rank'] = rank

    return results, total_combinations


def write_design_outputs(designs, genomes, outdir):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    primary_genome = genomes[0]
    summary_rows = []
    design_files = []

    for design in designs:
        combo = sorted(
            design['combo'],
            key=lambda ov: ov['coords'][primary_genome]['center'],
        )
        design_file = outdir / f"design_rank{design['rank']}_overlaps.tsv"
        design_files.append(design_file)
        rows = []

        for order, ov in enumerate(combo, start=1):
            row = {
                'junction_order': order,
                'cluster_region': ov['id'],
                'cluster': ov['cluster'],
                'region': ov['region'],
                'optimized_seq': ov['seq'],
                'optimized_length': ov['length'],
                'gc_content': ov['gc'],
                'overlap_penalty': ov['penalty'],
                'original_seq': ov['original_seq'],
                'original_length': ov['original_length'],
                'left_trim': ov['left_trim'],
                'right_trim': ov['right_trim'],
            }
            for genome in genomes:
                row[f'{genome}_start'] = ov['coords'][genome]['start']
                row[f'{genome}_stop'] = ov['coords'][genome]['stop']
                row[f'{genome}_center'] = ov['coords'][genome]['center']
            rows.append(row)

        pd.DataFrame(rows).to_csv(design_file, sep='\t', index=False)

        summary = {
            'rank': design['rank'],
            'score': design['score'],
            'group_fragment_size_std_dev': design['group_std_dev'],
            'average_overlap_penalty': design['average_overlap_penalty'],
            'large_fragment_penalty': design['large_fragment_penalty'],
            'junctions': ';'.join(ov['id'] for ov in combo),
        }
        for genome in genomes:
            summary[f'{genome}_fragment_sizes'] = ';'.join(
                f'{x:.1f}' for x in design['fragment_sizes'][genome]
            )
            summary[f'{genome}_fragment_size_std_dev'] = design['std_devs'][genome]
        summary_rows.append(summary)

    summary_file = outdir / 'design_summary.tsv'
    pd.DataFrame(summary_rows).to_csv(summary_file, sep='\t', index=False)
    return summary_file, design_files


def write_single_design(design, genomes, output_path):
    """Write one ranked or edited design as an overlap TSV."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    primary_genome = genomes[0]

    combo = sorted(
        design['combo'],
        key=lambda ov: ov['coords'][primary_genome]['center'],
    )

    rows = []
    for order, ov in enumerate(combo, start=1):
        row = {
            'junction_order': order,
            'cluster_region': ov['id'],
            'cluster': ov['cluster'],
            'region': ov['region'],
            'optimized_seq': ov['seq'],
            'optimized_length': ov['length'],
            'gc_content': ov['gc'],
            'overlap_penalty': ov['penalty'],
            'original_seq': ov['original_seq'],
            'original_length': ov['original_length'],
            'left_trim': ov['left_trim'],
            'right_trim': ov['right_trim'],
        }
        for genome in genomes:
            row[f'{genome}_start'] = ov['coords'][genome]['start']
            row[f'{genome}_stop'] = ov['coords'][genome]['stop']
            row[f'{genome}_center'] = ov['coords'][genome]['center']
        rows.append(row)

    pd.DataFrame(rows).to_csv(output_path, sep='\t', index=False)
    return output_path


def extract_sequence(sequence: str, start_1b: int, stop_1b_exclusive: int) -> str:
    start_idx = int(start_1b) - 1
    stop_idx = int(stop_1b_exclusive) - 1
    if stop_idx < start_idx:
        return sequence[start_idx:] + sequence[:stop_idx]
    return sequence[start_idx:stop_idx]


def wrap_coord(coord: int, seq_len: int) -> int:
    return ((int(coord) - 1) % seq_len) + 1


def build_fragment_bounds(design_overlaps, genome, seq_len, topology):
    num_overlaps = len(design_overlaps)
    num_fragments = num_overlaps if topology == 'circular' else num_overlaps + 1
    bounds = []

    if topology == 'circular':
        bounds.append((
            design_overlaps[-1]['coords'][genome]['start'],
            design_overlaps[0]['coords'][genome]['stop'],
        ))
        for i in range(1, num_fragments):
            bounds.append((
                design_overlaps[i - 1]['coords'][genome]['start'],
                design_overlaps[i]['coords'][genome]['stop'],
            ))
    else:
        bounds.append((1, design_overlaps[0]['coords'][genome]['stop']))
        for i in range(1, num_fragments - 1):
            bounds.append((
                design_overlaps[i - 1]['coords'][genome]['start'],
                design_overlaps[i]['coords'][genome]['stop'],
            ))
        bounds.append((
            design_overlaps[-1]['coords'][genome]['start'],
            seq_len + 1,
        ))

    return bounds


def write_fasta_record(handle, header, sequence):
    handle.write(f'>{header}\n')
    for i in range(0, len(sequence), 80):
        handle.write(sequence[i:i + 80] + '\n')


def materialize_design(
    design,
    genome_sequences,
    genomes,
    topology,
    outdir,
    chimera_prefix='Chimera',
    generate_chimeras=True,
    max_chimeras=DEFAULT_MAX_CHIMERAS,
    clean_outdir=False,
):
    outdir = Path(outdir)

    if clean_outdir and outdir.exists():
        shutil.rmtree(outdir)

    outdir.mkdir(parents=True, exist_ok=True)
    fragments_dir = outdir / 'fragments'
    fragments_dir.mkdir(parents=True, exist_ok=True)
    chimeras_dir = outdir / 'chimeras'

    primary_genome = genomes[0]
    overlaps = sorted(
        design['combo'],
        key=lambda ov: ov['coords'][primary_genome]['center'],
    )
    num_fragments = len(overlaps) if topology == 'circular' else len(overlaps) + 1
    bounds_by_genome = {}
    fragment_rows = []

    handles = {
        i: (fragments_dir / f'Fragment_{i}.fasta').open('w')
        for i in range(1, num_fragments + 1)
    }

    try:
        for genome in genomes:
            seq = genome_sequences[genome]
            bounds = build_fragment_bounds(overlaps, genome, len(seq), topology)
            bounds_by_genome[genome] = bounds

            for fragment_number, (start, stop) in enumerate(bounds, start=1):
                frag_seq = extract_sequence(seq, start, stop)
                header = f'{genome}_Fragment_{fragment_number}_{start}-{stop - 1}'
                write_fasta_record(handles[fragment_number], header, frag_seq)
                fragment_rows.append({
                    'fragment_number': fragment_number,
                    'source_genome': genome,
                    'parent_start': start,
                    'parent_stop': stop - 1,
                    'fragment_length': len(frag_seq),
                    'fasta_record': header,
                })
    finally:
        for handle in handles.values():
            handle.close()

    fragment_manifest = outdir / 'fragment_manifest.tsv'
    pd.DataFrame(fragment_rows).to_csv(fragment_manifest, sep='\t', index=False)

    chimera_count = 0
    chimera_log = None

    if generate_chimeras:
        total_chimeras = len(genomes) ** num_fragments - len(genomes)
        if total_chimeras > max_chimeras:
            raise ValueError(
                f'This design would generate {total_chimeras:,} non-parental '
                f'chimeras, exceeding max_chimeras={max_chimeras:,}.'
            )

        chimeras_dir.mkdir(parents=True, exist_ok=True)
        log_rows = []

        for chimera_index, combo in enumerate(
            (c for c in itertools.product(genomes, repeat=num_fragments) if len(set(c)) > 1),
            start=1,
        ):
            chimera_count += 1
            chimera_name = f'{chimera_prefix}_{chimera_index:04d}'
            chimera_seq = ''
            current_len = 0

            for frag_index, genome in enumerate(combo):
                seq = genome_sequences[genome]
                seq_len = len(seq)
                start, stop = bounds_by_genome[genome][frag_index]
                frag_seq = extract_sequence(seq, start, stop)

                left_ov_len = 0
                if topology == 'circular':
                    left_ov = overlaps[frag_index - 1]
                    left_ov_len = interval_length(
                        left_ov['coords'][genome]['start'],
                        left_ov['coords'][genome]['stop'],
                        seq_len,
                    )
                elif frag_index > 0:
                    left_ov = overlaps[frag_index - 1]
                    left_ov_len = interval_length(
                        left_ov['coords'][genome]['start'],
                        left_ov['coords'][genome]['stop'],
                        seq_len,
                    )

                trimmed = frag_seq[left_ov_len:]
                actual_parent_start = wrap_coord(start + left_ov_len, seq_len)
                chimera_start = current_len + 1
                chimera_stop = current_len + len(trimmed)
                chimera_seq += trimmed
                current_len = chimera_stop

                if topology == 'linear' and frag_index == num_fragments - 1:
                    right_overlap = None
                else:
                    right_overlap = overlaps[frag_index % len(overlaps)]['id']

                log_rows.append({
                    'chimera_id': chimera_name,
                    'fragment_number': frag_index + 1,
                    'source_genome': genome,
                    'actual_parent_start': actual_parent_start,
                    'actual_parent_stop': stop - 1,
                    'chimera_start': chimera_start,
                    'chimera_stop': chimera_stop,
                    'right_overlap': right_overlap,
                })

            with (chimeras_dir / f'{chimera_name}.fasta').open('w') as handle:
                write_fasta_record(
                    handle,
                    f"{chimera_name} fragments={'-'.join(combo)}",
                    chimera_seq,
                )

        chimera_log = chimeras_dir / 'chimera_assembly_log.tsv'
        pd.DataFrame(log_rows).to_csv(chimera_log, sep='\t', index=False)

    metadata_file = outdir / 'assembly_design.json'
    metadata_file.write_text(json.dumps({
        'rank': design['rank'],
        'score': design['score'],
        'topology': topology,
        'num_fragments': num_fragments,
        'genomes': genomes,
        'junctions': [ov['id'] for ov in overlaps],
        'chimera_count': chimera_count,
    }, indent=2))

    return {
        'fragments_dir': fragments_dir,
        'fragment_manifest': fragment_manifest,
        'chimeras_dir': chimeras_dir if generate_chimeras else None,
        'chimera_log': chimera_log,
        'metadata_file': metadata_file,
        'chimera_count': chimera_count,
    }


def design_assemblies(
    overlaps: str | Path,
    fasta: str | Path,
    outdir: str | Path,
    num_fragments: int,
    topology: str,
    top_n: int = DEFAULT_TOP_N,
    selected_rank: int | None = None,
    ideal_overlap_length: int = DEFAULT_IDEAL_OVERLAP_LENGTH,
    ideal_gc: float = DEFAULT_IDEAL_GC,
    max_fragment_size: int = DEFAULT_MAX_FRAGMENT_SIZE,
    max_combinations: int = DEFAULT_MAX_COMBINATIONS,
    weight_group_std_dev: float = DEFAULT_WEIGHT_GROUP_STD_DEV,
    weight_individual_penalty: float = DEFAULT_WEIGHT_INDIVIDUAL_PENALTY,
    short_overlap_multiplier: float = DEFAULT_SHORT_OVERLAP_MULTIPLIER,
    long_overlap_multiplier: float = DEFAULT_LONG_OVERLAP_MULTIPLIER,
    gc_multiplier: float = DEFAULT_GC_MULTIPLIER,
    large_fragment_multiplier: float = DEFAULT_LARGE_FRAGMENT_MULTIPLIER,
    generate_chimeras: bool = True,
    max_chimeras: int = DEFAULT_MAX_CHIMERAS,
    chimera_prefix: str = 'Chimera',
    analyze: bool = False,
    analysis_prepped_dir: str | Path | None = None,
    analysis_cluster_membership: str | Path | None = None,
    save_analysis_alignments: bool = False,
    top_variable_cds: int = 5,
    max_pairwise_assemblies: int = 250,
):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    topology = topology.lower()

    genome_sequences = read_genomes_fasta(fasta)
    candidates, genomes = load_overlap_candidates(
        overlaps,
        genome_sequences,
        ideal_overlap_length,
        ideal_gc,
        short_overlap_multiplier,
        long_overlap_multiplier,
        gc_multiplier,
    )
    genome_lengths = {g: len(genome_sequences[g]) for g in genomes}

    designs, total_combinations = find_top_designs(
        candidates,
        genomes,
        genome_lengths,
        num_fragments,
        topology,
        top_n=top_n,
        max_fragment_size=max_fragment_size,
        max_combinations=max_combinations,
        weight_group_std_dev=weight_group_std_dev,
        weight_individual_penalty=weight_individual_penalty,
        large_fragment_multiplier=large_fragment_multiplier,
    )

    summary_file, design_files = write_design_outputs(
        designs,
        genomes,
        outdir,
    )

    materialized = None
    analysis = None
    if selected_rank is not None:
        if not 1 <= selected_rank <= len(designs):
            raise ValueError(
                f'selected_rank must be between 1 and {len(designs)}'
            )
        materialized = materialize_design(
            designs[selected_rank - 1],
            genome_sequences,
            genomes,
            topology,
            outdir / f'design_rank{selected_rank}_assembly',
            chimera_prefix=chimera_prefix,
            generate_chimeras=generate_chimeras,
            max_chimeras=max_chimeras,
        )

        if analyze:
            if analysis_prepped_dir is None:
                raise ValueError(
                    'analysis_prepped_dir is required when analyze=True.'
                )
            if analysis_cluster_membership is None:
                raise ValueError(
                    'analysis_cluster_membership is required when analyze=True.'
                )

            # Imported lazily so design ranking/materialization remains usable
            # without loading matplotlib or the analysis module.
            from analyze_design import analyze_materialized_design

            assembly_dir = Path(materialized['metadata_file']).parent
            analysis = analyze_materialized_design(
                materialized_dir=assembly_dir,
                prepped_dir=analysis_prepped_dir,
                cluster_membership=analysis_cluster_membership,
                outdir=assembly_dir / 'analysis',
                save_alignments=save_analysis_alignments,
                top_variable_cds=top_variable_cds,
                max_pairwise_assemblies=max_pairwise_assemblies,
            )
    elif analyze:
        raise ValueError(
            'Analysis requires a materialized design. Supply selected_rank.'
        )

    return {
        'genomes': genomes,
        'candidate_count': len(candidates),
        'total_combinations': total_combinations,
        'designs': designs,
        'summary_file': summary_file,
        'design_files': design_files,
        'materialized': materialized,
        'analysis': analysis,
    }


def build_parser():
    parser = argparse.ArgumentParser(
        description=(
            'Rank CombPhage overlap sets, extract corresponding fragments, '
            'and optionally generate all non-parental chimeric assemblies.'
        )
    )
    parser.add_argument('--overlaps', required=True)
    parser.add_argument('--fasta', required=True)
    parser.add_argument('--outdir', required=True)
    parser.add_argument('--num-fragments', required=True, type=int)
    parser.add_argument('--topology', required=True, choices=['linear', 'circular'])
    parser.add_argument('--top-n', type=int, default=DEFAULT_TOP_N)
    parser.add_argument('--selected-rank', type=int, default=None)
    parser.add_argument('--ideal-overlap-length', type=int, default=DEFAULT_IDEAL_OVERLAP_LENGTH)
    parser.add_argument('--ideal-gc', type=float, default=DEFAULT_IDEAL_GC)
    parser.add_argument('--max-fragment-size', type=int, default=DEFAULT_MAX_FRAGMENT_SIZE)
    parser.add_argument('--max-combinations', type=int, default=DEFAULT_MAX_COMBINATIONS)
    parser.add_argument('--fragment-balance-weight', type=float, default=DEFAULT_WEIGHT_GROUP_STD_DEV)
    parser.add_argument('--overlap-quality-weight', type=float, default=DEFAULT_WEIGHT_INDIVIDUAL_PENALTY)
    parser.add_argument('--short-overlap-multiplier', type=float, default=DEFAULT_SHORT_OVERLAP_MULTIPLIER)
    parser.add_argument('--long-overlap-multiplier', type=float, default=DEFAULT_LONG_OVERLAP_MULTIPLIER)
    parser.add_argument('--gc-multiplier', type=float, default=DEFAULT_GC_MULTIPLIER)
    parser.add_argument('--large-fragment-multiplier', type=float, default=DEFAULT_LARGE_FRAGMENT_MULTIPLIER)
    parser.add_argument('--max-chimeras', type=int, default=DEFAULT_MAX_CHIMERAS)
    parser.add_argument('--chimera-prefix', default='Chimera')
    parser.add_argument('--fragments-only', action='store_true')
    parser.add_argument(
        '--analyze',
        action='store_true',
        help=(
            'After materializing --selected-rank, analyze fragment/CDS/assembly '
            'sequence diversity using analyze_design.py.'
        ),
    )
    parser.add_argument(
        '--prepped-dir',
        default=None,
        help=(
            'Prepared genome directory containing genomes.fasta and genbanks/. '
            'Required with --analyze.'
        ),
    )
    parser.add_argument(
        '--cluster-membership',
        default=None,
        help=(
            'Path to cluster_membership.tsv from find_overlaps.py. '
            'Required with --analyze.'
        ),
    )
    parser.add_argument(
        '--save-analysis-alignments',
        action='store_true',
        help='Save pairwise fragment alignments during design analysis.',
    )
    parser.add_argument(
        '--top-variable-cds',
        type=int,
        default=5,
        help='Number of most variable CDS clusters to plot per fragment.',
    )
    parser.add_argument(
        '--max-pairwise-assemblies',
        type=int,
        default=250,
        help=(
            'Maximum total assemblies for a full all-vs-all assembly identity '
            'matrix during analysis.'
        ),
    )
    return parser


def main():
    args = build_parser().parse_args()
    result = design_assemblies(
        overlaps=args.overlaps,
        fasta=args.fasta,
        outdir=args.outdir,
        num_fragments=args.num_fragments,
        topology=args.topology,
        top_n=args.top_n,
        selected_rank=args.selected_rank,
        ideal_overlap_length=args.ideal_overlap_length,
        ideal_gc=args.ideal_gc,
        max_fragment_size=args.max_fragment_size,
        max_combinations=args.max_combinations,
        weight_group_std_dev=args.fragment_balance_weight,
        weight_individual_penalty=args.overlap_quality_weight,
        short_overlap_multiplier=args.short_overlap_multiplier,
        long_overlap_multiplier=args.long_overlap_multiplier,
        gc_multiplier=args.gc_multiplier,
        large_fragment_multiplier=args.large_fragment_multiplier,
        generate_chimeras=not args.fragments_only,
        max_chimeras=args.max_chimeras,
        chimera_prefix=args.chimera_prefix,
        analyze=args.analyze,
        analysis_prepped_dir=args.prepped_dir,
        analysis_cluster_membership=args.cluster_membership,
        save_analysis_alignments=args.save_analysis_alignments,
        top_variable_cds=args.top_variable_cds,
        max_pairwise_assemblies=args.max_pairwise_assemblies,
    )

    print(
        f"Evaluated {result['total_combinations']:,} overlap-set combinations "
        f"from {result['candidate_count']} candidates."
    )
    print(f"Saved design summary: {result['summary_file']}")
    for path in result['design_files']:
        print(f'Saved: {path}')

    if result['materialized'] is not None:
        print(f"Fragments: {result['materialized']['fragments_dir']}")
        if result['materialized']['chimeras_dir'] is not None:
            print(
                f"Generated {result['materialized']['chimera_count']:,} "
                f"non-parental chimeras in "
                f"{result['materialized']['chimeras_dir']}"
            )

    if result.get('analysis') is not None:
        print(f"Analysis: {result['analysis']['outdir']}")
        print(json.dumps(result['analysis']['overview'], indent=2))


if __name__ == '__main__':
    main()
