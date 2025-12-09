#!/usr/bin/env python3
"""
Script to validate evidence files for OGRDB upload.
Checks sequence evidence CSV files against genomic sequences from GenBank or local FASTA files.
"""

import argparse
import csv
import os
import sys
from typing import Dict, List, Tuple

from receptor_utils import simple_bio_seq as simple
from Bio import Entrez, SeqIO


def basic_checks(evidence_file) -> Tuple[bool, List[str]]:
    """
    Perform basic checks:
     - Check that expected columns are present in the evidence file
     - Check that no cells contain leading or trailing spaces
    
    Returns:
        Tuple of (success, list_of_errors)
    """
    errors = []

    expected_columns = [
        'gene_label', 'sequence', 'sequence_type', 'repository', 'accession', 'patch',
        'start', 'end', 'sense', 'species_subgroup', 'subgroup_type', 'notes', 'gene_start', 'gene_end',
        'utr_5_prime_start', 'utr_5_prime_end',
        'leader_1_start', 'leader_1_end',
        'leader_2_start', 'leader_2_end',
        'v_rs_start', 'v_rs_end',
        'd_rs_3_prime_start', 'd_rs_3_prime_end',
        'd_rs_5_prime_start', 'd_rs_5_prime_end',
        'j_rs_start', 'j_rs_end',
        'j_codon_frame', 'j_cdr3_end',
        "c_exon_1_start", "c_exon_1_end", "c_exon_2_start", "c_exon_2_end", "c_exon_3_start",
        "c_exon_3_end", "c_exon_4_start", "c_exon_4_end", "c_exon_5_start", "c_exon_5_end",
        "c_exon_6_start", "c_exon_6_end", "c_exon_7_start", "c_exon_7_end", "c_exon_8_start", "c_exon_8_end",
        "c_exon_9_start", "c_exon_9_end", "utr_3_prime_start", "utr_3_prime_end", "c_tm_sequence", "c_sc_sequence"
    ]

    # Check if file exists
    if not os.path.exists(evidence_file):
        errors.append(f"Evidence file not found: {evidence_file}")
        return False, errors

    with open(evidence_file, 'r') as f:
        reader = csv.DictReader(f)
        # Check for expected columns
        missing_columns = [col for col in expected_columns if col not in reader.fieldnames]
        if missing_columns:
            errors.append(f"Missing columns in evidence file: {', '.join(missing_columns)}")

        # Check for leading/trailing spaces
        for row in reader:
            for col in expected_columns:
                if col in row:
                    if row[col].strip() != row[col]:
                        errors.append(f"Leading/trailing spaces found in row {row['gene_label']} column '{col}'")

    return len(errors) == 0, errors
    

def get_genbank_sequence(accession: str, email: str, save_to_file: bool = False) -> str:
    """
    Fetch a sequence from GenBank given the accession number.
    
    Args:
        accession: GenBank accession ID (with or without version)
        email: Email address for Entrez queries
        save_to_file: If True, save the sequence to a FASTA file
    
    Returns:
        The sequence string
    """
    print(f'Fetching from GenBank: {accession}')
    Entrez.email = email
    
    try:
        handle = Entrez.efetch(db='nucleotide', id=accession, rettype='fasta', retmode='text')
        seq_record = SeqIO.read(handle, "fasta")
        handle.close()
        
        sequence = str(seq_record.seq).upper()
        
        if save_to_file:
            filename = f"{accession}.fasta"
            simple.write_fasta(filename, {accession: sequence})
            print(f"Saved sequence to {filename}")
        
        return sequence
        
    except Exception as e:
        raise ValueError(f"Error fetching sequence {accession} from GenBank: {e}")


def load_sequence_from_file(accession: str) -> str:
    """
    Load sequence from local FASTA file.
    
    Args:
        accession: Accession ID used as filename (with .fasta extension)
    
    Returns:
        The sequence string
    """
    filename = f"{accession}.fasta"
    
    if not os.path.exists(filename):
        raise FileNotFoundError(f"FASTA file not found: {filename}")
    
    try:
        sequences = simple.read_fasta(filename)
        if len(sequences) == 0:
            raise ValueError(f"No sequences found in {filename}")
        
        # Return the first sequence (assuming single sequence per file)
        return list(sequences.values())[0].upper()
        
    except Exception as e:
        raise ValueError(f"Error reading sequence from {filename}: {e}")


def validate_required_fields(row: Dict, row_num: int) -> List[str]:
    """
    Validate that required fields are present and non-empty.
    
    Args:
        row: CSV row as dictionary
        row_num: Row number for error reporting
    
    Returns:
        List of error messages
    """
    errors = []
    required_fields = ['gene_label', 'sequence', 'sequence_type', 'repository',
                       'accession', 'start', 'end', 'sense', 'gene_start', 'gene_end']
    
    for field in required_fields:
        value = row[field]
        if not value:
            errors.append(f"Row {row_num}: Required field '{field}' is missing or empty")
    
    return errors


def validate_enumerated_fields(row: Dict, row_num: int) -> List[str]:
    """
    Validate fields with enumerated values.
    
    Args:
        row: CSV row as dictionary
        row_num: Row number for error reporting
    
    Returns:
        List of error messages
    """
    errors = []
    
    # Validate sequence_type
    valid_sequence_types = ['V', 'D', 'J', 'C']
    sequence_type = row['sequence_type']
    if sequence_type and sequence_type not in valid_sequence_types:
        errors.append(f"Row {row_num}: Invalid sequence_type '{sequence_type}'. Must be one of: {valid_sequence_types}")
    
    # Validate sense
    valid_senses = ['+', '-']
    sense = row['sense']
    if sense and sense not in valid_senses:
        errors.append(f"Row {row_num}: Invalid sense '{sense}'. Must be one of: {valid_senses}")
    
    return errors


def validate_coordinates(row: Dict, row_num: int) -> List[str]:
    """
    Validate coordinate fields are positive integers with start < end.
    
    Args:
        row: CSV row as dictionary
        row_num: Row number for error reporting
    
    Returns:
        List of error messages
    """
    errors = []
    
    # Coordinate pairs to validate
    coordinate_pairs = [
        ('start', 'end'),
        ('gene_start', 'gene_end'),
        ('utr_5_prime_start', 'utr_5_prime_end'),
        ('leader_1_start', 'leader_1_end'),
        ('leader_2_start', 'leader_2_end'),
        ('v_rs_start', 'v_rs_end'),
        ('d_rs_3_prime_start', 'd_rs_3_prime_end'),
        ('d_rs_5_prime_start', 'd_rs_5_prime_end'),
        ('j_rs_start', 'j_rs_end'),
        ('c_exon_1_start', 'c_exon_1_end'),
        ('c_exon_2_start', 'c_exon_2_end'),
        ('c_exon_3_start', 'c_exon_3_end'),
        ('c_exon_4_start', 'c_exon_4_end'),
        ('c_exon_5_start', 'c_exon_5_end'),
        ('c_exon_6_start', 'c_exon_6_end'),
        ('c_exon_7_start', 'c_exon_7_end'),
        ('c_exon_8_start', 'c_exon_8_end'),
        ('c_exon_9_start', 'c_exon_9_end'),
        ('utr_3_prime_start', 'utr_3_prime_end')
    ]
    
    for start_field, end_field in coordinate_pairs:
        start_val = row[start_field]
        end_val = row[end_field]
        
        # Skip if both are empty
        if not start_val and not end_val:
            continue
            
        # If one is provided, both must be provided
        if bool(start_val) != bool(end_val):
            errors.append(f"Row {row_num}: Both {start_field} and {end_field} must be provided together")
            continue
        
        # Validate they are positive integers
        try:
            start_int = int(start_val)
            end_int = int(end_val)
            
            if start_int < 1:
                errors.append(f"Row {row_num}: {start_field} must be >= 1 (1-based coordinates)")
            if end_int < 1:
                errors.append(f"Row {row_num}: {end_field} must be >= 1 (1-based coordinates)")
            if start_int >= end_int:
                errors.append(f"Row {row_num}: {start_field} ({start_int}) must be < {end_field} ({end_int})")
                
        except ValueError:
            errors.append(f"Row {row_num}: {start_field} and {end_field} must be integers")
    
    # Special validation for single coordinate fields
    single_coordinates = ['j_cdr3_end']
    for field in single_coordinates:
        val = row[field]
        if val:
            try:
                coord_int = int(val)
                if coord_int < 1:
                    errors.append(f"Row {row_num}: {field} must be >= 1 (1-based coordinates)")
            except ValueError:
                errors.append(f"Row {row_num}: {field} must be an integer")
    
    return errors


def validate_adjacency(row: Dict, row_num: int) -> List[str]:
    """
    Validate that features are adjacent according to OGRDB requirements.
    
    For + sense:
      V genes: utr_5 -> leader_1 (adjacent), leader_1 -> leader_2 (gap allowed),
               leader_1/2 -> gene -> v_rs (all adjacent)
      D genes: d_rs_5 -> gene -> d_rs_3 (all adjacent)
      J genes: j_rs -> gene (adjacent)
      C genes: c_exon_1 -> c_exon_2 -> ... -> c_exon_9 (gap allowed, null coordinates allowed)
    
    For - sense: functional order is reversed but genomic coordinates still start < end
      V genes: v_rs -> gene -> leader_2/1 -> leader_1 -> utr_5 (all adjacent except gap allowed)
      D genes: d_rs_3 -> gene -> d_rs_5 (all adjacent)
      J genes: gene -> j_rs (adjacent)
      C genes: c_exon_9 -> ... -> c_exon_2 -> c_exon_1 (gap allowed, null coordinates allowed)
    
    Also validates that start/end coordinates span the entire annotated sequence.
    
    Args:
        row: CSV row as dictionary
        row_num: Row number for error reporting
    
    Returns:
        List of error messages
    """
    errors = []
    sequence_type = row['sequence_type']
    sense = row['sense']
    
    if not sequence_type:
        return errors  # Skip if sequence_type is missing (caught by other validation)
    
    is_minus_sense = sense in ['-', 'minus']
    
    # Helper function to get coordinate value or None
    def get_coord(field_name):
        val = row[field_name]
        if not val or val == 'None':
            return None
        try:
            return int(val)
        except ValueError:
            return None
    
    # Get all coordinates
    start = get_coord('start')
    end = get_coord('end')
    gene_start = get_coord('gene_start')
    gene_end = get_coord('gene_end')
    
    if sequence_type == 'V':
        # V gene adjacency validation
        utr_5_start = get_coord('utr_5_prime_start')
        utr_5_end = get_coord('utr_5_prime_end')
        leader_1_start = get_coord('leader_1_start')
        leader_1_end = get_coord('leader_1_end')
        leader_2_start = get_coord('leader_2_start')
        leader_2_end = get_coord('leader_2_end')
        v_rs_start = get_coord('v_rs_start')
        v_rs_end = get_coord('v_rs_end')
        
        # Track expected features in genomic coordinate order
        features = []
        
        if not is_minus_sense:
            # Plus sense: utr_5 -> leader_1 -> leader_2 -> gene -> v_rs
            
            # UTR 5' prime
            if utr_5_start is not None and utr_5_end is not None:
                features.append(('utr_5_prime', utr_5_start, utr_5_end))
            
            # Leader 1
            if leader_1_start is not None and leader_1_end is not None:
                features.append(('leader_1', leader_1_start, leader_1_end))
                
                # Check adjacency with UTR 5'
                if utr_5_end is not None and leader_1_start != utr_5_end + 1:
                    errors.append(f"Row {row_num}: leader_1_start ({leader_1_start}) should be adjacent to utr_5_prime_end + 1 ({utr_5_end + 1})")
            
            # Leader 2 (gap allowed between leader_1 and leader_2)
            if leader_2_start is not None and leader_2_end is not None:
                features.append(('leader_2', leader_2_start, leader_2_end))
            
            # Gene (must be adjacent to last leader)
            if gene_start is not None and gene_end is not None:
                features.append(('gene', gene_start, gene_end))
                
                # Check adjacency with leader
                last_leader_end = None
                if leader_2_end is not None:
                    last_leader_end = leader_2_end
                elif leader_1_end is not None:
                    last_leader_end = leader_1_end
                
                if last_leader_end is not None and gene_start != last_leader_end + 1:
                    errors.append(f"Row {row_num}: gene_start ({gene_start}) should be adjacent to last leader end + 1 ({last_leader_end + 1})")
            
            # V recombination signal (must be adjacent to gene)
            if v_rs_start is not None and v_rs_end is not None:
                features.append(('v_rs', v_rs_start, v_rs_end))
                
                if gene_end is not None and v_rs_start != gene_end + 1:
                    errors.append(f"Row {row_num}: v_rs_start ({v_rs_start}) should be adjacent to gene_end + 1 ({gene_end + 1})")
        
        else:
            # Minus sense: v_rs -> gene -> leader_2 -> leader_1 -> utr_5 (in genomic coordinates)
            
            # V recombination signal comes first in genomic coordinates for minus sense
            if v_rs_start is not None and v_rs_end is not None:
                features.append(('v_rs', v_rs_start, v_rs_end))
            
            # Gene (must be adjacent to v_rs)
            if gene_start is not None and gene_end is not None:
                features.append(('gene', gene_start, gene_end))
                
                if v_rs_end is not None and gene_start != v_rs_end + 1:
                    errors.append(f"Row {row_num}: gene_start ({gene_start}) should be adjacent to v_rs_end + 1 ({v_rs_end + 1})")
            
            # Leader 2 (gap allowed between gene and leader_2)
            if leader_2_start is not None and leader_2_end is not None:
                features.append(('leader_2', leader_2_start, leader_2_end))
            
            # Leader 1 (must be adjacent to leader_2, or to gene if no leader_2)
            if leader_1_start is not None and leader_1_end is not None:
                features.append(('leader_1', leader_1_start, leader_1_end))
                
                # Check adjacency with previous feature
                prev_end = None
                if leader_2_end is not None:
                    prev_end = leader_2_end
                elif gene_end is not None:
                    prev_end = gene_end
                
                if prev_end is not None and leader_1_start != prev_end + 1:
                    prev_feature = 'leader_2' if leader_2_end is not None else 'gene'
                    errors.append(f"Row {row_num}: leader_1_start ({leader_1_start}) should be adjacent to {prev_feature}_end + 1 ({prev_end + 1})")
            
            # UTR 5' prime (must be adjacent to leader_1)
            if utr_5_start is not None and utr_5_end is not None:
                features.append(('utr_5_prime', utr_5_start, utr_5_end))
                
                if leader_1_end is not None and utr_5_start != leader_1_end + 1:
                    errors.append(f"Row {row_num}: utr_5_prime_start ({utr_5_start}) should be adjacent to leader_1_end + 1 ({leader_1_end + 1})")
    
    elif sequence_type == 'D':
        # D gene adjacency validation
        d_rs_5_start = get_coord('d_rs_5_prime_start')
        d_rs_5_end = get_coord('d_rs_5_prime_end')
        d_rs_3_start = get_coord('d_rs_3_prime_start')
        d_rs_3_end = get_coord('d_rs_3_prime_end')
        
        features = []
        
        if not is_minus_sense:
            # Plus sense: d_rs_5 -> gene -> d_rs_3
            
            # 5' recombination signal
            if d_rs_5_start is not None and d_rs_5_end is not None:
                features.append(('d_rs_5_prime', d_rs_5_start, d_rs_5_end))
                
                # Check adjacency with gene
                if gene_start is not None and gene_start != d_rs_5_end + 1:
                    errors.append(f"Row {row_num}: gene_start ({gene_start}) should be adjacent to d_rs_5_prime_end + 1 ({d_rs_5_end + 1})")
            
            # Gene
            if gene_start is not None and gene_end is not None:
                features.append(('gene', gene_start, gene_end))
            
            # 3' recombination signal
            if d_rs_3_start is not None and d_rs_3_end is not None:
                features.append(('d_rs_3_prime', d_rs_3_start, d_rs_3_end))
                
                # Check adjacency with gene
                if gene_end is not None and d_rs_3_start != gene_end + 1:
                    errors.append(f"Row {row_num}: d_rs_3_prime_start ({d_rs_3_start}) should be adjacent to gene_end + 1 ({gene_end + 1})")
        
        else:
            # Minus sense: d_rs_3 -> gene -> d_rs_5 (in genomic coordinates)
            
            # 3' recombination signal comes first in genomic coordinates for minus sense
            if d_rs_3_start is not None and d_rs_3_end is not None:
                features.append(('d_rs_3_prime', d_rs_3_start, d_rs_3_end))
            
            # Gene (must be adjacent to d_rs_3)
            if gene_start is not None and gene_end is not None:
                features.append(('gene', gene_start, gene_end))
                
                if d_rs_3_end is not None and gene_start != d_rs_3_end + 1:
                    errors.append(f"Row {row_num}: gene_start ({gene_start}) should be adjacent to d_rs_3_prime_end + 1 ({d_rs_3_end + 1})")
            
            # 5' recombination signal (must be adjacent to gene)
            if d_rs_5_start is not None and d_rs_5_end is not None:
                features.append(('d_rs_5_prime', d_rs_5_start, d_rs_5_end))
                
                if gene_end is not None and d_rs_5_start != gene_end + 1:
                    errors.append(f"Row {row_num}: d_rs_5_prime_start ({d_rs_5_start}) should be adjacent to gene_end + 1 ({gene_end + 1})")
    
    elif sequence_type == 'J':
        # J gene adjacency validation
        j_rs_start = get_coord('j_rs_start')
        j_rs_end = get_coord('j_rs_end')
        
        features = []
        
        if not is_minus_sense:
            # Plus sense: j_rs -> gene
            
            # J recombination signal
            if j_rs_start is not None and j_rs_end is not None:
                features.append(('j_rs', j_rs_start, j_rs_end))
                
                # Check adjacency with gene
                if gene_start is not None and gene_start != j_rs_end + 1:
                    errors.append(f"Row {row_num}: gene_start ({gene_start}) should be adjacent to j_rs_end + 1 ({j_rs_end + 1})")
            
            # Gene
            if gene_start is not None and gene_end is not None:
                features.append(('gene', gene_start, gene_end))
        
        else:
            # Minus sense: gene -> j_rs (in genomic coordinates)
            
            # Gene comes first in genomic coordinates for minus sense
            if gene_start is not None and gene_end is not None:
                features.append(('gene', gene_start, gene_end))
            
            # J recombination signal (must be adjacent to gene)
            if j_rs_start is not None and j_rs_end is not None:
                features.append(('j_rs', j_rs_start, j_rs_end))
                
                if gene_end is not None and j_rs_start != gene_end + 1:
                    errors.append(f"Row {row_num}: j_rs_start ({j_rs_start}) should be adjacent to gene_end + 1 ({gene_end + 1})")

    elif sequence_type == 'C':
        # C gene exon adjacency validation
        c_exon_errors = check_c_exons(row, sense)
        for field, message in c_exon_errors:
            errors.append(f"Row {row_num}: {message}")
        
        # Collect features for span validation
        features = []
        for i in range(1, 10):
            start_field = f'c_exon_{i}_start'
            end_field = f'c_exon_{i}_end'
            start_coord = get_coord(start_field)
            end_coord = get_coord(end_field)
            if start_coord is not None and end_coord is not None:
                features.append((f'c_exon_{i}', start_coord, end_coord))
    
    # Validate that start/end span the entire annotated sequence
    if features and start is not None and end is not None:
        expected_start = min(feat[1] for feat in features)  # Minimum start coordinate
        expected_end = max(feat[2] for feat in features)    # Maximum end coordinate
        
        if start != expected_start:
            feature_names = [feat[0] for feat in features]
            errors.append(f"Row {row_num}: start coordinate ({start}) should be {expected_start} (start of first feature: {feature_names[0]})")
        
        if end != expected_end:
            feature_names = [feat[0] for feat in features]
            errors.append(f"Row {row_num}: end coordinate ({end}) should be {expected_end} (end of last feature: {feature_names[-1]})")
    
    return errors


def validate_sequence_match(row: Dict, row_num: int, genomic_sequence: str) -> List[str]:
    """
    Validate that the sequence in the evidence file matches the genomic sequence.
    
    Args:
        row: CSV row as dictionary
        row_num: Row number for error reporting
        genomic_sequence: The full genomic sequence from GenBank/file
    
    Returns:
        List of error messages
    """
    errors = []
    
    try:
        # Get coordinates
        start = int(row['start'])
        end = int(row['end'])
        sense = row['sense']
        expected_sequence = row['sequence'].upper()
        
        # Extract sequence from genomic sequence (convert to 0-based for Python)
        extracted_sequence = genomic_sequence[start-1:end]

        # Check sequences are the same length
        if len(expected_sequence) != len(extracted_sequence):
            errors.append(
                f"Row {row_num}: Length mismatch for {row['gene_label']}. "
                f"Length of sequence specified in file: {len(expected_sequence)}, "
                f"Extracted length from coordinates: {len(extracted_sequence)}, "
                f"Coordinates: {start}-{end} ({sense} sense)"
            )
            return errors
        
        # Handle sense orientation
        if sense in ['-', 'minus']:
            extracted_sequence = simple.reverse_complement(extracted_sequence)
        
        # Compare sequences
        if expected_sequence != extracted_sequence:
            errors.append(
                f"Row {row_num}: Sequence mismatch for {row['gene_label']}. "
                f"Expected length: {len(expected_sequence)}, "
                f"Extracted length: {len(extracted_sequence)}, "
                f"Coordinates: {start}-{end} ({sense} sense)"
            )
            
            # Show first few differences for debugging
            min_len = min(len(expected_sequence), len(extracted_sequence))
            differences = []
            for i in range(min_len):
                if expected_sequence[i] != extracted_sequence[i]:
                    differences.append(f"pos {i+1}: expected '{expected_sequence[i]}', got '{extracted_sequence[i]}'")
                    if len(differences) >= 5:  # Limit to first 5 differences
                        differences.append("...")
                        break
            
            if differences:
                errors.append(f"Row {row_num}: First differences: {', '.join(differences)}")
    
    except Exception as e:
        errors.append(f"Row {row_num}: Error validating sequence match: {e}")
    
    return errors


def check_c_exons(gene_description, sense):
    """
    Validate C-gene exon coordinates.
    
    Returns a list of error messages. An empty list indicates no errors.
    
    Rules:
    - For each exon (c_exon_1 to c_exon_8), start and end should either both be null or both be positive integers
    - If both are provided, end must be > start
    - If any exon coordinates are null, all successive exon coordinates must also be null
    - Successive exons must have coordinates higher than the previous exon (+ sense)
    - c_exon_1 should never be null
    """
    errors = []
    
    # Check c_exon_1 is not null
    if gene_description['c_exon_1_start'] is None or gene_description['c_exon_1_end'] is None:
        errors.append(("c_exon_1_start", "C-gene exon 1 coordinates cannot be null"))
        return errors  # No point checking further if exon 1 is null
    
    # Check all 8 exons
    found_null = False
    if sense == '+':
        last_end = -1
    else:
        last_end = sys.maxsize
    
    for i in range(1, 9):
        start_attr = f'c_exon_{i}_start'
        end_attr = f'c_exon_{i}_end'
        
        start = gene_description.get(start_attr, None)
        end = gene_description.get(end_attr, None)
        
        # Check that both are null or both are not null
        if (not start) != (not end):
            errors.append((f"c_exon_{i}_start", f"Exon c_exon_{i}: start and end must both be null or both be provided"))
            continue
        
        # If both are null
        if (not start) and (not end):
            found_null = True
            continue
        
        # If we previously found a null exon, current exon should also be null
        if found_null:
            errors.append((f"c_exon_{i}_start", f"Exon c_exon_{i}: coordinates must be null because a previous exon has null coordinates"))
            continue

        try:
            start = int(start)
            end = int(end)
        except ValueError:
            errors.append((f"c_exon_{i}_start", f"Exon c_exon_{i}: start and end must be integers"))
            continue
        
        # Check that both are positive integers
        if start <= 0:
            errors.append((f"c_exon_{i}_start", f"Exon c_exon_{i}: start coordinate must be a positive integer"))
        if end <= 0:
            errors.append((f"c_exon_{i}_end", f"Exon c_exon_{i}: end coordinate must be a positive integer"))
        
        # Check that end > start
        if end <= start:
            errors.append((f"c_exon_{i}_end", f"Exon c_exon_{i}: end coordinate must be greater than start coordinate"))
        
        # Check that this exon's coordinates are higher than the previous exon
        if sense == '+':
            if start <= last_end:
                errors.append((f"c_exon_{i}_start", f"Exon c_exon_{i}: start coordinate must be greater than the previous exon's end coordinate"))
        else:
            if start >= last_end:
                errors.append((f"c_exon_{i}_start", f"Exon c_exon_{i}: start coordinate must be less than the previous exon's end coordinate (-ve sense sequence)"))
        
        last_end = end
    
    return errors


def check_evidence_file(filepath: str, use_genbank: bool, email: str,
                        save_sequences: bool, skip_sequence_validation: bool) -> Tuple[bool, List[str]]:
    """
    Check the evidence file for errors.
    
    Args:
        filepath: Path to the evidence CSV file
        use_genbank: If True, download sequences from GenBank; if False, read from local files
        email: Email address for GenBank queries
        save_sequences: If True and using GenBank, save sequences to FASTA files
        skip_sequence_validation: If True, skip sequence validation step
    
    Returns:
        Tuple of (success, list_of_errors)
    """
    errors = []
    sequence_cache = {}  # Cache sequences to avoid repeated downloads
    
    if not os.path.exists(filepath):
        return False, [f"Evidence file not found: {filepath}"]
    
    try:
        with open(filepath, 'r', newline='', encoding='utf-8') as csvfile:
            reader = csv.DictReader(csvfile)
            
            for row_num, row in enumerate(reader, start=2):  # Start at 2 (header is row 1)
                # Skip empty rows
                if not any(row.values()):
                    continue
                
                print(f"Checking row {row_num}: {row['gene_label']}")
                
                # Validate required fields
                errors.extend(validate_required_fields(row, row_num))
                
                # Validate enumerated fields
                errors.extend(validate_enumerated_fields(row, row_num))
                
                # Validate coordinates
                errors.extend(validate_coordinates(row, row_num))
                
                # Validate adjacency
                errors.extend(validate_adjacency(row, row_num))
                
                # Skip sequence validation if basic validation failed
                if skip_sequence_validation or any(f"Row {row_num}:" in error for error in errors[-10:]):  # Check recent errors
                    continue
                
                # Get genomic sequence
                accession = row['accession']
                patch = row['patch']
                if patch:
                    accession_with_patch = f"{accession}.{patch}"
                else:
                    accession_with_patch = accession
                
                try:
                    if accession_with_patch in sequence_cache:
                        genomic_sequence = sequence_cache[accession_with_patch]
                    elif use_genbank:
                        genomic_sequence = get_genbank_sequence(accession_with_patch, email, save_sequences)
                        sequence_cache[accession_with_patch] = genomic_sequence
                    else:
                        genomic_sequence = load_sequence_from_file(accession_with_patch)
                        sequence_cache[accession_with_patch] = genomic_sequence
                    
                    # Validate sequence match
                    errors.extend(validate_sequence_match(row, row_num, genomic_sequence))
                    
                except Exception as e:
                    errors.append(f"Row {row_num}: {str(e)}")
    
    except Exception as e:
        return False, [f"Error reading evidence file: {e}"]
    
    return len(errors) == 0, errors


def main():
    parser = argparse.ArgumentParser(
        description="Validate evidence file for OGRDB upload",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Check evidence file using GenBank (download sequences)
  python check_evidence_file.py evidence.csv --email your@email.com
  
  # Check evidence file using local FASTA files
  python check_evidence_file.py evidence.csv --use-local
  
  # Download and save sequences from GenBank for future use
  python check_evidence_file.py evidence.csv --email your@email.com --save-sequences
        """
    )
    
    parser.add_argument('evidence_file',
                        help='Path to the evidence CSV file to validate')
    
    parser.add_argument('--use-local', action='store_true',
                       help='Read sequences from local FASTA files instead of downloading from GenBank')
    
    parser.add_argument('--email', default='user@example.com',
                       help='Email address for GenBank queries (required when using GenBank)')
    
    parser.add_argument('--save-sequences', action='store_true',
                        help='Save downloaded GenBank sequences to FASTA files in current directory')
    
    parser.add_argument('--skip-sequence-validation', action='store_true',
                        help='Skip GenBank sequence validation')
    
    args = parser.parse_args()
    
    # Validate arguments
    if not args.use_local and args.email == 'user@example.com':
        print("Warning: Using default email for GenBank queries. Please provide your email with --email")
        print("This is required by NCBI's usage policies.")
    
    print(f"Checking evidence file: {args.evidence_file}")
    print(f"Sequence source: {'Local FASTA files' if args.use_local else 'GenBank'}")
    if not args.use_local and args.save_sequences:
        print("Will save downloaded sequences to FASTA files")
    print()

    if args.skip_sequence_validation and args.save_sequences:
        print("Error: saving sequences is incompatible with skipping sequence validation.")
        return 1

    success, errors = basic_checks(args.evidence_file)

    if not success:
        print("❌ Basic checks failed:")
        for error in errors:
            print(f"  • {error}")
        return 1
    
    success, errors = check_evidence_file(
        args.evidence_file,
        not args.use_local,
        args.email,
        args.save_sequences,
        args.skip_sequence_validation
    )
    
    if success:
        print("✅ Evidence file validation PASSED - no errors found!")
        return 0
    else:
        print("❌ Evidence file validation FAILED")
        print(f"\nFound {len(errors)} error(s):")
        for error in errors:
            print(f"  • {error}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
