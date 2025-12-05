#!/usr/bin/env python3
"""
Script to validate OGRDB upload files for completeness and consistency.
Checks upload CSV files against evidence files for proper formatting and data integrity.
"""

import argparse
import csv
import os
import sys
from typing import Dict, List, Tuple

from receptor_utils import simple_bio_seq as simple

from populate_upload import find_best_evidence_match


def basic_checks(upload_file) -> Tuple[bool, List[str]]:
    """
    Perform basic checks:
     - Check that expected columns are present in the evidence file
     - Check that no cells contain leading or trailing spaces
    
    Returns:
        Tuple of (success, list_of_errors)
    """
    errors = []

    expected_columns = [
        'gene_label', 'imgt', 'functionality', 'type', 'inference_type', 'sequence', 'sequence_gapped', 'species_subgroup', 'subgroup_type', 'alt_names',
        'affirmation', 'chromosome', 'paralogs', 'varb_rep', 'notes', 'inferred_extension', 'ext_3_prime', 'ext_5_prime', 'curational_tags', 'mapped',
        'gene_start', 'gene_end', 'utr_5_prime_start', 'utr_5_prime_end', 'leader_1_start', 'leader_1_end', 'leader_2_start', 'leader_2_end', 
        'v_rs_start', 'v_rs_end', 'd_rs_3_prime_start', 'd_rs_3_prime_end', 'd_rs_5_prime_start', 'd_rs_5_prime_end', 'j_codon_frame', 'j_rs_start', 'j_rs_end', 'j_cdr3_end'
    ]

    # Check if file exists
    if not os.path.exists(upload_file):
        errors.append(f"Upload file not found: {upload_file}")
        return False, errors

    with open(upload_file, 'r') as f:
        reader = csv.DictReader(f)
        # Check for expected columns
        missing_columns = [col for col in expected_columns if col not in reader.fieldnames]
        if missing_columns:
            errors.append(f"Missing columns in upload file: {', '.join(missing_columns)}")

        # Check for leading/trailing spaces
        for row in reader:
            for col in expected_columns:
                if col in row:
                    if row[col].strip() != row[col]:
                        errors.append(f"Leading/trailing spaces found in row {row['gene_label']} column '{col}'")

    return len(errors) == 0, errors


def validate_required_fields_upload(row: Dict, row_num: int) -> List[str]:
    """
    Validate that required fields are present and have appropriate types/enum values.
    
    Args:
        row: CSV row as dictionary
        row_num: Row number for error reporting
    
    Returns:
        List of error messages
    """
    errors = []
    
    # Required fields with their validation
    required_fields = [
        ('gene_label', str, None),
        ('functionality', str, ['F', 'ORF', 'P']),
        ('inference_type', str, ['Unrearranged', 'Unrearranged and rearranged', 'Rearranged only']),
        ('sequence', str, None),
        ('sequence_gapped', str, None),  # Will be validated separately for V genes
        ('affirmation', int, None),
        ('gene_start', int, None),
        ('gene_end', int, None)
    ]
    
    for field_name, field_type, enum_values in required_fields:
        value = row[field_name]
        
        # Check if field is empty
        if not value:
            # mapped can be empty for rearranged only
            if field_name == 'mapped' and row['inference_type'] == 'Rearranged only':
                continue
            # sequence_gapped can be empty for non-V genes
            elif field_name == 'sequence_gapped' and row['type'][:-1] != 'V':
                continue
            else:
                errors.append(f"Upload file row {row_num}: Required field '{field_name}' is missing or empty")
                continue
        
        # Type validation
        if field_type == int:
            try:
                int(value)
            except ValueError:
                errors.append(f"Row {row_num}: Field '{field_name}' must be an integer, got '{value}'")
                continue
        
        # Enum validation
        if enum_values and value not in enum_values:
            errors.append(f"Upload file row {row_num}: Invalid {field_name} '{value}'. Must be one of: {enum_values}")

        if len(row['type']) != 4:
            errors.append(f"Upload file row {row_num}: Invalid type '{row['type']}'. Must be 4 characters like 'IGHV', 'IGHD', 'IGHJ', 'IGHC', etc.")
        elif row['type'][3] not in ['V', 'D', 'J', 'C']:
            errors.append(f"Upload file row {row_num}: Invalid type '{row['type']}'. Final character must be one of V, D, J, C.")
        elif row['type'][:3] not in ['IGH', 'IGK', 'IGL', 'TRA', 'TRB', 'TRG', 'TRD']:
            errors.append(f"Upload file row {row_num}: Invalid type '{row['type']}'. Must start with one of IGH, IGK, IGL, TRA, TRB, TRG, TRD.")

    return errors


def validate_v_gene_gapped_sequence(row: Dict, row_num: int) -> List[str]:
    """
    Validate V gene gapped sequence requirements.
    
    Args:
        row: CSV row as dictionary
        row_num: Row number for error reporting
    
    Returns:
        List of error messages
    """
    errors = []
    
    if len(row['type']) != 4 or row['type'][-1] != 'V':
        return errors
        
    sequence = row['sequence']
    sequence_gapped = row['sequence_gapped']
    
    if not sequence_gapped:
        errors.append(f"Row {row_num}: V gene must have sequence_gapped field completed")
        return errors
    
    # Check that sequence_gapped contains periods
    if '.' not in sequence_gapped:
        errors.append(f"Row {row_num}: V gene sequence_gapped must contain periods (.)")
        return errors
    
    # Check that removing periods gives the same sequence as 'sequence'
    ungapped_sequence = sequence_gapped.replace('.', '').upper()
    start = int(row['gene_start'])
    end = int(row['gene_end'])
    if sequence and ungapped_sequence != sequence[start-1:end].upper():
        errors.append(f"Row {row_num}: V gene sequence_gapped without periods does not match sequence field")
        errors.append(f"  Expected: {sequence[start-1:end]}")
        errors.append(f"  Got: {ungapped_sequence}")
    
    return errors


def validate_coordinate_pairs_upload(row: Dict, row_num: int) -> List[str]:
    """
    Validate coordinate pairs in upload file.
    
    Args:
        row: CSV row as dictionary
        row_num: Row number for error reporting
    
    Returns:
        List of error messages
    """
    errors = []
    
    # Coordinate pairs to validate
    coordinate_pairs = [
        ('gene_start', 'gene_end'),
        ('utr_5_prime_start', 'utr_5_prime_end'),
        ('leader_1_start', 'leader_1_end'),
        ('leader_2_start', 'leader_2_end'),
        ('v_rs_start', 'v_rs_end'),
        ('d_rs_3_prime_start', 'd_rs_3_prime_end'),
        ('d_rs_5_prime_start', 'd_rs_5_prime_end'),
        ('j_rs_start', 'j_rs_end')
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
        
        # Validate they are positive integers with start < end
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


def validate_feature_adjacency_upload(row: Dict, row_num: int) -> List[str]:
    """
    Validate feature adjacency in upload file using same rules as check_evidence_file.py.
    Note: Upload coordinates are always in + sense (sequence-relative).
    
    Args:
        row: CSV row as dictionary
        row_num: Row number for error reporting
    
    Returns:
        List of error messages
    """
    errors = []
    sequence_type = row['type']
    
    if not sequence_type:
        return errors  # Skip if sequence_type is missing (caught by other validation)
    
    # Helper function to get coordinate value or None
    def get_coord(field_name):
        val = row[field_name]
        if not val or val == 'None':
            return None
        try:
            return int(val)
        except ValueError:
            return None
    
    # Get coordinates
    gene_start = get_coord('gene_start')
    gene_end = get_coord('gene_end')
    features = []

    if sequence_type == 'V':
        # V gene adjacency validation (always + sense in upload file)
        utr_5_start = get_coord('utr_5_prime_start')
        utr_5_end = get_coord('utr_5_prime_end')
        leader_1_start = get_coord('leader_1_start')
        leader_1_end = get_coord('leader_1_end')
        leader_2_start = get_coord('leader_2_start')
        leader_2_end = get_coord('leader_2_end')
        v_rs_start = get_coord('v_rs_start')
        v_rs_end = get_coord('v_rs_end')
            
        # Order: utr_5 -> leader_1 -> leader_2 -> gene -> v_rs
        if utr_5_start is not None and utr_5_end is not None:
            features.append(('utr_5_prime', utr_5_start, utr_5_end))
        
        if leader_1_start is not None and leader_1_end is not None:
            features.append(('leader_1', leader_1_start, leader_1_end))
            # Check adjacency with UTR 5'
            if utr_5_end is not None and leader_1_start != utr_5_end + 1:
                errors.append(f"Row {row_num}: leader_1_start ({leader_1_start}) should be adjacent to utr_5_prime_end + 1 ({utr_5_end + 1})")
        
        if leader_2_start is not None and leader_2_end is not None:
            features.append(('leader_2', leader_2_start, leader_2_end))
        
        if gene_start is not None and gene_end is not None:
            features.append(('gene', gene_start, gene_end))
            # Check adjacency with last leader
            last_leader_end = None
            if leader_2_end is not None:
                last_leader_end = leader_2_end
            elif leader_1_end is not None:
                last_leader_end = leader_1_end
            
            if last_leader_end is not None and gene_start != last_leader_end + 1:
                errors.append(f"Row {row_num}: gene_start ({gene_start}) should be adjacent to last leader end + 1 ({last_leader_end + 1})")
        
        if v_rs_start is not None and v_rs_end is not None:
            features.append(('v_rs', v_rs_start, v_rs_end))
            if gene_end is not None and v_rs_start != gene_end + 1:
                errors.append(f"Row {row_num}: v_rs_start ({v_rs_start}) should be adjacent to gene_end + 1 ({gene_end + 1})")
    
    elif sequence_type == 'D':
        # D gene adjacency validation
        d_rs_5_start = get_coord('d_rs_5_prime_start')
        d_rs_5_end = get_coord('d_rs_5_prime_end')
        d_rs_3_start = get_coord('d_rs_3_prime_start')
        d_rs_3_end = get_coord('d_rs_3_prime_end')
        
        # Order: d_rs_5 -> gene -> d_rs_3
        if d_rs_5_start is not None and d_rs_5_end is not None:
            features.append(('d_rs_5_prime', d_rs_5_start, d_rs_5_end))
            if gene_start is not None and gene_start != d_rs_5_end + 1:
                errors.append(f"Row {row_num}: gene_start ({gene_start}) should be adjacent to d_rs_5_prime_end + 1 ({d_rs_5_end + 1})")
        
        if gene_start is not None and gene_end is not None:
            features.append(('gene', gene_start, gene_end))
        
        if d_rs_3_start is not None and d_rs_3_end is not None:
            features.append(('d_rs_3_prime', d_rs_3_start, d_rs_3_end))
            if gene_end is not None and d_rs_3_start != gene_end + 1:
                errors.append(f"Row {row_num}: d_rs_3_prime_start ({d_rs_3_start}) should be adjacent to gene_end + 1 ({gene_end + 1})")
    
    elif sequence_type == 'J':
        # J gene adjacency validation
        j_rs_start = get_coord('j_rs_start')
        j_rs_end = get_coord('j_rs_end')
        
        # Order: j_rs -> gene
        if j_rs_start is not None and j_rs_end is not None:
            features.append(('j_rs', j_rs_start, j_rs_end))
            if gene_start is not None and gene_start != j_rs_end + 1:
                errors.append(f"Row {row_num}: gene_start ({gene_start}) should be adjacent to j_rs_end + 1 ({j_rs_end + 1})")
        
        if gene_start is not None and gene_end is not None:
            features.append(('gene', gene_start, gene_end))
    
    return errors


def validate_rearranged_only(row: Dict, row_num: int) -> List[str]:
    """
    Validate requirements for rearranged only inference type.
    
    Args:
        row: CSV row as dictionary
        row_num: Row number for error reporting
    
    Returns:
        List of error messages
    """
    errors = []
    
    inference_type = row['inference_type']
    if inference_type != 'Rearranged only':
        return errors
    
    gene_type = row['type'][3]
    
    # Check that gene_start, gene_end are completed
    gene_start = row['gene_start']
    gene_end = row['gene_end']
    
    if not gene_start or not gene_end:
        errors.append(f"Row {row_num}: Rearranged only entries must have gene_start and gene_end completed")
        return errors
    
    start = row['start']
    end = row['end']

    if not start or not end:
        errors.append(f"Row {row_num}: Rearranged only entries must have gene_start and gene_end completed")
        return errors
    
    # Check that no other coordinates have values, except j_cdr3_end for J genes
    forbidden_coords = [
        'utr_5_prime_start', 'utr_5_prime_end',
        'leader_1_start', 'leader_1_end', 'leader_2_start', 'leader_2_end',
        'v_rs_start', 'v_rs_end', 'd_rs_3_prime_start', 'd_rs_3_prime_end',
        'd_rs_5_prime_start', 'd_rs_5_prime_end', 'j_rs_start', 'j_rs_end'
    ]
    
    # For J genes, allow j_cdr3_end
    if gene_type != 'J':
        forbidden_coords.append('j_cdr3_end')
    
    for coord_field in forbidden_coords:
        value = row[coord_field]
        if value:
            errors.append(f"Row {row_num}: Rearranged only entries should not have {coord_field} coordinate")
    
    # For J genes, check that j_cdr3_end is present
    if gene_type == 'J':
        j_cdr3_end = row['j_cdr3_end']
        if not j_cdr3_end:
            errors.append(f"Row {row_num}: J gene rearranged only entries must have j_cdr3_end completed")
    
    return errors


def validate_against_evidence(row: Dict, row_num: int, evidence_rows: List[Dict]) -> List[str]:
    """
    Validate upload row against corresponding evidence file entry.
    
    Args:
        row: Upload CSV row as dictionary
        row_num: Row number for error reporting
        evidence_dict: Evidence file lookup dictionary
    
    Returns:
        List of error messages
    """
    errors = []
    
    inference_type = row['inference_type']
    if inference_type == 'Rearranged only':
        return errors  # Skip for rearranged only   
    
    # Look up evidence entry

    evidence_row = find_best_evidence_match(row['gene_label'], evidence_rows)

    if not evidence_row:
        errors.append(f"Row {row_num}: No evidence entry found for gene_label '{row['gene_label']}'")
        return errors

    # Check that same features have coordinates and have equal lengths
    coordinate_pairs = [
        ('gene_start', 'gene_end'),
        ('utr_5_prime_start', 'utr_5_prime_end'),
        ('leader_1_start', 'leader_1_end'),
        ('leader_2_start', 'leader_2_end'),
        ('v_rs_start', 'v_rs_end'),
        ('d_rs_3_prime_start', 'd_rs_3_prime_end'),
        ('d_rs_5_prime_start', 'd_rs_5_prime_end'),
        ('j_rs_start', 'j_rs_end')
    ]
    
    for start_field, end_field in coordinate_pairs:
        upload_start = row[start_field]
        upload_end = row[end_field]
        evidence_start = evidence_row[start_field]
        evidence_end = evidence_row[end_field]
        
        # Check if both files have or don't have this coordinate pair
        upload_has_coords = bool(upload_start and upload_end)
        evidence_has_coords = bool(evidence_start and evidence_end)
        
        if upload_has_coords != evidence_has_coords:
            status = "present" if upload_has_coords else "missing"
            evidence_status = "present" if evidence_has_coords else "missing"
            errors.append(f"Row {row_num}: Feature {start_field[:-6]} coordinates {status} in upload but {evidence_status} in evidence")
            continue
        
        # If both have coordinates, check lengths match
        if upload_has_coords and evidence_has_coords:
            try:
                upload_length = int(upload_end) - int(upload_start) + 1
                evidence_length = int(evidence_end) - int(evidence_start) + 1
                
                if upload_length != evidence_length:
                    errors.append(f"Row {row_num}: Feature {start_field[:-6]} length mismatch - upload: {upload_length}, evidence: {evidence_length}")
            except ValueError:
                # Coordinate validation errors will be caught elsewhere
                pass
    
    # For J genes, check j_cdr3_end matches
    gene_type = evidence_row['sequence_type']
    if gene_type == 'J':
        upload_cdr3_end = row['j_cdr3_end']
        evidence_cdr3_end = evidence_row['j_cdr3_end']
        
        # Both should have j_cdr3_end if it's a J gene
        upload_has_cdr3 = bool(upload_cdr3_end)
        evidence_has_cdr3 = bool(evidence_cdr3_end)
        
        if upload_has_cdr3 != evidence_has_cdr3:
            status = "present" if upload_has_cdr3 else "missing"
            evidence_status = "present" if evidence_has_cdr3 else "missing"
            errors.append(f"Row {row_num}: j_cdr3_end {status} in upload but {evidence_status} in evidence")
        elif upload_has_cdr3 and evidence_has_cdr3:
            if upload_cdr3_end != evidence_cdr3_end:
                errors.append(f"Row {row_num}: j_cdr3_end mismatch - upload: {upload_cdr3_end}, evidence: {evidence_cdr3_end}")
    
    return errors


def check_upload_file(upload_path: str, evidence_path: str) -> Tuple[bool, List[str]]:
    """
    Check upload file for completeness and consistency.
    
    Args:
        upload_path: Path to upload CSV file
        evidence_path: Path to evidence CSV file
    
    Returns:
        Tuple of (success, list_of_errors)
    """
    errors = []
    
    # Read evidence file
    evidence_rows = simple.read_csv(evidence_path)
    print(f"Loaded {len(evidence_rows)} evidence entries for validation")
    
    # Read upload file
    if not os.path.exists(upload_path):
        return False, [f"Upload file not found: {upload_path}"]
    
    try:
        with open(upload_path, 'r', newline='', encoding='utf-8') as csvfile:
            reader = csv.DictReader(csvfile)
            
            for row_num, row in enumerate(reader, start=2):
                # Skip empty rows
                if not any(row.values()):
                    continue
                
                gene_label = row['gene_label']
                print(f"Checking row {row_num}: {gene_label}")
                
                # Basic field validation
                errors.extend(validate_required_fields_upload(row, row_num))
                
                # Coordinate pair validation
                errors.extend(validate_coordinate_pairs_upload(row, row_num))
                
                # Feature adjacency validation
                errors.extend(validate_feature_adjacency_upload(row, row_num))
                
                # Rearranged only specific validation
                errors.extend(validate_rearranged_only(row, row_num))
                
                # V gene gapped sequence validation
                errors.extend(validate_v_gene_gapped_sequence(row, row_num))
                
                # Evidence file consistency validation
                errors.extend(validate_against_evidence(row, row_num, evidence_rows))
    
    except Exception as e:
        return False, [f"Error reading upload file: {e}"]
    
    return len(errors) == 0, errors


def main():
    parser = argparse.ArgumentParser(
        description="Validate upload file for OGRDB submission",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Check upload file against evidence file
  python check_upload_file.py evidence.csv upload.csv
  
  # The script will validate:
  # 1. Required field completion and types
  # 2. V gene gapped sequence consistency
  # 3. Coordinate pair integrity and adjacency
  # 4. Rearranged only specific requirements
  # 5. Consistency with evidence file data
        """
    )
    parser.add_argument('evidence_file',
                        help='Path to the evidence CSV file')    
    parser.add_argument('upload_file',
                        help='Path to the upload CSV file to validate')
    args = parser.parse_args()
    
    print(f"Evidence file: {args.evidence_file}")
    print(f"Upload file: {args.upload_file}")
    print()

    success, errors = basic_checks(args.upload_file)
    if not success:
        print("❌ Basic checks failed:")
        for error in errors:
            print(f"  • {error}")
        return 1
    
    success, errors = check_upload_file(args.upload_file, args.evidence_file)
    
    if success:
        print("✅ Upload file validation PASSED - no errors found!")
        return 0
    else:
        print("❌ Upload file validation FAILED")
        print(f"\nFound error(s):")
        for error in errors:
            print(f"  • {error}")
        return 1


if __name__ == "__main__":
    sys.exit(main())