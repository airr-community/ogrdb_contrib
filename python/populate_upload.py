#!/usr/bin/env python3
"""
Script to transfer genomic evidence data to OGRDB upload files.
Transfers sequence and coordinate data from evidence files to upload files
for entries with genomic inference types.
"""

import argparse
import os
import sys
import csv
from typing import Dict, List, Tuple

from receptor_utils import simple_bio_seq as simple


def basic_checks_upload(upload_file) -> Tuple[bool, List[str]]:
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


def basic_checks_evidence(evidence_file) -> Tuple[bool, List[str]]:
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
        'j_codon_frame', 'j_cdr3_end'
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
 
 
def count_annotated_features(evidence_row: Dict) -> int:
    """
    Count the number of annotated features in an evidence row.
    
    Args:
        evidence_row: Evidence file row dictionary
    
    Returns:
        Number of annotated features
    """
    feature_count = 0
    
    # Define all possible coordinate pair features
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
    
    # Count coordinate pairs that are both non-empty
    for start_field, end_field in coordinate_pairs:
        start_val = evidence_row[start_field]
        end_val = evidence_row[end_field]
        
        if start_val and end_val and start_val != 'None' and end_val != 'None':
            feature_count += 1
    
    # Count single coordinate features
    single_coordinates = ['j_cdr3_end']
    for field in single_coordinates:
        val = evidence_row[field]
        
        if val and val != 'None':
            feature_count += 1
    
    return feature_count


def find_best_evidence_match(gene_label: str, evidence_dict: Dict) -> Dict:
    """
    Find the best evidence record for a gene_label when repository/accession not specified.
    
    Args:
        gene_label: Gene label to search for
        evidence_dict: Dictionary of all evidence records
    
    Returns:
        Best matching evidence record, or None if not found
    """
    # Find all evidence records for this gene_label
    matching_records = []
    for (label, repo, acc, patch), evidence_row in evidence_dict.items():
        if label == gene_label:
            matching_records.append(evidence_row)
    
    if len(matching_records) == 0:
        return None
    elif len(matching_records) == 1:
        return matching_records[0]
    else:
        # Multiple records found - select the one with most annotated features
        best_record = None
        max_features = -1
        
        for record in matching_records:
            feature_count = count_annotated_features(record)
            if feature_count > max_features:
                max_features = feature_count
                best_record = record
        
        return best_record


def read_evidence_file(evidence_path: str) -> Dict[str, Dict]:
    """
    Read evidence file and create lookup dictionary.
    
    Args:
        evidence_path: Path to evidence CSV file
    
    Returns:
        Dictionary with keys as (gene_label, repository, accession, patch) tuples
        and values as evidence row dictionaries
    """
    evidence_dict = {}
    
    if not os.path.exists(evidence_path):
        raise FileNotFoundError(f"Evidence file not found: {evidence_path}")
    
    try:
        evidence_rows = simple.read_csv(evidence_path)
        
        for row in evidence_rows:
            # Skip empty rows
            if not row['gene_label']:
                continue
                
            gene_label = row['gene_label']
            repository = row['repository']
            accession = row['accession']
            patch = row['patch']
            
            key = (gene_label, repository, accession, patch)
            
            if key in evidence_dict:
                print(f"Warning: Duplicate evidence entry found for {key}")
            
            evidence_dict[key] = row
            
    except Exception as e:
        raise ValueError(f"Error reading evidence file: {e}")
    
    return evidence_dict


def transform_coordinates(evidence_row: Dict, sequence_type: str) -> Dict[str, str]:
    """
    Transform coordinates from genomic assembly coordinates to sequence-relative coordinates.
    
    Args:
        evidence_row: Evidence file row dictionary
        sequence_type: V, D, J, or C
    
    Returns:
        Dictionary of transformed coordinate fields
    """
    transformed = {}
    
    # Get the sense to determine if we need to handle reverse complement
    sense = evidence_row['sense']
    
    # Helper function to get coordinate value
    def get_coord(field_name):
        val = evidence_row[field_name]
        if not val or val in ['', 'None']:
            return None
        try:
            return int(val)
        except ValueError:
            return None
    
    # For coordinate transformation, we need to know the order of features
    # in the evidence sequence and transform them to be relative to position 1
    
    # Collect all features with their coordinates
    features = []
    
    # Common coordinates
    gene_start = get_coord('gene_start')
    gene_end = get_coord('gene_end')
    if gene_start is not None and gene_end is not None:
        features.append(('gene', gene_start, gene_end))
    
    if sequence_type == 'V':
        # V-specific features
        utr_5_start = get_coord('utr_5_prime_start')
        utr_5_end = get_coord('utr_5_prime_end')
        if utr_5_start is not None and utr_5_end is not None:
            features.append(('utr_5_prime', utr_5_start, utr_5_end))
        
        leader_1_start = get_coord('leader_1_start')
        leader_1_end = get_coord('leader_1_end')
        if leader_1_start is not None and leader_1_end is not None:
            features.append(('leader_1', leader_1_start, leader_1_end))
        
        leader_2_start = get_coord('leader_2_start')
        leader_2_end = get_coord('leader_2_end')
        if leader_2_start is not None and leader_2_end is not None:
            features.append(('leader_2', leader_2_start, leader_2_end))
        
        v_rs_start = get_coord('v_rs_start')
        v_rs_end = get_coord('v_rs_end')
        if v_rs_start is not None and v_rs_end is not None:
            features.append(('v_rs', v_rs_start, v_rs_end))
    
    elif sequence_type == 'D':
        # D-specific features
        d_rs_5_start = get_coord('d_rs_5_prime_start')
        d_rs_5_end = get_coord('d_rs_5_prime_end')
        if d_rs_5_start is not None and d_rs_5_end is not None:
            features.append(('d_rs_5_prime', d_rs_5_start, d_rs_5_end))
        
        d_rs_3_start = get_coord('d_rs_3_prime_start')
        d_rs_3_end = get_coord('d_rs_3_prime_end')
        if d_rs_3_start is not None and d_rs_3_end is not None:
            features.append(('d_rs_3_prime', d_rs_3_start, d_rs_3_end))
    
    elif sequence_type == 'J':
        # J-specific features
        j_rs_start = get_coord('j_rs_start')
        j_rs_end = get_coord('j_rs_end')
        if j_rs_start is not None and j_rs_end is not None:
            features.append(('j_rs', j_rs_start, j_rs_end))
        
        j_cdr3_end = get_coord('j_cdr3_end')
        if j_cdr3_end is not None:
            transformed['j_cdr3_end'] = str(j_cdr3_end)
    
    if not features:
        return transformed
    
    # Sort features by start coordinate (handle both + and - sense)
    if sense in ['-', 'minus']:
        # For minus sense, features are in reverse order in the genomic sequence
        # but should be in forward order in the final sequence
        features.sort(key=lambda x: x[1], reverse=True)
    else:
        features.sort(key=lambda x: x[1])
    
    # Get the start of the evidence sequence for coordinate adjustment
    evidence_start = get_coord('start')
    if evidence_start is None:
        return transformed
    
    # Transform coordinates using the approach from make_ogrdb_upload.py
    p = 1  # Current position in the transformed sequence
    
    for feature_name, feat_start, feat_end in features:
        # Calculate the feature length
        feature_length = feat_end - feat_start + 1
        
        # Set the transformed start and end coordinates
        if feature_name == 'gene':
            transformed['gene_start'] = str(p)
            transformed['gene_end'] = str(p + feature_length - 1)
        elif feature_name == 'utr_5_prime':
            transformed['utr_5_prime_start'] = str(p)
            transformed['utr_5_prime_end'] = str(p + feature_length - 1)
        elif feature_name == 'leader_1':
            transformed['leader_1_start'] = str(p)
            transformed['leader_1_end'] = str(p + feature_length - 1)
        elif feature_name == 'leader_2':
            transformed['leader_2_start'] = str(p)
            transformed['leader_2_end'] = str(p + feature_length - 1)
        elif feature_name == 'v_rs':
            transformed['v_rs_start'] = str(p)
            transformed['v_rs_end'] = str(p + feature_length - 1)
        elif feature_name == 'd_rs_5_prime':
            transformed['d_rs_5_prime_start'] = str(p)
            transformed['d_rs_5_prime_end'] = str(p + feature_length - 1)
        elif feature_name == 'd_rs_3_prime':
            transformed['d_rs_3_prime_start'] = str(p)
            transformed['d_rs_3_prime_end'] = str(p + feature_length - 1)
        elif feature_name == 'j_rs':
            transformed['j_rs_start'] = str(p)
            transformed['j_rs_end'] = str(p + feature_length - 1)
        
        # Move to next position
        p += feature_length
    
    # Set overall sequence coordinates
    transformed['start'] = '1'
    transformed['end'] = str(p - 1)
    
    return transformed


def rearranged_fixups(upload_row):
    """
    Apply any necessary fixups for rearranged sequences.
    
    Args:
        upload_row: Upload file row dictionary
    
    Returns:
        None (modifies upload_row in place)
    """
    # set the type based on gene_label if not already set
    if not upload_row['type'] and len(upload_row['gene_label']) >= 4:
        locus = upload_row['gene_label'][:3]
        gtype = upload_row['gene_label'][3]
        if gtype in ['V', 'D', 'J', 'C'] and locus in ['IGH', 'IGK', 'IGL', 'TRA', 'TRB', 'TRG', 'TRD']:
            upload_row['type'] = upload_row['gene_label'][:4]
        else:
            print(f"Warning: Could not determine type from gene_label '{upload_row['gene_label']}'")

    upload_row['mapped'] = 'N'
    upload_row['start'] = upload_row['gene_start'] = 1
    upload_row['end'] = upload_row['gene_end'] = len(upload_row['sequence']) if upload_row['sequence'] else 0


def transfer_evidence_to_upload(evidence_path: str, upload_path: str) -> Tuple[bool, List[str]]:
    """
    Transfer genomic evidence data to upload file.
    
    Args:
        evidence_path: Path to evidence CSV file
        upload_path: Path to upload CSV file
    
    Returns:
        Tuple of (success, list_of_errors)
    """
    errors = []
    
    # Read evidence file
    try:
        evidence_dict = read_evidence_file(evidence_path)
        print(f"Loaded {len(evidence_dict)} evidence entries")
    except Exception as e:
        return False, [str(e)]
    
    # Read upload file
    if not os.path.exists(upload_path):
        return False, [f"Upload file not found: {upload_path}"]
    
    try:
        upload_rows = simple.read_csv(upload_path)
        print(f"Processing {len(upload_rows)} upload entries")
    except Exception as e:
        return False, [f"Error reading upload file: {e}"]
    
    # Process each upload row
    updated_rows = []
    updates_made = 0
    
    for row_num, upload_row in enumerate(upload_rows, start=2):
        # Skip empty rows
        if not upload_row['gene_label']:
            updated_rows.append(upload_row)
            continue
        
        gene_label = upload_row['gene_label']
        inference_type = upload_row['inference_type']
        
        print(f"Processing row {row_num}: {gene_label} ({inference_type})")
        
        # Only process genomic inference types
        if inference_type not in ['Unrearranged', 'Unrearranged and rearranged']:
            rearranged_fixups(upload_row)
            updated_rows.append(upload_row)
            continue
        
        # Check required fields - but allow automatic lookup if missing
        repository = upload_row['repository']
        accession = upload_row['accession']
        
        evidence_row = None
        evidence_source = "exact match"
        
        if repository and accession:
            # Repository and accession specified - do exact lookup
            patch = upload_row['patch']
            evidence_key = (gene_label, repository, accession, patch)
            
            if evidence_key not in evidence_dict:
                errors.append(f"Row {row_num}: No matching evidence entry found for {evidence_key}")
                updated_rows.append(upload_row)
                continue
            
            evidence_row = evidence_dict[evidence_key]
        
        elif not repository and not accession:
            # Neither specified - find best match by gene_label
            evidence_row = find_best_evidence_match(gene_label, evidence_dict)
            evidence_source = "automatic best match"
            
            if evidence_row is None:
                errors.append(f"Row {row_num}: No evidence records found for gene_label '{gene_label}'")
                updated_rows.append(upload_row)
                continue
            
            # Update upload row with found repository/accession info
            upload_row['repository'] = evidence_row['repository']
            upload_row['accession'] = evidence_row['accession']
            upload_row['patch'] = evidence_row['patch']
            
        else:
            # Only one of repository/accession specified - this is an error
            missing_field = 'accession' if repository else 'repository'
            errors.append(f"Row {row_num}: {missing_field} field is blank but repository/accession should be specified together for genomic inference type")
            updated_rows.append(upload_row)
            continue
        
        # Get sequence type
        sequence_type = evidence_row['sequence_type']
        if not sequence_type:
            errors.append(f"Row {row_num}: sequence_type not specified in evidence file")
            updated_rows.append(upload_row)
            continue
        
        # Transfer sequence
        evidence_sequence = evidence_row['sequence']
        if evidence_sequence:
            upload_row['sequence'] = evidence_sequence
        
        # Transform and transfer coordinates
        try:
            transformed_coords = transform_coordinates(evidence_row, sequence_type)
            
            # Update upload row with transformed coordinates
            for coord_field, coord_value in transformed_coords.items():
                if coord_field in upload_row:  # Only update fields that exist in original CSV
                    upload_row[coord_field] = coord_value
            
            updates_made += 1
            print(f"  ✅ Updated {gene_label} with genomic evidence ({evidence_source})")
            
        except Exception as e:
            errors.append(f"Row {row_num}: Error transforming coordinates for {gene_label}: {e}")
        
        updated_rows.append(upload_row)
    
    # Write updated upload file
    try:
        # Create backup
        backup_path = upload_path + ".backup"
        simple.write_csv(backup_path, simple.read_csv(upload_path))
        print(f"Created backup: {backup_path}")
        
        # Write updated file
        simple.write_csv(upload_path, updated_rows)
        print(f"Updated upload file: {upload_path}")
        print(f"Made {updates_made} updates")
        
    except Exception as e:
        errors.append(f"Error writing updated upload file: {e}")
        return False, errors
    
    return len(errors) == 0, errors


def main():
    parser = argparse.ArgumentParser(
        description="Transfer genomic evidence data to OGRDB upload file",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Transfer evidence data to upload file
  python transfer_evidence_to_upload.py evidence.csv upload.csv
  
  # The script will:
  # 1. Read the evidence file and create a lookup dictionary
  # 2. Process each upload file row with genomic inference type
  # 3. Find matching evidence entries and transfer sequence/coordinates
  # 4. Transform coordinates from genomic to sequence-relative
  # 5. Create backup of original upload file before updating
        """
    )
    
    parser.add_argument('evidence_file',
                        help='Path to the evidence CSV file')
    
    parser.add_argument('upload_file',
                        help='Path to the upload CSV file to update')
    
    args = parser.parse_args()
    
    print(f"Evidence file: {args.evidence_file}")
    print(f"Upload file: {args.upload_file}")
    print()

    # Basic checks
    success, errors = basic_checks_evidence(args.evidence_file)
    if not success:
        print("❌ Evidence file checks failed:")
        for error in errors:
            print(f"  • {error}")
        return 1
    
    success, errors = basic_checks_upload(args.upload_file)
    if not success:
        print("❌ Upload file checks failed:")
        for error in errors:
            print(f"  • {error}")
        return 1
    
    success, errors = transfer_evidence_to_upload(args.evidence_file, args.upload_file)
    
    if success and len(errors) == 0:
        print("✅ Evidence transfer completed successfully!")
        return 0
    else:
        if success:
            print("⚠️  Evidence transfer completed with warnings")
        else:
            print("❌ Evidence transfer failed")
        
        if errors:
            print(f"\nFound {len(errors)} issue(s):")
            for error in errors:
                print(f"  • {error}")
        
        return 1 if not success else 0


if __name__ == "__main__":
    sys.exit(main())