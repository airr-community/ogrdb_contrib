#!/usr/bin/env python3
"""
Debug version to check coordinate transformation for J gene
"""

def transform_coordinates_debug():
    # Simulate J gene data from evidence file
    evidence_row = {
        'gene_label': 'IGHJ1*01',
        'sequence_type': 'J',
        'sense': '-',
        'start': '3000',
        'end': '3100', 
        'gene_start': '3050',
        'gene_end': '3100',
        'j_rs_start': '3000',
        'j_rs_end': '3049',
        'j_cdr3_end': '3075'
    }
    
    print("Evidence row data:")
    for key, value in evidence_row.items():
        print(f"  {key}: {value}")
    
    # Features in order for J genes: j_rs -> gene
    ordered_features = []
    
    # Add j_rs feature if present
    if evidence_row.get('j_rs_start') and evidence_row.get('j_rs_end'):
        j_rs_start = int(evidence_row['j_rs_start'])
        j_rs_end = int(evidence_row['j_rs_end'])
        ordered_features.append(('j_rs', j_rs_start, j_rs_end))
        print(f"Added j_rs feature: {j_rs_start}-{j_rs_end}")
    
    # Add gene feature if present  
    if evidence_row.get('gene_start') and evidence_row.get('gene_end'):
        gene_start = int(evidence_row['gene_start'])
        gene_end = int(evidence_row['gene_end'])
        ordered_features.append(('gene', gene_start, gene_end))
        print(f"Added gene feature: {gene_start}-{gene_end}")
    
    print(f"Ordered features: {ordered_features}")
    
    # Transform coordinates
    transformed = {}
    p = 1  # Position counter for upload file coordinates
    
    for feature_name, orig_start, orig_end in ordered_features:
        feature_length = orig_end - orig_start + 1
        print(f"Processing {feature_name}: length={feature_length}, position={p}")
        
        if feature_name == 'gene':
            transformed['gene_start'] = str(p)
            transformed['gene_end'] = str(p + feature_length - 1)
            print(f"  Set gene_start={p}, gene_end={p + feature_length - 1}")
        elif feature_name == 'j_rs':
            transformed['j_rs_start'] = str(p)
            transformed['j_rs_end'] = str(p + feature_length - 1)
            print(f"  Set j_rs_start={p}, j_rs_end={p + feature_length - 1}")
        
        p += feature_length
    
    transformed['start'] = '1'
    transformed['end'] = str(p - 1)
    
    print("Transformed coordinates:")
    for key, value in transformed.items():
        print(f"  {key}: {value}")

if __name__ == "__main__":
    transform_coordinates_debug()