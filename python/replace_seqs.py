# replace sequences in evidence file with sequences downloaded from NCBI

import argparse
import csv
import os
import sys
from typing import Dict, List, Tuple

from receptor_utils import simple_bio_seq as simple
from Bio import Entrez, SeqIO

def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Replace sequences in evidence file with sequences downloaded from NCBI.")
    parser.add_argument("evidence_file", type=str, help="Path to the evidence TSV file.")
    parser.add_argument("output_file", type=str, help="Path to the output TSV file with replaced sequences.")
    return parser.parse_args()


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


sequence_cache = {}

args = parse_arguments()
evidence_file = args.evidence_file
output_file = args.output_file
errors = []

with open(evidence_file, "r") as ef, open(output_file, "w", newline='') as of:
    reader = csv.DictReader(ef, delimiter=",")
    writer = csv.DictWriter(of, fieldnames=reader.fieldnames, delimiter=",")
    writer.writeheader()

    for row_num, row in enumerate(reader, start=2):  # Start at 2 (header is row 1)
        # Skip empty rows
        if not any(row.values()):
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
            else:
                genomic_sequence = load_sequence_from_file(accession_with_patch)
                sequence_cache[accession_with_patch] = genomic_sequence
            
        except Exception as e:
            errors.append(f"Row {row_num}: {str(e)}")
            continue

        start = int(row['start'])
        end = int(row['end'])
        sense = row['sense']
        required_sequence = genomic_sequence[start-1:end]

        if sense == '-':
            required_sequence = simple.reverse_complement(required_sequence)

        row['sequence'] = required_sequence.upper()
        writer.writerow(row)

if errors:
    sys.stderr.write("Errors encountered during processing:\n")
    for error in errors:
        sys.stderr.write(f"{error}\n")
    sys.exit(1)

sys.exit(0)