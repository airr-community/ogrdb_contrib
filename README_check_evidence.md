# Evidence File Validation Script

## Overview

`check_evidence_file.py` validates evidence CSV files for OGRDB upload by checking:

1. **Required field completion** - ensures all mandatory fields are filled
2. **Field value validation** - checks enumerated values and data types
3. **Coordinate validation** - verifies coordinates are positive integers with start < end
4. **Sequence matching** - compares sequences against genomic references

## Prerequisites

- Python 3.6+
- `receptor_utils` library
- `biopython` library

Install dependencies:
```bash
pip install biopython receptor-utils
```

## Usage

### Basic Usage

Check evidence file using GenBank (downloads sequences):
```bash
python check_evidence_file.py evidence.csv --email your@email.com
```

Check evidence file using local FASTA files:
```bash
python check_evidence_file.py evidence.csv --use-local
```

Download and save sequences for future use:
```bash
python check_evidence_file.py evidence.csv --email your@email.com --save-sequences
```

### Arguments

- `evidence_file`: Path to the evidence CSV file to validate
- `--use-local`: Read sequences from local FASTA files instead of GenBank
- `--email`: Email address for GenBank queries (required by NCBI)
- `--save-sequences`: Save downloaded GenBank sequences to FASTA files

### Local FASTA Files

When using `--use-local`, the script expects FASTA files named `<accession>.fasta` in the current directory.

For example, if your evidence file references accession `NC_059450.1`, the script will look for `NC_059450.1.fasta`.

## Validation Rules

### Required Fields
- gene_label
- sequence  
- sequence_type
- repository
- accession
- start, end
- sense
- gene_start, gene_end

### Field Validation
- `sequence_type`: Must be V, D, J, or C
- `sense`: Must be +, -, plus, or minus
- Coordinates: Must be positive integers (1-based) with start < end
- `j_cdr3_end`: Must be positive integer for J genes

### Adjacency Validation
The script validates that features are properly adjacent according to OGRDB requirements:

**V genes:**
- `utr_5_prime` → `leader_1` (must be adjacent)
- `leader_1` → `leader_2` (gap allowed)  
- `leader_1/leader_2` → `gene` → `v_rs` (must be adjacent)
- Overall `start` coordinate must point to first feature
- Overall `end` coordinate must point to last feature

**D genes:**
- `d_rs_5_prime` → `gene` → `d_rs_3_prime` (all must be adjacent)
- Overall `start`/`end` must span the entire sequence

**J genes:**
- `j_rs` → `gene` (must be adjacent)
- Overall `start`/`end` must span the entire sequence

### Sequence Validation
The script extracts the sequence from the genomic reference using the provided coordinates and compares it with the sequence field. For minus sense, it uses reverse complement.

## Output

The script provides detailed error reporting:
- ✅ Success: No errors found
- ❌ Failure: Lists all validation errors with row numbers and descriptions

## Example Error Output

```
❌ Evidence file validation FAILED

Found 5 error(s):
  • Row 2: Required field 'gene_start' is missing or empty
  • Row 3: Invalid sequence_type 'X'. Must be one of: ['V', 'D', 'J', 'C']
  • Row 4: start (100) must be < end (50)
  • Row 5: leader_1_start (125) should be adjacent to utr_5_prime_end + 1 (121)
  • Row 6: Sequence mismatch for IGHV1-1. Expected length: 295, Extracted length: 290
```