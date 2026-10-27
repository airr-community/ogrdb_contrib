# An easy way to replace aligned V sequences in an upload file with revised sequences from an OGRDB alignment file
# Text layout of the alignment file:
# <gene label> whitespace <aligned V sequence using - or . for gaps> whitespace> rest of line ignored

import argparse
from receptor_utils import simple_bio_seq as simple
from receptor_utils import number_v


def parse_alignment_file(alignment_filepath):
    alignment_dict = {}
    with open(alignment_filepath, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            gene_label = parts[0]
            aligned_sequence = parts[1]

            if '*' not in gene_label or ('TR' not in gene_label and 'IG' not in gene_label):
                print(f"Skipping line: {line}")
                continue

            if len(aligned_sequence) < 100 or len(aligned_sequence) > 120:
                print(f"Skipping line with unexpected sequence length: {gene_label} {aligned_sequence}")
                continue 

            alignment_dict[gene_label] = aligned_sequence

    return alignment_dict


def replace_sequences_in_upload(upload_filepath, alignment_dict, output_filepath):
    upload_recs = simple.read_csv(upload_filepath)
    upload_dict = {rec['gene_label']: rec for rec in upload_recs}

    for gene_label, revised_alignment in alignment_dict.items():
        if gene_label not in upload_dict:
            print(f"Gene label {gene_label} not found in upload file, skipping.")
            continue

        nt_seq = upload_dict[gene_label]['sequence_gapped'].replace('.', '')
        aa_seq = simple.translate(nt_seq)

        if aa_seq != alignment_dict[gene_label].replace('-', '').replace('.', ''):
            print(f"Amino acid sequence mismatch for {gene_label}\n{aa_seq}\n{alignment_dict[gene_label].replace('-', '').replace('.', '')}, skipping.")

        upload_dict[gene_label]['sequence_gapped'] = number_v.gap_nt_from_aa(nt_seq, revised_alignment)

    simple.write_csv(output_filepath, upload_recs)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Replace aligned V sequences in an upload file with revised sequences from an OGRDB alignment file.")
    parser.add_argument("alignment_filepath", help="Path to the alignment file")
    parser.add_argument("upload_filepath", help="Path to the upload CSV file")
    parser.add_argument("output_filepath", help="Path to the output CSV file with updated sequences")

    args = parser.parse_args()

    replace_sequences_in_upload(args.upload_filepath, parse_alignment_file(args.alignment_filepath), args.output_filepath)
