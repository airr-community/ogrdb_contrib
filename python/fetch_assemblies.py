from receptor_utils import simple_bio_seq as simple
from Bio import Entrez, SeqIO


# Fetch a sequence from genbank given the accession number
def get_genbank_sequence(acc):
    print(f'fetching: {acc}')
    Entrez.email = 'william@lees.org.uk'
    handle = Entrez.efetch(db='nucleotide', id=acc, rettype='fasta', retmode='text')
    patch = None

    seq_record = SeqIO.read(handle, "fasta")
    if '.' in seq_record.id:
        patch = seq_record.id.split('.')[-1:][0]
        accession_id = '.'.join(seq_record.id.split('.')[:-1])
    else:
        patch = ''
        accession_id = seq_record.id

    simple.write_fasta(f'assemblies/{accession_id}.{patch}.fasta', {f'{accession_id}.{patch}': str(seq_record.seq)})


loci = ('IGH', 'IGK', 'IGL')
for locus in loci:
    recs = simple.read_csv(f'for_pub/original/gene_usage_summary_{locus}.csv')
    fetched_sequences = []

    for rec in recs:
        if 'Contig' in rec['loc']:
            contig = rec['loc'].replace('Contig:', '')
            contig = contig.split('|')[0]
            if contig not in fetched_sequences:
                get_genbank_sequence(contig)
                fetched_sequences.append(contig)
