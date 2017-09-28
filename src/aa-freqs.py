import gzip
import os.path
from Bio import SeqIO


UNIPROT = "/nfs/ftp/pub/databases/uniprot/current_release/"
PROTEOMES = "knowledgebase/reference_proteomes/"
HUMAN_PROTEOME = "Eukaryota/UP000005640_9606.fasta.gz"


def calc_aa_freqs(proteome):
    aa_freqs = dict()
    aa_count = 0
    for seqrec in proteome:
        for aa in seqrec.seq:
            if aa not in aa_freqs:
                aa_freqs[aa] = 1
            else:
                aa_freqs[aa] = aa_freqs[aa] + 1
            aa_count = aa_count + 1
    for aa in aa_freqs:
        aa_freqs[aa] = float(aa_freqs[aa])/float(aa_count)
    return aa_freqs


if __name__ == "__main__":
    prot_file = UNIPROT + PROTEOMES + HUMAN_PROTEOME
    with gzip.open(prot_file, 'rt') as h:
        proteome = list(SeqIO.parse(h, "fasta"))
    aa_freqs = calc_aa_freqs(proteome)
    for aa in aa_freqs:
        print("{0}\t{1}".format(aa, aa_freqs[aa]))
