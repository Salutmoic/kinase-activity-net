import sys
import re
import csv
from Bio import SeqIO


def read_proteome(proteome_file):
    proteome = dict()
    proteome_gen = SeqIO.parse(proteome_file, "fasta")
    gene_name_re = re.compile(r"GN=(\S+)")
    for protein in proteome_gen:
        prot_acc = protein.id.split("|")[1]
        prot_gn_match = gene_name_re.search(protein.description)
        if prot_gn_match is None:
            continue
        prot_gn = prot_gn_match.group(1)
        proteome[prot_gn] = {"acc": prot_acc, "seq": str(protein.seq)}
    return proteome


def read_phosphoproteome(phosphoproteome_file):
    # phosphoproteome_file is almost certainly the
    # "phosphoproteome_annotated" file from David.
    phosphoproteome = dict()
    with open(phosphoproteome_file) as h:
        table = csv.DictReader(h, delimiter=",", quotechar="\"")
        for row in table:
            if row["acc"] not in phosphoproteome:
                phosphoproteome[row["acc"]] = [(int(row["position"]),
                                                row["residue"])]
            else:
                phosphoproteome[row["acc"]].append((int(row["position"]),
                                                    row["residue"]))
    return phosphoproteome


def get_motif_seq(seq, pos, res, tail_len):
    i = pos - 1
    if i < 0 or i >= len(seq):
        raise ValueError
    if i - tail_len < 0:
        extra = -1 * (i - tail_len)
        bottom = extra*'_' + seq[0:i]
    else:
        bottom = seq[i-tail_len:i]
    if i + tail_len >= len(seq):
        extra = i + tail_len - len(seq) + 1
        top = seq[i:len(seq)] + extra*'_'
    else:
        top = seq[i:i+tail_len+1]
    return bottom + top


def get_phosphosites(proteome, phosphoproteome, prot_list):
    psites = []
    for prot in prot_list:
        if prot not in proteome:
            continue
        prot_acc = proteome[prot]["acc"]
        if prot_acc not in phosphoproteome:
            continue
        prot_phospho = phosphoproteome[prot_acc]
        for pos, res in prot_phospho:
            motif = get_motif_seq(proteome[prot]["seq"], pos, res, 7)
            psites.append((prot, str(pos), res, motif))
    return psites


if __name__ == "__main__":
    if len(sys.argv) not in [3, 4]:
        sys.exit("USAGE: <script> PROTEOME PHOSPHOPROTEOME [PROT_LIST]")
    proteome_file = sys.argv[1]
    phosphoproteome_file = sys.argv[2]
    proteome = read_proteome(proteome_file)
    if len(sys.argv) == 4:
        prot_list_file = sys.argv[3]
        prot_list = set()
        with open(prot_list_file) as h:
            for line in h:
                prot_list.add(line.strip())
    else:
        prot_list = set(proteome.keys())
    phosphoproteome = read_phosphoproteome(phosphoproteome_file)
    psites = get_phosphosites(proteome, phosphoproteome, prot_list)
    for psite in psites:
        print("\t".join(psite))
