import gzip
import numpy as np
from math import log
from statistics import median
from dbm import gnu as gdbm
from Bio import SeqIO
import re


AMINO_ACIDS = {"G": 0, "A": 1, "V": 2, "L": 3, "I": 4, "P": 5, "F": 6, "Y": 7,
               "W": 8, "S": 9, "T": 10, "C": 11, "M": 12, "N": 13, "Q": 14,
               "K": 15, "R": 16, "H": 17, "D": 18, "E": 19}


def read_proteome():
    uniprot = "/nfs/ftp/pub/databases/uniprot/current_release/"
    proteomes = "knowledgebase/reference_proteomes/"
    human_proteome = "Eukaryota/UP000005640_9606.fasta.gz"
    human_add_proteome = "Eukaryota/UP000005640_9606_additional.fasta.gz"
    prot_file = uniprot + proteomes + human_proteome
    prot_add_file = uniprot + proteomes + human_add_proteome
    id_map = dict()
    with open("data/uniprot-id-map.tsv") as h:
        for line in h:
            uniprot, hgnc = line.strip().split()
            id_map[uniprot] = hgnc
    proteome = dict()
    for proteome_file in [prot_file, prot_add_file]:
        with gzip.open(proteome_file, 'rt') as h:
            for seqrec in SeqIO.parse(h, "fasta"):
                uniprot = seqrec.id.split("|")[1]
                if uniprot not in id_map:
                    continue
                hgnc = id_map[uniprot]
                if hgnc not in proteome:
                    proteome[hgnc] = [(uniprot, seqrec.seq)]
                else:
                    proteome[hgnc].append((uniprot, seqrec.seq))
    return proteome


def get_median_disopred(uniprot, start, end, diso_db):
    disos = []
    for i in range(start+1, end+2):
        db_key = "-".join([uniprot, str(i)])
        diso = diso_db.get(db_key)
        if diso is not None:
            disos.append(float(diso.decode()))
    if not disos:
        return None
    return median(disos)


def get_sty_motif_seqs(uniprot, seq, disordered=False):
    if disordered:
        diso_db = gdbm.open("data/external/disopred.db", 'ru')
    seq_mod = "_"*7 + seq + "_"*7
    motif_seqs_all = list(re.finditer(r"(?=(.{7}[STY].{7}))", str(seq_mod)))
    motif_seqs = []
    for motif_match in motif_seqs_all:
        motif_str = motif_match.group(1)
        # if motif_str.count("_") > int(0.15*len(motif_str)):
        #     continue
        if disordered:
            # +7 for the leading "_" that we added, +7 to reach the
            # middle of the motif
            median_diso = get_median_disopred(uniprot, motif_match.start()+7,
                                              motif_match.end()+7, diso_db)
            if median_diso is None or median_diso < 0.5:
                continue
        motif_seqs.append(motif_str)
    if disordered:
        diso_db.close()
    return motif_seqs


def read_kin_sub_data(kin_sub_file):
    # Reads kinase-substrate table where substrate and kinase are
    # found in

    ktable = {}
    with open(kin_sub_file, "r") as f:
        columns = f.readline().split("\t")
        seq_col = columns.index("SITE_...7_AA")
        sub_col = columns.index("SUB_GENE")
        for line in f:
            line_spl = line.split("\t")
            kinase = line_spl[0]
            if kinase not in ktable:
                ktable[kinase] = []
            motif_seq = line_spl[seq_col].upper()
            substrate = line_spl[sub_col]
            if kinase == substrate:
                continue
            ktable[kinase].append((substrate, motif_seq))
    return ktable


def read_aa_freqs(aa_freqs_file):
    aa_freqs = {}
    with open(aa_freqs_file, "r") as h:
        for line in h:
            aa, freq = line.strip().split()
            aa_freqs[aa] = float(freq)
    return aa_freqs


def read_all_motif_seqs(full_kin_sub_file):
    # This function extracts peptides surrounding phospho sites found
    # in the phosphosite plus human_kinase_table

    with open(full_kin_sub_file, "r") as f:
        columns = f.readline().split("\t")
        seq_col = columns.index("SITE_...7_AA")
        sub_col = columns.index("SUB_GENE")
        motif_seqs = []
        for line in f:
            line_spl = line.split("\t")
            motif_seq = line_spl[seq_col].upper()
            substrate = line_spl[sub_col]
            motif_seqs.append((substrate, motif_seq))
    return motif_seqs


def weight_motif_seqs(motif_seqs):
    # Weight motif sequences according to the position-based method of
    # Henikoff & Henikoff J. Mol. Biol. 1994

    seq_lens = set([len(s) for s in motif_seqs])
    if len(seq_lens) > 1:
        raise ValueError("Inconsistent motif lengths")
    seq_len = list(seq_lens)[0]
    res_counts = np.zeros((20, seq_len))
    res_pres = np.zeros((20, seq_len))
    # In each column, count the occurrences of each residue and also
    # keep track of the number of unique residues
    for seq in motif_seqs:
        for i, aa in enumerate(seq):
            if aa not in AMINO_ACIDS:
                continue
            res_counts[AMINO_ACIDS[aa], i] = (
                res_counts[AMINO_ACIDS[aa], i] + 1.0)
            res_pres[AMINO_ACIDS[aa], i] = 1.0
    pos_res_counts = np.sum(res_pres, axis=0)
    res_weights = np.zeros((20, seq_len))
    # Calculate residue weights
    for seq in motif_seqs:
        for i, aa in enumerate(seq):
            if aa not in AMINO_ACIDS:
                continue
            res_weights[AMINO_ACIDS[aa], i] = (
                1.0/(pos_res_counts[i]*res_counts[AMINO_ACIDS[aa], i]))
    seq_weights = []
    # Calculate sequence weights from residue weights
    for seq in motif_seqs:
        weight = 0.0
        for i, aa in enumerate(seq):
            if aa not in AMINO_ACIDS:
                continue
            weight = weight + res_weights[AMINO_ACIDS[aa], i]
        seq_weights.append(weight)
    # Normalize
    norm_seq_weights = [weight/sum(seq_weights) for weight in seq_weights]
    return norm_seq_weights


def calc_num_pseudocounts(motif_seqs, m):
    # Calculate the number of pseudocounts according to the
    # position-based method of Henikoff & Henikoff CABIOS 1996
    # m is a tunable parameter

    seq_lens = set([len(s) for s in motif_seqs])
    if len(seq_lens) > 1:
        raise ValueError("Inconsistent motif lengths")
    seq_len = list(seq_lens)[0]
    # Determine the number of unique residues present at each column.
    res_pres = np.zeros((20, seq_len))
    for seq in motif_seqs:
        for i, aa in enumerate(seq):
            if aa not in AMINO_ACIDS:
                continue
            res_pres[AMINO_ACIDS[aa], i] = 1.0
    pos_res_counts = np.sum(res_pres, axis=0)
    # N = m*R_c.  Just in case no residues were present (rare case of
    # all gaps), add one pseudocount.
    num_pc = [m*max(pos_res_count, 1) for pos_res_count in pos_res_counts]
    return num_pc


def calc_pssm(motif_seqs, proteome_aa_freqs):
    # returns a pssm matrix made from a list of known substrate
    # sequences for a kinase

    if len(motif_seqs) < 5:
        return None
    seq_lens = set([len(s) for s in motif_seqs])
    if len(seq_lens) > 1:
        raise ValueError("Inconsistent motif lengths")
    seq_len = list(seq_lens)[0]
    aa_counts = np.zeros((20, seq_len))
    seq_weights = weight_motif_seqs(motif_seqs)
    # Count residue occurrences using the sequence weights
    for n, seq in enumerate(motif_seqs):
        seq = seq.upper()
        seq_weight = seq_weights[n]
        for i, aa in enumerate(seq):
            if aa not in AMINO_ACIDS:
                continue
            aa_counts[AMINO_ACIDS[aa], i] = (
                aa_counts[AMINO_ACIDS[aa], i] + seq_weight*len(motif_seqs))
    num_nc = np.sum(aa_counts, axis=0)
    # Add in pseudocounts and calculate residue frequencies
    motif_aa_freqs = np.zeros((20, seq_len))
    # The second argument, m, is a tunable parameter
    num_pc = calc_num_pseudocounts(motif_seqs, 1)
    for i in range(seq_len):
        for aa in AMINO_ACIDS:
            # Use the background AA frequency to calculate
            # pseudocounts (see Henikoff & Henikoff 1996, eq. 7)
            pc = proteome_aa_freqs[aa] * num_pc[i]
            # Weighted average of normal counts and pseudocounts
            # (see Henikoff & Henikoff 1996, eq. 5)
            motif_aa_freqs[AMINO_ACIDS[aa], i] = (
                (aa_counts[AMINO_ACIDS[aa], i]+pc)/(num_nc[i]+num_pc[i]))
    # Calculate the final PSSM as log2 fold-change versus the
    # background proteome residue frequencies.
    pssm = np.zeros((20, seq_len))
    for i in range(seq_len):
        for aa in AMINO_ACIDS:
            pssm[AMINO_ACIDS[aa], i] = (
                log(motif_aa_freqs[AMINO_ACIDS[aa], i]/proteome_aa_freqs[aa])
                / log(2.0))
    return pssm


def is_tyr_kinase(pssm):
    return(pssm[AMINO_ACIDS['Y'], 7] > pssm[AMINO_ACIDS['S'], 7]
           and pssm[AMINO_ACIDS['Y'], 7] > pssm[AMINO_ACIDS['T'], 7])


def score_motif_seq(motif_seq, pssm):
    # scores peptites with pssm and normalizes the score by the
    # minimum and maximum scores possible with that PSSM.  Since we're
    # only scoring ST or Y sites, we just ignore the middle column.

    pssm_min = np.min(pssm, axis=0)
    pssm_max = np.max(pssm, axis=0)
    min_score = 0
    max_score = 0
    cur_score = 0
    for i, aa in enumerate(motif_seq):
        if i == 7 or aa not in AMINO_ACIDS:
            continue
        min_score = min_score + pssm_min[i]
        max_score = max_score + pssm_max[i]
        cur_score = cur_score + pssm[AMINO_ACIDS[aa], i]
    score = (cur_score - min_score) / (max_score - min_score)
    return score
