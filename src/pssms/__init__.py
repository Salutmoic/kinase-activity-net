import numpy as np
from math import log


AMINO_ACIDS = {"G": 0, "A": 1, "V": 2, "L": 3, "I": 4, "P": 5, "F": 6, "Y": 7,
              "W": 8, "S": 9, "T": 10, "C": 11, "M": 12, "N": 13, "Q": 14,
              "K": 15, "R": 16, "H": 17, "D": 18, "E": 19}


def read_kin_sub_data(kin_sub_file):
    #Reads kinase-substrate table where substrate and kinase are found in 

    ktable = {}
    with open(kin_sub_file, "r") as f:
        columns = f.readline().split("\t")
        seq_col = columns.index("SITE_...7_AA")
        sub_col = columns.index("SUB_GENE")
        dom_col = columns.index("DOMAIN")
        for line in f:
            line_spl = line.split("\t")
            kinase = line_spl[0]
            if kinase not in ktable:
                ktable[kinase] = []
            motif_seq = line_spl[seq_col].upper()
            substrate = line_spl[sub_col]
            # Really we want to omit autophosphorylation sites that
            # fall in the kinases' activation loops, but finding the
            # activation loops requires a bit of work so this is just
            # a quick cheat for now:
            # domain = line_spl[dom_col]
            # if substrate == kinase and domain in ["Pkinase", "Pkinase_Tyr"]:
            #     continue
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

    with open(full_kin_sub_file,"r") as f:
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

    res_counts = np.zeros((20, 15))
    res_pres = np.zeros((20, 15))
    for seq in motif_seqs:
        for i, aa in enumerate(seq):
            if aa not in AMINO_ACIDS:
                continue
            res_counts[AMINO_ACIDS[aa]][i] = res_counts[AMINO_ACIDS[aa]][i] + 1.0
            res_pres[AMINO_ACIDS[aa]][i] = 1.0
    pos_res_counts = np.sum(res_pres, axis=0)
    res_weights = np.zeros((20, 15))
    for seq in motif_seqs:
        for i, aa in enumerate(seq):
            if aa not in AMINO_ACIDS:
                continue
            res_weights[AMINO_ACIDS[aa]][i] = 1.0/(pos_res_counts[i]*res_counts[AMINO_ACIDS[aa]][i])
    seq_weights = []
    for seq in motif_seqs:
        weight = 0.0
        for i, aa in enumerate(seq):
            if aa not in AMINO_ACIDS:
                continue
            weight = weight + res_weights[AMINO_ACIDS[aa]][i]
        seq_weights.append(weight)
    norm_seq_weights = [weight/sum(seq_weights) for weight in seq_weights]
    return norm_seq_weights
    

def calc_num_pseudocounts(motif_seqs, m):
    # Calculate the number of pseudocounts according to the
    # position-based method of Henikoff & Henikoff CABIOS 1996
    # m is a tunable parameter

    res_counts = np.zeros((20, 15))
    res_pres = np.zeros((20, 15))
    for seq in motif_seqs:
        for i, aa in enumerate(seq):
            if aa not in AMINO_ACIDS:
                continue
            res_pres[AMINO_ACIDS[aa]][i] = 1.0
    pos_res_counts = np.sum(res_pres, axis=0)
    num_pc = m * pos_res_counts
    return num_pc


def calc_pssm(motif_seqs, aa_freqs):
    # returns a pssm matrix made from a list of known substrate
    # sequences for a kinase

    pssms = {}
    min_seqs = 10
    if len(motif_seqs) < min_seqs:
        return None
    cmat = np.zeros((20, 15))
    seq_weights = weight_motif_seqs(motif_seqs)
    for n, seq in enumerate(motif_seqs):
        seq = seq.upper()
        seq_weight = seq_weights[n]
        for i, aa in enumerate(seq):
            if aa not in AMINO_ACIDS:
                continue
            cmat[AMINO_ACIDS[aa]][i] = cmat[AMINO_ACIDS[aa]][i] + seq_weight*len(motif_seqs)

    # position-based pseudocounts 
    num_nc= len(motif_seqs)
    # The second argument, m, is a tunable parameter
    num_pc = calc_num_pseudocounts(motif_seqs, 8)
    for i in range(15):
        for aa in aa_freqs:
            # Use the background AA frequency to calculate
            # pseudocounts (see Henikoff & Henikoff 1996, eq. 7)
            pc = aa_freqs[aa] * num_pc[i]
            # Weighted average of normal counts and pseudocounts
            # (see Henikoff & Henikoff 1996, eq. 5)
            cmat[AMINO_ACIDS[aa]][i] = (cmat[AMINO_ACIDS[aa]][i]+pc)/(num_nc+num_pc[i])

    for i in range(15):
        for aa in aa_freqs:
            cmat[AMINO_ACIDS[aa]][i] = log(cmat[AMINO_ACIDS[aa]][i]/aa_freqs[aa])/log(2.0)
    return cmat


def score_motif_seq(motif_seq, pssm):
    #scores peptites with pssm

    score = 1
    for i, aa in enumerate(motif_seq):
        if aa != "_":
            score = score + pssm[AMINO_ACIDS[aa]][i]
        else:
            score = score + np.min(pssm, axis=0)[i]
    return score
