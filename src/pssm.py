import numpy as np
import random
from os import listdir
from os.path import isfile, join, walk, basename, splitext
from math import log, ceil, sqrt
import sys
import itertools


def read_kin_sub_data():
    #Reads kinase-substrate table where substrate and kinase are found in 

    ktable = {}
    with open("data/reduced_kinase_table.tsv","r") as f:
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
            domain = line_spl[dom_col]
            # Really we want to omit autophosphorylation sites that
            # fall in the kinases' activation loops, but finding the
            # activation loops requires a bit of work so this is just
            # a quick cheat for now:
            # if substrate == kinase and domain in ["Pkinase", "Pkinase_Tyr"]:
            #     continue
            ktable[kinase].append(motif_seq)
    return ktable


def amino_acids():
    #returns a dictionary containing a list of the common 20 amino acids 

    aminos = {"G": 0, "A": 1, "V": 2, "L": 3, "I": 4, "P": 5, "F": 6, "Y": 7,
              "W": 8, "S": 9, "T": 10, "C": 11, "M": 12, "N": 13, "Q": 14,
              "K": 15, "R": 16, "H": 17, "D": 18, "E": 19}
    return aminos


def read_aa_freqs():
    aa_freqs = {}
    with open("data/aa-freqs.tsv", "r") as h:
        for line in h:
            aa, freq = line.strip().split()
            aa_freqs[aa] = float(freq)
    return aa_freqs


def weight_motif_seqs(motif_seqs):
    # Weight motif sequences according to the position-based method of
    # Henikoff & Henikoff J. Mol. Biol. 1994

    aminos = amino_acids()
    res_counts = np.zeros((20, 15))
    res_pres = np.zeros((20, 15))
    for seq in motif_seqs:
        for i, aa in enumerate(seq):
            if aa not in aminos:
                continue
            res_counts[aminos[aa]][i] = res_counts[aminos[aa]][i] + 1.0
            res_pres[aminos[aa]][i] = 1.0
    pos_res_counts = np.sum(res_pres, axis=0)
    res_weights = np.zeros((20, 15))
    for seq in motif_seqs:
        for i, aa in enumerate(seq):
            if aa not in aminos:
                continue
            res_weights[aminos[aa]][i] = 1.0/(pos_res_counts[i]*res_counts[aminos[aa]][i])
    seq_weights = []
    for seq in motif_seqs:
        weight = 0.0
        for i, aa in enumerate(seq):
            if aa not in aminos:
                continue
            weight = weight + res_weights[aminos[aa]][i]
        seq_weights.append(weight)
    norm_seq_weights = [weight/sum(seq_weights) for weight in seq_weights]
    return norm_seq_weights
    

def calc_num_pseudocounts(motif_seqs, m):
    # Calculate the number of pseudocounts according to the
    # position-based method of Henikoff & Henikoff CABIOS 1996

    aminos = amino_acids()
    res_counts = np.zeros((20, 15))
    res_pres = np.zeros((20, 15))
    for seq in motif_seqs:
        for i, aa in enumerate(seq):
            if aa not in aminos:
                continue
            res_pres[aminos[aa]][i] = 1.0
    pos_res_counts = np.sum(res_pres, axis=0)
    num_pc = m * pos_res_counts
    return num_pc


def calc_pssms(ktable):
    # returns a dictionary where keys are kinases and values are pmms
    # matrices made from the known substrates for each kinase

    aminos = amino_acids()
    pssms = {}
    aa_freqs = read_aa_freqs()
    min_seqs = 10
    for kinase in ktable:
        motif_seqs = ktable[kinase]
        if len(motif_seqs) < min_seqs:
            pssms[kinase] = None
            continue
        cmat = np.zeros((20, 15))
        seq_weights = weight_motif_seqs(motif_seqs)
        for n, seq in enumerate(motif_seqs):
            seq = seq.upper()
            seq_weight = seq_weights[n]
            for i, aa in enumerate(seq):
                if aa not in aminos:
                    continue
                cmat[aminos[aa]][i] = cmat[aminos[aa]][i] + seq_weight*len(motif_seqs)

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
                cmat[aminos[aa]][i] = (cmat[aminos[aa]][i]+pc)/(num_nc+num_pc[i])

        for i in range(15):
            for aa in aa_freqs:
                cmat[aminos[aa]][i] = log(cmat[aminos[aa]][i]/aa_freqs[aa])/log(2.0)
        pssms[kinase] = cmat
    return pssms


def get_all_motif_seqs():
    # This function extracts peptides surrounding phospho sites found
    # in the phosphosite plus human_kinase_table

    with open("data/human_kinase_table.tsv","r") as f:
        columns = f.readline()
        columns =columns.split("\t")
        seqind = columns.index("SITE_...7_AA")
        motif_seqs = []
        for line in f:
            line_spl = line.split("\t")
            motif_seqs.append(line_spl[seqind].upper())
    return motif_seqs


def score_motif_seq(motif_seq, pssm):
    #scores peptites with pssm

    aminos = amino_acids()
    score = 1
    for i, aa in enumerate(motif_seq):
        if aa != "_":
            score = score + pssm[aminos[aa]][i]
        else:
            score = score + np.min(pssm, axis=0)[i]
    return score
        
        
def calc_dist(motif_seqs, pssms):
    #creates distribution of pssm scores
    #scores are calculated by scoring a random phosphosite relative to to kinase i

    dist = {}
    for kinase in pssms:
        if pssms[kinase] is None:
            continue
        d = []
        for j in range(1000):
            motif_seq = random.choice(motif_seqs)
            d.append(score_motif_seq(motif_seq, pssms[kinase]))
        dist[kinase] = d
    return dist
        

def kin_phosphosites(ktable):
    #returns dictionary with every known phosphosite found on each kinase in the
    #Ochoa et al. dataset
    
    psites = {}
    with open("data/phosphosites_reduced.tsv") as f:
        columns = f.readline()
        columns =columns.split("\t")
        seqind = columns.index("SITE_...7_AA")
        for line in f:
            line_spl = line.split("\t")
            protein = line_spl[0]
            motif_seq = line_spl[seqind].upper()
            if protein not in ktable:
                continue
            if protein not in psites:
                psites[protein] = [motif_seq]
            else:
                psites[protein].append(motif_seq)
    return psites


def score_kin_pairs(ktable, pssms, kin_act_file):
    # loads a netwoek of kinase-> kinase interactions and scores them with the pssm
    # it returns a dictionary with kinasa A and kinase B and list of scores for each phosphosite found
    # on B

    prots = set()
    with open(kin_act_file) as h:
        for line in h:
            prot, _, _ = line.split()
            prots.add(prot)
    # Just get interaction pairs for kinases that are in our kinase activity table
    pairs = itertools.product(prots, repeat=2)
    psites = kin_phosphosites(ktable)
    scores = {}
    for pair in pairs:
        A, B = pair
        if A == B:
            continue
        if A not in pssms or pssms[A] is None or B not in psites:
            scores[(A, B)] = []
            continue
        motif_seqs = psites[B]
        for motif_seq in motif_seqs:
            score = score_motif_seq(motif_seq, pssms[A])
            if (A, B) in scores:
                scores[(A, B)].append(score)
            else:
                scores[(A, B)] = [score]
    return scores


def assess_edges(scores, dist):
    #Finds the highest scoreing phosphosite for ech A,B pair and compares it to the
    #distribution generated for B (A-> B) from the dist() function
    zscores = {}

    for pair in scores:
        if not scores[pair]:
            zscores[pair] = None
            continue
        max_score = max(scores[pair])
        A, B = pair
        z = (max_score-np.mean(dist[A]))/np.std(dist[A])
        zscores[pair] = z
    return zscores


def score_network(kin_act_file):
    #This takes network files and scores all edges according to pssm

    ktable = read_kin_sub_data()
    pssms = calc_pssms(ktable)
    motif_seqs = get_all_motif_seqs()
    dist = calc_dist(motif_seqs, pssms)
    out_file_base = splitext(basename(kin_act_file))[0]
    dist_out_file = "out/" + out_file_base + "-pssm-dists.tsv"
    with open(dist_out_file, "w") as v:
        for kinase in dist:
            line = [kinase]
            for j in dist[kinase]:
                line.append(str(j))
            v.write("\t".join(line) + "\n")     
    scores = score_kin_pairs(ktable, pssms, kin_act_file)
    out_file = "out/" + out_file_base + "-pssm.tsv"
    with open(out_file, "w") as v:
        v.write("\t".join(["node1", "node2", "pssm.score"])+"\n")
        for pair in scores:
            if not scores[pair]:
                v.write("\t".join([pair[0], pair[1], "NA"])+"\n")
                continue
            max_score = max(scores[pair])
            v.write("\t".join([pair[0], pair[1], str(max_score)])+"\n")
    
    
def known_scores(pssms, ktable):
    v = open("out/known_kinase_psite-score","w")
    with open("data/human_kinase_table","r") as f:
        columns = f.readline()
        columns =columns.split("\t")
        seqind = columns.index("SITE_...7_AA")
        for line in f:
            line_spl = line.split("\t")
            kinase = line_spl[0] 
            motif_seq = line_spl[seqind]
            score = score_motif_seq(motif_seq, pssms[kinase])
            v.write("\t".join([kinase, motif_seq, str(score)])+"\n")
    v.close()
    return    


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("USAGE: <script> DATA_FILE")
    kin_act_file = sys.argv[1]
    score_network(kin_act_file)
