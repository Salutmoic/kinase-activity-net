import sys
import os.path
import random
import numpy as np
import pssms


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


def calc_dist(motif_seqs, pssm):
    #creates distribution of pssm scores
    #scores are calculated by scoring a random phosphosite relative to to kinase i

    if pssm is None:
        return []
    dist = []
    for j in range(1000):
        motif_seq = random.choice(motif_seqs)
        dist.append(pssms.score_motif_seq(motif_seq, pssm))
    return dist
        

def calc_distr(kin_act_file):
    ktable = pssms.read_kin_sub_data("data/reduced_kinase_table.tsv")
    aa_freqs = pssms.read_aa_freqs("data/aa-freqs.tsv")
    all_motif_seqs = pssms.read_all_motif_seqs("data/human_kinase_table.tsv")
    out_file_base = os.path.splitext(os.path.basename(kin_act_file))[0]
    dist_out_file = "out/" + out_file_base + "-pssm-dists.tsv"
    with open(dist_out_file, "w") as v:
        for kinase in ktable:
            kin_motif_seqs = ktable[kinase]
            motif_seqs = [seq for (substrate, seq) in kin_motif_seqs]
            pssm = pssms.calc_pssm(motif_seqs, aa_freqs)
            sub_motif_seqs = [seq for (substrate, seq) in all_motif_seqs
                              if (substrate, seq) not in kin_motif_seqs]
            dist = calc_dist(sub_motif_seqs, pssm)
            line = [kinase]
            line.extend([str(x) for x in dist])
            v.write("\t".join(line) + "\n")     


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("USAGE: <script> DATA_FILE")
    kin_act_file = sys.argv[1]
    calc_distr(kin_act_file)
