import os.path
import sys
import itertools
import pssms


def read_phosphosites(psites_file):
    #returns dictionary with every known phosphosite found on each kinase in the
    #Ochoa et al. dataset
    
    psites = {}
    with open(psites_file) as f:
        columns = f.readline()
        columns =columns.split("\t")
        seqind = columns.index("SITE_...7_AA")
        for line in f:
            line_spl = line.split("\t")
            protein = line_spl[0]
            motif_seq = line_spl[seqind].upper()
            if protein not in psites:
                psites[protein] = [motif_seq]
            else:
                psites[protein].append(motif_seq)
    return psites


def score_kin_pairs(ktable, kin_act_file, aa_freqs, psites):
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
    scores = {}
    for pair in pairs:
        A, B = pair
        if A == B:
            continue
        if A not in ktable:
            scores[(A, B)] = []
            continue
        if B not in psites:
            scores[(A, B)] = []
            continue
        # Don't use any of the substrate sequences in calculating the
        # PSSM to avoid circularity
        kin_motif_seqs = [seq for (substrate, seq) in ktable[A]
                          if substrate != B]
        pssm = pssms.calc_pssm(kin_motif_seqs, aa_freqs)
        if pssm is None:
            scores[(A, B)] = []
            continue
        sub_motif_seqs = psites[B]
        for motif_seq in sub_motif_seqs:
            score = pssms.score_motif_seq(motif_seq, pssm)
            if (A, B) in scores:
                scores[(A, B)].append(score)
            else:
                scores[(A, B)] = [score]
    return scores


def score_network(kin_act_file):
    #This takes network files and scores all edges according to pssm

    ktable = pssms.read_kin_sub_data("data/reduced_kinase_table.tsv")
    aa_freqs = pssms.read_aa_freqs("data/aa-freqs.tsv")
    psites = read_phosphosites("data/phosphosites_reduced.tsv")
    scores = score_kin_pairs(ktable, kin_act_file, aa_freqs, psites)
    out_file_base = os.path.splitext(os.path.basename(kin_act_file))[0]
    out_file = "out/" + out_file_base + "-pssm.tsv"
    with open(out_file, "w") as v:
        v.write("\t".join(["node1", "node2", "pssm.score"])+"\n")
        for pair in scores:
            if not scores[pair]:
                v.write("\t".join([pair[0], pair[1], "NA"])+"\n")
                continue
            max_score = max(scores[pair])
            v.write("\t".join([pair[0], pair[1], str(max_score)])+"\n")
    
    
if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("USAGE: <script> DATA_FILE")
    kin_act_file = sys.argv[1]
    score_network(kin_act_file)
