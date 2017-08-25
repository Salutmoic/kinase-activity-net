import sys
import os.path
import random
import numpy as np
import pssms


def get_random_sty_motifs(proteome, omit, n):
    rand_motif_seqs = []
    for _ in range(n):
        substrate = random.choice(list(proteome.keys()))
        seqs = proteome[substrate]
        sub_motif_seqs = set()
        for uniprot, isoform in seqs:
            sub_motif_seqs.update(pssms.get_sty_motif_seqs(uniprot, isoform))
        if sub_motif_seqs:
            motif_seq = random.choice(list(sub_motif_seqs))
        else:
            motif_seq = None
        while (motif_seq is None
               or 'X' in motif_seq
               or 'U' in motif_seq
               or (substrate, motif_seq) in omit):
            substrate = random.choice(list(proteome.keys()))
            seqs = proteome[substrate]
            sub_motif_seqs = set()
            for uniprot, isoform in seqs:
                sub_motif_seqs.update(
                    pssms.get_sty_motif_seqs(uniprot, isoform))
            if sub_motif_seqs:
                motif_seq = random.choice(list(sub_motif_seqs))
            else:
                motif_seq = None
        rand_motif_seqs.append(motif_seq)
    return rand_motif_seqs


def calc_dist(motif_seqs, pssm):
    #creates distribution of pssm scores scores are calculated by
    #scoring a random phosphosite relative to to kinase i

    if pssm is None:
        return []
    dist = []
    for motif_seq in motif_seqs:
        dist.append(pssms.score_motif_seq(motif_seq, pssm))
    return dist


def calc_distr(kin_act_file):
    ktable = pssms.read_kin_sub_data("data/reduced_kinase_table.tsv")
    aa_freqs = pssms.read_aa_freqs("data/aa-freqs.tsv")
    proteome = pssms.read_proteome()
    out_file_base = os.path.splitext(os.path.basename(kin_act_file))[0]
    dist_out_file = "out/" + out_file_base + "-pssm-dists.tsv"
    with open(dist_out_file, "w") as v:
        for kinase in ktable:
            kin_motif_seqs = ktable[kinase]
            motif_seqs = [seq for (substrate, seq) in kin_motif_seqs]
            pssm = pssms.calc_pssm(motif_seqs, aa_freqs)
            sub_motif_seqs = get_random_sty_motifs(proteome, kin_motif_seqs,
                                                   1000)
            dist = calc_dist(sub_motif_seqs, pssm)
            if not dist:
                continue
            line = [kinase]
            line.extend([str(x) for x in dist])
            v.write("\t".join(line) + "\n")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("USAGE: <script> DATA_FILE")
    kin_act_file_arg = sys.argv[1]
    calc_distr(kin_act_file_arg)
