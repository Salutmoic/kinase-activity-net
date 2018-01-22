import os.path
import sys
import itertools
import pssms
from math import log2, log, exp


def read_phosphosites(psites_file):
    # returns dictionary with every known phosphosite found on each
    # kinase in the Ochoa et al. dataset

    psites = {}
    with open(psites_file) as f:
        columns = f.readline()
        columns = columns.split("\t")
        seqind = columns.index("SITE_...7_AA")
        posind = columns.index("MOD_RSD")
        for line in f:
            line_spl = line.split("\t")
            protein = line_spl[0]
            motif_seq = line_spl[seqind].upper()
            pos = line_spl[posind][1:-2]
            if protein not in psites:
                psites[protein] = [(pos, motif_seq)]
            else:
                psites[protein].append((pos, motif_seq))
    return psites


def read_distrs(pssm_distr_file):
    distrs = dict()
    with open(pssm_distr_file) as h:
        for line in h:
            line_spl = line.strip().split()
            kinase = line_spl[0]
            distr = [float(x) for x in line_spl[1:]]
            distrs[kinase] = distr
    return distrs


def calc_zscore(pssm_score, distr_mean, distr_stdev):
    z = (pssm_score-distr_mean)/distr_stdev
    return z


def read_kinase_file(kinase_file):
    kinases = set()
    with open(kinase_file) as h:
        for line in h:
            kinases.add(line.strip())
    return(kinases)


def read_phosfun_file(phosfun_file):
    phosfun = dict()
    with open(phosfun_file) as h:
        for line in h:
            prot, pos, score = line.strip().split()
            if prot not in phosfun:
                phosfun[prot] = {pos: float(score)}
            else:
                phosfun[prot][pos] = float(score)
    return phosfun


def read_reg_sites_file(reg_sites_file):
    reg_sites = dict()
    with open(reg_sites_file) as h:
        header = h.readline().split("\t")
        geneind = header.index("GENE")
        posind = header.index("MOD_RSD")
        for line in h:
            line_spl = line.strip().split("\t")
            protein = line_spl[0]
            pos = line_spl[posind][1:]
            if protein not in reg_sites:
                reg_sites[protein] = [pos]
            else:
                reg_sites[protein].append(pos)
    return reg_sites


def score_kin_pairs(ktable, kinases, aa_freqs, psites):
    # loads a netwoek of kinase-> kinase interactions and scores them
    # with the pssm it returns a dictionary with kinasa A and kinase B
    # and list of scores for each phosphosite found on B

    # Just get interaction pairs for kinases that are in our kinase
    # activity table
    pairs = itertools.product(kinases, repeat=2)
    scores = {}
    full_pssms = {}
    for kin_a, kin_b in pairs:
        if (kin_a == kin_b
                or kin_a not in ktable
                or kin_b not in psites):
            scores[(kin_a, kin_b)] = []
            continue

        # # Don't use any of the substrate sequences in calculating the
        # # PSSM to avoid circularity
        kin_motif_seqs = [seq for (substrate, seq) in ktable[kin_a]
                          if substrate != kin_b]
        # If kin_b isn't a known substrate, check if we've already
        # calculated the full PSSM for kin_a and, if so, take it
        # rather than wasting time recalculating it.
        if len(kin_motif_seqs) == len(ktable[kin_a]) and kin_a in full_pssms:
            pssm = full_pssms[kin_a]
        else:
            # It's either the first time calculating a PSSM for kin_a
            # or we need to recalculate it for a subset of its known
            # substrates.
            pssm = pssms.calc_pssm(kin_motif_seqs, aa_freqs)
            if pssm is None:
                scores[(kin_a, kin_b)] = []
                continue
            if len(kin_motif_seqs) == len(ktable[kin_a]):
                full_pssms[kin_a] = pssm

        # Take only tyrosine or serine/threonine motif sequences from
        # kin_b, depending on the specificity of kin_a.
        kin_a_is_tyr_kin = pssms.is_tyr_kinase(pssm)
        sub_motif_seqs = [(pos, seq) for pos, seq in psites[kin_b]
                          if ((kin_a_is_tyr_kin and seq[7] == 'Y') or
                              (not kin_a_is_tyr_kin and seq[7] in ['S', 'T']))]
        if not sub_motif_seqs:
            scores[(kin_a, kin_b)] = []
            continue

        for pos, motif_seq in sub_motif_seqs:
            score = pssms.score_motif_seq(motif_seq, pssm)
            if (kin_a, kin_b) in scores:
                scores[(kin_a, kin_b)].append((pos, score))
            else:
                scores[(kin_a, kin_b)] = [(pos, score)]
    return scores


def calc_dcg(scores):
    return sum([score/log2(n+2) for n, score in enumerate(scores)])


def regulation_score(pair, scores, phosfun, reg_sites):
    # Rank the PSSM scores
    scores.sort(key=lambda t: t[1], reverse=True)
    kin_b_reg_sites = reg_sites.get(pair[1])
    kin_b_phosfun = phosfun.get(pair[1])
    hs_scores = []
    for pos, score in scores:
        if kin_b_reg_sites and pos in kin_b_reg_sites:
            hs_score = 1.0
        elif kin_b_phosfun is None or pos not in kin_b_phosfun:
            hs_score = 0.001
        else:
            hs_score = kin_b_phosfun[pos]
        hs_scores.append(hs_score)
    # Discounted Cumulative Gain
    dcg = calc_dcg(hs_scores)
    # Ideal Discounted Cumulative Gain
    hs_scores.sort(reverse=True)
    idcg = calc_dcg(hs_scores)
    return (dcg/idcg)*max(hs_scores)


def score_network(kinase_file, out_file):
    # This takes network files and scores all edges according to pssm

    ktable = pssms.read_kin_sub_data("data/reduced_kinase_table.tsv")
    aa_freqs = pssms.read_aa_freqs("data/aa-freqs.tsv")
    psites = read_phosphosites("data/phosphosites_reduced.tsv")
    kinases = read_kinase_file(kinase_file)
    phosfun = read_phosfun_file("data/phosfun.tsv")
    reg_sites = read_reg_sites_file("data/reg_sites.tsv")
    scores = score_kin_pairs(ktable, kinases, aa_freqs, psites)
    with open(out_file, "w") as v:
        v.write("\t".join(["node1", "node2", "pssm.score"])+"\n")
        for pair in scores:
            if not scores[pair]:
                v.write("\t".join([pair[0], pair[1], "NA"])+"\n")
                continue
            max_score = max([score for pos, score in scores[pair]])
            dcg = regulation_score(pair, scores[pair], phosfun, reg_sites)
            v.write("\t".join([pair[0], pair[1], str(dcg*max_score)])+"\n")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("USAGE: <script> DATA_FILE OUT_FILE")
    kin_act_file_arg = sys.argv[1]
    out_file = sys.argv[2]
    score_network(kin_act_file_arg, out_file)
