import os.path
import sys
import itertools
import pssms
from math import log2, log10, log, exp, inf, isinf
import random
from statistics import mean

def read_phosphosites(psites_file):
    # returns dictionary with every known phosphosite found on each
    # kinase in the Ochoa et al. dataset

    psites = {}
    with open(psites_file) as f:
        for line in f:
            protein, pos, res, motif_seq = line.strip().split("\t")
            if protein not in psites:
                psites[protein] = [(pos, motif_seq, res)]
            else:
                psites[protein].append((pos, motif_seq, res))
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
    from statistics import median
    phosfun = dict()
    scores = []
    with open(phosfun_file) as h:
        for line in h:
            prot, pos, score = line.strip().split()
            if prot not in phosfun:
                phosfun[prot] = {pos: float(score)}
            else:
                phosfun[prot][pos] = float(score)
            scores.append(float(score))
    return (phosfun, median(scores))


def read_reg_sites_file(reg_sites_file):
    reg_sites = dict()
    with open(reg_sites_file) as h:
        header = h.readline().split("\t")
        geneind = header.index("GENE")
        posind = header.index("MOD_RSD")
        for line in h:
            line_spl = line.strip().split("\t")
            protein = line_spl[0]
            pos = line_spl[posind].replace("-p", "")[1:]
            if protein not in reg_sites:
                reg_sites[protein] = [pos]
            else:
                reg_sites[protein].append(pos)
    return reg_sites


def read_pfam_file(pfam_file, kinases):
    kin_types = dict()
    with open(pfam_file) as h:
        for line in h:
            prot, dom, start, end = line.strip().split()
            if (prot not in kinases or dom not in
                ["Pkinase", "Pkinase_Tyr", "PI3_PI4_kinase", "Alpha_kinase"]):
                continue
            if dom in ["Pkinase", "PI3_PI4_kinase", "Alpha_kinase"]:
                kin_types[prot] = "ST"
            else:
                kin_types[prot] = "Y"
    return kin_types


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
                or len(ktable[kin_a]) < 10
                or kin_b not in psites):
            # if kin_a not in ktable:
            #     print("Skipping because no substrates for {0}".format(kin_a))
            # if kin_b not in psites:
            #     print("Skipping because no sites for {0}".format(kin_b))
            scores[(kin_a, kin_b)] = []
            continue

        # # Don't use any of the substrate sequences in calculating the
        # # PSSM to avoid circularity
        kin_motif_seqs = []
        for substrate, pos, res, seq in ktable[kin_a]:
            if substrate == kin_b:
                continue
            if substrate not in psites:
                continue
            for pos_b, seq_b, res_b in psites[substrate]:
                if pos == pos_b:
                    kin_motif_seqs.append(seq_b)
                    break
        # kin_motif_seqs = [seq for substrate, pos, res, seq in ktable[kin_a]
        #                   if substrate != kin_b and
        #                   substrate in psites and
        #                   (pos, seq, res) in psites[substrate]]
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
                # print("Skipping because no PSSM calculated for {0} - {1} motifs".format(kin_a, len(kin_motif_seqs)))
                scores[(kin_a, kin_b)] = []
                continue
            if len(kin_motif_seqs) == len(ktable[kin_a]):
                full_pssms[kin_a] = pssm

        # Take only tyrosine or serine/threonine motif sequences from
        # kin_b, depending on the specificity of kin_a.
        kin_a_is_tyr_kin = pssms.is_tyr_kinase(pssm)
        sub_motif_seqs = [(pos, seq, res) for pos, seq, res in psites[kin_b]
                          if ((kin_a_is_tyr_kin and res == 'Y') or
                              (not kin_a_is_tyr_kin and res in ['S', 'T']))]
        if not sub_motif_seqs:
            scores[(kin_a, kin_b)] = []
            continue

        for pos, motif_seq, res in sub_motif_seqs:
            score = pssms.score_motif_seq(motif_seq, pssm)
            if (kin_a, kin_b) in scores:
                scores[(kin_a, kin_b)].append((pos, res, score))
            else:
                scores[(kin_a, kin_b)] = [(pos, res, score)]
    return scores


def normalize_scores(scores):
    scores_minmax = dict()
    for pair in scores:
        kin = pair[0]
        if kin not in scores_minmax:
            scores_minmax[kin] = [inf, -inf]
        if not scores[pair]:
            continue
        for pos, res, score in scores[pair]:
            if score < scores_minmax[kin][0]:
                scores_minmax[kin][0] = score
            if score > scores_minmax[kin][1]:
                scores_minmax[kin][1] = score
    norm_scores = dict()
    for pair in scores:
        if not scores[pair]:
            continue
        kin = pair[0]
        kin_min = scores_minmax[kin][0]
        kin_max = scores_minmax[kin][1]
        if isinf(kin_min) or isinf(kin_max):
            continue
        pair_norm_scores = []
        for pos, res, score in scores[pair]:
            if kin_min == kin_max:
                norm_score = 1.0
            else:
                norm_score = (score-kin_min)/(kin_max-kin_min)
            pair_norm_scores.append((pos, res, norm_score))
        norm_scores[pair] = pair_norm_scores
    return norm_scores


def calc_dcg(scores):
    return sum([score/log2(n+2) for n, score in enumerate(scores)])


def regulation_score(pair, scores, phosfun, med_phosfun, reg_sites):
    # Rank the PSSM scores
    scores.sort(key=lambda t: t[2], reverse=True)
    kin_b_reg_sites = reg_sites.get(pair[1])
    kin_b_phosfun = phosfun.get(pair[1])
    hs_scores = []
    for pos, res, score in scores:
        if kin_b_phosfun is None or pos not in kin_b_phosfun:
            continue
        hs_score = kin_b_phosfun[pos]
        hs_scores.append(hs_score)
    if not hs_scores:
        return ("NA", "NA")
    if len(hs_scores) < 2:
        return ("NA", max(hs_scores))
    # Discounted Cumulative Gain
    dcg = calc_dcg(hs_scores)
    # Ideal Discounted Cumulative Gain
    hs_scores.sort(reverse=True)
    max_dcg = calc_dcg(hs_scores)
    hs_scores.sort(reverse=False)
    min_dcg = calc_dcg(hs_scores)
    if max_dcg == min_dcg:
        return ("NA", max(hs_scores))
    return ((dcg-min_dcg)/(max_dcg-min_dcg), max(hs_scores))


def score_network(kinase_file, out_file):
    # This takes network files and scores all edges according to pssm

    ktable = pssms.read_kin_sub_data("data/psiteplus-kinase-substrates.tsv")
    aa_freqs = pssms.read_aa_freqs("data/aa-freqs.tsv")
    psites = read_phosphosites("data/psiteplus-phosphosites.tsv")
    kinases = read_kinase_file(kinase_file)
    # phosfun, med_phosfun = read_phosfun_file("data/phosfun-alt.tsv")
    phosfun_st, med_phosfun_st = read_phosfun_file("data/phosfun-ST.tsv")
    phosfun_y, med_phosfun_y = read_phosfun_file("data/phosfun-Y.tsv")
    reg_sites = read_reg_sites_file("data/psiteplus-reg-sites.tsv")
    kin_types = read_pfam_file("data/human-pfam.tsv", kinases)
    scores = score_kin_pairs(ktable, kinases, aa_freqs, psites)
    norm_scores = normalize_scores(scores)
    # scores = norm_scores
    with open(out_file, "w") as v:
        v.write("\t".join(["node1", "node2", "kinase.type", "sub.kinase.type",
                           "max.pssm.score", "dcg", "max.func.score"])+"\n")
        for pair in scores:
            if pair[0] not in kin_types:
                if pair[0] == "MTOR":
                    kin_type = "ST"
                else:
                    kin_type = "NA"
            else:
                kin_type = kin_types[pair[0]]
            if pair[1] not in kin_types:
                if pair[1] == "MTOR":
                    sub_type = "ST"
                else:
                    sub_type = "NA"
            else:
                sub_type = kin_types[pair[1]]
            if not scores[pair]:
                v.write("\t".join([pair[0], pair[1], kin_type, sub_type,
                                   "NA", "NA", "NA"])+"\n")
                continue
            residues = set([res for pos, res, score in scores[pair]])
            max_score = max([score for pos, res, score in scores[pair]])
            if residues == set(['Y']):
                dcg, max_func_score = regulation_score(
                    pair, scores[pair], phosfun_y, med_phosfun_y, reg_sites)
            else:
                dcg, max_func_score  = regulation_score(
                    pair, scores[pair], phosfun_st, med_phosfun_st, reg_sites)
            v.write("\t".join([pair[0], pair[1], kin_type, sub_type,
                               str(max_score), str(dcg), str(max_func_score)])+"\n")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("USAGE: <script> DATA_FILE OUT_FILE")
    kin_act_file_arg = sys.argv[1]
    out_file = sys.argv[2]
    score_network(kin_act_file_arg, out_file)
