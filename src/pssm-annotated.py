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
        posind = header.index("position")
        resind = header.index("residue")
        for line in h:
            line_spl = line.strip().split("\t")
            protein = line_spl[0]
            pos = line_spl[posind]
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


def read_kinase_families(kin_fam_file):
    kinase_fams = dict()
    family_kins = dict()
    with open(kin_fam_file) as h:
        for line in h:
            kinase, family = line.strip().split()
            kinase_fams[kinase] = family
            if family not in family_kins:
                family_kins[family] = set([kinase])
            else:
                family_kins[family].add(kinase)
    return (kinase_fams, family_kins)


def read_sdr_assignments(sdr_assign_file):
    sdr_assigns = dict()
    with open(sdr_assign_file) as h:
        h.readline()
        for line in h:
            kinase, assign, score = line.strip().split()
            if float(score) > 0.8:
                sdr_assigns[kinase] = assign
            else:
                sdr_assigns[kinase] = None
    return sdr_assigns


def assign_pssms_by_sdr(kinases, kinase_pssms, sdr_assigns):
    for kinase in kinases:
        if kinase in kinase_pssms:
            continue
        if kinase not in sdr_assigns:
            continue
        if sdr_assigns[kinase] is None:
            continue
        assign = sdr_assigns[kinase]
        assign_pssm = kinase_pssms[assign][0]
        kinase_pssms[kinase] = (assign_pssm, "sdr")
    return kinase_pssms


def score_kin_pairs(kinases, kinase_pssms, psites, kinase_fams):
    # loads a netwoek of kinase-> kinase interactions and scores them
    # with the pssm it returns a dictionary with kinasa A and kinase B
    # and list of scores for each phosphosite found on B

    # Just get interaction pairs for kinases that are in our kinase
    # activity table
    pairs = itertools.product(kinases, repeat=2)
    scores = {}
    for kin_a, kin_b in pairs:
        if kin_a == kin_b or kin_b not in psites:
            scores[(kin_a, kin_b)] = []
            continue
        if kin_a not in kinase_pssms:
            scores[(kin_a, kin_b)] = []
            continue
        pssm, _ = kinase_pssms[kin_a]
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
    kin_b_reg_sites = reg_sites.get(pair[1])
    kin_b_phosfun = phosfun.get(pair[1])
    hs_scores = []
    for pos, res, score in scores:
        if kin_b_phosfun is None or pos not in kin_b_phosfun:
            hs_score = "NA"
        else:
            hs_score = kin_b_phosfun[pos]
        hs_scores.append(hs_score)
    all_scores = ",".join([str(score) for score in hs_scores])
    hs_scores = []
    # Rank the PSSM scores
    scores.sort(key=lambda t: t[2], reverse=True)
    for pos, res, score in scores:
        if kin_b_phosfun is None or pos not in kin_b_phosfun:
            continue
        hs_score = kin_b_phosfun[pos]
        hs_scores.append(hs_score)
    if not hs_scores:
        return ("NA", "NA", "NA", "NA")
    sites = ["{0}{1}".format(res, pos) for pos, res, score in scores]
    max_score = max(hs_scores)
    max_score_i = hs_scores.index(max_score)
    if len(hs_scores) < 2:
        return ("NA", max_score, sites[max_score_i], all_scores)
    # Discounted Cumulative Gain
    dcg = calc_dcg(hs_scores)
    # Ideal Discounted Cumulative Gain
    hs_scores.sort(reverse=True)
    max_dcg = calc_dcg(hs_scores)
    hs_scores.sort(reverse=False)
    min_dcg = calc_dcg(hs_scores)
    if max_dcg == min_dcg:
        return ("NA", max(hs_scores))
    return ((dcg-min_dcg)/(max_dcg-min_dcg), max_score, sites[max_score_i],
            all_scores)


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
    kinase_fams, family_kins = read_kinase_families("data/kinase-families.tsv")
    sdr_assigns = read_sdr_assignments("data/sdr-pssm-assignments.tsv")
    family_pssms = pssms.build_family_pssms(family_kins, ktable, aa_freqs)
    kinase_pssms = pssms.build_pssms(kinases, kinase_fams, family_pssms, ktable,
                                     aa_freqs)
    kinase_pssms = assign_pssms_by_sdr(kinases, kinase_pssms, sdr_assigns)
    scores = score_kin_pairs(kinases, kinase_pssms, psites, family_pssms)
    norm_scores = normalize_scores(scores)
    # scores = norm_scores
    kinases_sorted = list(kinases)
    kinases_sorted.sort()
    pairs = itertools.product(kinases_sorted, repeat=2)
    with open(out_file, "w") as v:
        v.write("\t".join(["node1", "node2", "kinase.type", "sub.kinase.type",
                           "pssm.source", "max.pssm.score", "dcg",
                           "max.func.score", "max.pssm.score.res",
                           "max.func.score.res", "all.pssm.scores",
                           "all.func.scores", "all.sites"])+"\n")
        for pair in pairs:
            if pair[0] == pair[1]:
                continue
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
            if pair[0] in kinase_pssms:
                source = kinase_pssms[pair[0]][1]
            else:
                source = "NA"
            if not scores[pair]:
                v.write("\t".join([pair[0], pair[1], kin_type, sub_type, source,
                                   "NA", "NA", "NA", "NA", "NA", "NA", "NA",
                                   "NA"])+"\n")
                continue
            own_pssm = str(kinase_pssms[pair[0]][1]).upper()
            residues = set([res for pos, res, score in scores[pair]])
            pssm_scores = ([score for pos, res, score in scores[pair]])
            sites = (["{0}{1}".format(res, pos) for pos, res, score in scores[pair]])
            max_pssm_score = max(pssm_scores)
            max_pssm_score_i = pssm_scores.index(max_pssm_score)
            max_pssm_score_site = sites[max_pssm_score_i]
            all_pssm_scores = ",".join([str(score) for score in pssm_scores])
            all_sites = ",".join(sites)
            if residues == set(['Y']):
                (dcg, max_func_score, max_func_score_site,
                 all_func_scores) = regulation_score(
                     pair, scores[pair], phosfun_y, med_phosfun_y, reg_sites)
            else:
                (dcg, max_func_score, max_func_score_site,
                 all_func_scores) = regulation_score(
                    pair, scores[pair], phosfun_st, med_phosfun_st, reg_sites)
            v.write("\t".join([pair[0], pair[1], kin_type, sub_type, source,
                               str(max_pssm_score), str(dcg), str(max_func_score),
                               max_pssm_score_site, max_func_score_site,
                               all_pssm_scores, all_func_scores,
                               all_sites])+"\n")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("USAGE: <script> DATA_FILE OUT_FILE")
    kin_act_file_arg = sys.argv[1]
    out_file = sys.argv[2]
    score_network(kin_act_file_arg, out_file)
