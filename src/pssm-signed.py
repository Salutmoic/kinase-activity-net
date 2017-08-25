import os.path
import sys
import itertools
import pssms
from math import log2, log, exp, copysign


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


def read_sign_file(sign_file):
    signs = dict()
    with open(sign_file) as h:
        h.readline()
        for line in h:
            prot, pos, score = line.strip().split()
            if prot not in signs:
                signs[prot] = {pos: float(score)}
            else:
                signs[prot][pos] = float(score)
    return signs


def read_reg_sites_file(reg_sites_file):
    reg_sites = dict()
    with open(reg_sites_file) as h:
        header = h.readline().split("\t")
        geneind = header.index("GENE")
        posind = header.index("position")
        funcind = header.index("ON_FUNCTION")
        for line in h:
            line_spl = line.strip().split("\t")
            protein = line_spl[0]
            pos = line_spl[posind].replace("-p", "")[1:]
            func = line_spl[funcind]
            if "activity, induced" in func:
                sign = 1.0
            elif "activity, inhibited" in func:
                sign = -1.0
            else:
                sign = 0.0
            if protein not in reg_sites:
                reg_sites[protein] = {pos: sign}
            else:
                reg_sites[protein][pos] = sign
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
        kinase_pssms[kinase] = (assign_pssm, False)
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


def calc_dcg(scores):
    # return sum([score/log2(n+2) for n, score in enumerate(scores)])
    dcg = list(itertools.accumulate([score/log2(n+2) for n, score in enumerate(scores)]))
    abs_dcg = [abs(x) for x in dcg]
    max_dcg = max(abs_dcg)
    max_dcg_i = abs_dcg.index(max_dcg)
    return(dcg[max_dcg_i])


def signed_regulation_score(pair, scores, phosfun, med_phosfun, reg_sites, signs):
    # Rank the PSSM scores
    scores.sort(key=lambda t: t[2], reverse=True)
    kin_b_reg_sites = reg_sites.get(pair[1])
    kin_b_phosfun = phosfun.get(pair[1])
    kin_b_signs = signs.get(pair[1])
    hs_scores = []
    hs_signs = []
    for pos, res, score in scores:
        if kin_b_phosfun is None or pos not in kin_b_phosfun:
            continue
        hs_score = kin_b_phosfun[pos]
        hs_sign = kin_b_signs[pos]
        hs_scores.append(hs_score)
        hs_signs.append(hs_sign)
    if not hs_scores or not hs_signs:
        return("NA", "NA", "NA")
    if abs(min(hs_signs)) > max(hs_signs):
        site_score = min(hs_signs)
    else:
        site_score = max(hs_signs)
    # Discounted Cumulative Gain
    # signs = [copysign(1, x) for x in hs_signs]
    # signed_func_scores = [sign * score for (sign, score) in zip(signs, hs_scores)]
    signed_func_scores = [copysign(1.0, sign) * score for (sign, score) in
                          zip(hs_signs, hs_scores)]
    # signed_func_scores = [sign * score for (sign, score) in zip(hs_signs, hs_scores)]
    if abs(min(signed_func_scores)) > max(signed_func_scores):
        signed_func_score = min(signed_func_scores)
    else:
        signed_func_score = max(signed_func_scores)
    dcg = calc_dcg(signed_func_scores)
    if dcg < 0:
        signed_func_scores.sort()
        best_dcg = calc_dcg(signed_func_scores)
        signed_func_scores.sort(reverse=True)
        # worst_dcg = calc_dcg(signed_func_scores)
        worst_dcg = 0
        sign = -1
    else:
        signed_func_scores.sort()
        # worst_dcg = calc_dcg(signed_func_scores)
        worst_dcg = 0
        signed_func_scores.sort(reverse=True)
        best_dcg = calc_dcg(signed_func_scores)
        sign = 1
    if best_dcg == worst_dcg or len(hs_scores) < 2:
        return ("NA", site_score, signed_func_score)
    idcg = sign*(dcg-worst_dcg)/(best_dcg-worst_dcg)
    return (idcg, site_score, signed_func_score)


def score_network(kinase_file, out_file):
    # This takes network files and scores all edges according to pssm

    ktable = pssms.read_kin_sub_data("data/psiteplus-kinase-substrates.tsv")
    aa_freqs = pssms.read_aa_freqs("data/aa-freqs.tsv")
    psites = read_phosphosites("data/psiteplus-phosphosites.tsv")
    kinases = read_kinase_file(kinase_file)
    phosfun_st, med_phosfun_st = read_phosfun_file("data/phosfun-ST.tsv")
    phosfun_y, med_phosfun_y = read_phosfun_file("data/phosfun-Y.tsv")
    signs = read_sign_file("out/reg-site-bart-sign-preds.tsv")
    reg_sites = read_reg_sites_file("data/psiteplus-reg-sites.tsv")
    kin_types = read_pfam_file("data/human-pfam.tsv", kinases)
    kinase_fams, family_kins = read_kinase_families("data/kinase-families.tsv")
    sdr_assigns = read_sdr_assignments("data/sdr-pssm-assignments.tsv")
    family_pssms = pssms.build_family_pssms(family_kins, ktable, aa_freqs)
    kinase_pssms = pssms.build_pssms(kinases, kinase_fams, family_pssms, ktable,
                                     aa_freqs)
    kinase_pssms = assign_pssms_by_sdr(kinases, kinase_pssms, sdr_assigns)
    scores = score_kin_pairs(kinases, kinase_pssms, psites, family_pssms)
    # scores = norm_scores
    with open(out_file, "w") as v:
        v.write("\t".join(["node1", "node2", "kinase.type", "sub.kinase.type",
                           "max.pssm.score", "signed.dcg", "site.sign.score",
                           "signed.func.score"])+"\n")
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
                v.write("\t".join([pair[0], pair[1], kin_type, sub_type, "NA",
                                   "NA", "NA", "NA"])+"\n")
                continue
            own_pssm = str(kinase_pssms[pair[0]][1]).upper()
            residues = set([res for pos, res, score in scores[pair]])
            max_score = max([score for pos, res, score in scores[pair]])
            if residues == set(['Y']):
                dcg, site_score, signed_func_score = signed_regulation_score(
                    pair, scores[pair], phosfun_y, med_phosfun_y, reg_sites, signs)
            else:
                dcg, site_score, signed_func_score = signed_regulation_score(
                    pair, scores[pair], phosfun_st, med_phosfun_st, reg_sites, signs)
            v.write("\t".join([pair[0], pair[1], kin_type, sub_type, str(max_score),
                               str(dcg), str(site_score), str(signed_func_score)])+"\n")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("USAGE: <script> DATA_FILE OUT_FILE")
    kin_act_file_arg = sys.argv[1]
    out_file = sys.argv[2]
    score_network(kin_act_file_arg, out_file)
