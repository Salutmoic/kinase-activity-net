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
        posind = header.index("MOD_RSD")
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
                ["Pkinase", "Pkinase_Tyr", "PI3_PI4_kinase", "Alpha_kinase"]:
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
                or kin_b not in psites):
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


def calc_dcg(scores):
    return sum([score/log2(n+2) for n, score in enumerate(scores)])


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
            sys.exit("\t".join([pair[1], pos]))
        hs_score = kin_b_phosfun[pos]
        hs_sign = kin_b_signs[pos]
        hs_scores.append(hs_score)
        hs_signs.append(hs_sign)
    if abs(min(hs_signs)) > max(hs_signs):
        site_score = min(hs_signs)
    else:
        site_score = max(hs_signs)
    # Discounted Cumulative Gain
    # signs = [copysign(1, x) for x in hs_signs]
    # signed_func_scores = [sign * score for (sign, score) in zip(signs, hs_scores)]
    signed_func_scores = [sign * score for (sign, score) in zip(hs_signs, hs_scores)]
    if abs(min(signed_func_scores)) > max(signed_func_scores):
        signed_func_score = min(signed_func_scores)
    else:
        signed_func_score = max(signed_func_scores)
    dcg = calc_dcg(signed_func_scores)
    if dcg < 0:
        signed_func_scores.sort()
        best_dcg = calc_dcg(signed_func_scores)
        signed_func_scores.sort(reverse=True)
        worst_dcg = calc_dcg(signed_func_scores)
        sign = -1
    else:
        signed_func_scores.sort()
        worst_dcg = calc_dcg(signed_func_scores)
        signed_func_scores.sort(reverse=True)
        best_dcg = calc_dcg(signed_func_scores)
        sign = 1
    if best_dcg == worst_dcg:
        return (0.0, site_score, signed_func_score)
    idcg = sign*(dcg-worst_dcg)/(best_dcg-worst_dcg)
    return (idcg, site_score, signed_func_score)


def score_network(kinase_file, out_file):
    # This takes network files and scores all edges according to pssm

    ktable = pssms.read_kin_sub_data("data/psiteplus-kinase-substrates.tsv")
    aa_freqs = pssms.read_aa_freqs("data/aa-freqs.tsv")
    psites = read_phosphosites("data/pride-phosphosites.tsv")
    kinases = read_kinase_file(kinase_file)
    phosfun_st, med_phosfun_st = read_phosfun_file("data/phosfun-ST.tsv")
    phosfun_y, med_phosfun_y = read_phosfun_file("data/phosfun-Y.tsv")
    signs = read_sign_file("out/reg-site-bart-sign-preds.tsv")
    reg_sites = read_reg_sites_file("data/psiteplus-reg-sites.tsv")
    kin_types = read_pfam_file("data/human-pfam.tsv", kinases)
    scores = score_kin_pairs(ktable, kinases, aa_freqs, psites)
    with open(out_file, "w") as v:
        v.write("\t".join(["node1", "node2", "kinase.type", "sub.kinase.type",
                           "max.pssm.score", "signed.dcg", "site.sign.score",
                           "signed.func.score"])+"\n")
        for pair in scores:
            if pair[0] not in kin_types:
                kin_type = "NA"
            else:
                kin_type = kin_types[pair[0]]
            if pair[1] not in kin_types:
                sub_type = "NA"
            else:
                sub_type = kin_types[pair[1]]
            if not scores[pair]:
                v.write("\t".join([pair[0], pair[1], kin_type, sub_type, "NA",
                                   "NA", "NA", "NA"])+"\n")
                continue
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
