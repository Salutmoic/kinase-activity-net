import os.path
import sys
import itertools
import pssms


def read_phosphosites(psites_file):
    # returns dictionary with every known phosphosite found on each
    # kinase in the Ochoa et al. dataset

    psites = {}
    with open(psites_file) as f:
        columns = f.readline()
        columns = columns.split("\t")
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
        sub_motif_seqs = [seq for seq in psites[kin_b]
                          if ((kin_a_is_tyr_kin and seq[7] == 'Y') or
                              (not kin_a_is_tyr_kin and seq[7] in ['S', 'T']))]
        if not sub_motif_seqs:
            scores[(kin_a, kin_b)] = []
            continue

        for motif_seq in sub_motif_seqs:
            score = pssms.score_motif_seq(motif_seq, pssm)
            if (kin_a, kin_b) in scores:
                scores[(kin_a, kin_b)].append(score)
            else:
                scores[(kin_a, kin_b)] = [score]
    return scores


def score_network(kinase_file):
    # This takes network files and scores all edges according to pssm

    ktable = pssms.read_kin_sub_data("data/reduced_kinase_table.tsv")
    aa_freqs = pssms.read_aa_freqs("data/aa-freqs.tsv")
    psites = read_phosphosites("data/phosphosites_reduced.tsv")
    kinases = read_kinase_file(kinase_file)
    scores = score_kin_pairs(ktable, kinases, aa_freqs, psites)
    out_file_base = os.path.splitext(os.path.basename(kinase_file))[0]
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
    kin_act_file_arg = sys.argv[1]
    score_network(kin_act_file_arg)
