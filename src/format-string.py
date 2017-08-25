import sys
from itertools import product


def read_id_map(id_map_h):
    id_map = dict()
    for line in id_map_h:
        from_id, to_id = line.strip().split()
        id_map[from_id] = to_id
    return(id_map)


def read_kinact_tbl(h):
    kinases = set()
    for line in h:
        prot, cond, score = line.strip().split()
        kinases.add(prot)
    return kinases


def reformat_string(string_h, string_id_map, kinases, kinase_pairs, score_col):
    intxns = dict()
    score_title = None
    if kinases is None:
        prots = set()
    else:
        prots = None
        prot_pairs = kinase_pairs
    for line in string_h:
        line_spl = line.strip().split()
        ensp1 = line_spl[0]
        ensp2 = line_spl[1]
        score = line_spl[score_col]
        if not score_title:
            # Basically, are we on the first line?  If so, grab the
            # name of our score column
            score_title = "string." + score
            continue
        if ensp1 not in string_id_map or ensp2 not in string_id_map:
            continue
        prot1 = string_id_map[ensp1]
        prot2 = string_id_map[ensp2]
        if kinases is not None:
            if prot1 not in kinases or prot2 not in kinases:
                continue
        else:
            prots.update([prot1, prot2])
        score_mod = float(score)/1000
        if (prot1, prot2) not in intxns:
            intxns[(prot1, prot2)] = score_mod
        else:
            if score_mod > intxns[(prot1, prot2)]:
                intxns[(prot1, prot2)] = score_mod
    print("\t".join(["node1", "node2", score_title]))
    for intxn in intxns:
        prot1, prot2 = intxn
        score_mod = intxns[intxn]
        print("\t".join([prot1, prot2, str(score_mod)]))
    if kinases is None:
        prot_pairs = product(prots, repeat=2)
    if prot_pairs is not None:
        for prot_pair in prot_pairs:
            if prot_pair in intxns:
                continue
            prot1, prot2 = prot_pair
            print("\t".join([prot1, prot2, "NA"]))


if __name__ == "__main__":
    if len(sys.argv) not in [5, 6]:
        sys.exit("USAGE: <script> STRING STRING_ID_MAP KINASE_LIST COLUMN [--all]")
    string_file = sys.argv[1]
    string_id_map_file = sys.argv[2]
    kinase_file = sys.argv[3]
    score_col = int(sys.argv[4]) - 1
    if len(sys.argv) == 6 and sys.argv[5] == "--all":
        keep_all = True
    else:
        keep_all = False
    with open(string_id_map_file) as string_h:
        string_id_map = read_id_map(string_h)
    if not keep_all:
        kinases = set()
        with open(kinase_file) as h:
            for line in h:
                kinases.add(line.strip())
        kinase_pairs = product(kinases, repeat=2)
    else:
        kinases = None
        kinase_pairs = None
    with open(string_file) as string_h:
        reformat_string(string_h, string_id_map, kinases, kinase_pairs,
                        score_col)
