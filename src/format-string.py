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


def reformat_string(string_h, string_id_map, kinases, kinact_intxns, score_col):
    intxns = dict()
    score_title = None
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
        if prot1 not in kinases or prot2 not in kinases:
            continue
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
    for kinact_intxn in kinact_intxns:
        if kinact_intxn in intxns:
            continue
        prot1, prot2 = kinact_intxn
        print("\t".join([prot1, prot2, str(0.0)]))


if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.exit("USAGE: <script> STRING STRING_ID_MAP KINACT_DATA COLUMN")
    string_file = sys.argv[1]
    string_id_map_file = sys.argv[2]
    kinact_file = sys.argv[3]
    score_col = int(sys.argv[4]) - 1
    with open(string_id_map_file) as string_h:
        string_id_map = read_id_map(string_h)
    with open(kinact_file) as h:
        kinases = read_kinact_tbl(h)
    kinact_intxns = product(kinases, repeat=2)
    with open(string_file) as string_h:
        reformat_string(string_h, string_id_map, kinases, kinact_intxns, score_col)
