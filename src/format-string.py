import sys


def read_id_map(id_map_h):
    id_map = dict()
    for line in id_map_h:
        from_id, to_id = line.strip().split()
        id_map[from_id] = to_id
    return(id_map)


def read_kinact_tbl(h):
    kinases = set()
    for line in h:
        prot1, prot2, score = line.strip().split()
        kinases.add(prot1)
        kinases.add(prot2)
    return kinases


def reformat_string(string_h, string_id_map, kinases):
    intxns = dict()
    for line in string_h:
        ensp1, ensp2, score = line.strip().split()
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
    for intxn in intxns:
        prot1, prot2 = intxn
        score_mod = intxns[intxn]
        print("\t".join([prot1, prot2, str(score_mod)]))


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit("USAGE: <script> STRING STRING_ID_MAP KINACT_DATA")
    string_file = sys.argv[1]
    string_id_map_file = sys.argv[2]
    kinact_file = sys.argv[3]
    with open(string_id_map_file) as string_h:
        string_id_map = read_id_map(string_h)
    with open(kinact_file) as h:
        kinases = read_kinact_tbl(h)
    with open(string_file) as string_h:
        reformat_string(string_h, string_id_map, kinases)
