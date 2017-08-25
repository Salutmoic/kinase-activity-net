import sys
import os.path
import pypath


def write_pathways(pa, prots, h):
    paths = ["kegg", "netpath", "signalink"]
    for prot in prots:
        if pa.gs(prot) is None:
            continue
        for path in paths:
            key = path + "_pathways"
            if key not in pa.gs(prot).attribute_names():
                continue
            prot_paths = pa.gs(prot)[path+"_pathways"]
            for prot_path in prot_paths:
                h.write("\t".join([path+"--"+prot_path.replace(" ", "-"), prot])+"\n")


if __name__ == "__main__":
    if not os.path.exists("cache/default_network.pickle"):
        sys.exit("run init-omnipath.py first")
    if len(sys.argv) != 2:
        sys.exit("USAGE: <script> KINS_TO_USE")
    kin_list_file = sys.argv[1]
    kinases = set()
    with open(kin_list_file) as h:
        for line in h:
            kinases.add(line.strip())
    pa = pypath.PyPath()
    pa.init_network(pfile="cache/default_network.pickle")
    pa.kegg_pathways()
    # pa.signor_pathways()
    pa.pathway_attributes()
    with open("data/omnipath-path-ref.tsv", 'w') as h:
        write_pathways(pa, kinases, h)
