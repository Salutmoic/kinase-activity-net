#!/bin/env python2

import sys
import os.path
import pypath


def get_reg_rels(pa, prot_list):
    reg_rels = set()
    for prot1 in prot_list:
        for prot2 in prot_list:
            if prot1 == prot2:
                reg_rels.add((prot1, prot2, "NA"))
            elif prot1 in list(pa.gs_stimulated_by(prot2).gs()):
                reg_rels.add((prot1, prot2, "activ"))
            elif prot1 in list(pa.gs_inhibited_by(prot2).gs()):
                reg_rels.add((prot1, prot2, "inhib"))
            else:
                reg_rels.add((prot1, prot2, "NA"))
    return reg_rels

if __name__ == "__main__":
    if not os.path.exists("cache/default_network.pickle"):
        sys.exit("You must first run init-omnipath.py")
    pa = pypath.PyPath()
    pa.init_network(pfile="cache/default_network.pickle")
    pa.get_directed()
    kinases = set(["AKT1", "BRAF", "MAP2K1", "MAPK1", "PDPK1", "RPS6KB1"])
    reg_rels = get_reg_rels(pa, kinases)
    with open("data/egf-validation-set-omnipath.tsv", 'w') as h:
        h.write("\t".join(["node1", "node2", "val"])+"\n")
        for rel in reg_rels:
            h.write("\t".join(rel)+"\n")
