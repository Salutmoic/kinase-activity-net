#!/bin/env python2

import sys
import os.path
import pypath


def get_assocs(pa, prot_list):
    assocs = set()
    for prot in prot_list:
        neighbs = list(pa.gs_neighbors(prot).gs())
        for neighb in neighbs:
            if neighb == prot or neighb not in prot_list:
                continue
            srcs = pa.get_edge(prot, neighb)["sources"]
            if len(srcs) < 2:
                continue
            assocs.add((prot, neighb))
    return assocs


def get_reg_rels(pa, prot_list):
    reg_rels = set()
    for prot in prot_list:
        # if pa.dgenesymbol(prot) is None:
        #     continue
        # print(prot)
        # regs = list(pa.gs_stimulated_by(prot).gs())
        # regs.extend(list(pa.gs_inhibited_by(prot).gs()))
        regs = list(pa.gs_affected_by(prot).gs())
        for reg in regs:
            if reg == prot or reg not in prot_list:
                continue
            srcs = pa.get_edge(reg, prot)["sources"]
            if len(srcs) < 2:
                continue
            reg_rels.add((reg, prot))
    return reg_rels


if __name__ == "__main__":
    if not os.path.exists("cache/default_network.pickle"):
        sys.exit("You must first run init-omnipath.py")
    if len(sys.argv) != 2:
        sys.exit("USAGE: <script> KINS_TO_USE")
    kin_list_file = sys.argv[1]
    pa = pypath.PyPath()
    pa.init_network(pfile="cache/default_network.pickle")
    kinases = set()
    with open(kin_list_file) as h:
        for line in h:
            kinases.add(line.strip())
    assocs = get_assocs(pa, kinases)
    with open("data/assoc-validation-set-omnipath.tsv", 'w') as h:
        for assoc in assocs:
            h.write("\t".join(assoc)+"\n")
    pa.get_directed()
    reg_rels = get_reg_rels(pa, kinases)
    with open("data/direct-validation-set-omnipath.tsv", 'w') as h:
        for rel in reg_rels:
            h.write("\t".join(rel)+"\n")
