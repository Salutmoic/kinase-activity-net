#!/bin/env python2

import sys
import os.path
import pypath
from pypath import mapping


def get_reg_rels(pa, prot_list):
    reg_rels = set()
    m = mapping.Mapper(9606)
    for prot in prot_list:
        # if pa.dgenesymbol(prot) is None:
        #     continue
        # print(prot)
        pos_regs = set(pa.gs_stimulated_by(prot).gs())
        neg_regs = set(pa.gs_inhibited_by(prot).gs())
        ambig_regs = pos_regs.intersection(neg_regs)
        pos_regs = pos_regs - ambig_regs
        neg_regs = neg_regs - ambig_regs
        for reg in pos_regs:
            if reg == prot or reg not in prot_list:
                continue
            srcs = pa.get_edge(reg, prot)["sources"]
            if len(srcs) < 2:
                continue
            reg_rels.add((reg, prot, "activates"))
        for reg in neg_regs:
            if reg == prot or reg not in prot_list:
                continue
            srcs = pa.get_edge(reg, prot)["sources"]
            if len(srcs) < 2:
                continue
            reg_rels.add((reg, prot, "inhibits"))
        for reg in ambig_regs:
            if reg == prot or reg not in prot_list:
                continue
            srcs = pa.get_edge(reg, prot)["sources"]
            if len(srcs) < 2:
                continue
            dirs = pa.get_edge(reg, prot)["dirs"]
            signs_dict = dirs.majority_sign()
            signs = None
            for sign_pair in signs_dict:
                sign_reg = m.map_name(sign_pair[0], "uniprot", "genesymbol")[0]
                sign_sub = m.map_name(sign_pair[1], "uniprot", "genesymbol")[0]
                if sign_reg == reg and sign_sub == prot:
                    sign = signs_dict[sign_pair]
                    break
            if sign is None:
                print("crap", reg, prot)
            if sign is None or sign[0] and sign[1]:
                # ambiguous sign
                continue
            elif sign[0]:
                reg_rels.add((reg, prot, "activates"))
            elif sign[1]:
                reg_rels.add((reg, prot, "inhibits"))
            else:
                continue
    return reg_rels


if __name__ == "__main__":
    if not os.path.exists("cache/default_network.pickle"):
        sys.exit("You must first run init-omnipath.py")
    if len(sys.argv) != 2:
        sys.exit("USAGE: <script> KINS_TO_USE")
    kin_list_file = sys.argv[1]
    pa = pypath.PyPath()
    pa.init_network(pfile="cache/default_network.pickle")
    pa.get_directed()
    kinases = set()
    with open(kin_list_file) as h:
        for line in h:
            kinases.add(line.strip())
    reg_rels = get_reg_rels(pa, kinases)
    with open("data/sign-validation-set-omnipath.tsv", 'w') as h:
        for rel in reg_rels:
            h.write("\t".join(rel)+"\n")
