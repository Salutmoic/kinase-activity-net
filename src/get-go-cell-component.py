########################################################################
#                                                                      #
#            script to take the cellular complexes 		       #
#	     from GO-cellular component       			       #    
#                                                                      #
# Copyright (c) 2016 - 2019, EMBL - European Bioinformatics Institute  #
#       Author: Bruno Ariano (ariano@ebi.ac.uk)                        #
#       Distributed under the GPLv3 License.                           #
########################################################################

import re
from goatools.obo_parser import GODag
from goatools.associations import read_associations
import sys


def get_really_all_parents(term_rec):
    all_parents = set()
    for p in term_rec.parents:
        all_parents.add(p.id)
        all_parents |= get_really_all_parents(p)
    if hasattr(term_rec, "relationship") and "part_of" in term_rec.relationship:
        for p in term_rec.relationship["part_of"]:
            all_parents.add(p.id)
            all_parents |= get_really_all_parents(p)
    if hasattr(term_rec, "relationship") and "occurs_in" in term_rec.relationship:
        for p in term_rec.relationship["part_of"]:
            all_parents.add(p.id)
            all_parents |= get_really_all_parents(p)
    return all_parents


def get_associations(assoc):
    GO_term = {}
    for prot in assoc:
        terms = assoc[prot]
        for term in terms:
            if term in GO_term:
                GO_term[term].append(prot)
            else:
                GO_term[term] = [prot]
    return GO_term


def print_locations(GO_term, cytoplasm, plasma_membrane, nucleus):
    for term in cytoplasm:
        if term not in GO_term:
            continue
        for prot_id in GO_term[term]:
            print("GO:0005737\t" + prot_id)
    for term in plasma_membrane:
        if term not in GO_term:
            continue
        for prot_id in GO_term[term]:
            print("GO:0005886\t" + prot_id)
    for term in nucleus:
        if term not in GO_term:
            continue
        for prot_id in GO_term[term]:
            print("GO:0005634\t" + prot_id)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("USAGE: <script> ASSOCS OBO")
    assoc_file = sys.argv[1]
    obo_file = sys.argv[2]
    obo_dag = GODag(obo_file=obo_file, optional_attrs=["relationship"])
    cytoplasm = set()
    plasma_membrane = set()
    nucleus = set()
    for term in obo_dag:
        term_rec = obo_dag[term]
        parents = get_really_all_parents(term_rec)
        if "GO:0005737" in parents:
            cytoplasm.add(term)
        if "GO:0005886" in parents:
            plasma_membrane.add(term)
        if "GO:0005634" in parents:
            nucleus.add(term)
    assoc = read_associations(assoc_file)
    GO_term = get_associations(assoc)
    print_locations(GO_term, cytoplasm, plasma_membrane, nucleus)
