########################################################################
#                                                                      #
#            script to take the cellular complexes 		       #
#	     from GO-cellular component       			       #    
#                                                                      #
# Copyright (c) 2016 - 2019, 2017, EMBL - European Bioinformatics Institute  #
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


def print_locations(assoc, cytoplasm, plasma_membrane, nucleus):
    for prot in assoc:
        in_cyto = False
        in_memb = False
        in_nucl = False
        for term in assoc[prot]:
            in_cyto = in_cyto or term in cytoplasm
            in_memb = in_memb or term in plasma_membrane
            in_nucl = in_nucl or term in nucleus
        if in_cyto:
            print("GO:0005737\t" + prot)
        # if in_memb:
        #     print("GO:0005886\t" + prot)
        if in_nucl:
            print("GO:0005634\t" + prot)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("USAGE: <script> ASSOCS OBO")
    assoc_file = sys.argv[1]
    obo_file = sys.argv[2]
    obo_dag = GODag(obo_file=obo_file, optional_attrs=["relationship"])
    cytoplasm = set(["GO:0005737"])
    plasma_membrane = set(["GO:0005886"])
    nucleus = set(["GO:0005634"])
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
    print_locations(assoc, cytoplasm, plasma_membrane, nucleus)
