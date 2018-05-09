import os.path
from Bio import Entrez
import time

from urllib.error import HTTPError

def get_gene_ids(kinome):
    kinome_entrez = dict()
    genename_uniprot = dict()
    uniprot_entrez = dict()
    with open("data/hgnc-to-entrez.tsv") as h:
        for line in h:
            hgnc, entrez = line.strip().split("\t")
            if hgnc in kinome:
                kinome_entrez[entrez] = hgnc
    kinome_entrez["83903"] = "GSG2"
    kinome_entrez["124923"] = "SGK494"
    return kinome_entrez


def read_entrez_data():
    kinome_citations = dict()
    kinome_citations_filt = dict()
    with open("data/kinome-entrez-data.tsv") as h:
        for line in h:
            kinase, num_cites, num_cites_filt = line.strip().split()
            kinome_citations[kinase] = num_cites
            kinome_citations_filt[kinase] = num_cites_filt
    return (kinome_citations, kinome_citations_filt)


def get_citations(kinome_entrez):
    if os.path.exists("data/kinome-entrez-data.tsv"):
        kinome_citations, kinome_citations_filt = read_entrez_data()
        return (kinome_citations, kinome_citations_filt)
    Entrez.email = "invergo@ebi.ac.uk"
    kinome_citations = dict()
    kinome_num_citations = dict()
    kinome_num_citations_filt = dict()
    pubs = dict()
    for entrez in kinome_entrez.keys():
        kinase = kinome_entrez[entrez]
        print(kinase)
        attempt = 0
        while attempt < 3:
            attempt += 1
            try:
                handle = Entrez.elink(dbfrom="gene", id=entrez, db="pubmed",
                                      retmode="xml")
            except HTTPError as err:
                if 500 <= err.code <= 599:
                    print("Received error from server {0}".format(err))
                    print("Attempt {0} of 3".format(attempt))
                    time.sleep(15)
                else:
                    raise
            record = Entrez.read(handle)
            handle.close()
            kinome_num_citations[kinase] = len(record[0]["LinkSetDb"][0]["Link"])
            kinome_citations[kinase] = set([pub["Id"] for pub in
                                            record[0]["LinkSetDb"][0]["Link"]])
            for pub_rec in record[0]["LinkSetDb"][0]["Link"]:
                pub = pub_rec["Id"]
                if pub not in pubs:
                    pubs[pub] = 1
                else:
                    pubs[pub] = pubs[pub] + 1
    bad_pubs = set()
    for pub in pubs:
        if pubs[pub] <= 10:
            continue
        bad_pubs.add(pub)
    for kinase in kinome_citations:
        kin_pubs = kinome_citations[kinase]
        good_pubs = set(kin_pubs) - bad_pubs
        kinome_num_citations_filt[kinase] = len(good_pubs)
    return (kinome_num_citations, kinome_num_citations_filt)


def write_num_citations(kinome_citations, kinome_citations_filt):
    with open("data/kinome-entrez-data.tsv", 'w') as h:
        for kinase in kinome_citations:
            num_citations = kinome_citations[kinase]
            num_citations_filt = kinome_citations_filt[kinase]
            h.write("\t".join([kinase, str(num_citations), str(num_citations_filt)])+"\n")


def get_num_substrates(kinome):
    kinase_subs = dict()
    kinase_num_subs = dict()
    with open("data/psiteplus-kinase-substrates-unfiltered.tsv") as h:
        for line in h:
            line_spl = line.strip().split("\t")
            kinase = line_spl[0]
            sub = line_spl[7]
            site = line_spl[10]
            if kinase not in kinase_subs:
                kinase_subs[kinase] = set([(sub, site)])
            else:
                kinase_subs[kinase].add((sub, site))
    for kinase in kinome:
        if kinase not in kinase_subs:
            kinase_num_subs[kinase] = 0
        else:
            kinase_num_subs[kinase] = len(kinase_subs[kinase])
    return kinase_num_subs


def get_num_kin_substrates(kinome):
    kinase_subs = dict()
    kinase_num_subs = dict()
    with open("data/psiteplus-kinase-substrates-unfiltered.tsv") as h:
        for line in h:
            line_spl = line.strip().split("\t")
            kinase = line_spl[0]
            sub = line_spl[7]
            if sub not in kinome:
                continue
            site = line_spl[10]
            if kinase not in kinase_subs:
                kinase_subs[kinase] = set([(sub, site)])
            else:
                kinase_subs[kinase].add((sub, site))
    for kinase in kinome:
        if kinase not in kinase_subs:
            kinase_num_subs[kinase] = 0
        else:
            kinase_num_subs[kinase] = len(kinase_subs[kinase])
    return kinase_num_subs


def get_kin_phos_rels(kinome):
    kinase_rels = dict()
    with open("data/psiteplus-kinase-substrates-unfiltered.tsv") as h:
        for line in h:
            line_spl = line.strip().split("\t")
            kinase = line_spl[0]
            sub = line_spl[7]
            if sub not in kinome:
                continue
            if kinase not in kinase_rels:
                kinase_rels[kinase] = set([sub])
            else:
                kinase_rels[kinase].add(sub)
    return kinase_rels


def get_kin_num_phos_kins(kinome):
    kinase_regs = dict()
    with open("data/psiteplus-kinase-substrates-unfiltered.tsv") as h:
        for line in h:
            line_spl = line.strip().split("\t")
            sub = line_spl[7]
            if sub not in kinome:
                continue
            if sub not in kinase_regs:
                kinase_regs[sub] = 1
            else:
                kinase_regs[sub] = kinase_regs[sub] + 1
    return kinase_regs


def get_num_reg_sites(kinome):
    reg_sites = dict()
    kinase_num_reg_sites = dict()
    with open("data/psiteplus-reg-sites-unfiltered.tsv") as h:
        for line in h:
            line_spl = line.strip().split("\t")
            gene = line_spl[0]
            site = line_spl[8]
            func = line_spl[11]
            if "activity, induced" not in func and  "activity, inhibited" not in func:
                continue
            if gene not in reg_sites:
                reg_sites[gene] = set([site])
            else:
                reg_sites[gene].add(site)
    for kinase in kinome:
        if kinase not in reg_sites:
            kinase_num_reg_sites[kinase] = 0
        else:
            kinase_num_reg_sites[kinase] = len(reg_sites[kinase])
    return kinase_num_reg_sites


def get_kin_reg_rels(kinome):
    kinase_regs = dict()
    with open("data/validation-set-omnipath.tsv") as h:
        for line in h:
            kinase, sub = line.strip().split()
            if kinase not in kinase_regs:
                kinase_regs[kinase] = set([sub])
            else:
                kinase_regs[kinase].add(sub)
    return kinase_regs


if __name__ == "__main__":
    kinome = set()
    with open("data/human-kinome.txt") as h:
        for line in h:
            kinome.add(line.strip())
    kinome_entrez = get_gene_ids(kinome)
    kinome_citations, kinome_citations_filt = get_citations(kinome_entrez)
    write_num_citations(kinome_citations, kinome_citations_filt)
    kinome_num_subs = get_num_substrates(kinome)
    kinome_num_kin_subs = get_num_kin_substrates(kinome)
    kinome_num_reg_sites = get_num_reg_sites(kinome)
    kinome_num_phos_kins = get_kin_num_phos_kins(kinome)
    kinase_rels = get_kin_phos_rels(kinome)
    kinase_regs = get_kin_reg_rels(kinome)
    with open("out/kinase-citations.tsv", 'w') as h:
        h.write("\t".join(["kinase", "count.type", "num.cites", "num.cites.filt",
                           "count"])+"\n")
        # h.write("\t".join(["kinase", "num.cites", "num.subs", "num.kin.subs",
        #                    "num.reg.sites"])+"\n")
        for kinase in kinome:
            if kinase not in kinome_citations:
                num_cites = 0
                num_cites_filt = 0
            else:
                num_cites = kinome_citations[kinase]
                num_cites_filt = kinome_citations_filt[kinase]
            for dat_type in ["num.subs", "num.kin.subs", "num.reg.sites",
                             "num.phos.kins"]:
                if dat_type == "num.subs":
                    dat = kinome_num_subs[kinase]
                elif dat_type == "num.kin.subs":
                    dat = kinome_num_kin_subs[kinase]
                elif dat_type == "num.phos.kins":
                    if kinase not in kinome_num_phos_kins:
                        dat = 0
                    else:
                        dat = kinome_num_phos_kins[kinase]
                else:
                    dat = kinome_num_reg_sites[kinase]
                h.write("\t".join([kinase, dat_type, str(num_cites),
                                   str(num_cites_filt), str(dat)])+"\n")
            # num_subs = kinome_num_subs[kinase]
            # num_kin_subs = kinome_num_kin_subs[kinase]
            # num_reg_sites = kinome_num_reg_sites[kinase]
            # h.write("\t".join([kinase, str(num_cites), str(num_subs),
            #                    str(num_kin_subs), str(num_reg_sites)])+"\n")
    with open("out/kinase-kinase-rels.tsv", 'w') as h:
        h.write("\t".join(["kinase1", "num.cites1", "num.cites.filt1",
                           "kinase2", "num.cites2", "num.cites.filt2",
                           "rel"])+"\n")
        for kinase in kinome:
            if kinase not in kinome_citations:
                num_cites1 = 0
                num_cites_filt1 = 0
            else:
                num_cites1 = kinome_citations[kinase]
                num_cites_filt1 = kinome_citations_filt[kinase]
            if kinase in kinase_rels:
                for substrate in kinase_rels[kinase]:
                    if substrate not in kinome_citations:
                        num_cites2 = 0
                        num_cites_filt2 = 0
                    else:
                        num_cites2 = kinome_citations[substrate]
                        num_cites_filt2 = kinome_citations_filt[substrate]
                    h.write("\t".join([kinase, str(num_cites1),
                                       str(num_cites_filt1), substrate,
                                       str(num_cites2), str(num_cites_filt2),
                                       "phos"])+"\n")
            if kinase in kinase_regs:
                for substrate in kinase_regs[kinase]:
                    if substrate not in kinome_citations:
                        num_cites2 = 0
                        num_cites_filt2 = 0
                    else:
                        num_cites2 = kinome_citations[substrate]
                        num_cites_filt2 = kinome_citations_filt[substrate]
                    h.write("\t".join([kinase, str(num_cites1),
                                       str(num_cites_filt1),
                                       substrate, str(num_cites2),
                                       str(num_cites_filt2),
                                       "reg"])+"\n")
