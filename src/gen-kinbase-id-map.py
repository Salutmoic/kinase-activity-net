import re


# ID mappings missing from KinHub
EXTRA_ID_MAP = {"MLK4": "MAP3K21",
                "ZAK": "MAP3K20",
                "SgK494": "SGK494",
                "SgK223": "PRAG1"}


def read_uniprot_id_map(id_map_handle):
    gene_map = dict()
    syn_map = dict()
    id_map = dict()
    for line in id_map_handle:
        uniprot, id_type, other_id = line.strip().split("\t")
        if id_type == "Gene_Name":
            gene_map[uniprot] = other_id
        elif id_type == "Gene_Synonym":
            syn_map[other_id] = uniprot
    for syn_id in syn_map:
        uniprot = syn_map[syn_id]
        gene_id = gene_map[uniprot]
        id_map[syn_id] = gene_id
    return id_map


def read_kinhub_id_map(id_map_handle):
    gene_map = dict()
    syn_map = dict()
    id_map = dict()
    for line in id_map_handle:
        xname_id, _, hgnc_id = line.strip().split("\t")[:3]
        id_map[xname_id] = hgnc_id
    return id_map


def translate_id(prot_id, uniprot_id_map, kinhub_id_map):
    if prot_id in kinhub_id_map:
        new_id = kinhub_id_map[prot_id]
        if not new_id and prot_id in EXTRA_ID_MAP:
            new_id = EXTRA_ID_MAP[prot_id]
        return new_id
    prot_id = prot_id.upper()
    if prot_id in uniprot_id_map:
        return uniprot_id_map[prot_id]
    return prot_id


def print_id_map(kb_fasta_h, uniprot_id_map, kinhub_id_map):
    gene_re = re.compile(r"gene=(\S+)")
    for line in h:
        if not line.startswith(">"):
            continue
        m = gene_re.search(line)
        if m is None:
            continue
        gene = m.group(1)
        proper_id = translate_id(gene, uniprot_id_map, kinhub_id_map)
        if not proper_id:
            continue
        print(gene, proper_id, sep="\t")


if __name__ == "__main__":
    uniprot_id_map_file = "data/external/HUMAN_9606_idmapping.dat"
    kinhub_id_map_file = "data/external/kinhub-families.tsv"
    kinbase_seq_file = "data/external/kinbase-sequences.fasta"
    with open(uniprot_id_map_file) as h:
        uniprot_id_map = read_uniprot_id_map(h)
    with open(kinhub_id_map_file) as h:
        kinhub_id_map = read_kinhub_id_map(h)
    with open(kinbase_seq_file) as h:
        print_id_map(h, uniprot_id_map, kinhub_id_map)
