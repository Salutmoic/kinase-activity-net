import pssms


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


def read_sdr_file(sdr_handle, ktable, uniprot_id_map, kinhub_id_map):
    sdr_handle.readline()
    for line in sdr_handle:
        kinase, sim_kinases_str, score = line.strip().split()
        kinase = translate_id(kinase, uniprot_id_map, kinhub_id_map)
        if kinase in ktable and len(ktable[kinase]) > 10:
            print(kinase, kinase, 1.0)
            continue
        sim_kinases = sim_kinases_str.split("/")
        if len(sim_kinases) == 1:
            sim_kin = translate_id(sim_kinases[0], uniprot_id_map, kinhub_id_map)
            if sim_kin and sim_kin in ktable and len(ktable[sim_kin])>10:
                print(kinase, sim_kin, float(score))
            else:
                print(kinase, "None", "NA")
            continue
        best_kinase = None
        best_n = 0
        for sim_kin in sim_kinases:
            sim_kin = translate_id(sim_kin, uniprot_id_map, kinhub_id_map)
            if not sim_kin:
                continue
            if sim_kin not in ktable or len(ktable[sim_kin])<=10:
                continue
            n = len(ktable[sim_kin])
            if n > best_n:
                best_n = n
                best_kinase = sim_kin
        if best_kinase is None:
            print(kinase, "None", "NA")
            continue
        print(kinase, best_kinase, float(score))


if __name__ == "__main__":
    sdr_file = "data/external/sdr-kinase-specificity-predictions.txt"
    uniprot_id_map_file = "data/external/HUMAN_9606_idmapping.dat"
    kinhub_id_map_file = "data/external/kinhub-families.tsv"
    kinome_file = "data/human-kinome.txt"
    ktable = pssms.read_kin_sub_data("data/psiteplus-kinase-substrates.tsv")
    with open(uniprot_id_map_file) as h:
        uniprot_id_map = read_uniprot_id_map(h)
    with open(kinhub_id_map_file) as h:
        kinhub_id_map = read_kinhub_id_map(h)
    with open(sdr_file) as h:
        sdr_matches = read_sdr_file(h, ktable, uniprot_id_map, kinhub_id_map)
