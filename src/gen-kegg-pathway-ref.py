import os


def parse_id_map():
    id_map_file = "data/uniprot-id-map.tsv"
    id_map = dict()
    with open(id_map_file, 'r') as h:
        for line in h:
            uniprot, gene_name = line.strip().split()
            id_map[uniprot] = gene_name
    return id_map


def get_path_members(path_rel_handle, path_name, id_map):
    path_members = set()
    for line in path_rel_handle:
        line_split = line.strip().split()
        uniprot1 = line_split[0]
        if uniprot1 not in id_map:
            continue
        uniprot2 = line_split[1]
        if uniprot2 not in id_map:
            continue
        path_members.add((path_name, id_map[uniprot1]))
        path_members.add((path_name, id_map[uniprot2]))
    return path_members


def parse_kegg_paths(id_map):
    kegg_rel_path = "data/external/kegg-relationships/"
    path_files = os.listdir(kegg_rel_path)
    path_prots = set()
    for f in path_files:
        path_name = os.path.splitext(f)[0]
        with open(kegg_rel_path+f, 'r') as h:
            path_members = get_path_members(h, path_name, id_map)
            path_prots.update(path_members)
    return(path_prots)


def print_paths(path_prots):
    for path_prot in path_prots:
        print("\t".join(path_prot))


if __name__ == "__main__":
    id_map = parse_id_map()
    path_prots = parse_kegg_paths(id_map)
    print_paths(path_prots)
    
