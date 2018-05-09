import networkx as nx
import sys
import numpy as np
import os.path
from sklearn.neighbors import LocalOutlierFactor
from sklearn.cluster import KMeans

def read_network(path,t):
    network = nx.Graph()
    
    with open(path,"r") as f:
        head = f.readline()
        head = head.split("\t")
        head[len(head)-1] = head[len(head)-1].replace("\n","")
        print(head)
        inda = head.index('SymbolA')
        indb = head.index('SymbolB')
        
        i = 0        
        for line in f:          
            edge = line.split("\t") 
            
            if len(edge) < indb:
                continue
            if len(edge) > 2:
                pind = head.index('p(Interaction)')
                if float(edge[pind]) > t:
                    network.add_edge(edge[inda],edge[indb].replace("\n",""))
            else:
                
                network.add_edge(edge[inda],edge[indb].replace("\n",""))
            i += 1
    
    return network

def read_pred(path,prior,threshold):
    network = {}
    with open(path,"r") as f:
        head = f.readline()
        for line in f:
                        
            edges = line.split("\t")
            if edges[2].replace("\n","") == "NA":
                continue
            if float(edges[2].replace("\n","")) > threshold and edges[0] in prior and edges[1] in prior:
                if edges[0] in network:
                    network[edges[0]] =  network[edges[0]] + [edges[1]]
                else:
                    network[edges[0]] = [edges[1]]
    return network

def calc_dists(subs,prior):
    
    if len(subs) <= 2:
        return({'sub':0 })
    avg_dists = {}
    dist_matrix = make_dist_matrix(prior,subs)
    
    for i in range(len(dist_matrix)):
        avg_dists[subs[i]] = np.mean(dist_matrix[i])
    return avg_dists

def calc_bg_dist(prior):
	
    nodes = prior.nodes()
    sources = np.random.choice(nodes, 100,replace = False)
    targets = np.random.choice(nodes, 100,replace = False)
    dists = []

    for node1 in sources:
        for node2 in targets:
            if node1 == node2:
                continue  
            if not nx.has_path(prior, node1,node2):
                
                continue
            d= len(nx.shortest_path(prior,node1,node2))-1
            dists.append(d)
    return dists

def calc_cl_dist(prior,kinases):
    
    cl_dists = []
    for kinase in kinases:
        if kinase not in prior:
            continue
        subs = list(prior[kinase])
       
        if len(subs) <=1:
            continue 
        dist_matrix = make_dist_matrix(prior,subs)    
        kmeans = KMeans(n_clusters=1, random_state=0).fit(dist_matrix)        
        inertia = kmeans.inertia_
        
        for row in dist_matrix:
            d_new= dist_matrix[:]
            d_new.remove(row)
            knew = KMeans(n_clusters=1, random_state=0).fit(d_new)
            n_inert = inertia-knew.inertia_
            cl_dists.append(inertia-n_inert)
    return cl_dists
            
def check_lof(subs,kinase,prior):
    
    if len(subs) <= 2:
        return( [0])
    
    dist_matrix = make_dist_matrix(prior,subs)
    
    lof = LocalOutlierFactor(n_neighbors=int(len(subs) * 0.75))  
    y = lof.fit_predict(dist_matrix)
    print(y)
    results = {}
    for i in range(len(y)):
        if y[i] == 1:
            results[subs[i]] = 0
        else:
            results[subs[i]] = 1   
    return(results)

def make_dist_matrix(prior,subs):
    
    dist_matrix = []

    for sub1 in subs:
        row = []
        if sub1 not in prior:
            row = [10]*len(subs)
            continue
        for sub2 in subs:
            if sub2 not in prior or not nx.has_path(prior, sub1,sub2):
                row.append(10)
            elif sub2 == sub1:
                continue
            else:
                row.append(len(nx.shortest_path(prior,sub1,sub2))-1)
        dist_matrix.append(row)

    return(dist_matrix)
	

def check_cluster(subs,prior, kinase):
    
    
    if len(subs) <= 2:
        return {"substrate":0}
    dist_matrix = make_dist_matrix(prior,subs)
    kmeans = KMeans(n_clusters=1, random_state=0).fit(dist_matrix)
    inertia = kmeans.inertia_
    
    results = {}
    for sub in subs:
        
        d_new = dist_matrix[:]
        del d_new[subs.index(sub)]
        new_kmeans = KMeans(n_clusters=1, random_state=0).fit(d_new)
        new_inertia = new_kmeans.inertia_
        results[sub] = inertia - new_inertia
    return(results) 

def read_kin_file(f):
    kins = []
    
    with open(f,"r") as f:
        for line in f:
            kin = line.replace("\n","")
            kins.append(kin)
    return kins
 
         
def prune_network(pred_network,avg_dist,bg,prior,method):
    
    if method == 'bg_dist':
        allowed_dist = np.percentile(bg,50)

    elif method == 'lof':
        allowed_dist = 0

    elif method == "cluster":
        kinases = read_kin_file("data/human-kinome.txt")
        cl_dist = calc_cl_dist(prior,kinases)
        allowed_dist = np.percentile(cl_dist,80)
 
    for kinase in pred_network:
        if method!= 'lof' and method != "cluster":
            dists =  avg_dist[kinase]  
            print(dists)        
        if method == 'sub_dist': 
            allowed_dist = np.mean(list(dists.values())) + np.std(list(dists.values()))
        
        if method== 'lof':
            dists = check_lof(pred_network[kinase],kinase,prior) 

        if method == 'cluster':
            dists = check_cluster(pred_network[kinase],prior, kinase)    
    
        for substrates in dists:
            if dists[substrates] > allowed_dist:
                  
                pred_network[kinase].remove(substrates)
                
    return pred_network     
    
if __name__ == "__main__":
    
    threshold = float(sys.argv[3])
    dset = sys.argv[2]
    prior_network = read_network(dset,float(sys.argv[5]))
    bg_dists = calc_bg_dist(prior_network)
    method = sys.argv[4]
    pred_network = read_pred(sys.argv[1],prior_network,threshold)    
    clust_dist = calc_cl_dist(prior_network,list(pred_network.keys()))    
	
    avg_dist = {}	
    
    for kinase in pred_network:
        pred_subs = pred_network[kinase]
        
        if len(pred_subs) == 0:
            del pred_network[kinase]
            continue 
                   
        dists = calc_dists(pred_subs,prior_network)
        avg_dist[kinase] = dists 
    pruned_network = prune_network(pred_network,avg_dist,bg_dists,prior_network,method)             

    with open(sys.argv[6], "w") as f:
        for kinase in pruned_network:
            for substrate in pruned_network[kinase]:
                f.write(kinase+"\t"+substrate + "\n")                  	

