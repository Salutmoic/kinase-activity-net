import networkx as nx
import sys
import numpy as np
import os.path

def read_network(path):
    network = nx.Graph()
    
    with open(path,"r") as f:
        head = f.readline()
        head = head.split("\t")
        inda = head.index('SymbolA')
        indb = head.index('SymbolB')
        i = 0        
        for line in f:          
            edge = line.split("\t") 
            
            if(len(edge) < indb):
                continue
            network.add_edge(edge[inda],edge[indb])
            i += 1
    print(i)
    return network

def read_pred(path,prior,threshold):
    network = {}
    with open(path,"r") as f:
        for line in f:
                        
            edges = line.split("\t")
            
            if float(edges[2]) > threshold and edges[0] in prior and edges[1] in prior:
                if edges[0] in network:
                    network[edges[0]] =  network[edges[0]] + [edges[1]]
                else:
                    network[edges[0]] = [edges[1]]
    return network

def calc_dists(subs,prior):
    
    avg_dists = {}
    for i in subs:
        if i not in prior:
            avg_dist[i] = 10
            continue
        dist = []
        for j in subs:
            if j not in prior or not nx.has_path(prior, i,j):
                dist.append(numpy.nan)
            elif i == j:
               
                continue                 
            else:
                d= len(nx.shortest_path(prior,i,j))-1
                dist.append(d)
        
        if len(dist) == 0:
            dist = [0]       
        avg_dists[i] = np.mean(dist)
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

def prune_network(pred_network,avg_dist,bg):
    
    allowed_dist = np.percentile(bg,15)
    for kinase in pred_network:
        dists =  avg_dist[kinase]  
               
        for substrates in dists:
            if dists[substrates] > allowed_dist:
                print(pred_network[kinase])   
                pred_network[kinase].remove(substrates)
                print(pred_network[kinase])
    return pred_network     
    
if __name__ == "__main__":
    
    threshold = float(sys.argv[3])
    dset = sys.argv[2]
    prior_network = read_network(dset)
    bg_dists = calc_bg_dist(prior_network)
   
    pred_network = read_pred(sys.argv[1],prior_network,threshold)    
	
    avg_dist = {}
    print(max(bg_dists))	
  	  
    for kinase in pred_network:
        pred_subs = pred_network[kinase]
        
        if len(pred_subs) == 0:
            del pred_network[kinase]
            continue 
                   
        dists = calc_dists(pred_subs,prior_network)
        avg_dist[kinase] = dists 
    pruned_network = prune_network(pred_network,avg_dist,bg_dists)             
    
    with open(sys.argv[3], "w") as f:
        for kinase in pruned_network:
            for substrate in pruned_network[kinase]:
                f.write(kinase+"\t"+substrate + "\n")                  	

