import numpy as np
import random
from os import listdir
from os.path import isfile, join, walk, basename, splitext
from math import log
import sys
import itertools


def read_data():
    #Reads kinase-substrate table where substrate and kinase are found in 
    f = open("data/reduced_kinase_table.tsv","r")

    

    columns = f.readline()

    columns =columns.split("\t")
    tside = columns.index("SITE_...7_AA")
    ktable = {}
    for line in f:
        
        line = line.split("\t")
        if line[0] in ktable:
            ktable[line[0]].append(line[tside])
        else:
            ktable[line[0]] = [line[tside]]
       
    return tside,ktable

def amino_acids():
    #returns a dictionary containing a list of the common 20 amino acids 
    aminos = {"G":0,"A":1,"V":2,"L":3,"I":4,"P":5,"F":6,"Y":7,"W":8,"S":9,"T":10,"C":11,"M":12,"N":13,"Q":14,"K":15,"R":16,"H":17,"D":18,"E":19}

    return aminos


def read_aa_freqs():
    aa_freqs = dict()
    with open("data/aa-freqs.tsv", "r") as h:
        for line in h:
            aa, freq = line.strip().split()
            aa_freqs[aa] = float(freq)
    return aa_freqs


def pmms(ks):
    # returns a dictionary where keys are kinases and values are pmms matrices made from the known substrates
    # for each kinase
    
    aminos = amino_acids()
    pmms = {}

    aa_freqs = read_aa_freqs()
    kinases = ks.keys()

    for i in kinases:
        seqs = ks[i]
        cmat = np.zeros((15,20))
        for j in seqs:
            j = j.upper()
            for k in range(0,len(j)):
                if j[k] in aminos:
                    cmat[k][aminos[j[k]]] = cmat[k][aminos[j[k]]] + 1.0
            
        pc = 0.1
        cmat = cmat + pc #(float(1)/(len(seqs)+2))*len(seqs)
        cmat = cmat/(len(seqs)+pc*20.0)

        for k in range(15):
            for aa in aa_freqs:
                cmat[k][aminos[aa]] = log(cmat[k][aminos[aa]]/aa_freqs[aa])/log(2.0)
        
        pmms[i] = cmat

    return pmms



def subseqs():
    #This function extracts peptides surrounding phospho sites found in the phosphosite plus
    # human_kinase_table
    f = open("data/human_kinase_table.tsv","r")
    
    columns = f.readline()

    columns =columns.split("\t")
    seqind = columns.index("SITE_...7_AA")
    seqs = []
    
    for line in f:
        
        line = line.split("\t")
        seqs.append(line[seqind])


    return seqs


def score(seq,pmms,ks):
    #scores peptites with pssm
    aminos = amino_acids()
    score = 1
    for i in range(0,len(seq)):
        if seq[i] != "_":
            score = score + pmms[i][aminos[seq[i].upper()]]
        else:
            score = score + np.min(pmms[i]) #(float(1)/(len(ks)+2))

    return score
        
        
def dist(seqs,pmms,ks):
    #creates distribution of pmms scores
    #scores are calculated by scoring a random phosphosite relative to to kinase i
    dist = {}

    for i in pmms:
        d = []
        for j in range(0,1000):

            rand =  random.randint(0,(len(seqs)-1))
            seq = seqs[rand]

            d.append(score(seq,pmms[i],ks[i]))
        dist[i] = d

    return dist
        
def kin_phosphosites(ks):
    #returns dictionary with every known phosphosite found on each kinase in the
    #Ochoa et al. dataset
    
    f =  open("data/phosphosites_reduced.tsv")
    psites = {}
    columns = f.readline()

    columns =columns.split("\t")
    seqind = columns.index("SITE_...7_AA")
    
    for line in f:
        
        line = line.split("\t")
        if line[0] in ks.keys():
            if line[0] not in psites:
                psites[line[0]] = [line[seqind]]
            else:
                psites[line[0]].append(line[seqind])

    return psites

def score_kin_pairs(ks,pmms,netdir):
    # loads a netwoek of kinase-> kinase interactions and scores them with the pssm
    # it returns a dictionary with kinasa A and kinase B and list of scores for each phosphosite found
    # on B

    prots = set()
    with open(netdir) as h:
        for line in h:
            prot, _, _ = line.split()
            prots.add(prot)
    pairs = itertools.product(prots, repeat=2)
    from pprint import pprint

    psites = kin_phosphosites(ks)
    
    scores = {}
    for pair in pairs:
        A, B = pair
        if A == B:
            continue
        sites = psites[B]
        for i in sites:
            sc = score(i,pmms[A],ks[A])
            if (A,B) in scores.keys():
                scores[(A,B)].append(sc)
               
            else:
                scores[(A,B)] = [sc]
    return scores

def assess_edges(scores,dist):
    #Finds the highest scoreing phosphosite for ech A,B pair and compares it to the
    #distribution generated for B (A-> B) from the dist() function
    zscores = {}

    for i in scores.keys():
        m = max(scores[i])
        z = (m-np.mean(dist[i[0]]))/np.std(dist[i[0]])
        zscores[i] = z

    return zscores

def score_network(filename):
    #This takes network files and scores all edges according to pssm

    t,ks = read_data()
    pmm = pmms(ks)

    seqs = subseqs()

    d = dist(seqs,pmm,ks)
    

    v = open("out/kinase_distributions.tsv","w")
    for i in d:

        line = i + "\t"
        for j in d[i]:
            line = line + str(j) + "\t"

        v.write(line + "\n")     
    v.close()
    
    scores = score_kin_pairs(ks,pmm,filename)
    out_file = "out/" + splitext(basename(filename))[0] + "-pssm.tsv"
    v = open(out_file, "w")
    for i in scores:
        
        m_score = max(scores[i])
        v.write(str(i[0])+"\t"+str(i[1])+ "\t"+ str(m_score)+"\n")

    v.close()

    
    
def known_scores(pmms,ks):

    
    v = open("out/known_kinase_psite-score","w")
    f = open("data/human_kinase_table","r")
    columns = f.readline()
    columns =columns.split("\t")
    seqind = columns.index("SITE_...7_AA")

    for line in f:
        line = line.split("\t")
        kinase = line[0] 
        psite = line[seqind]
        sc = score(psite,pmms[kinase],ks[kinase])
        v.write(kinase+"\t"+psite+"\t"+str(sc)+ "\n")

    v.close()
    f.close()
    return    
   

if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("USAGE: <script> DATA_FILE")
    filename = sys.argv[1]
    score_network(filename)
