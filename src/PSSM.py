import numpy as np
import random
from os import listdir
from os.path import isfile, join, walk


def read_data():
    #Reads kinase-substrate table where substrate and kinase are found in 
    f = open("data/reduced_kinase_table","r")

    

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


def pmms(ks):
    # returns a dictionary where keys are kinases and values are pmms matrices made from the known substrates
    # for each kinase
    
    aminos = amino_acids()
    pmms = {}

    kinases = ks.keys()

    for i in kinases:
        seqs = ks[i]
        cmat = np.zeros((15,20))
        for j in seqs:
            j = j.upper()
            for k in range(0,len(j)):
                if j[k] in aminos:
                    cmat[k][aminos[j[k]]] = cmat[k][aminos[j[k]]] +1
                
            
        #I use 0.01 as pseudocount 
        cmat = cmat + (float(1)/(len(seqs)+2))*len(seqs)
        cmat = cmat/len(seqs)
        pmms[i] = cmat

    return pmms



def subseqs():
    #This function extracts peptides surrounding phospho sites found in the phosphosite plus
    # human_kinase_table
    f = open("data/human_kinase_table","r")
    
    columns = f.readline()

    columns =columns.split("\t")
    seqind = columns.index("SITE_...7_AA")
    seqs = []
    
    for line in f:
        
        line = line.split("\t")
        seqs.append(line[seqind])


    return seqs


def score(seq,pmms):
    #scores peptites with pssm
    aminos = amino_acids()
    score = 1
    for i in range(0,len(seq)):
        if seq[i] != "_":
            score = score*pmms[i][aminos[seq[i].upper()]]
        else:
           score = score*(float(1)/(len(seqs)+2))

    return score
        
        
def dist(seqs,pmms):
    #creates distribution of pmms scores
    #scores are calculated by scoring a random phosphosite relative to to kinase i
    dist = {}

    for i in pmms:
        d = []
        for j in range(0,1000):

            rand =  random.randint(0,(len(seqs)-1))
            seq = seqs[rand]

            d.append(score(seq,pmms[i]))
        dist[i] = d

    return dist
        
def kin_phosphosites(ks):
    #returns dictionary with every known phosphosite found on each kinase in the
    #Ochoa et al. dataset
    
    f =  open("data/phosphosites_reduced.txt")
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

    #f is some imaginary network
    f =  open(netdir)
    psites = kin_phosphosites(ks)
    
    scores = {}
    for line in f:
        line = line.split("\t")
        A = line[0]
        B= line[1][0:(len(line[1])-1)]
        
        sites = psites[B]
        for i in sites:
            sc = score(i,pmms[A])
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

    d = dist(seqs,pmm)
    files = listdir(path)

    v = open("kinase_distributions.tsv","w")
    for i in dist:

        line = i + "\t"
        for j in dist[i]:
            line = line + str(j) + "\t"

        distf.write(line + "\n")     
    v.close()
    
    scores = score_kin_pairs(ks,pmm,filename)
    v = open("kinase_kinase_scores.tsv","w")
    for i in scores:
        
        m_score = max(scores[i])
        n.write(str(i[0])+"\t"+str(i[1])+ "\t"+ str(m_score)+"\n")

    v.close()
    
    
def known_scores(pmms,ks):

    
    v = open("known_kinase_psite-score","w")
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
   


score_network(path)

    

        








    
