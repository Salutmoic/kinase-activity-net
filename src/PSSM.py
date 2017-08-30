import numpy as np
import random

def read_data():
    #Reads kinase-substrate table where substrate and kinase are found in 
    f = open("/home/borgthor/Kinase_activities/data/reduced_kinase_table","r")

    

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
        cmat = cmat + 0.01*len(seqs)
        cmat = cmat/len(seqs)
        pmms[i] = cmat

    return pmms



def subseqs():
    #This function extracts peptides surrounding phospho sites found in the phosphosite plus
    # human_kinase_table
    f = open("/home/borgthor/Kinase_activities/data/human_kinase_table","r")
    
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
           score = score*0.01 

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
    
    f =  open("/home/borgthor/phosphosites_reduced.txt")
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

def score_kin_pairs(ks,pmms):
    # loads a netwoek of kinase-> kinase interactions and scores them with the pssm
    # it returns a dictionary with kinasa A and kinase B and list of scores for each phosphosite found
    # on B

    #f is some imaginary network
    f =  open("/home/borgthor/network.txt")
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
    
    
               
t,ks = read_data()

pmms = pmms(ks)

seqs = subseqs()

dist = dist(seqs,pmms)

psites = kin_phosphosites(ks)

scores = score_kin_pairs(ks,pmms)

zscores = assess_edges(scores,dist)

np.save('/home/borgthor/Kinase_activities/cons_matrices', pmms) 
np.save('/home/borgthor/Kinase_activities/score_distributions', dist)
np.save('/home/borgthor/Kinase_activities/kinase pairs score', score) 
np.save('/home/borgthor/Kinase_activities/phosphosites per kinases', psites) 
np.save('/home/borgthor/Kinase_activities/pcals for each kinase pair', pvals) 

  
    

        








    
