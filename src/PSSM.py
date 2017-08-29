
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
    #
    #
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
                
            

        cmat = cmat + 0.01*len(seqs)
        cmat = cmat/len(seqs)
        pmms[i] = cmat

    return pmms



def subseqs():
    #This function extracts peptides surrounding phospho sites
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
    #creates distribution of 
    dist = {}

    for i in pmms:
        d = []
        for j in range(0,1000):

            rand =  random.randint(0,(len(seqs)-1))
            seq = seqs[rand]

            d.append(score(seq,pmms[i]))
        dist[i] = d

    return dist
        

    
    
    


t,ks = read_data()

pmms = pmms(ks)

seqs = subseqs()

dist = dist(seqs,pmms)
np.save('/home/borgthor/Kinase_activities/cons_matrices', pmms) 
np.save('/home/borgthor/Kinase_activities/score_distributions', dist) 
    
read_dictionary = np.load('/home/borgthor/Kinase_activities/cons_matrices.npy').item()    
    

        








    
