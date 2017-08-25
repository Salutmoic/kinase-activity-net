import pssms
from pssm import *
import numpy as np
from math import log

AMINO_ACIDS = {"G": 0, "A": 1, "V": 2, "L": 3, "I": 4, "P": 5, "F": 6, "Y": 7, "W": 8, "S": 9, "T": 10, "C": 11, "M": 12, "N": 13, "Q": 14,"K": 15, "R": 16, "H": 17, "D": 18, "E": 19}

def shannon_entropy(pssm):
    
    dim = pssm.shape
    
    entropy = [0]*dim[1]
    for i in range(dim[1]):
        for j in range(dim[0]):
            if(pssm[j,i]) == 0:
                continue
            if i == 7:
                continue
            print(pssm[j,i])
            print(pssm[j,i])
            entropy[i] = entropy[i] + 2**pssm[j,i]*pssm[j,i]
    entropy[i] = (-1)*entropy[i]
    return entropy

def column_entropy(pssm):
    if pssm is None:
        return(-100)
    if np.all(pssm==0):
        return(-100)
    entropies = shannon_entropy(pssm)
    print(entropies)
    return np.mean(entropies)

def pssm_prop(seqs):
    if len(seqs)== 0:
        return(np.zeros((20,2)))
    #this function returns simple proportinal pssm
    pssm = np.zeros((20, len(seqs[0])))
    for seq in seqs:
        for i,aa in enumerate(seq):
            if aa not in AMINO_ACIDS:
                continue
            pssm[AMINO_ACIDS[aa],i] += 1/len(seqs)
    return(pssm) 	
			
		

if __name__ == "__main__":
    ktable = pssms.read_kin_sub_data("data/psiteplus-kinase-substrates.tsv")
    aa_freqs = pssms.read_aa_freqs("data/aa-freqs.tsv")
    psites = read_phosphosites("data/pride-phosphosites.tsv")
    kinases = read_kinase_file("data/human-kinome.txt")
    # phosfun, med_phosfun = read_phosfun_file("data/phosfun-alt.tsv")
    phosfun_st, med_phosfun_st = read_phosfun_file("data/phosfun-ST.tsv")
    phosfun_y, med_phosfun_y = read_phosfun_file("data/phosfun-Y.tsv")
    reg_sites = read_reg_sites_file("data/psiteplus-reg-sites.tsv")
    entropies = {}
    pssm_div = {}    
    for kinase in kinases:
        if kinase not in ktable:
            continue    
        kin_motif_seqs = []
        for substrate, pos, res, seq in ktable[kinase]:
            
            if substrate not in psites:
                continue
            for pos_b, seq_b, res_b in psites[substrate]:
                if pos == pos_b:
                    kin_motif_seqs.append(seq_b)
 
        pssm = pssms.calc_pssm(kin_motif_seqs, aa_freqs)
        naivepssm = pssm_prop(kin_motif_seqs)
        entropies[kinase] = column_entropy(pssm)
		
        if pssm is not None:
            
            print(kinase)
           
            pssm_min = np.min(pssm, axis=0)
            pssm_max = np.max(pssm, axis=0)
            pssm_div[kinase] = sum( pssm_max)-sum( pssm_min)
    
    with open("out/pssm-motifness-score.tsv","w") as f:
        for kinase in pssm_div:
            
            f.write(kinase + "\t"+ str(pssm_div[kinase]) + "\n")
    f.close()        
	
    with open("out/kinase-entropy-score.tsv","w") as f:
        for kinase in entropies:
            f.write(kinase +"\t"+ str(entropies[kinase]) + "\n")
    f.close()


	
      
