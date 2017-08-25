#Test the Imputation method described in Macau: Scalable Bayesian Multi-relational
#Factorization with Side Information using MCMC
# Might also want check out: https://github.com/jaak-s/macau

import macau
from scipy import sparse
import mmap
from openpyxl import load_workbook
import numpy as np
import pandas as pd


# Here I read in the data
def read_kin():

    wb = load_workbook("/home/borgthor/Downloads/Supplementary_Table_2.xlsx", read_only=True)

    data= wb['KSEA activities']

    ksea = []
    i = 0
    for row in data.rows:
        
        if i > 1:
            j =0
            kinase = []
            for cell in row:
                if j > 0:
                    if cell.value == 'NA':
                        kinase.append(np.nan)

                    else:
                        kinase.append(float(cell.value))
                j = j+1
            
            ksea.append(kinase)
        i = i+1

    return ksea


kin = read_kin()
n = np.matrix(kin)
#Here I factorize the kinase activity values with the Macau method
#not that this needs some tinkering, for instance how to best incorporate
# prior knowledge
result = macau.macau(Y = sparse.csr_matrix(n),
                     Ytest = 0.1,
                     num_latent = 32,
                     precision = "adaptive",
                     burnin = 500,
                     nsamples = 3500,
                     save_prefix = 'Ochao')
#Here I recunstruct the matrix from the model
meanvalue = np.loadtxt("Ochao-meanvalue.csv").tolist()
U = np.loadtxt("Ochao-sample%d-U1-latents.csv" % 1, delimiter=",")
V = np.loadtxt("Ochao-sample%d-U2-latents.csv" % 1, delimiter=",")

kinhat = U.transpose().dot(V) + meanvalue

#The model is then tested... check the script imputation_test

f = open("/home/borgthor/impute_kinases.txt","w")

for i in kinhat:
    s = ""
    for j in i:
        s = s+ str(j)+"\t"

    f.write(s[0:(len(s)-1)] + "\n")

f.close()
    
        

    
