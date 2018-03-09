#code that writes training data set
import scipy
import numpy as np
import sys
from random import *


def pos_set(f,data):
    #returns a set of known positives
    pos = {}
    ypos = {}
    datadict = prediction_matrix(data)
    
    with open(f,"r") as f1:
        first_line = f1.readline()
        for line in f1:
            if 'NA' not in line:            
                line = line.split("\t")
                line[1] = line[1].replace("\n","")
                if (line[0],line[1]) in datadict:
                    pos[line[0],line[1]] = datadict[line[0],line[1]]
                    ypos[line[0],line[1]] = 1
    return pos,ypos

def train_random(f1,data):
    #uses random set of edges as a negative set for training
    train, y = pos_set(f1,data)
    datadict = prediction_matrix(data)
    n = len(y)

    samp =  sample(range(len(datadict)), n)   
    keys = list(datadict.keys())

    for i in samp:
        if keys[i] in train:
            o = 0
            while o == 0:
                j = int(random() *len(datadict))
                if keys[j] not in train:
                    o = 1
                    train[keys[j]] = datadict[keys[j]]
                    y[keys[j]] = 0
        else:
            train[keys[i]] = datadict[keys[i]]
            y[keys[i]] = 0 

    return train, y

def fix(x):
    if x == "NA":
        return np.nan
    else:
        return float(x) 

def train_prior(f1,f2,data):
    #uses known negatives as a negative set for training
    train, y = pos_set(f1,data)
    datadict = prediction_matrix(data)
    print(len(datadict))
    n = len(y)
    print(n)
    truenegs = []

    with open(f2,"r") as negs:
        first_line = negs.readline()
        for line in negs:
            
            line = line.split("\t")
            if (line[0],line[1].replace("\n","")) in datadict:
                truenegs.append((line[0],line[1].replace("\n","")))
     
     
    samp = sample(range(len(truenegs)), n)
    print(len(samp))    
    for i in samp:     
        train[truenegs[i]] = datadict[truenegs[i]]
        y[truenegs[i]] = 0
    print(len(train))
    return train, y
    
def prediction_matrix(data):
    #makes matrix with pssm values string score and corr
    v = {}
    with open(data,"r") as d:
        first_line = d.readline()
        for line in d:
            if 'NA' not in line:
                line = line.split("\t")
                line[len(line)-1] = line[len(line)-1].replace("\n","")
              
                v[line[0],line[1]] = [fix(x) for x in line[2:]]
                
    return v
