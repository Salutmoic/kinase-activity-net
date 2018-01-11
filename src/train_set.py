#code that writes training data set
import scipy
import numpy as np
import sys
from random import *


def pos_set(f,data):
    pos = {}
    ypos = []
    datadict = prediction_matrix(data)
    
    with open(f,"r") as f1:
        first_line = f1.readline()
        for line in f1:
            if 'NA' not in line:
                line = line.split("\t")
                line[1] = line[1].replace("\n","")
                if (line[0],line[1]) in datadict:
                    pos[line[0],line[1]] = datadict[line[0],line[1]]
    ypos = [1]*len(pos)
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
        else:
            train[keys[i]] = datadict[keys[i]]
    y = y + [0]*n

    return train, y

def train_prior(f1,f2,data):
    #uses known negatives as a negative set for training
    train, y = pos_set(f1,data)
    datadict = prediction_matrix(data)
    n = len(y)

    truenegs = []

    with open(f2,"r") as negs:
        first_line = negs.readline()
        for line in negs:
            if 'NA' not in line:
                line = line.split("\t")
                truenegs.append((line[0],line[1].replace("\n","")))

    samp = sample(range(len(truenegs)), n)

    for i in samp:
        train[truenegs[i]] = datadict[truenegs[i]]
    y = y + [0]*n
    return train, y
    
def prediction_matrix(data):
    #makes matrix with pssm values string score and corr
    v = {}
    with open(data,"r") as d:
        first_line = d.readline()
        for line in d:
            if 'NA' not in line:
                line = line.split("\t")
                v[line[0],line[1]] = [float(x) for x in line[2:]]
    return v
