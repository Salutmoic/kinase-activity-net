import scipy
import numpy as np
import sklearn
from sklearn import svm, linear_model
import sys
from sklearn.naive_bayes import GaussianNB
from sklearn.neighbors.nearest_centroid import NearestCentroid
from random import *
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score

def traindata(file1,file2,data,threshold):
# Creates matrix (list of list) containing true pos and true negatives
# for training and predictions, f1 is true positives and f2 true negatives
    # x is predictive variables (ppsm, string,cor) and y outcomes
    pos = {}
    neg={}
    ypos = []
    yneg = []
    datadict = prediction_matrix(data)[0]
    with open(file2,"r") as f2:
        first_line = f2.readline()
        for line in f2:
            if 'NA' not in line:
                line.replace("\n","")
                line = line.split("\t")
                line[1] = line[1].replace('\n','')
                if (line[0],line[1]) in datadict:
                    neg[line[0],line[1]] = datadict[line[0],line[1]]
    yneg = [0] *len(neg)

    with open(file1,"r") as f1:
        first_line = f1.readline()
        for line in f1:
            if 'NA' not in line:
                line = line.split("\t")
                line[1] = line[1].replace("\n","")
                if (line[0],line[1]) in datadict:
                    
                    pos[line[0],line[1]] = datadict[line[0],line[1]]
    ypos = [1] *len(pos) 

    whole = pos
    y = ypos
    k = list(pos.keys())
    
    for i in k:
        rand = random()
        if rand > threshold:
            if (i[1],i[0]) not in neg and (i[1],i[0]) not in pos:
                whole[i[1],i[0]] = datadict[i[1],i[0]]
                y.append(0)

    
    m = shuffle(neg.keys())

    while len(whole) < 2*len(pos):
        key = m.pop()
        whole[key] = neg[key]
        y.append(0)

    
    return whole,y
    


def prediction_matrix(data):
    #makes matrix with pssm values string score and corr
    v = {}
    with open(data,"r") as d:
        first_line = d.readline()
        for line in d:
            if 'NA' not in line:
                line = line.split("\t")

                ppairs.append([line[0],line[1]])
                v[line[0],line[1]] = [float(x) for x in line[2:]]
    return v

  

def find_centers(test,results):
#works with nearestneighbour this function finds the centers of
#the two clusters
    center1 = [[],[],[],[]]
    center2 = [[],[],[],[]]

    for i in range(0,len(test)):
        if results[i] == 0:
            center1[0].append(test[i][0])
            center1[1].append(test[i][1])
            center1[2].append(test[i][2])
            center1[3].append(test[i][3])
        else:
            center2[0].append(test[i][0])
            center2[1].append(test[i][1])
            center2[2].append(test[i][2])
            center2[3].append(test[i][3])
    if len(center1[0]) == 0 or len(center2[0]) == 0:
        return [[],[]]
    c1 = [np.mean(center1[0]),np.mean(center1[1]),np.mean(center1[2]),np.mean(center1[3])]
    c2 = [np.mean(center2[0]),np.mean(center2[1]),np.mean(center2[2]),np.mean(center2[3])]

    return [c1,c2]

def dist(t,c):
#calculates eucledian distance between point x and y
    dist = float(abs(t[0]-c[0])**2+abs(t[1]-c[1])**2+abs(t[2]-c[2])**2 + abs(t[3]-c[3])**2)**0.5

    return dist

def probability(c1,c2,test,results):
#lukewarm function that calculates the probability that given value is in 
#cluster x rather than y        
    probability = []
    for i in range(0,len(results)):
        if results[i] == 0:
            p = dist(test[i],c2)/(dist(test[i],c2)+dist(test[i],c1))
            probability.append([p,1-p])
        else:
            p = dist(test[i],c1)/(dist(test[i],c2)+dist(test[i],c1))
            probability.append([p,1-p])
    return probability

 
if __name__ == "__main__":

    truePos = sys.argv[1]
    trueNeg = sys.argv[2]
    data = sys.argv[3]
    method = sys.argv[4]
    threshold = float(sys.argv[5])

    testdict = prediction_matrix(data)
    test = list(testdict.values())
  
    traindict,outcomes = traindata(truePos,trueNeg,data,threshold)
    train = list(traindict.values())
    
    if method == "NB":
        model = GaussianNB()
        
    if method == "SVM":
        model = svm.SVC(probability = True)
       
    if method == "kmeans":
        model = NearestCentroid()
        model.fit(train, outcomes)
        prediction = model.predict(test)

        centers = find_centers(test,results)
        c1 = centers[0]
        c2 = centers[1]

        if len(c1) > 0 and len(c2) > 0:
        
            probability = probability(c1,c2,test,results)

    if method == "logreg":
        model = linear_model.LogisticRegression()
        

    if method == "randfor":
        model = RandomForestClassifier()

    if model != "kmeans":
        
        prediction = model.fit(train,outcomes).predict(test)
        probability = model.fit(train,outcomes).predict_proba(test)
    
    with open(sys.argv[7],"w") as crossval:t 
        # repeated cross validation, is this something along the lines of what you were thinking 
        for i in range(50):
            scores = cross_val_score(model, test, prediction, scoring = 'roc_auc', cv=5)
            for score in scores:
                crossval.write(score + "\t")
            crossval.write("\n")
        
    with open(sys.argv[6],"w") as output:
        for i in range(len(testdict.keys())):
            s = str.join("\t",[str(x) for x in list(testdict.keys())[i]])
            s = s + "\t"+ str(prediction[i]) +"\t"+ str(probability[i][1])

            output.write(s+"\n")
                                                                      
                                                                                                                                                                                          1,1           Top
