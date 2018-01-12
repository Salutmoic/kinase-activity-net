import scipy
import numpy as np
import sklearn
from sklearn import svm, linear_model, preprocessing
import sys
from sklearn.naive_bayes import GaussianNB
from sklearn.neighbors.nearest_centroid import NearestCentroid
from random import *
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
from train_set import *
  
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
    train_method = sys.argv[5]
    
    min_max_scaler = preprocessing.MinMaxScaler()

    testdict = prediction_matrix(data)
    test = list(testdict.values())
    min_max_scaler.fit(test)

    test_scale = min_max_scaler.transform(test)

    if train_method == "random":
        traindict,outcomes = train_random(truePos,data)
    else:
        traindict,outcomes = train_prior(truePos,trueNeg,data)
    train = list(traindict.values())
    min_max_scaler.fit(train)
    train_scale = min_max_scaler.transform(train)

    
    if method == "NB":
        model = GaussianNB()
        
    if method == "SVM":
        model = svm.SVC(probability = True)
       
    if method == "kmeans":
        model = NearestCentroid()
        model.fit(train_scale, outcomes)
        prediction = model.predict(test_scale)

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
        
        prediction = model.fit(train_scale,outcomes).predict(test_scale)
        probability = model.fit(train_scale,outcomes).predict_proba(test_scale)
    
    with open(sys.argv[7],"w") as crossval: 
        # repeated cross validation, is this something along the lines of what you were thinking 
        for i in range(50):
            if train_method == "random":
                traindict,outcomes = train_random(truePos,data)
            else:
               traindict,outcomes = train_prior(truePos,trueNeg,data)
            train = list(traindict.values())
            min_max_scaler = preprocessing.MinMaxScaler()
            min_max_scaler.fit(train)
            train_scale = min_max_scaler.transform(train)
          
            scores = cross_val_score(model, train_scale, outcomes, scoring = 'roc_auc', cv=5)
            for score in scores:
                crossval.write(str(score) + "\t")
            crossval.write("\n")
        
    with open(sys.argv[6],"w") as output:
        for i in range(len(testdict.keys())):
            s = str.join("\t",[str(x) for x in list(testdict.keys())[i]])
            s = s + "\t"+ str(prediction[i]) +"\t"+ str(probability[i][1])

            output.write(s+"\n")
