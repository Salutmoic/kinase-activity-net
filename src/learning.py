import scipy
import numpy
import sklearn
from sklearn import svm, linear_model
import sys
from sklearn.naive_bayes import GaussianNB
from sklearn.neighbors.nearest_centroid import NearestCentroid


def traindata(file1):
# Creates matrix (list of list) containing true pos and true negatives
# for training and predictions, f1 is true positives and true negatives
# x is predictive variables (ppsm, string,cor) and y outcomes
    x = []
    y=[]
    ppairs = []

    with open(file1,"r") as f1:
        for line in f1:
            line = line.split("\t")
            if 'NA' not in line:
               x.append([float(v) for v in line[2:6]])
               ppairs.append([line[0],line[1]])
               y.append(line[6])

    return x, y, ppairs

def prediction_matrix(data):
    #makes matrix with pssm values string score and corr
    ppairs = []
    v = []
    with open(data,"r") as d:
        d.readline()
        for line in d:
            line = line.split("\t")
            if 'NA' not in line:
                v.append([float(x) for x in line[2:]])
                ppairs.append([line[0],line[1]])
    return v, ppairs

def SuppVectMach(data,f1):
# uses Support Vector Machines to classify data into 1 (interaction) and 0 (no interaction) and #returns probabilities
    test,ppairs = prediction_matrix(data)
    train,outcomes,ppairs_train = traindata(f1)
    clf = svm.SVC(probability = True)
    clf.fit(train, outcomes)

    results = clf.predict(test)
    probability = clf.predict_proba(test)
    return ppairs,results,probability

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
    c1 = [numpy.mean(center1[0]),numpy.mean(center1[1]),numpy.mean(center1[2]),numpy.mean(center1[3])]
    c2 = [numpy.mean(center2[0]),numpy.mean(center2[1]),numpy.mean(center2[2]),numpy.mean(center2[3])]

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

def nearestneighbour(data,f1):
# clusters test data into 2 clusters (0 and 1)
    test,ppairs = prediction_matrix(data)
    train,outcomes,ppairs_train = traindata(f1)
    clf = NearestCentroid()
    clf.fit(train, outcomes)
    results = clf.predict(test)

    centers = find_centers(test,results)
    c1 = centers[0]
    c2 = centers[1]

    if len(c1) > 0 and len(c2) > 0:

        probs = probability(c1,c2,test,results)

        return ppairs,results,probs
    else:
        return ppairs,results

def naive_bayes(data,f1):
# calculates naive bayes 
    model = GaussianNB()
    test,ppairs = prediction_matrix(data)
    train,outcomes,ppairs_train = traindata(f1)

    prediction = model.fit(train,outcomes).predict(test)
    probability = model.fit(train,outcomes).predict_proba(test)

    return ppairs,prediction,probability

def log_reg(data,f1):
    model = linear_model.LogisticRegression()
    test,ppairs = prediction_matrix(data)
    train,outcomes,ppairs_train = traindata(f1)

    prediction = model.fit(train,outcomes).predict(test)
    probability = model.fit(train,outcomes).predict_proba(test)

    return ppairs,prediction,probability

if __name__ == "__main__":

    trues = sys.argv[1]
    data = sys.argv[2]
    method = sys.argv[3]

    if method == "NB":
        model = naive_bayes(data,trues)

    if method == "SVM":
        model = SuppVectMach(data,trues)

    if method == "kmeans":
        model = nearestneighbour(data,trues)

    if method == "logreg":
        model = log_reg(data,trues)

    ppairs = model[0]
    prediction = model[1]
    if len(model) > 2:
        probability = model[2]

    with open(sys.argv[4],"w") as output:
        for i in range(0,len(ppairs)):
            if len(model) > 2:
                output.write(str(ppairs[i][0]) + "\t"+ str(ppairs[i][1]) \
                + "\t"+str(prediction[i]) +"\t" + str(max(probability[i])) + "\n" )
            else:
                output.write(str(ppairs[i][0]) + "\t"+ str(ppairs[i][1])\
                 + "\t"+str(prediction[i]) +"\n" )