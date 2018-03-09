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
from scipy import interp
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt
from itertools import cycle
from sklearn.metrics import roc_curve, auc, average_precision_score, precision_recall_curve
from sklearn.model_selection import StratifiedKFold
from sklearn.ensemble import AdaBoostClassifier
from sklearn.neural_network import MLPClassifier
import random
from sklearn.preprocessing import Imputer
  
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

def plot_ROC(model, truePos,trueNeg,data,train_method):
    
    splits = StratifiedKFold(n_splits=5)
    mean_fpr = np.linspace(0, 1, 100)
    aucs= []
    tprs = []
    rcls = []
    for k in range(50):
        if train_method == "random":
            traindict,outcomesdict = train_random(truePos,data)
        else:
            traindict,outcomesdict = train_prior(truePos,trueNeg,data)
        keys = list(traindict.keys())
        random.shuffle(keys)
        train = []
        outcomes = []
        for key in keys:

            train.append(traindict[key])
            outcomes.append(outcomesdict[key])
        min_max_scaler = preprocessing.MinMaxScaler()
        train = np.asarray(train)
        outcomes = np.asarray(outcomes)
        min_max_scaler.fit(train)
        train_scale = min_max_scaler.transform(train)
        



        i = 0
        for train,test in splits.split(train_scale,outcomes):
            probs = model.fit(train_scale[train],outcomes[train]).predict_proba(train_scale[test])
            fpr, tpr, thresholds = roc_curve(outcomes[test],probs[:,1])
            tprs.append(interp(mean_fpr, fpr, tpr))
            tprs[-1][0] = 0.0
            roc_auc = auc(fpr,tpr)
            aucs.append(roc_auc)
            plt.plot(fpr,tpr,lw =1, alpha = 0.3)
            i +=1
        plt.draw()
    plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r',
         label='random', alpha=.8)


    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    std_auc = np.std(aucs)
    plt.plot(mean_fpr, mean_tpr, color='b',
         label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc),
         lw=2, alpha=.8)

    std_tpr = np.std(tprs, axis=0)
    tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
    tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
    plt.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,
                 label=r'$\pm$ 1 std. dev.')
    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.legend(loc="lower right")
    
    
    plt.savefig("img/mldata/"+ method + "_unscale_" + train_method+'_'+'ROC2'+'.pdf')
    plt.gcf().clear()     

def plot_PR(model, truePos,trueNeg,data,train_method):
    rcls = []
    prcs = []
    aprcs = []
    aucs = []
    splits = StratifiedKFold(n_splits=5)

    for k in range(50):
        if train_method == "random":
            traindict,outcomesdict = train_random(truePos,data)
        else:
            traindict,outcomesdict = train_prior(truePos,trueNeg,data)
        keys = list(traindict.keys())
        random.shuffle(keys)
        train = []
        outcomes = []
        for key in keys:

            train.append(traindict[key])
            outcomes.append(outcomesdict[key])

        min_max_scaler = preprocessing.MinMaxScaler()
        train = np.asarray(train)
        outcomes = np.asarray(outcomes)
        min_max_scaler.fit(train)
        train_scale = min_max_scaler.transform(train)
        finval = []
     
        for train,test in splits.split(train_scale,outcomes):
            y_test = outcomes[test]
        
            y_scores = model.fit(train_scale[train],outcomes[train]).predict_proba(train_scale[test])
            # y_scores = [x[v] for x,v in zip(y_scores,y_test)]
            y_scores = [x[1] for x in y_scores]
         
            precision, recall, _ = precision_recall_curve(y_test, y_scores)
            prcs.append(precision)
            print(precision)
            rcls.append(recall)
            finval.append(precision[0])   
            aucs.append(auc( recall,precision)) 
            plt.step(recall, precision, alpha=0.2,where='post')
            aprcs.append(average_precision_score(y_test, y_scores))
            plt.draw()
    ap = np.median(aprcs)
    aucm = np.mean(aucs)
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.ylim([0.0, 1.05])
    plt.xlim([0.0, 1.0])
    plt.plot(label='mean AUC =' + str(aucm) )
    plt.legend([ 'average precision = ' + str(round(ap,2))],loc="lower right")
    print(np.median(finval))
    plt.savefig("img/mldata/"+ method + "_unscale_"+"_" + train_method + '_Precision_Recall_curve2'+'.pdf')
    
def predict(model,truePos,trueNeg,data,train_method):
    
    min_max_scaler = preprocessing.MinMaxScaler()
    testdict = prediction_matrix(data)
    test = list(testdict.values())
    pairs = list(testdict.keys())
    min_max_scaler.fit(test)

    test_scale = min_max_scaler.transform(test)

    predictions = {}    
   
    for k in range(50):
        if train_method == "random":
            traindict,outcomesdict = train_random(truePos,data)
        else:
            traindict,outcomesdict = train_prior(truePos,trueNeg,data)
        keys = list(traindict.keys())
        random.shuffle(keys)
        train = []
        outcomes = []
        for key in keys:

            train.append(traindict[key])
            outcomes.append(outcomesdict[key])

        min_max_scaler = preprocessing.MinMaxScaler()
        train = np.asarray(train)
        outcomes = np.asarray(outcomes)
        min_max_scaler.fit(train)
        train_scale = min_max_scaler.transform(train)
        
        probability = model.fit(train_scale,outcomes).predict_proba(test_scale)
        
        for i in range(len(pairs)):
            
            if pairs[i] in predictions:
                predictions[pairs[i]].append(probability[i][1])
            else:
                predictions[pairs[i]] = [probability[i][1]]     
            if k == 49:  
                predictions[pairs[i]] = np.mean(predictions[pairs[i]])
    return predictions

def cross_val(model,method,truePos,trueNeg,data,train_method):
    yvals = []
    pvals = []
      
    splits = StratifiedKFold(n_splits=5)
    for k in range(50):
        if train_method == "random":
            traindict,outcomesdict = train_random(truePos,data)
        else:
            traindict,outcomesdict = train_prior(truePos,trueNeg,data)

        keys = list(traindict.keys())
        random.shuffle(keys)
        train = []
        outcomes = []
        for key in keys:
       
            train.append(traindict[key])
            outcomes.append(outcomesdict[key])
        min_max_scaler = preprocessing.MinMaxScaler()
        train = np.asarray(train)
        outcomes = np.asarray(outcomes) 
        min_max_scaler.fit(train)
        train_scale = min_max_scaler.transform(train)
        for train,test in splits.split(train_scale,outcomes):
            probs = model.fit(train_scale[train],outcomes[train]).predict_proba(train_scale[test])
            y = outcomes[test]
            pvals.append([p[1] for p in probs])
            yvals.append(y)  
   
    print(len(pvals))
    with open("mlres/"+method+"prob.tsv","w") as pfile:
        for i in range(len(pvals)):
            print([str(x) for x in pvals[i]])    
            pfile.write("\t".join([str(x) for x in pvals[i]])+"\n")
    with open("mlres/"+method+"yvals.tsv","w") as yfile:
        for i in range(len(yvals)):
            yfile.write("\t".join([str(x) for x in yvals[i]])+"\n")     
    return
 
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
        traindict,outcomesdict = train_random(truePos,data)
    else:
        traindict,outcomesdict = train_prior(truePos,trueNeg,data)

    keys = list(traindict.keys())
    random.shuffle(keys)
    train = []
    outcomes = []
    for key in keys:
       print(key)
       train.append(traindict[key])
       outcomes.append(outcomesdict[key])
    min_max_scaler.fit(train)
    train_scale = min_max_scaler.transform(train) 

    
    if method == "NB":
        model = GaussianNB()
        
    if method == "SVM":
        model = svm.SVC(kernel = 'rbf',probability = True)
    
    if method == "boost":
        model = AdaBoostClassifier(n_estimators=150)     
   
    if method == "kmeans":
        model = NearestCentroid()
        model.fit(train_scale, outcomes)
        prediction = model.predict(test_scale)

        centers = find_centers(test,prediction)
        c1 = centers[0]
        c2 = centers[1]

        if len(c1) > 0 and len(c2) > 0:
        
            probability = probability(c1,c2,test,prediction)

    if method == "logreg":
        model = linear_model.LogisticRegression()
    
    if method == "neural":
       model =  MLPClassifier()       

    if method == "randfor":
        model = RandomForestClassifier(n_estimators = 150)

    if model != "kmeans":
        
        prediction = model.fit(train_scale,outcomes).predict(test_scale)
        probability = model.fit(train_scale,outcomes).predict_proba(test_scale)
    
    plot_ROC(model,truePos,trueNeg,data,train_method)
    plot_PR(model, truePos,trueNeg,data,train_method)

    probs = predict(model,truePos,trueNeg,data,train_method)
    cross_val(model,method,truePos,trueNeg,data,train_method) 
    # repeated cross validation, is this something along the lines of what you were thinkin
      
    with open(sys.argv[6],"w") as output:
        for i in range(len(testdict.keys())):
            s = str.join("\t",[str(x) for x in list(testdict.keys())[i]])
            s = s + "\t"+ str(probs[list(testdict.keys())[i]])

            output.write(s+"\n")


