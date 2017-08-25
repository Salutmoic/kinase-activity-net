import scipy
from train_set import *
import sklearn
import sys
from sklearn.model_selection import cross_val_score
from random import *
from sklearn import preprocessing

def cross_validation(model,file,truePos,trueNeg,data,j,method):
    # function to use with the machine learning methods, model is the model that is being cross-validated.
    # file is the output file. method is the methot used to provide true negatives (prior or random)
    # data is the kinase-pair feature matrix.
    with open("cross_res/"+ file,"w") as crossval:

        for i in range(j):
            if method == "random":
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
    return

