#code that uses support vector machines to classify kinase-kinase interactions
from train_set import * 
import scipy
import numpy as np
import sklearn
from sklearn import svm, preprocessing
from sklearn.model_selection import cross_val_score
from crossval import cross_validation
import sys

def suppvectm(train,outcomes,test,kernel,j):
    
   #note that the probability in SVM is based on cross validation, therefore the cross val assessment
   #might be
    if kernel == "rbf":
        model = svm.SVC(C = j,cache_size = 1000,probability = True)

    if kernel == "linear":
        model = svm.SVC(C = j,cache_size = 1000,kernel = 'linear', probability = True)
        print(outcomes)    
    prediction = model.fit(train,outcomes).predict(test_scale)
    probability = model.fit(train,outcomes).predict_proba(test)

    return model, prediction, probability
    
if __name__ == "__main__":

    truePos = sys.argv[1]
    trueNeg = sys.argv[2]
    data = sys.argv[3]
    kernel = sys.argv[4]
    train_method = sys.argv[5]
    out = "svmres/"+sys.argv[6]+"_" + kernel
    
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

    
    c = [float(i)/20 + 0.01 for i in range(20)]

    for j in c:
        model, prediction, probability = suppvectm(train_scale,outcomes,test_scale,kernel,j)

        with open(out+"_" + str(int(j*100))+".tsv","w") as output:
            for i in range(len(testdict.keys())):
                s = str.join("\t",[str(x) for x in list(testdict.keys())[i]])
                s = s + "\t"+ str(prediction[i]) +"\t"+ str(probability[i][1])
                output.write(s+"\n")

        
        cross_validation(model, out+"_"+ train_method + "_" +str(j),truePos,trueNeg,data,50,train_method)
       
            


    


            


            

    



