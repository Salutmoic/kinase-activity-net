from train_set import *
import scipy
import numpy as np
import sklearn
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
from crossval import cross_validation
from sklearn import preprocessing
#Script to test out different number of estimators (number of trees in the forest). The probabilities seem to behave strangely when the number is low.

    
if __name__ == "__main__":

    truePos = sys.argv[1]
    trueNeg = sys.argv[2]
    data = sys.argv[3]
    train_method = sys.argv[5]
    out = "randres/"+sys.argv[6] 
    
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

    #number of trees in the forest
    c = [float(i) for i in range(1,138,2)]

    for j in c:
        model = RandomForestClassifier(n_estimators = int(j))
        prediction = model.fit(train,outcomes).predict(test)
        probability = model.fit(train,outcomes).predict_proba(test)	

        with open(out +"_"+ str(j)+".tsv","w") as output:
            for i in range(len(testdict.keys())):
                s = str.join("\t",[str(x) for x in list(testdict.keys())[i]])
                s = s + "\t"+ str(prediction[i]) +"\t"+ str(probability[i][1])
                output.write(s+"\n")

        
            
        cross_validation(model, out + "_"+ train_method +  "_" +str(j),truePos,trueNeg,data, 50, train_method)
