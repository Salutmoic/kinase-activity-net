
#Use logistic regression to test and predict kinase pairs


library(glm)

argv <- commandArgs(TRUE)
if (length(argv) != 2){
  stop("USAGE: <script> DATA_FILE OUT_FILE METHOD")
}

#I envision these tables as follows:
#train contains: kinase1, kinase2 true.pos/true.neg[0/1] pssm, association (cor-score), string
#test: kin1, kin2, pssm, association (cor-score), string

train <- argv[1]
test <- argv[2]

train <- function(train){

  model = glm(train[,3]~train[,4]+train[,5]+train[,6]) 

  
  return(model)  
}


test <- function(test,train){
  
  model = train(train)

  predition = predict(model,test[,3:5], type = response)
  residuals = predict(model,test[,3:5], type = "deviance")

  test = cbind(test,prediction,residuals)
  return(test)
}

test.res = test(test,train)
write.table(test.res,"logistic_regression.tsv",quote = F, sep = "\t")
