#Use logistic regression to test and predict kinase pairs


library(glm2)

argv <- commandArgs(TRUE)
if (length(argv) != 2){
  stop("USAGE: <script> DATA_FILE OUT_FILE METHOD")
}

#I envision these tables as follows:
#train contains: kinase1, kinase2 true.pos/true.neg[0/1] pssm, association (cor-score), string
#test: kin1, kin2, pssm, association (cor-score), string

train <- argv[1]
test <- argv[2]

#take part of the training set (known outcomes) to cross validate the regression model 
train.valid = train[(round(0.75*nrow(train))):nrow(train),]
train = train[1:(round(0.75*nrow(train))),]

train <- function(train){

  model = glm(train[,3]~train[,4]+train[,5]+train[,6],family = gaussian()) 

  
  return(model)  
}


test <- function(test,train){
  
  model = train(train)

  predition = predict(model,test[,3:5], type = response)
  residuals = predict(model,test[,3:5], type = "deviance")

  test = cbind(test,prediction,residuals)
  return(test)
}

#test performance of the model on the train.valid data frame



test.res = test(test,train)
write.table(test.res,"logistic_regression.tsv",quote = F, sep = "\t")
