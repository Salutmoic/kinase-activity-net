make.set = function(data,pairs){
  #This function finds true neg and true pos 
  #kinase pairs in the merged data set and returns
  #The matrices.
  data.pairs = paste(data[,1],data[,2])
  true.pairs = paste(pairs[,1],pairs[,2])
  
  true.matrix = data[which(data.pairs %in% true.pairs),]
  
  return(true.matrix)
}

argv <- commandArgs(TRUE)
if(length(argv) != 4){
  stop("USAGE: <script> DATA_FILE OUT_FILE")
  
}

data = read.delim(argv[1])
true.pos =read.delim(argv[2])
true.neg = read.delim(argv[3])
out.file = argv[4]

pos.mat = make.set(data,true.pos)
neg.mat = make.set(data,true.neg)

out.mat = cbind(rbind(pos.mat,neg.mat),c(rep(1,nrow(pos.mat)),
                                         rep(0,nrow(neg.mat))))

write.table(out.mat,out.file,quote = F,sep = "\t",
            row.names = F,col.names = F)







