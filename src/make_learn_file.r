#This scripts writes files that the learning methods can use 

argv <- commandArgs(TRUE)
if (length(argv) != 4){
  stop("USAGE: <script> DATA_FILE OUT_FILE METHOD")
}


cor = read.table(argv[1],as.is = TRUE)
pssm = read.table(argv[2],as.is = TRUE)
string = read.table(argv[3],as.is = TRUE)
out.file = argv[4]

cor = cor[order(cor[,1],cor[2]),]
pssm = pssm[order(pssm[,1],pssm[2]),]
string = string[order(pssm[,1],pssm[2]),]


if((paste(cor[,1],cor[,2]) == paste(ppsm[,1],pssm[,2])) &  (paste(cor[,1],cor[,2]) == paste(string[,1],string[,2]))){
	out = bind(cor,pssm[,3],string[,3])

}



