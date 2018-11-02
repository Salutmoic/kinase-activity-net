source("phon-dat-functions.r")

argv <- commandArgs(TRUE)

t.a = argv[1]
t.d = argv[2]
status = argv[3]

if(status == "reg"){

	my.map = read.delim("my-map-reg.tsv")
	data = read.delim("log-subset-reg.csv")

}else if(status == "func"){
	my.map = read.delim("my-map-func.tsv")
	data = read.delim("log-subset-func.csv")

}else{

	my.map = read.delim("my-map.tsv")
	data = read.delim("log-subset.csv")


}

uniprot = read.delim("uniprot-id-map.tsv",header = F)
load("Data_AND_Background_Network/allD_HUMAN_Omnipath.Rdata")

d.rows = strsplit(as.character(data[,1])," ")
d.rows = lapply(d.rows,function(x) x[1])

netw = read.delim("separated-predictor-annot.tsv")
kins = unique(as.character(netw[,1]))

netw = netw[which(netw[,3] > t.a & netw$direct.pred.mean > t.d),]
netw = netw[which(netw[,2] %in% d.rows),]

prior = make.prior(netw,allD,uniprot,data,d.rows,my.map)
prior = as.data.frame(prior)

if(status == "reg"){
	save("prior",file= "prior-matrix-reg.Rdata")
}else if(status == "func"){
	save("prior",file= "prior-matrix-func.Rdata")
}else{
	save("prior",file= "prior-matrix.Rdata")
}


