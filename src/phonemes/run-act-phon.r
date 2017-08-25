#system("Rscript make-act-data.r")

targets <- as.data.frame(read.csv("Data_AND_Background_Network/targets.csv"))
targets  <- unique(as.character(targets$condition))
#add sensible pairs
m = seq(1,length(targets),2)
all.ts = ""

for(j in m){
	if(j > 13){targets = c(targets,paste(targets[j-1],targets[j],sep = "-"))
			all.ts = paste(all.ts,paste(targets[j-1],targets[j],sep = "-"),sep = "-")}
	if(j < 13){targets = c(targets,paste(targets[j],targets[j+1],sep = "-"))
	all.ts = paste(all.ts,paste(targets[j],targets[j+1],sep = "-"),sep = "-")}
}
all.ts = substr(all.ts,2,nchar(all.ts))
targets = c(targets,all.ts)


####First of all make the network to use for this


ts =seq(0.5,0.7,0.025)
system("Rscript make-act-map.r validation-set-omnipath.tsv")


for(t in ts){
	system(paste("Rscript make-act-map.r merged-predictor-hinegs-annot.tsv",t))
}

for(target in targets){

system(paste("Rscript Scripts/buildTT.R act-tbl-norm.tsv",target))
	system(paste("Rscript Scripts/executionScript.R my-act-networks/prior-matrix-act-omni.Rdata map-act.tsv",sep = ""))
	load(paste("my-act-networks-hinegs/prior-matrix-act-omni.Rdata",sep = ""))	
	system(paste("mv Scripts/resultsSIF.txt my-results-act-hinegs/resultsSIF-omni-",as.character(nrow(prior)),"-",target,"-act.txt",sep = ""))
	system(paste("mv act-species/species.tsv act-species/act-species-",target,"-omni.tsv",sep = ""))



for(t in ts){
	n = length(list.files("my-results-act-hinegs/"))
	system(paste("Rscript Scripts/buildTT.R act-tbl-norm.tsv",target))
	system(paste("Rscript Scripts/executionScript.R my-act-networks-hinegs/prior-matrix-act-",t,".Rdata map-act.tsv",sep = ""))
	load(paste("my-act-networks-hinegs/prior-matrix-act-",t,".Rdata",sep = ""))	
	system(paste("mv Scripts/resultsSIF.txt my-results-act-hinegs/resultsSIF-",as.character(t),"-",as.character(nrow(prior)),"-",target,"-act.txt",sep = ""))
	system(paste("mv act-species/species.tsv act-species/act-species-",target,"-",t,".tsv",sep = ""))
print(t)
	n2 = length(list.files("my-results-act-hinegs/"))
	#if(n == n2){break }
}
}



