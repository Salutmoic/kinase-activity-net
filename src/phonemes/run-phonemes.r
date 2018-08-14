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

for(target in targets){

#Omnipath
system("Rscript phonm-data.r 0.0 0.0 func")
system(paste("Rscript Scripts/buildTT.R log-subset-func.csv",target))
system("Rscript Scripts/executionScript.R omni-prior-func.Rdata my-map-func.tsv")
system(paste("mv Scripts/resultsSIF.txt my-results-2/resultsSIF-omni-prior-func",target,".txt",sep = "-"))

system("Rscript phonm-data.r 0.0 0.0 n")
system(paste("Rscript Scripts/buildTT.R log-subset.csv",target))
system("Rscript Scripts/executionScript.R omni-prior.Rdata my-map.tsv")
system(paste("mv Scripts/resultsSIF.txt my-results-2/resultsSIF-omni-prior",target,".txt",sep = "-"))

system("Rscript phonm-data.r 0.0 0.0 reg")
system(paste("Rscript Scripts/buildTT.R log-subset-reg.csv",target))
system("Rscript Scripts/executionScript.R omni-prior-reg.Rdata my-map-reg.tsv")
system(paste("mv Scripts/resultsSIF.txt my-results-2/resultsSIF-omni-prior-reg",target,".txt",sep = "-"))

#Bioplex
bio.ts = seq(0.75,0.95,0.05)

for(t in bio.ts){
system(paste("Rscript make-bioplex-prior.r",t,"log-subset.csv my-map.tsv"))
system(paste("Rscript Scripts/buildTT.R log-subset.csv",target))
system("Rscript Scripts/executionScript.R bioplex-prior-matrix.Rdata my-map.tsv")
load("bioplex-prior-matrix.Rdata")
system(paste("mv Scripts/resultsSIF.txt my-results-2/resultsSIF-bio-prior","-",nrow(prior),target,".txt",sep = "-"))


system(paste("Rscript make-bioplex-prior.r",t, "log-subset-func.csv my-map-func.tsv"))
system(paste("Rscript Scripts/buildTT.R log-subset-func.csv",target))
system("Rscript Scripts/executionScript.R bioplex-prior-matrix.Rdata my-map-func.tsv")
load("bioplex-prior-matrix.Rdata")
system(paste("mv Scripts/resultsSIF.txt my-results-2/resultsSIF-bio-prior",target,"-",nrow(prior),"-func.txt",sep = "-"))

system(paste("Rscript make-bioplex-prior.r",t,"log-subset-reg.csv my-map-reg.tsv"))
system(paste("Rscript Scripts/buildTT.R log-subset-reg.csv",target))
system("Rscript Scripts/executionScript.R bioplex-prior-matrix.Rdata my-map-reg.tsv")
load("bioplex-prior-matrix.Rdata")
system(paste("mv Scripts/resultsSIF.txt my-results-2/resultsSIF-bio-prior",target,"-",nrow(prior),"-reg.txt",sep = "-"))

}

#My-network
a.t = c(0.0,seq(0.6,0.9,0.05))
d.t = c(0.0,seq(0.5,0.7,0.05))

for(a in a.t){
	for(d in d.t){
		system(paste("Rscript phonm-data.r",a,d,"func",sep = " "))
# NB gave good results		system(paste("Rscript Scripts/buildTT.R Data_&_Background_Network/log-subset.csv",target))
		system(paste("Rscript Scripts/buildTT.R log-subset-func.csv",target))
		system("Rscript Scripts/executionScript.R prior-matrix-func.Rdata my-map-func.tsv")
		load("prior-matrix-func.Rdata")	
		system(paste("mv Scripts/resultsSIF.txt my-results-2/resultsSIF-",as.character(a),"-",as.character(d),"-",as.character(nrow(prior)),"-",target,"-func.txt",sep = ""))
		
		system(paste("Rscript phonm-data.r",a,d,"n",sep = " "))
		system(paste("Rscript Scripts/buildTT.R log-subset.csv",target))	
		system("Rscript Scripts/executionScript.R prior-matrix.Rdata my-map.tsv")	
		load("prior-matrix.Rdata")
		system(paste("mv Scripts/resultsSIF.txt my-results-2/resultsSIF-",as.character(a),"-",as.character(d),"-",as.character(nrow(prior)),"-",target,"-2.txt",sep = ""))

		system(paste("Rscript phonm-data.r",a,d,"reg",sep = " "))
		system(paste("Rscript Scripts/buildTT.R log-subset-reg.csv",target))	
		system("Rscript Scripts/executionScript.R prior-matrix-reg.Rdata my-map-reg.tsv")	
		load("prior-matrix-reg.Rdata")
		system(paste("mv Scripts/resultsSIF.txt my-results-2/resultsSIF-",as.character(a),"-",as.character(d),"-",as.character(nrow(prior)),"-",target,"-reg.txt",sep = ""))
	}
}

}


