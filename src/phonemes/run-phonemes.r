
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

#make data-sets and networks
system("Rscript phonm-data.r func")
system("Rscript phonm-data.r n")
system("Rscript phonm-data.r reg")

#bio-plex networks
if(length(list.files("bio-prior/")) ==0){
bio.ts = seq(0.75,0.95,0.05)
message("making-bioplex-priors this might take a while")
for(t in bio.ts){
	system(paste("Rscript make-bioplex-prior.r",t,"log-subset.csv my-map.tsv"))
	system(paste("Rscript make-bioplex-prior.r",t, "log-subset-func.csv my-map-func.tsv"))
	system(paste("Rscript make-bioplex-prior.r",t,"log-subset-reg.csv my-map-reg.tsv"))
	print(t)
}
}

if(length(list.files("my-networks")) == 0){
#My -netw
message("making my-priors this might take a while")
a.t = seq(0.6,0.9,0.05)
d.t = seq(0.5,0.7,0.05)
for(a in a.t){
	for(d in d.t){
		system(paste("Rscript make-netw-for-phonms.r",a,d,"func",sep = " "))
		system(paste("Rscript make-netw-for-phonms.r",a,d,"n",sep = " "))
		system(paste("Rscript make-netw-for-phonms.r",a,d,"reg",sep = " "))		
	}
	print(a)
}

}
#Run phonemes
for(target in targets){

#Omnipath
	
	system(paste("Rscript Scripts/buildTT.R Data_AND_Background_Network/log_measurments.csv",target))
	system("Rscript Scripts/executionScript.R omni-prior-func.Rdata my-map-func.tsv")
	system(paste("mv Scripts/resultsSIF.txt my-results-3/resultsSIF-omni-prior-func",target,".txt",sep = "-"))

	system(paste("Rscript Scripts/buildTT.R Data_AND_Background_Network/log_measurments.csv",target))
	system("Rscript Scripts/executionScript.R omni-prior.Rdata my-map.tsv")
	system(paste("mv Scripts/resultsSIF.txt my-results-3/resultsSIF-omni-prior",target,".txt",sep = "-"))

	system(paste("Rscript Scripts/buildTT.R Data_AND_Background_Network/log_measurments.csv",target))
	system("Rscript Scripts/executionScript.R omni-prior-reg.Rdata my-map-reg.tsv")
	system(paste("mv Scripts/resultsSIF.txt my-results-3/resultsSIF-omni-prior-reg",target,".txt",sep = "-"))

	#Bioplex
	bio.ts = seq(0.75,0.95,0.05)

	for(t in bio.ts){
		
		system(paste("Rscript Scripts/buildTT.R Data_AND_Background_Network/log_measurments.csv",target))
		system(paste("Rscript Scripts/executionScript.R bio-prior/bioplex-prior-matrix-",t,".Rdata my-map.tsv"))
		load(paste("bio-prior/bioplex-prior-matrix-",t,".Rdata",sep=""))
		system(paste("mv Scripts/resultsSIF.txt my-results-3/resultsSIF-bio-prior",t,nrow(prior),target,".txt",sep = ""))

		system(paste("Rscript Scripts/buildTT.R Data_AND_Background_Network/log_measurments.csv",target))
		system(paste("Rscript Scripts/executionScript.R bio-prior/bioplex-prior-matrix-",t,"-func.Rdata my-map-func.tsv",sep=""))
		load(paste("bio-prior/bioplex-prior-matrix-",t,"-func.Rdata",sep = ""))
		system(paste("mv Scripts/resultsSIF.txt my-results-3/resultsSIF-bio-prior",t,target,nrow(prior),"-func.txt",sep = ""))
		
		system(paste("Rscript Scripts/buildTT.R Data_AND_Background_Network/log_measurments.csv",target))
		system(paste("Rscript Scripts/executionScript.R bio-prior/bioplex-prior-matrix-",t,"-reg.Rdata my-map-reg.tsv",sep = ""))
		load(paste("bio-prior/bioplex-prior-matrix-",t,"-reg.Rdata",sep = ""))
		system(paste("mv Scripts/resultsSIF.txt my-results-3/resultsSIF-bio-prior",t,target,"-",nrow(prior),"-reg.txt",sep = ""))

	}

#My-network
	a.t = seq(0.6,0.9,0.05)
	d.t = seq(0.5,0.7,0.05)
	for(a in a.t){
		for(d in d.t){
			
# NB gave good results		system(paste("Rscript Scripts/buildTT.R Data_&_Background_Network/log-subset.csv",target))
			system(paste("Rscript Scripts/buildTT.R Data_AND_Background_Network/log_measurments.csv",target))
			system(paste("Rscript Scripts/executionScript.R my-networks/prior-matrix-func-",a,"-",d,".Rdata my-map-func.tsv my-map-func.tsv",sep = ""))
			load(paste("my-networks/prior-matrix-func-",a,"-",d,".Rdata",sep = ""))	
			system(paste("mv Scripts/resultsSIF.txt my-results-3/resultsSIF-",as.character(a),"-",as.character(d),"-",as.character(nrow(prior)),"-",target,"-func.txt",sep = ""))
		
			system(paste("Rscript Scripts/buildTT.R Data_AND_Background_Network/log_measurments.csv",target))	
			system(paste("Rscript Scripts/executionScript.R my-networks/prior-matrix-",a,"-",d,".Rdata my-map.tsv",sep =""))	
			load(paste("my-networks/prior-matrix-",a,"-",d,".Rdata",sep = ""))
			system(paste("mv Scripts/resultsSIF.txt my-results-3/resultsSIF-",as.character(a),"-",as.character(d),"-",as.character(nrow(prior)),"-",target,"-2.txt",sep = ""))

			system(paste("Rscript Scripts/buildTT.R Data_AND_Background_Network/log_measurments.csv",target))	
			system(paste("Rscript Scripts/executionScript.R my-networks/prior-matrix-reg-",a,"-",d,".Rdata my-map-reg.tsv",sep=""))	
			load(paste("my-networks/prior-matrix-reg-",a,"-",d,".Rdata",sep=""))
			system(paste("mv Scripts/resultsSIF.txt my-results-3/resultsSIF-",as.character(a),"-",as.character(d),"-",as.character(nrow(prior)),"-",target,"-reg.txt",sep = ""))
		}
	}

}


