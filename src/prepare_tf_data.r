
find_tfs = function(){

	tfs = read.delim("TFs/tfClass2ensembl.txt")
	tfs = tfs[tfs$specie == "Homo sapiens",]    

	map = read.delim("TFs/ensembl_map.txt")
	
	tfs.gene = unique(map[map[,1] %in% tfs[,2],3])
    write.table(tfs.gene,"TFs/TFs-human.txt",quote = F,sep = "\t", row.names = F)
}
find_tf_kinases = function(){

	ksub = read.delim("data/external/Kinase_Substrate_Dataset_2017-09-08",skip =2)
	ksub = ksub[ksub$KIN_ORGANISM == "human" & ksub$SUB_ORGANISM == "human",]

	tfs = read.delim("TFs/TFs-human.txt")

	kinases = read.delim("data/human-kinome.txt")

	kins.to.use = ksub[ksub$SUB_GENE %in% tfs[,1],1]

	kinases = kinases[kinases[,1] %in% kins.to.use,]
	ktable = ksub[ksub$GENE %in% kinases,]

	write.table(kinases,"TFs/TF-kinases.tsv",row.names = F, sep = "\t",quote = F)
	write.table(ktable,"TFs/tf_kinase-table.tsv",row.names = F, sep = "\t",quote = F)

}

find_tf_regsites = function(){
	regs = read.delim("data/external/Regulatory_sites_2017-09-08",skip =2 )
	regs = regs[regs$ORGANISM == "human",]

	tfs = read.delim("TFs/TFs-human.txt")

	tf.reg = regs[regs[,1] %in% tfs[,1],]

	write.table(tf.reg,"TFs/tf-reg-sites.tsv",row.names = F, sep = "\t",quote = F)
}


find_tf_psites = function(){

	psites = read.delim("data/external/Phosphorylation_site_dataset_2017-09-08",skip =2)

	psites = psites[psites$ORGANISM == "human",]
	tfs = read.delim("TFs/TFs-human.txt")

	tf.p = psites[psites$GENE %in% tfs[,1],]

	write.table(tf.p,"TFs/TF-sites.tsv",row.names = F, sep = "\t",quote = F)

}

find_phosfun = function(){
	st.func = read.delim("data/external/phosfun-predictions-ST.tsv", sep = ",")
	y.func  = read.delim("data/external/phosfun-predictions-Y.tsv", sep = ",")

	ps = read.delim("TFs/TF-sites.tsv")

	res = gsub("-p","",ps$MOD_RSD)
	res = substr(res,2,nchar(res))
	ps = cbind(ps,res)
	psites = paste(ps[,3],res,sep = "_")

	func = rbind(st.func,y.func)

	d = ps[which(psites %in% func[,1]),]
	psites.func = psites[which(psites %in% func[,1])]

	func.tf = func[which(func[,1] %in% psites),]

	func.tf = func.tf[order(as.character(func.tf[,1])),]

	d[,16] = paste(d[,3],d[,15],sep = "_")

	d = d[order(d[,16]), ]

	print(table(func.tf[,1] == d[,16]))

	d = cbind(as.character(d[,2]),as.character(d[,15]),func.tf[,2])
	colnames(d) = c("prot","pos","score")

	write.table(d,"TFs/TFs-functional_score.tsv",sep = "\t",quote = F, row.names = F)

}


find_tfs()
find_tf_kinases()
find_tf_psites()
find_tf_regsites()
find_phosfun()
