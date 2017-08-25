#The idea here is to trim the sequences that are used to make the pssms
#Not sure if this works for kinases with few substrates or diverse set of substrates. But, ideally,
#this should make the pssm predictions more conservative
library(odseq)
library(kebabs)

args <- commandArgs(TRUE)

kinases = read.delim(args[1],header = F)
ksub = read.delim("data/psiteplus-kinase-substrates.tsv")

ksub.new = c()

for(kinase in kinases){
	if(kinase %in% ksub$GENE){
		kin.m = ksub[ksub$GENE == kinase,]
		kin.seqs = ksub[ksub$GENE == kinase,]$SITE_...7_AA
        kin.seqs = toupper(kin.seqs)
		sp = spectrumKernel(k = 3)
		kin.seqs = AAVector(kin.seqs)
		mat = getKernelMatrix(sp, kin.seqs)
		v = odseq_unaligned(mat, B = 1000, threshold = 0.025, type = "similarity")
		kin.m = kin.m[-which(v==T),]
		ksub.new = rbind(ksub.new,kin.m)
	}

}

write.table(ksub.new,"kinase-substrate-filt.tsv",quote = F, sep = "\t",col.names = F)


