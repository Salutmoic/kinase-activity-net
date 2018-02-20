#made by acessed via Luz Garcia Alonso  https://github.com/saezlab/DoRothEA/blob/master/src/README.md
source('DoRothEA/src/lib_enrichment_scores.r')


#activity claculated from the tcga data set
make_regulon = function(tfs){

regulon = list()
regulon[[1]] = as.character(unique(tfs[,1]))


regulon[[2]] = list()
names(regulon) = c("NAME","GENES")


for(t in tfs[,1]){

tf = tfs[tfs[,1] == t,]

regulon[[2]][[t]] = rep(1,nrow(tf))
names(regulon[[2]][[t]]) = tf[,2]

}

return(regulon)
}


TFs = read.delim("TFs/tfinter.tsv",header = F)

expr = read.delim("TCGA/tcge_expr_brca_ov.txt")
rownames(expr) = expr[,1]
expr = data.matrix(expr[,2:ncol(expr)])
expr = gene_expression_statistic(expr,method = 'scale')

regulon = make_regulon(TFs)


activity = SLEA(E = expr, genesets = regulon, method = 'VIPER')$NES

write.table(activity,"TFs/tcga-TF-activities.tsv",quote = F, sep = "\t")







