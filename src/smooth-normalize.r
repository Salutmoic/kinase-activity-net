source("qsmooth/R/AllClasses.R")

suppressMessages(library(Biobase))
suppressMessages(library(reshape2))
load("data/external/esetNR.Rdata")

print("here")
cond.anno <- pData(esetNR)
phospho.anno <- fData(esetNR)
phospho.vals.full <- exprs(esetNR)

o = grep("Ovarian_Cancer_TCGA",cond.anno$description)
b1 = grep("C8",cond.anno$description)
b2 = grep("E2",cond.anno$description)

cptac = phospho.vals.full[,c(o,b1,b2)]

nas = apply(cptac,1,function(x) sum(is.na(x)))

cptac = cptac[nas < 100,]

cptac.norm = qsmooth(cptac,c(rep(1,length(o)),rep(2,length(b1)),rep(3,length(b2))))

norm.data = cptac.norm@qsmoothData

write.table(norm.data,"data/smooth-norm-cptac-data.tsv",quote = F, row.names = F,sep = "\t")

