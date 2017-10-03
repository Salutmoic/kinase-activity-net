load("data/external/log.wGSEA.kinase_condition.clean.Rdata")

argv <- commandArgs(TRUE)
if (length(argv) != 2){
    stop("USAGE: <script> KIN_SUB_TBL OUT_FILE")
}

psite.plus.file <- argv[1]
out.file <- argv[2]

psite.plus.tbl <- read.delim(psite.plus.file, as.is=TRUE, sep="\t")
kin.act <- log.wGSEA.kinase_condition.clean

psite.plus.filt <- subset(psite.plus.tbl, GENE %in% rownames(kin.act))

write.table(psite.plus.filt, out.file, sep="\t", quote=FALSE,
            row.names=FALSE, col.names=TRUE)
