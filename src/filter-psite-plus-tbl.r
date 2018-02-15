load("data/log.wKSEA.kinase_condition.clean.Rdata")

argv <- commandArgs(TRUE)
if (length(argv) != 2){
    stop("USAGE: <script> KIN_SUB_TBL OUT_FILE")
}

psite.plus.file <- argv[1]
out.file <- argv[2]

kinases <- read.table("data/human-kinome.txt", as.is=TRUE)[,1]

psite.plus.tbl <- read.delim(psite.plus.file, as.is=TRUE, sep="\t")

psite.plus.filt <- subset(psite.plus.tbl, GENE %in% kinases)

if ("LT_LIT" %in% colnames(psite.plus.filt)){
    psite.plus.filt$LT_LIT[which(is.na(psite.plus.filt$LT_LIT))] <- 0
    psite.plus.filt$MS_LIT[which(is.na(psite.plus.filt$MS_LIT))] <- 0
    psite.plus.filt <- subset(psite.plus.filt, LT_LIT+MS_LIT >= 5)
}

write.table(psite.plus.filt, out.file, sep="\t", quote=FALSE,
            row.names=FALSE, col.names=TRUE)
