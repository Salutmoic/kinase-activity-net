library(reshape2)

argv <- commandArgs(TRUE)
if (length(argv) != 2){
    stop("USAGE: <script> DATA_FILE OUT_FILE")
}

kin.act.file <- argv[1]
out.file <- argv[2]

kin.act.tbl <- read.table(kin.act.file, as.is=TRUE)
names(kin.act.tbl) <- c("kinase", "condition", "activity")
kin.act.m <- as.data.frame(acast(kin.act.tbl, kinase ~ condition))

## kin.cond.tbl <- read.delim("data/external/kinase_invivoconditions.csv",
##                            as.is=TRUE, sep=",")
kin.cond.tbl <- read.delim("data/external/kinase-condition-pairs.tsv", as.is=TRUE)

intxns <- NULL
kinases1 <- c()
kinases2 <- c()
cors <- c()

for (kinase1 in rownames(kin.act.m)){
    message(paste(kinase1, "-", paste0(match(kinase1, rownames(kin.act.m)), "/", nrow(kin.act.m))))
    for (kinase2 in rownames(kin.act.m)){
        if (kinase2 %in% kin.cond.tbl$Kinase){
            kin2.conds <- subset(kin.cond.tbl, Kinase==kinase2)
            unperturbed2 <- setdiff(colnames(kin.act.m), kin2.conds$Condition)
        }else{
            unperturbed2 <- colnames(kin.act.m)
        }
        kin.act.sub <- subset(kin.act.m, select=unperturbed2)
        pair.cor <- cor(t(kin.act.sub[kinase1,]), t(kin.act.sub[kinase2,]), method="spearman")
        kinases1 <- c(kinases1, kinase1)
        kinases2 <- c(kinases2, kinase2)
        cors <- c(cors, pair.cor)
    }
}
intxns <- cbind(kinases1, kinases2, cors)
intxns <- as.data.frame(intxns)
names(intxns) <- c("node1", "node2", "pair.cor")
intxns <- intxns[order(intxns$node1, intxns$node2),]
write.table(intxns, out.file, row.names=FALSE, col.names=TRUE, sep="\t",
            quote=FALSE)
