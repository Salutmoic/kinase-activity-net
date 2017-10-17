library(reshape2)

argv <- commandArgs(TRUE)
if (length(argv) != 1){
    stop("USAGE: <script> KSEA_DATA")
}
ksea.data.file <- argv[1]
ksea.data <- read.table(ksea.data.file, sep="\t", as.is=TRUE)
names(ksea.data) <- c("kinase", "condition", "activity")

log.wKSEA.kinase_condition.clean <- acast(ksea.data, kinase ~ condition)

save(log.wKSEA.kinase_condition.clean,
     file="data/log.wKSEA.kinase_condition.clean.Rdata")
