suppressPackageStartupMessages(library(data.table))
options(java.parameters = "-Xmx16g")
suppressPackageStartupMessages(library(bartMachine))
set_bart_machine_num_cores(24)
suppressPackageStartupMessages(library(ROCR))

argv <- commandArgs(TRUE)
if (length(argv) != 2){
    stop("USAGE: <script> MERGED_PREDS RAND_NEGS")
}

merged_pred_file <- argv[1]
use_rand_negs <- as.logical(argv[2])

merged_pred <- read.delim(merged_pred_file, as.is=TRUE)
names(merged_pred)[1:2] <- c("prot1", "prot2")
## merged_pred <- subset(merged_pred, prot1 != prot2)
rownames(merged_pred) <- paste(merged_pred$prot1, merged_pred$prot2, sep="-")

file_base <- strsplit(merged_pred_file, split="\\.")[[1]][1]
if (use_rand_negs){
    load(paste0(file_base, "-rand-negs-bart.Rdata"))
}else{
    load(paste0(file_base, "-bart.Rdata"))
}

preds <- predict(bart, merged_pred[3:ncol(merged_pred)])

out_tbl <- data.frame(prot1=merged_pred$prot1,
                      prot2=merged_pred$prot2,
                      bart.pred=preds)

if (use_rand_negs){
    write.table(out_tbl, paste0(file_base, "-rand-negs-bart.tsv"), sep="\t",
                quote=FALSE, col.names=TRUE, row.names=FALSE)
}else{
    write.table(out_tbl, paste0(file_base, "-bart.tsv"), sep="\t",
                quote=FALSE, col.names=TRUE, row.names=FALSE)
}
