suppressPackageStartupMessages(library(data.table))
options(java.parameters = "-Xmx16g")
suppressPackageStartupMessages(library(bartMachine))
set_bart_machine_num_cores(24)
suppressPackageStartupMessages(library(ROCR))

argv <- commandArgs(TRUE)
if (length(argv) != 2){
    stop("USAGE: <script> MERGED_PREDS BART_MODEL")
}

merged_pred_file <- argv[1]
bart_model_file <- argv[2]

merged_pred <- read.delim(merged_pred_file, as.is=TRUE)
names(merged_pred)[1:2] <- c("prot1", "prot2")
## Remove missing data
merged_pred <- merged_pred[complete.cases(merged_pred),]
rownames(merged_pred) <- paste(merged_pred$prot1, merged_pred$prot2, sep="-")

load(bart_model_file)

## preds <- 1.0-predict(bart, merged_pred[3:ncol(merged_pred)])
preds <- 1.0-bart_machine_get_posterior(bart, merged_pred[3:ncol(merged_pred)])$y_hat

out_tbl <- data.frame(prot1=merged_pred$prot1,
                      prot2=merged_pred$prot2,
                      bart.pred=preds)

file_base <- strsplit(bart_model_file, split="\\.")[[1]][1]
write.table(out_tbl, paste0(file_base, "-preds.tsv"), sep="\t",
            quote=FALSE, col.names=TRUE, row.names=FALSE)
