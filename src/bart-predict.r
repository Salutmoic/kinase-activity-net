suppressPackageStartupMessages(library(data.table))
options(java.parameters = "-Xmx16g")
suppressPackageStartupMessages(library(bartMachine))
set_bart_machine_num_cores(24)
suppressPackageStartupMessages(library(ROCR))

argv <- commandArgs(TRUE)
if (length(argv) != 3){
    stop("USAGE: <script> MERGED_PREDS BART_MODEL_DIR N")
}

merged_pred_file <- argv[1]
bart_model_dir <- argv[2]
n <- as.integer(sub("^\\\\", "", argv[3]))

merged_pred <- read.delim(merged_pred_file, as.is=TRUE)
names(merged_pred)[1:2] <- c("prot1", "prot2")
## ## Remove missing data
## merged_pred <- merged_pred[complete.cases(merged_pred),]
rownames(merged_pred) <- paste(merged_pred$prot1, merged_pred$prot2, sep="-")

file_base <- sub("out/", "", bart_model_dir)

bart_model_file <- paste0(bart_model_dir, "/", file_base, "-bart-", n, ".Rdata")
load(bart_model_file)
preds <- 1.0-bart_machine_get_posterior(bart,
                                        merged_pred[3:ncol(merged_pred)])$y_hat

out_tbl <- data.frame(prot1=merged_pred$prot1,
                      prot2=merged_pred$prot2,
                      bart.pred=preds)

bad_rows <- which(apply(merged_pred[,3:ncol(merged_pred)], 1,
                        function(row){
                            all(is.na(row))
                        }))
out_tbl[bad_rows, "bart.pred"] <- NA

colnames(out_tbl)[3] <- paste0("bart.pred.", n)

write.table(out_tbl,
            paste0(bart_model_dir, "/", file_base, "-bart-", n, "-preds.tsv"),
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
