suppressMessages(library(VIM))
suppressMessages(library(reshape2))
source("src/optimize-activity-table.r")

argv <- commandArgs(TRUE)
if (length(argv) != 2){
    stop("USAGE: <script> IN_RDS_FILE PROT_LIST OUT_FILE")
}
in_file <- argv[1]
prot_list_file <- argv[2]
out_file <- argv[3]

in_data <- as.matrix(readRDS(in_file))

in_data_imp <- impute.on.table(in_data)

in_data_df <- melt(in_data_imp)
colnames(in_data_df) <- c("gene", "condition", "rnai")

prot_list <- read.table(prot_list_file)
in_data_df <- in_data_df[which(in_data_df$gene %in% prot_list),]

write.table(in_data_df, out_file, sep="\t", quote=FALSE, row.names=FALSE,
            col.names=TRUE)
