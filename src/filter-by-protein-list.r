argv <- commandArgs(TRUE)
if (length(argv) != 3){
    stop("USAGE: <script> IN_FILE PROT_LIST OUT_FILE")
}

in_file <- argv[1]
prot_list_file <- argv[2]
out_file <- argv[3]

in_data <- read.delim(in_file, as.is=TRUE)
prot_list <- read.table(prot_list_file, as.is=TRUE)

out_data <- in_data[which(in_data[,1] %in% prot_list$V1 &
                          in_data[,2] %in% prot_list$V1),]

write.table(out_data, out_file, sep="\t", quote=FALSE,
            row.names=FALSE, col.names=TRUE)
