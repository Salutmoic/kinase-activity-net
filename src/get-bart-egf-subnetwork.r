argv <- commandArgs(TRUE)
if (length(argv) != 3){
    stop("USAGE: <script> EGF_DATA BART_NET OUT_FILE")
}
egf.data.file <- argv[1]
bart.net.file <- argv[2]
out.file <- argv[3]

egf.data <- read.table(egf.data.file, as.is=TRUE)
bart.net <- read.delim(bart.net.file, as.is=TRUE)

egf.kins <- unique(egf.data[,1])

bart.net.sub <- subset(bart.net, prot1 %in% egf.kins & prot2 %in% egf.kins)

write.table(bart.net.sub, out.file, quote=FALSE, sep="\t", row.names=FALSE,
            col.names=TRUE)
