argv <- commandArgs(TRUE)
if (length(argv) != 2){
    stop("USAGE: <script> NETWORK CLUSTS")
}

net_file <- argv[1]
clusts_file <- argv[2]

network <- read.delim(net_file, as.is=TRUE)
clusts <- read.table(clusts_file, as.is=TRUE)
names(clusts) <- c("prot", "clust")

file_base <- strsplit(net_file, split="\\.")[[1]][1]

for (i in 1:max(clusts$clust)){
    clust <- subset(clusts, clust==i)
    clust.net <- subset(network, prot1 %in% clust$prot & prot2 %in% clust$prot)
    write.table(clust.net, paste0(file_base, "-clust", i, ".tsv"), quote=FALSE,
                sep="\t", row.names=FALSE, col.names=TRUE)
}
