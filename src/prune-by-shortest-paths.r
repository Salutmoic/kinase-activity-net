library(reshape2)
library(igraph)
library(parallel)

argv <- commandArgs(TRUE)
if (length(argv) != 1){
    stop("USAGE: <script> NETWORK")
}

network_file <- argv[1]

network <- read.delim(network_file, as.is=TRUE)

network$edge.cost <- 1/network$assoc.pred.mean
network$edge.cost[is.na(network$edge.cost)] <- Inf
network.m <- acast(network, prot1 ~ prot2, value.var="edge.cost")
graph <- graph_from_adjacency_matrix(network.m, mode="undirected", weighted=TRUE,
                                     diag=FALSE)

is.direct.edge <- function(pair, graph){
    kin1 <- pair[1]
    kin2 <- pair[2]
    path <- shortest_paths(graph, from=kin1, to=kin2, mode="out")
    return (length(path$vpath[[1]]) == 2)
}

cl <- makeCluster(31)
clusterEvalQ(cl, library(igraph))
direct.edges <- parApply(cl, combn(unique(network$prot1), 2), 2, is.direct.edge, graph)
stopCluster(cl)
network.prune <- network[direct.edges,]

file_base <- strsplit(network_file, split="\\.")[[1]][1]

save(direct.edges, file=paste0("tmp/", basename(file_base), "-direct-edges.Rdata")
write.table(network.prune, paste0(file_base, "-pruned.tsv"),
            quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

