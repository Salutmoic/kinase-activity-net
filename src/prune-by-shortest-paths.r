library(reshape2)
library(igraph)
library(parallel)

network <- read.delim("out/full-annotated-predictor.tsv", as.is=TRUE)
cluster <- read.table("tmp/cluster1.txt", as.is=TRUE)$V1

network$direct.pred.alt <- network$direct.pred.mean
network$direct.pred.alt[is.na(network$direct.pred.mean)] <- 0.1 ## mean(network$direct.pred.mean, na.rm=TRUE)

## network$edge.pred <- apply(network[,c("assoc.pred.mean", "direct.pred.alt")], 1,
##                            weighted.mean, w=c(0.2, 0.8))
## network$edge.pred <- network$assoc.pred.mean * network$direct.pred.alt

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

save(direct.edges, file="tmp/direct-edges.Rdata")
write.table(network.prune, "out/full-annotated-predictor-pruned.tsv",
            quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

## (length(all.edges)-length(direct.edges))/length(all.edges)
## cluster.net.direct <- cluster.net[direct.edges,]
## table(cluster.net$assoc.is.true & cluster.net$direct.is.true)/nrow(cluster.net)
## table(cluster.net.direct$assoc.is.true & cluster.net.direct$direct.is.true)/nrow(cluster.net.direct)


## cluster.net <- subset(network, prot1 %in% cluster & prot2 %in% cluster)
## all.edges <- paste(cluster.net$prot1, cluster.net$prot2)
## rownames(cluster.net) <- all.edges
## cluster.net$assoc.pred.mean[cluster.net$assoc.pred.mean < 0.8] <- NA
## cluster.net.min <- min(cluster.net$assoc.pred.mean, na.rm=TRUE)
## cluster.net.max <- max(cluster.net$assoc.pred.mean, na.rm=TRUE)
## cluster.net$edge.pred <- (cluster.net$assoc.pred.mean-cluster.net.min + 1e-3)/(cluster.net.max-cluster.net.min + 1e-3)
## cluster.net$edge.cost <- 1/cluster.net$assoc.pred.mean
## cluster.net$edge.cost[is.na(cluster.net$edge.cost)] <- Inf
## cluster.net.m <- acast(cluster.net, prot1 ~ prot2, value.var="edge.cost")
## cluster.graph <- graph_from_adjacency_matrix(cluster.net.m, mode="undirected",
##                                              weighted=TRUE, diag=FALSE)
## direct.edges <- c()
## for (kin1 in cluster){
##     for (kin2 in cluster){
##         if (kin1==kin2)
##             next
##         path <- shortest_paths(cluster.graph, from=kin1, to=kin2, mode="out")
##         if (length(path$vpath[[1]]) == 2)
##             direct.edges <- c(direct.edges, paste(kin1, kin2))
##     }
## }
## (length(all.edges)-length(direct.edges))/length(all.edges)
## cluster.net.direct <- cluster.net[direct.edges,]
## table(cluster.net$assoc.is.true & cluster.net$direct.is.true)/nrow(cluster.net)
## table(cluster.net.direct$assoc.is.true & cluster.net.direct$direct.is.true)/nrow(cluster.net.direct)
