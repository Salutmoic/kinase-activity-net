library(dplyr)

argv <- commandArgs(TRUE)
if (length(argv) != 2){
    stop("USAGE: <script> NETWORK EXPR_SPECIFICITY")
}

network_file <- argv[1]
expr_spec_file <- argv[2]

network <- read.delim(network_file, as.is=TRUE)
expr_spec_full <- read.delim(expr_spec_file, as.is=TRUE)

expr_spec <- unique(expr_spec_full[,c(1,3)])

network_summ <- network %>%
    group_by(prot1) %>%
    summarise(med.assoc=median(assoc.pred.mean, na.rm=TRUE))

overrep_df <- merge(expr_spec, network_summ, by.x="node1", by.y="prot1")
overrep_df$tissue.broad <- 1.0 - overrep_df$tissue.sel1
## Is a kinase highly associated with most kinases, and is this most
## likely due to it being broadly expressed?  If so, it's
## overrepresented in the network.
overrep_df$overrep <- overrep_df$tissue.broad * overrep_df$med.assoc

## Remove the top 5% of the nodes based on their "overrepresentation".
## This is specifically for the purposes of clustering.
top0.05 <- tail(overrep_df[order(overrep_df$overrep),],
                n=round(0.05*nrow(overrep_df)))$node1
network_pruned <- subset(network, !(prot1 %in% top0.05 | prot2 %in% top0.05))

file_base <- strsplit(network_file, split="\\.")[[1]][1]

write.table(network_pruned, paste0(file_base, "-to-cluster.tsv"),
            quote=FALSE, sep="\t", col.names=TRUE, row.names=TRUE)
