argv <- commandArgs(TRUE)
if (length(argv) != 1){
    stop("USAGE: <script> NETWORK")
}

net_file <- argv[1]

network <- read.delim(net_file, as.is=TRUE)

egf_kins <- c("EGFR", "PDPK1", "ARAF", "BRAF", "RAF1", "MTOR",
              "MAP2K1", "MAP2K2", "MAP2K3", "AKT1", "AKT2", "AKT3",
              "MAPK1", "MAPK3", "GSK3A", "GSK3B", "RPS6KB1",
              "RPS6KB2", "RPS6KA3", "RPS6KA4")

cell_cyc_kins <- c("ATM", "ATR", "CHK1", "CHK2", "GSK3B", "CDK4",
                   "CDK6", "CDK2", "CDK7", "CDK1", "WEE1", "WEE2",
                   "CDC7", "BUB1", "BUB1B", "PLK1", "EGFR", "RAF1",
                   "MAP2K1", "MAP2K2", "MAPK1", "MAPK3", "PRKDC",
                   "CHEK1", "CHEK2", "TTK")

kinome <- read.table("data/human-kinome.txt", as.is=TRUE)$V1
mapks <- kinome[which(sapply(kinome, function(kin) grepl("MAPK[0-9]+", kin)))]
map2ks <- kinome[which(sapply(kinome, function(kin) grepl("MAP2K[0-9]+", kin)))]
map3ks <- kinome[which(sapply(kinome, function(kin) grepl("MAP3K[0-9]+", kin)))]
map4ks <- kinome[which(sapply(kinome, function(kin) grepl("MAP4K[0-9]+", kin)))]
map_kins <- c(mapks, map2ks, map3ks, map4ks)

egf_net <- subset(network, prot1 %in% egf_kins & prot2 %in% egf_kins)
map_net <- subset(network, prot1 %in% map_kins & prot2 %in% map_kins)
cell_cyc_net <- subset(network, prot1 %in% cell_cyc_kins & prot2 %in% cell_cyc_kins)


file_base <- strsplit(net_file, split="\\.")[[1]][1]

write.table(egf_net, paste0(file_base, "-egf.tsv"), quote=FALSE, sep="\t",
            row.names=FALSE, col.names=TRUE)
write.table(map_net, paste0(file_base, "-mapk.tsv"), quote=FALSE, sep="\t",
            row.names=FALSE, col.names=TRUE)
write.table(cell_cyc_net, paste0(file_base, "-cell-cyc.tsv"), quote=FALSE, sep="\t",
            row.names=FALSE, col.names=TRUE)
