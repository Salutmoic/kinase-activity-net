argv <- commandArgs(TRUE)
if (length(argv) != 3){
    stop("USAGE: <script> FINAL_NETWORK OUT_FILE ASSOC_CUTOFF")
}

network_file <- argv[1]
out_file <- argv[2]
cutoff <- as.numeric(argv[3])

network <- read.delim(network_file, as.is=TRUE)

network$direct.pred.mean[network$assoc.pred.mean < cutoff] <- NA
network$direct.pred.sd[network$assoc.pred.mean < cutoff] <- NA
network$prob.act.mean[network$assoc.pred.mean < cutoff] <- NA
network$prob.act.sd[network$assoc.pred.mean < cutoff] <- NA
network$assoc.pred.sd[network$assoc.pred.mean < cutoff] <- NA
network$assoc.pred.mean[network$assoc.pred.mean < cutoff] <- NA


write.table(network, out_file, quote=FALSE, sep="\t", row.names=FALSE,
            col.names=TRUE)
