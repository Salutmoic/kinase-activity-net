argv <- commandArgs(TRUE)
if (length(argv) != 1){
    stop("USAGE: <script> MERGED_BART_PREDS")
}

bart.preds.file <- argv[1]

bart.preds <- read.delim(bart.preds.file, as.is=TRUE)
rownames(bart.preds) <- paste(bart.preds$prot1, bart.preds$prot2, sep="-")

mean.preds <- apply(bart.preds[3:ncol(bart.preds)], 1, mean, na.rm=TRUE)
median.preds <- apply(bart.preds[3:ncol(bart.preds)], 1, median, na.rm=TRUE)
sd.preds <- apply(bart.preds[3:ncol(bart.preds)], 1, sd, na.rm=TRUE)
se.preds <- sd.preds/sqrt(ncol(bart.preds)-2)

bart.final <- data.frame(prot1=bart.preds$prot1,
                         prot2=bart.preds$prot2,
                         bart.pred.mean=mean.preds,
                         bart.pred.median=median.preds,
                         bart.pred.sd=sd.preds,
                         bart.pred.se=se.preds)


file.base <- strsplit(basename(bart.preds.file), split="\\.")[[1]][1]

pdf(paste0("img/", file.base, "-bart-pred-dists.pdf"))
## for (edge in c("BRAF-MAP2K1", "MAP2K1-MAPK1", "MAPK1-BRAF", "PDPK1-AKT1",
##                "MAPK1-AKT1", "AKT1-MAPK1", sample(rownames(bart.final), 10))){
for (edge in sample(rownames(bart.final), 10)){
    edge.preds <- t(bart.preds[edge, 3:ncol(bart.preds)])
    plot(density(edge.preds, from=0, to=1),
         xlim=c(0, 1), main=edge, xlab="BART edge probability")
}
dev.off()

write.table(bart.final,
            paste0("out/", file.base, "-final.tsv"),
            sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
