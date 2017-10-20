suppressMessages(library(ROCR))

gen.possible.false.intxns <- function(all.kinases, prot.groups){
    possible.false.intxns <- c()
    for (k1 in all.kinases){
        for (k2 in all.kinases){
            if (k1 == k2)
                next
            if (!(k1 %in% prot.groups$protein) || !(k2 %in% prot.groups$protein))
                next
            k1.groups <- unique(subset(prot.groups, protein == k1)$group)
            k2.groups <- unique(subset(prot.groups, protein == k2)$group)
            if (length(k1.groups)>0 && length(k2.groups)>0 && length(intersect(k1.groups, k2.groups)) == 0)
                possible.false.intxns <- c(possible.false.intxns, paste(k1, k2, sep="-"))
        }
    }
    return (possible.false.intxns)
}

argv <- commandArgs(TRUE)
if (length(argv) != 3){
    stop("USAGE: <script> DATA_FILE VAL_SET PROTEIN_GROUPING")
}

pred.score.file <- argv[1]
val.set.file <- argv[2]
prot.groups.file <- argv[3]

## Get information from the prediction file name.  Build up a name of
## the prediction method for adding a title to plots, and figure out
## what kind of method was in use in order to determine any kind of
## normalisation that should be done before making predictions.
pred.name <- strsplit(basename(pred.score.file), split="\\.")[[1]][1]
pred.name <- sub("kinact-", "", pred.name)
assoc.method <- "[undefined]"
method <- "other"
if (grepl("pcor", pred.score.file)){
    assoc.method <- "Pearson's Correlation"
    method <- "cor"
}else if (grepl("scor", pred.score.file)){
    assoc.method <- "Spearman's Correlation"
    method <- "cor"
}else if (grepl("paircor", pred.score.file)){
    assoc.method <- "Pairwise Spearman's Correlation"
    method <- "cor"
}else if (grepl("nfchisq", pred.score.file)){
    assoc.method <- "FunChiSq"
    method <- "nfchisq"
}else if (grepl("fnn_mut_info", pred.score.file)){
    assoc.method <- "FNN Mutual Information"
    method <- "mutinfo"
}else if (grepl("mut_info", pred.score.file)){
    assoc.method <- "Mutual Information"
    method <- "mutinfo"
}else if (grepl("partcor", pred.score.file)){
    assoc.method <- "Partial Correlation"
    method <- "cor"
}else if (grepl("fvalue", pred.score.file)){
    assoc.method <- "ANOVA F value."
    method <- "meansq"
}else if (grepl("pssm", pred.score.file)){
    assoc.method <- "PSSM Score"
}else if (grepl("string-coexp", pred.score.file)){
    assoc.method <- "STRING Coexpression Score"
}else if (grepl("string-exper", pred.score.file)){
    assoc.method <- "STRING Experimental Score"
}else if (grepl("all", pred.score.file)){
    assoc.method <- "Mean Score"
}else{
    assoc.method <- pred.score.file
}
## The next two append info to the plot title and override for
## normalisation purposes the method used, in order of increasing
## precedence.
if (grepl("filter", pred.score.file)){
    assoc.method <- paste0(assoc.method, " (filtered)")
    method <- "filtered"
}
if (grepl("final-predictor", pred.score.file)){
    assoc.method <- paste0(assoc.method, " (merged)")
    method <- "predictor"
}
if (grepl("balanced", pred.score.file)){
    table.method <- "balanced table"
}else if (grepl("max-rows", pred.score.file)){
    table.method <- "max kinases"
}else if (grepl("max-cols", pred.score.file)){
    table.method <- "max conditions"
}else{
    table.method <- ""
}

pred.score <- read.delim(pred.score.file, as.is=TRUE)
names(pred.score) <- c("prot1", "prot2", "pred.score")
pred.score <- subset(pred.score, prot1 != prot2)
rownames(pred.score) <- paste(pred.score$prot1, pred.score$prot2, sep="-")

if (method %in% c("cor", "nfchisq")){
    min.pred <- min(pred.score$pred.score)
    max.pred <- max(pred.score$pred.score)
    pred.score$pred.score <- (pred.score$pred.score-min.pred)/(max.pred-min.pred)
}else if (method %in% c("mutinfo", "fvalue")){
    pred.score$pred.score <- pred.score$pred.score/max(pred.score$pred.score)
}


true.intxns.tbl <- read.table(val.set.file, as.is=TRUE)
possible.true.intxns <- paste(true.intxns.tbl[,1], true.intxns.tbl[,2], sep="-")
possible.true.intxns <- intersect(possible.true.intxns, rownames(pred.score))

prot.groups.full <- read.table(prot.groups.file, as.is=TRUE)
names(prot.groups.full) <- c("group", "protein")

all.kinases <- unique(c(pred.score$prot1, pred.score$prot2))
prot.groups <- subset(prot.groups.full, protein %in% all.kinases)
possible.false.intxns <- gen.possible.false.intxns(all.kinases, prot.groups)
## Remove true interactions from the set of possible false interactions
possible.false.intxns <- setdiff(possible.false.intxns, possible.true.intxns)
## Keep only false interactions for which we have predictions
possible.false.intxns <- intersect(possible.false.intxns, rownames(pred.score))

rev.true.intxns <- paste(true.intxns.tbl[,2], true.intxns.tbl[,1], sep="-")
rev.true.intxns <- intersect(rev.true.intxns, rownames(pred.score))
rev.true.intxns <- setdiff(rev.true.intxns, possible.true.intxns)
rev.true.intxns <- setdiff(rev.true.intxns, possible.false.intxns)

pred.data <- NULL
label.data <- NULL

print(c(length(possible.false.intxns), length(possible.true.intxns), length(rev.true.intxns)))

n <- 100

## sample.size <- 0.5*min(length(possible.false.intxns), length(possible.true.intxns))
sample.size <- round(0.8*min(length(possible.false.intxns)+0.5*length(rev.true.intxns), length(possible.true.intxns)))

## This procedure is sufficient for training-free predictions.
for (i in 1:n){
    ## false.intxns <- sample(possible.false.intxns, size=sample.size)
    false.intxns <- sample(union(possible.false.intxns,
                                 sample(rev.true.intxns, size=as.integer(0.5*length(rev.true.intxns)))),
                                 size=sample.size)
    true.intxns <- sample(possible.true.intxns, size=sample.size)
    ## Put together the tables of predictions and TRUE/FALSE labels
    intxns <- c(true.intxns, false.intxns)
    preds <- pred.score[intxns, "pred.score"]
    labels <- c(rep(TRUE, length(true.intxns)), rep(FALSE, length(false.intxns)))
    if (is.null(pred.data)){
        pred.data <- preds
        label.data <- labels
    }else{
        pred.data <- cbind(pred.data, preds)
        label.data <- cbind(label.data, labels)
    }
    if (i==1){
        out.tbl <- cbind(intxns, pred.data, label.data)
        write.table(out.tbl, "tmp/validation-test-set.tsv", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
    }
}

pred <- prediction(pred.data, label.data)

data.basename <- strsplit(basename(pred.score.file), split="\\.")[[1]][1]
out.img <- paste0("img/", data.basename, "-val.pdf")

pdf(out.img)
par(cex=1.25, cex.main=0.8)

## ROC curve w/ AUC info
perf.roc <- performance(pred, measure="tpr", x.measure="fpr")
perf.auc <- performance(pred, "auc")
mean.auc <- mean(unlist(perf.auc@y.values))
se.auc <- sd(unlist(perf.auc@y.values))/sqrt(n)
message(assoc.method)
message(paste("Mean AUC =", format(mean.auc, digits=2)))
message(paste("S.E.M. =", format(se.auc, digits=2)))
message(paste("n =", n))
message(paste("sample size =", sample.size))
plot(perf.roc, avg="vertical", spread.estimate="boxplot",
     main=paste(assoc.method,
                table.method,
                paste0("Mean AUC=", format(mean.auc, digits=2)),
                paste0("S.E.M.=", format(se.auc, digits=2)),
                paste0("Sample size=", sample.size),
                sep="\n"))
## Precision-recall curve
## perf.pr <- performance(pred, measure="prec", x.measure="rec")
## plot(perf.pr, avg="vertical", spread.estimate="boxplot",
##      main=paste(assoc.method, table.method, sep="\n"))

dev.off()
