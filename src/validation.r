library(ROCR)

argv <- commandArgs(TRUE)
if (length(argv) != 2){
    stop("USAGE: <script> DATA_FILE VAL_SET")
}

pred.score.file <- argv[1]
val.set.file <- argv[2]

pred.score <- read.delim(pred.score.file, as.is=TRUE)
rownames(pred.score) <- paste(pred.score[,1], pred.score[,2], sep="-")

true.intxns.tbl <- read.table(val.set.file, as.is=TRUE)
true.intxns <- paste(true.intxns.tbl[,1], true.intxns.tbl[,2], sep="-")
true.intxns <- true.intxns[which(true.intxns %in% rownames(pred.score))]

all.kinases <- unique(c(true.intxns.tbl[,1], true.intxns.tbl[,2]))
possible.intxns <- as.vector(outer(all.kinases, all.kinases, paste, sep="-"))
possible.false.intxns <- setdiff(possible.intxns, true.intxns)

pred.data <- NULL
label.data <- NULL

n <- 50

## This procedure is sufficient for training-free predictions.
for (i in 1:n){
    false.intxns <- sample(possible.false.intxns, size=3*length(true.intxns))
    ## Put together the tables of predictions and TRUE/FALSE labels
    intxns <- c(true.intxns, false.intxns)
    preds <- pred.score[intxns, 3]
    labels <- c(rep(TRUE, length(true.intxns)), rep(FALSE, length(false.intxns)))
    if (is.null(pred.data)){
        pred.data <- preds
        label.data <- labels
    }else{
        pred.data <- cbind(pred.data, preds)
        label.data <- cbind(label.data, labels)
    }
}


pred <- prediction(pred.data, label.data)

data.basename <- strsplit(basename(pred.score.file), split="\\.")[[1]][1]
out.img <- paste("img/", data.basename, "-val.pdf")

pdf(out.img)
par(cex=1.25)

## ROC curve w/ AUC info
perf.roc <- performance(pred, measure="tpr", x.measure="fpr")
perf.auc <- performance(pred, "auc")
mean.auc <- mean(unlist(perf.auc@y.values))
se.auc <- sd(unlist(perf.auc@y.values))/sqrt(n)
message(paste("Mean AUC =", mean.auc, "S.E.M. =", se.auc))
plot(perf.roc, avg="vertical",
     main=paste(data.basename, paste0("Mean AUC=", mean.auc), sep="\n"))

## Precision-recall curve
perf.pr <- performance(pred, measure="prec", x.measure="rec")
plot(perf.pr, avg="vertical", main=data.basename)

dev.off()
