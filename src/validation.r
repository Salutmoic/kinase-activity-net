library(ROCR)

argv <- commandArgs(TRUE)
if (length(argv) != 2){
    stop("USAGE: <script> DATA_FILE VAL_SET")
}

pred.score.file <- argv[1]
val.set.file <- argv[2]

pred.score <- read.delim(pred.score.file, as.is=TRUE)
rownames(pred.score) <- paste(pred.score[,1], pred.score[,2], sep="-")

true.intxns.tbl <- read.table(val.set.file)
true.intxns <- paste(true.intxns.tbl[,1], true.intxns.tbl[,2], sep="-")
true.intxns <- true.intxns[which(true.intxns %in% rownames(pred.score))]

pred.data <- NULL
label.data <- NULL

## This procedure is sufficient for training-free predictions.
for (i in 1:100){
    false.intxns <- c()
    ## Build up a set of random interactions that are not in the
    ## true-interaction set and have a prediction score
    while (length(false.pos.intxns) < 3*length(true.pos.intxns)){
        false.intxn.pair <- c(true.intxns.tbl[sample(1:nrow(true.intxns.tbl)), 1],
                              true.intxns.tbl[sample(1:nrow(true.intxns.tbl)), 2])
        if (false.intxn.pair[1] == false.intxn.pair[2])
            next
        false.intxn <- paste(false.intxn.pair[1], false.intxn.pair[2], sep="-")
        if (false.intxn %in% false.intxns || !(false.intxn %in% rownames(pred.score)))
            next
        false.intxns <- c(false.intxns, false.intxn)
    }
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
plot(perf.roc, avg="vertical",
     main=paste(data.basename, paste0("Mean AUC=", mean.auc), sep="\n"))

## Precision-recall curve
perf.pr <- performance(pred, measure="prec", x.measure="rec")
plot(perf.pr, avg="vertical", main=data.basename)

dev.off()
