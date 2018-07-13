suppressMessages(library(ROCR))
library(ggplot2)
source("src/rocr-helpers.r")

argv <- commandArgs(TRUE)
if (length(argv) != 3){
    stop("USAGE: <script> DATA_FILE VAL_SET COL_NUM")
}

pred_score_file <- argv[1]
val_set_file <- argv[2]
col_num <- as.integer(argv[3])

sign_preds <- read.delim(pred_score_file, as.is=TRUE)
names(sign_preds)[1:2] <- c("node1", "node2")
rownames(sign_preds) <- paste(sign_preds$node1, sign_preds$node2, sep="_")
val_set <- read.table(val_set_file)
names(val_set) <- c("node1",  "node2", "label")
rownames(val_set) <- paste(val_set$node1, val_set$node2, sep="_")

col_name <- colnames(sign_preds)[col_num]

test_data <- merge(sign_preds, val_set)
test_data <- test_data[which(!is.na(test_data[col_name])),]
test_data_act <- subset(test_data, label=="activates")
test_data_inhib <- subset(test_data, label=="inhibits")

test_data_rand <- NULL
test_labels_rand <- NULL
n <- 0.9*min(nrow(test_data_act), nrow(test_data_inhib))
reps <- 1000
for (i in 1:reps){
    rand_data <- c(sample(test_data_act[,col_name], n),
                   sample(test_data_inhib[,col_name], n))
    rand_labels <- c(rep(TRUE, n), rep(FALSE, n))
    if (is.null(test_data_rand)){
        test_data_rand <- rand_data
        test_labels_rand <- rand_labels
    }else{
        test_data_rand <- cbind(test_data_rand, rand_data)
        test_labels_rand <- cbind(test_labels_rand, rand_labels)
    }
}

score_name <- gsub("\\.", ". ", col_name)

rocr_pred <- prediction(test_data_rand, test_labels_rand)
rocr_mcc <- performance(rocr_pred, "mat")
rocr_balance <- performance(rocr_pred, measure="tpr", x.measure="tnr")
avg_mccs <- perf_vert_avg(rocr_mcc)
max_mcc_cutoff <- avg_mccs@x.values[[1]][which.max(avg_mccs@y.values[[1]])]
max_mcc <- max(avg_mccs@y.values[[1]])
message(score_name)
message(paste("Max MCC =", format(max_mcc, digits=3)))
avg_balance <- perf_thresh_avg(rocr_balance)
approx_max_mcc_cutoff_i <- which.min(abs(avg_balance@alpha.values[[1]]-max_mcc_cutoff))
max_mcc_tnr <- avg_balance@x.values[[1]][approx_max_mcc_cutoff_i]
max_mcc_tpr <- avg_balance@y.values[[1]][approx_max_mcc_cutoff_i]

data_basename <- strsplit(basename(pred_score_file), split = "\\.")[[1]][1]
if (ncol(sign_preds) > 3){
    data_basename <- paste(data_basename,
                           gsub("\\.", "-", col_name),
                           sep="-")
}
data_basename <- paste(data_basename, "val", sep="-")
out_img <- paste0("img/", data_basename, ".pdf")

pdf(out_img)
par(cex=1.25, cex.main=0.8)
plot(density(test_data_act[,col_name]))
plot(rocr_mcc, spread.estimate="stderror", avg="vertical",
     main=score_name)
abline(v=max_mcc_cutoff, col="blue", lty=2)
points(max_mcc_cutoff, max_mcc, pch=19, col="blue")
text(max_mcc_cutoff, max_mcc*1.05, pos=4, col="blue",
     labels=paste("max MCC =", format(max_mcc, digits=3)))
legend("bottomleft",
       legend=c(paste0("max-MCC cutoff: ", format(max_mcc_cutoff, digits=2))),
       lty=2, col=c("blue", "red"), pch=19)
plot(rocr_balance, spread.estimate="stderror", avg="threshold",
     main=paste(score_name, "BART"), xlab="Average true inhibitory rate",
     ylab="Average true activating rate")
lines(c(max_mcc_tnr, max_mcc_tnr),
      c(-0.03, max_mcc_tpr),
      lty=2, col="blue")
lines(c(-0.03, max_mcc_tnr),
       c(max_mcc_tpr, max_mcc_tpr),
       lty=2, col="blue")
points(max_mcc_tnr, max_mcc_tpr, col="blue", pch=19)
legend("bottomleft",
       legend=c(paste0("max-MCC cutoff: ", format(max_mcc_cutoff, digits=2))),
       lty=2, col=c("blue", "red"), pch=19)
dev.off()

