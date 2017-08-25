suppressMessages(library(ROCR))
library(ggplot2)
library(reshape2)
library(dplyr)
library(cowplot)
source("src/rocr-helpers.r")
source("src/rocr-ggplot2.r")

argv <- commandArgs(TRUE)
if (length(argv) != 2){
    stop("USAGE: <script> DATA_FILE VAL_SET")
}

pred_score_file <- argv[1]
val_set_file <- argv[2]

sign_preds <- read.delim(pred_score_file, as.is=TRUE)
names(sign_preds)[1:2] <- c("node1", "node2")
rownames(sign_preds) <- paste(sign_preds$node1, sign_preds$node2, sep="_")
val_set <- read.table(val_set_file)
names(val_set) <- c("node1",  "node2", "label")
rownames(val_set) <- paste(val_set$node1, val_set$node2, sep="_")

pred_cols <- 3:ncol(sign_preds)

ns <- c()
for (i in 3:ncol(sign_preds)){
    col_name <- colnames(sign_preds)[i]
    pred_score <- sign_preds[, c(1, 2, i)]
    pred_score <- pred_score[complete.cases(pred_score),]
    if (!all(is.numeric(pred_score[,3]))){
        next
    }
    test_data <- merge(pred_score, val_set)
    test_data_act <- subset(test_data, label=="activates")
    test_data_inhib <- subset(test_data, label=="inhibits")

    n <- 0.9*min(nrow(test_data_act), nrow(test_data_inhib))
    ns <- c(ns, n)
}
n <- min(ns)

reps <- 100
preds <- vector("list", length(pred_cols))
for (i in pred_cols){
    col_name <- colnames(sign_preds)[i]
    pred_score <- sign_preds[, c(1, 2, i)]
    pred_score <- pred_score[complete.cases(pred_score),]
    if (!all(is.numeric(pred_score[,3]))){
        preds[[i-2]] <- NA
        next
    }
    test_data <- merge(pred_score, val_set)
    test_data_act <- subset(test_data, label=="activates")
    test_data_inhib <- subset(test_data, label=="inhibits")

    test_data_rand <- NULL
    test_labels_rand <- NULL
    ## n <- 0.9*min(nrow(test_data_act), nrow(test_data_inhib))
    for (j in 1:reps){
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
    pred <- prediction(test_data_rand, test_labels_rand)
    preds[[i-2]] <- pred
}

print(length(preds))

if (all(is.na(preds))){
    stop("No numeric predictors present")
}

data_basename <- strsplit(basename(pred_score_file), split = "\\.")[[1]][1]
data_basename <- paste(data_basename, "val", sep="-")
out_img <- paste0("img/", data_basename, ".pdf")
rds_basename <- paste0("img/rds/", data_basename)

head(melt(sign_preds[, c(1,2,pred_cols)],
          id.vars=c("node1", "node2")))

pred_score_melt <- melt(sign_preds[,c(1, 2, pred_cols)],
                        id.vars=c("node1", "node2"))[c("node1", "node2", "variable", "value")]
names(pred_score_melt)[3] <- "feature"

feature_names <- list(log10_hotspot_pval_min="log10(hotspot p-val.)",
                      disopred_score="DISOPRED score",
                      pos_in_domain="% pos. in kinase domain",
                      pos_perc="% pos. in protein",
                      cor.est.272="Wilkes signed coreg. score",
                      cor.est.293="CPTAC signed coreg. score",
                      signed.func.score="signed functional score",
                      signed.dcg="signed.dcg")

pred_score_melt <- subset(pred_score_melt, feature %in% names(feature_names))

pred_score_melt$feature <- factor(pred_score_melt$feature,
                                  levels=levels(pred_score_melt$feature),
                                  labels=feature_names[levels(pred_score_melt$feature)])

pdf(out_img, title=data_basename, width=7, height=7)

mccs <- build.perf1d.df(preds, "mat",
                        colnames(sign_preds)[pred_cols],
                        reps, TRUE)
mccs$feature <- factor(mccs$feature,
                       levels=levels(mccs$feature),
                       labels=feature_names[levels(mccs$feature)])
mcc_plot <- rocr2d.ggplot2(mccs, "\"activation\" cutoff",
                           "Matthews correlation coefficient", reps, NULL,
                           FALSE) +
    facet_wrap(~feature, scales="free_x") +
    guides(colour=FALSE)
saveRDS(mcc_plot, file=paste0(rds_basename, "-mcc.rds"))
write.table(mccs, paste0(rds_basename, "-mcc-data.tsv"), col.names=TRUE,
            row.names=FALSE, quote=FALSE, sep="\t")
print(mcc_plot)

dev.off()

## score_name <- gsub("\\.", ". ", col_name)

## rocr_pred <- prediction(test_data_rand, test_labels_rand)
## rocr_mcc <- performance(rocr_pred, "mat")
## rocr_balance <- performance(rocr_pred, measure="tpr", x.measure="tnr")
## avg_mccs <- perf_vert_avg(rocr_mcc)
## max_mcc_cutoff <- avg_mccs@x.values[[1]][which.max(avg_mccs@y.values[[1]])]
## max_mcc <- max(avg_mccs@y.values[[1]])
## message(score_name)
## message(paste("Max MCC =", format(max_mcc, digits=3)))
## avg_balance <- perf_thresh_avg(rocr_balance)
## approx_max_mcc_cutoff_i <- which.min(abs(avg_balance@alpha.values[[1]]-max_mcc_cutoff))
## max_mcc_tnr <- avg_balance@x.values[[1]][approx_max_mcc_cutoff_i]
## max_mcc_tpr <- avg_balance@y.values[[1]][approx_max_mcc_cutoff_i]

## data_basename <- strsplit(basename(pred_score_file), split = "\\.")[[1]][1]
## if (ncol(sign_preds) > 3){
##     data_basename <- paste(data_basename,
##                            gsub("\\.", "-", col_name),
##                            sep="-")
## }
## data_basename <- paste(data_basename, "val", sep="-")
## out_img <- paste0("img/", data_basename, ".pdf")

## pdf(out_img)
## par(cex=1.25, cex.main=0.8)
## plot(density(test_data_act[,col_name]))
## plot(rocr_mcc, spread.estimate="stderror", avg="vertical",
##      main=score_name)
## abline(v=max_mcc_cutoff, col="blue", lty=2)
## points(max_mcc_cutoff, max_mcc, pch=19, col="blue")
## text(max_mcc_cutoff, max_mcc*1.05, pos=4, col="blue",
##      labels=paste("max MCC =", format(max_mcc, digits=3)))
## legend("bottomleft",
##        legend=c(paste0("max-MCC cutoff: ", format(max_mcc_cutoff, digits=2))),
##        lty=2, col=c("blue", "red"), pch=19)
## plot(rocr_balance, spread.estimate="stderror", avg="threshold",
##      main=paste(score_name, "BART"), xlab="Average true inhibitory rate",
##      ylab="Average true activating rate")
## lines(c(max_mcc_tnr, max_mcc_tnr),
##       c(-0.03, max_mcc_tpr),
##       lty=2, col="blue")
## lines(c(-0.03, max_mcc_tnr),
##        c(max_mcc_tpr, max_mcc_tpr),
##        lty=2, col="blue")
## points(max_mcc_tnr, max_mcc_tpr, col="blue", pch=19)
## legend("bottomleft",
##        legend=c(paste0("max-MCC cutoff: ", format(max_mcc_cutoff, digits=2))),
##        lty=2, col=c("blue", "red"), pch=19)
## dev.off()

