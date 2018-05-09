options(java.parameters = "-Xmx16g")
suppressPackageStartupMessages(library(bartMachine))
set_bart_machine_num_cores(1)
suppressPackageStartupMessages(library(ROCR))
library(reshape2)
source("src/rocr-helpers.r")

argv <- commandArgs(TRUE)
if (length(argv) != 4){
    stop("USAGE: <script> MERGED_SIGN_PREDS EDGE_PREDS SIGN_GUESSES VAL_SET")
}

merged_pred_file <- argv[1]
edge_pred_file <- argv[2]
sign_guess_file <- argv[3]
val_set_file <- argv[4]

merged_pred <- read.delim(merged_pred_file, as.is=TRUE)
names(merged_pred)[1:2] <- c("prot1", "prot2")

edge_pred <- read.delim(edge_pred_file, as.is=TRUE)

sign_guess <- read.delim(sign_guess_file, as.is=TRUE)
rownames(sign_guess) <- sign_guess$kinase

good_rows <- which(apply(merged_pred[,3:ncol(merged_pred)], 1,
                        function(row){
                            any(!is.na(row))
                        }))
merged_pred <- merged_pred[good_rows,]
merged_pred <- subset(merged_pred, prot1 != prot2)

rownames(merged_pred) <- paste(merged_pred$prot1, merged_pred$prot2, sep="-")

## Temporary fix:
if (any(is.infinite(merged_pred$reg.site.cor.est))){
    inf.reg.site.cor <- which(is.infinite(merged_pred$reg.site.cor.est))
    max.noninf.reg.site.cor <- max(subset(merged_pred, !is.infinite(reg.site.cor.est))$reg.site.cor.est, na.rm=TRUE)
    merged_pred$reg.site.cor.est[inf.reg.site.cor] <- max.noninf.reg.site.cor*sign(merged_pred$reg.site.cor.est[inf.reg.site.cor])
}

edge_pred_m <- acast(edge_pred, prot1 ~ prot2, value.var="bart.pred.mean")
edge_pred_m <- edge_pred_m[,rownames(edge_pred_m)]
reg_site_cor_m <- acast(merged_pred, prot1 ~ prot2, value.var="cor.est.293")
good_kins <- intersect(rownames(reg_site_cor_m), rownames(edge_pred_m))
good_kins <- intersect(colnames(reg_site_cor_m), good_kins)
reg_site_cor_m <- reg_site_cor_m[good_kins, good_kins]

reg_site_cor_med <- sapply(rownames(reg_site_cor_m),
                           function(kinase){
                               if (!(kinase %in% rownames(edge_pred_m)))
                                   return(NA)
                               kinase_cors <- reg_site_cor_m[kinase,]
                               kinase_edges <- edge_pred_m[kinase, colnames(reg_site_cor_m)]
                               return(weighted.mean(x=kinase_cors,
                                                    w=kinase_edges,
                                                    na.rm=TRUE))
                           })

reg_site_cor_med[is.nan(reg_site_cor_med)] <- NA

## reg_site_new_pred_m <- t(apply(edge_pred_m, 1,
##                                function(row) return(row*reg_site_cor_med)))
## reg_site_new_pred <- melt(reg_site_new_pred_m)
## names(reg_site_new_pred) <- c("prot1", "prot2", "reg.site.sign")
## merged_pred <- merge(merged_pred, reg_site_new_pred)
merged_pred$two.step.cors <- reg_site_cor_med[merged_pred$prot2]

## merged_pred$sign.guess <- sapply(rownames(merged_pred),
##                                      function(pair){
##                                          prot2 <- new_merged_pred[pair, "prot2"]
##                                          guess <- sign_guess[prot2, "sign.guess"]
##                                          edge <- edge_pred[pair, "bart.pred"]
##                                          return(guess*edge)
##                                      })
merged_pred$sign.guess <- sign_guess[merged_pred$prot2, "sign.guess"]
rownames(merged_pred) <- paste(merged_pred$prot1, merged_pred$prot2,
                                   sep="-")

val_set_tbl <- read.table(val_set_file, as.is = TRUE)
val_set_tbl[,3] <- as.factor(val_set_tbl[,3])
rownames(val_set_tbl) <- paste(val_set_tbl[,1], val_set_tbl[,2], sep="-")

good_pairs <- intersect(rownames(merged_pred), rownames(val_set_tbl))

val_pred_full <- merged_pred[good_pairs,]
val_set_tbl <- val_set_tbl[good_pairs,]

## ## Remove missing data
## new_merged_pred <- new_merged_pred[complete.cases(new_merged_pred),]

samp_size = min(nrow(subset(val_set_tbl, V3=="inhibits")),
                nrow(subset(val_set_tbl, V3=="activates")))
random_pairs <- c(sample(rownames(subset(val_set_tbl, V3=="inhibits")),
                         samp_size),
                  sample(rownames(subset(val_set_tbl, V3=="activates")),
                         samp_size))
random_pairs <- sample(random_pairs)

## init_feats <- c("signed.dcg", "signed.func.score")## , "kinact.cor.est",
##                 ## "two.step.cors", "sign.guess")
init_feats <- colnames(val_pred_full)[3:ncol(val_pred_full)]

preds <- val_pred_full[random_pairs, init_feats]
labels <- val_set_tbl[random_pairs, 3]

sink("/dev/null")
bart_init <- bartMachineCV(preds, labels,
                           use_missing_data=TRUE,
                           use_missing_data_dummies_as_covars=TRUE,
                           replace_missing_data_with_x_j_bar=FALSE,
                           impute_missingness_with_x_j_bar_for_lm=FALSE,
                           prob_rule_class=0.5,
                           verbose=FALSE,
                           serialize=FALSE)
sink()

var_select <- var_selection_by_permute(bart_init)

var_select$important_vars_local_names

## final_feats <- colnames(preds)
## final_feats <- c("signed.dcg", "signed.func.score", "kinase.type",
##                  "sub.kinase.type")
## final_feats <- c("signed.dcg", "signed.func.score", "kinact.cor.est.293",
##                  "kinase.type", "sub.kinase.type", "cor.est.293")
final_feats <- c("signed.dcg", "signed.func.score",
                 "kinase.type", "sub.kinase.type", "cor.est.293")
## final_feats <- init_feats

k <- 3
n <- 20
all_probs <- NULL
all_labels <- NULL
for (j in 1:n){
    message(paste("Run:", j))
    random_pairs <- c(sample(rownames(subset(val_set_tbl, V3=="inhibits")),
                             samp_size),
                      sample(rownames(subset(val_set_tbl, V3=="activates")),
                             samp_size))
    random_pairs <- sample(random_pairs)
    preds <- val_pred_full[random_pairs, final_feats]
    labels <- val_set_tbl[random_pairs, 3]
    val_n <- ceiling(nrow(preds)/k)
    train_n <- nrow(preds) - val_n
    for (i in 1:k){
        val_i <- (val_n*(i-1))+1
        if (i < k){
            val_rows <- val_i:(val_i + val_n - 1)
        }else{
            val_rows <- val_i:nrow(preds)
        }
        train_rows <- setdiff(1:nrow(preds), val_rows)
        val_set <- preds[val_rows,]
        val_labels <- labels[val_rows]
        train_set <- preds[train_rows,]
        train_labels <- labels[train_rows]
        sink("/dev/null")
        bart <- bartMachineCV(train_set,
                              train_labels,
                              use_missing_data=TRUE,
                              use_missing_data_dummies_as_covars=TRUE,
                              replace_missing_data_with_x_j_bar=FALSE,
                              impute_missingness_with_x_j_bar_for_lm=FALSE,
                              prob_rule_class=0.5,
                              verbose=FALSE,
                              serialize=FALSE)
        ## ROCR assumes that the first factor level corresponds to
        ## lower scores, but BART assigns interactions higher
        ## probabilities to the first level.  So, to accommodate ROCR,
        ## here we're taking 1.0-prob
        val_preds <- 1.0-bart_machine_get_posterior(bart, val_set)$y_hat
        sink()
        if (is.null(all_probs)){
            all_probs <- val_preds
            all_labels <- val_labels
        }else{
            if (length(val_preds) < val_n){
                num.na <- val_n - length(val_preds)
                val_preds <- c(val_preds, rep(NA, num.na))
                val_labels <- c(val_labels, rep(NA, num.na))
                val_labels <- as.factor(val_labels)
                levels(val_labels) <- levels(labels)
            }
            all_probs <- cbind.data.frame(all_probs, val_preds)
            all_labels <- cbind.data.frame(all_labels, val_labels)
        }
    }
}

rocr_pred <- prediction(all_probs, all_labels)
rocr_roc <- performance(rocr_pred, measure="tpr", x.measure="fpr")
rocr_roc2 <- performance(rocr_pred, measure="tnr", x.measure="fnr")
rocr_prec <- performance(rocr_pred, measure="prec", x.measure="rec")
rocr_auc <- performance(rocr_pred, "auc")
rocr_mcc <- performance(rocr_pred, "mat")
mean_auc <- mean(unlist(rocr_auc@y.values), na.rm=TRUE)
se_auc <- sd(unlist(rocr_auc@y.values), na.rm=TRUE)/sqrt(k)
act.fprs <- rocr_roc@x.values
act.fpr0.1s <- unlist(lapply(act.fprs,
                        function(f){
                            fpr0.1 <-max(which(f<0.1))
                            rocr_roc@alpha.values[[1]][fpr0.1]
                        }))
act.cutoff <- mean(act.fpr0.1s[!is.infinite(act.fpr0.1s)])
inhib.fprs <- rocr_roc2@x.values
inhib.fpr0.1s <- unlist(lapply(inhib.fprs,
                        function(f){
                            fpr0.1 <-max(which(f<0.1))
                            rocr_roc2@alpha.values[[1]][fpr0.1]
                        }))
inhib.cutoff <- mean(inhib.fpr0.1s[!is.infinite(inhib.fpr0.1s)], na.rm=TRUE)
rocr_balance <- performance(rocr_pred, measure="tpr", x.measure="tnr")
avg_mccs <- perf_vert_avg(rocr_mcc)
max_mcc_cutoff <- avg_mccs@x.values[[1]][which.max(avg_mccs@y.values[[1]])]
max_mcc <- max(avg_mccs@y.values[[1]])
print(max_mcc)
avg_balance <- perf_thresh_avg(rocr_balance)
approx_max_mcc_cutoff_i <- which.min(abs(avg_balance@alpha.values[[1]]-max_mcc_cutoff))
max_mcc_tnr <- avg_balance@x.values[[1]][approx_max_mcc_cutoff_i]
max_mcc_tpr <- avg_balance@y.values[[1]][approx_max_mcc_cutoff_i]

file_base <- strsplit(basename(merged_pred_file), split="\\.")[[1]][1]
img_file <- paste0(file_base, "-bart-sign")
img_file <- paste0("img/", img_file, "-roc.pdf")
pdf(img_file)
plot(rocr_roc, spread.estimate="stderror", avg="vertical",
     xlab="False 'activation' rate", ylab="Avg. true 'activation' rate",
     main=paste("FPR=0.1 cutoff =", format(act.cutoff, digits=3)))
abline(v=0.1, lty=2, col="blue")
plot(rocr_roc2, spread.estimate="stderror", avg="vertical",
     xlab="False 'inhibition' rate", ylab="Avg. true 'inhibition' rate",
     main=paste("FPR=0.1 cutoff =", format(inhib.cutoff, digits=3)))
abline(v=0.1, lty=2, col="blue")
plot(rocr_mcc, spread.estimate="stderror", avg="vertical")
abline(v=max_mcc_cutoff, col="blue", lty=2)
points(max_mcc_cutoff, max_mcc, pch=19, col="blue")
text(max_mcc_cutoff, max_mcc*1.05, pos=4, col="blue",
     labels=paste("max MCC =", format(max_mcc, digits=3)))
legend("bottomleft",
       legend=c(paste0("max-MCC cutoff: ", format(max_mcc_cutoff, digits=2))),
       lty=2, col=c("blue", "red"), pch=19)
plot(rocr_balance, spread.estimate="stderror", avg="threshold", colorize=TRUE,
     main=paste0("BART"), xlab="Average true inhibitory rate",
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
investigate_var_importance(bart)
dev.off()

samp_size = min(nrow(subset(val_set_tbl, V3=="inhibits")),
                nrow(subset(val_set_tbl, V3=="activates")))

n <- 20
prob_tbl <- NULL
for (i in 1:n){
    print(i)
    random_pairs <- c(sample(rownames(subset(val_set_tbl, V3=="inhibits")),
                             samp_size),
                      sample(rownames(subset(val_set_tbl, V3=="activates")),
                             samp_size))
    random_pairs <- sample(random_pairs)
    preds <- merged_pred[random_pairs, init_feats]
    labels <- val_set_tbl[random_pairs, 3]
    bart <- bartMachineCV(preds, labels,
                          use_missing_data=TRUE,
                          use_missing_data_dummies_as_covars=TRUE,
                          replace_missing_data_with_x_j_bar=FALSE,
                          impute_missingness_with_x_j_bar_for_lm=FALSE,
                          prob_rule_class=0.5,
                          verbose=FALSE,
                          serialize=FALSE)
    ## new_preds <- predict(bart, merged_pred[colnames(preds)],
    ##                      type="class", prob_rule_class=0.5)
    new_probs <- bart_machine_get_posterior(bart, merged_pred[colnames(preds)])$y_hat
    if (is.null(prob_tbl)){
        prob_tbl <- new_probs
    }else{
        prob_tbl <- cbind(prob_tbl, new_probs)
    }
}

prob_means <- apply(prob_tbl, 1, mean, na.rm=TRUE)
prob_sd <- apply(prob_tbl, 1, sd, na.rm=TRUE)

out_tbl <- data.frame(prot1=merged_pred$prot1,
                      prot2=merged_pred$prot2,
                      prob.act.mean=prob_means
                      prob.act.sd=prob_sd)
out_tbl <- out_tbl[order(out_tbl$prot1, out_tbl$prot2),]

out_file <- paste0(file_base, "-bart-sign-preds")
out_file <- paste0("out/", out_file, ".tsv")

write.table(out_tbl, out_file, quote=FALSE,
            row.names=FALSE, col.names=TRUE, sep="\t")
