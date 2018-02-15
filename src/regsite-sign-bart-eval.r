suppressPackageStartupMessages(library(data.table))
options(java.parameters = "-Xmx16g")
suppressPackageStartupMessages(library(bartMachine))
set_bart_machine_num_cores(24)
suppressPackageStartupMessages(library(ROCR))
source("src/rocr-helpers.r")

get_pos_in_domain <- function(pos, dom_start, dom_end){
    if (pos < dom_start || pos > dom_end)
        return(NA)
    return((pos-dom_start)/(dom_end-dom_start))
}

reg_sites <- read.delim("data/reg_sites.tsv")
reg_sites$MOD_RSD <- sub("^[STY]", "", reg_sites$MOD_RSD)
reg_sites$MOD_RSD <- sub("-p$", "", reg_sites$MOD_RSD)
## reg_sites$ACC_ID <- sub("-[0-9]+", "", reg_sites$ACC_ID)
rownames(reg_sites) <- paste(reg_sites$ACC_ID, reg_sites$MOD_RSD, sep="_")

phosfun_feat <- read.table("data/external/phosphoproteome_annotated",
                           header=TRUE, row.names=1, sep=",")

kin.domains <- read.delim("data/human-pfam.tsv", as.is=TRUE, sep="\t")

phosfun <- read.table("data/phosfun.tsv", sep="\t")
names(phosfun) <- c("kinase", "position", "phosfun_score")

sign_feats <- merge(reg_sites, phosfun_feat,
                    by.x=c("ACC_ID", "MOD_RSD"),
                    by.y=c("acc", "position"))
rownames(sign_feats) <- paste(sign_feats$GENE, sign_feats$MOD_RSD, sep="_")
sign_feats$pos_perc <- sapply(1:nrow(sign_feats),
                              function(i){
                                  as.integer(sign_feats[i, "MOD_RSD"])/sign_feats[i, "prot_length"]
                              })
sign_feats$pos_in_domain <- apply(sign_feats, 1,
                                  function(row){
                                      kinase <- row[["GENE"]]
                                      pos <- as.integer(row[["MOD_RSD"]])
                                      if (!(kinase %in% kin.domains$protein)){
                                          return(NA)
                                      }
                                      domains <- subset(kin.domains,
                                                        protein==kinase)
                                      pos_in_domain <-
                                          apply(domains, 1,
                                                function(domain){
                                                    dom_start <- as.integer(domain[["start"]])
                                                    dom_end <- as.integer(domain[["end"]])
                                                    get_pos_in_domain(pos, dom_start, dom_end)
                                                })
                                      in_dom <- which(!is.na(pos_in_domain) &
                                                      pos_in_domain >= 0 &
                                                      pos_in_domain <= 1)
                                      if (length(in_dom)>0)
                                          return(pos_in_domain[in_dom])
                                      dist_from_dom <- sapply(pos_in_domain,
                                                              function(x){
                                                                  if (is.na(x))
                                                                      return(NA)
                                                                  if (x < 0){
                                                                      return(abs(x))
                                                                  }else{
                                                                      return(x-1)
                                                                  }
                                                              })
                                      if (all(is.na(dist_from_dom)))
                                          return(NA)
                                      return(pos_in_domain[which.min(dist_from_dom)])
                                  })
sign_feats$pfam <- as.factor(apply(sign_feats, 1,
                                  function(row){
                                      kinase <- row[["GENE"]]
                                      pos <- as.integer(row[["MOD_RSD"]])
                                      if (!(kinase %in% kin.domains$protein)){
                                          return(NA)
                                      }
                                      domains <- subset(kin.domains,
                                                        protein==kinase)
                                      pos_in_domain <-
                                          apply(domains, 1,
                                                function(domain){
                                                    dom_start <- as.integer(domain[["start"]])
                                                    dom_end <- as.integer(domain[["end"]])
                                                    return(pos >= dom_start && pos <= dom_end)
                                                })
                                      if (any(pos_in_domain)){
                                          dom <- domains[which(pos_in_domain), "domain"]
                                          if (dom %in% c("Pkinase", "Pkinase_Tyr"))
                                              return(dom)
                                      }
                                      return(NA)
                                  }))
sign_feats$is_tyr_kin <- apply(sign_feats, 1,
                               function(row){
                                   kinase <- row[["GENE"]]
                                   if (!(kinase %in% kin.domains$protein)){
                                       return(NA)
                                   }
                                   domains <- subset(kin.domains,
                                                     protein==kinase)
                                   return ("Pkinase_Tyr" %in% domains$domain)
                               })

sign_feats$DOMAIN[which(sign_feats$DOMAIN=="")] <- NA
sign_feats$comp_bias_uniprot[which(sign_feats$comp_bias_uniprot=="")] <- NA
sign_feats$domain_uniprot[which(sign_feats$domain_uniprot=="")] <- NA
sign_feats$motif_uniprot[which(sign_feats$motif_uniprot=="")] <- NA
sign_feats$region_uniprot[which(sign_feats$region_uniprot=="")] <- NA
sign_feats$repeat_uniprot[which(sign_feats$repeat_uniprot=="")] <- NA
sign_feats$zn_finger_uniprot[which(sign_feats$zn_finger_uniprot=="")] <- NA
sign_feats <- merge(sign_feats, phosfun, by.x=c("GENE", "MOD_RSD"),
                    by.y=c("kinase", "position"))
sign_feats <- sign_feats[which(grepl("activity, induced", sign_feats$ON_FUNCTION) |
                               grepl("activity, inhibited", sign_feats$ON_FUNCTION)),]

rownames(sign_feats) <- paste(sign_feats$GENE, sign_feats$MOD_RSD, sep="_")
rownames(phosfun) <- paste(phosfun$kinase, phosfun$position,  sep="_")
phosfun_nonreg <- phosfun[which(!(rownames(phosfun) %in% rownames(sign_feats))),]

non_feats <- c("ACC_ID", "MOD_RSD", "GENE", "PROTEIN", "PROT_TYPE", "GENE_ID",
               "HU_CHR_LOC", "ORGANISM", "SITE_GRP_ID", "SITE_...7_AA",
               "ON_FUNCTION", "ON_PROCESS", "ON_PROT_INTERACT",
               "ON_OTHER_INTERACT", "PMIDs", "LT_LIT", "MS_LIT", "MS_CST",
               "NOTES", "X", "mq_siteid", "isHotspot", "phosfun_score")
phospho_feats <- c("isELMkinaseMotif", "netpho_max_all",
                   "netpho_max_KIN", "netpho_max_PTP",
                   "netpho_max_STdomain", "netpho_max_Ydomain", "ptmdb_coreg_kinases",
                   "ptmdb_coreg_kinases_len", "ptmdb_maxcoreg_kinase", "PWM_max_mss",
                   "PWM_nkinTop005", "PWM_nkinTop01", "PWM_nkinTop02")
init_feats_full <- setdiff(colnames(sign_feats), c(non_feats, phospho_feats))

## init_feats_full <- c("disopred_score", "comp_bias_uniprot",
##                      "domain_uniprot", "motif_uniprot",
##                      "region_uniprot", "repeat_uniprot",
##                      "exp3d_acid_ddG_prob", "exp3d_ala_ddG_prob",
##                      "exp3d_max_ddG_prob", "exp3d_mean_ddG_prob",
##                      "log10_hotspot_pval_min", "isInterface",
##                      "intfc_max_mean_ddG_prob",
##                      "intfc_max_max_ddG_prob",
##                      "intfc_max_ala_ddG_prob",
##                      "intfc_max_acid_ddG_prob", "prot_length",
##                      "ACCpro", "SSpro", "SSpro8", "sift_min_score",
##                      "sift_mean_score", "sift_ala_score",
##                      "sift_acid_score", "topology_uniprot",
##                      "transitpep_uniprot", "pfam", "pos_in_domain")
init_feats <- names(which(apply(sign_feats[init_feats_full], 2,
                                function(col){
                                    length(which(is.na(col)))/length(col)
                                }) < 0.75))
preds <- sign_feats[init_feats]
rownames(preds) <- rownames(sign_feats)
## preds$ACCpro <- as.factor(preds$ACCpro)
## preds$SSpro <- as.factor(preds$SSpro)
## preds$SSpro8 <- as.factor(preds$SSpro8)
## preds$domain_uniprot <- as.factor(preds$domain_uniprot)
## preds$pfam <- as.factor(preds$pfam)
labels <- as.factor(sapply(1:nrow(sign_feats),
                                  function(n){
                                      act <- sign_feats$ON_FUNCTION[n]
                                      if (grepl("activity, induced", act)){
                                          return ("up")
                                      }else if (grepl("activity, inhibited", act)){
                                          return ("down")
                                      }else{
                                          return (NA)
                                      }
                                  }))
preds <- preds[which(!is.na(labels)),]
labels <- labels[which(!is.na(labels))]

eval_rows <- c(which(labels=="down"), sample(which(labels=="up"), 56))
labels_sub <- labels[eval_rows]
preds_sub <- preds[eval_rows,]

sink("/dev/null")
bart_init <- bartMachineCV(preds_sub, labels_sub,
                           use_missing_data=TRUE,
                           use_missing_data_dummies_as_covars=TRUE,
                           replace_missing_data_with_x_j_bar=FALSE,
                           impute_missingness_with_x_j_bar_for_lm=FALSE,
                           prob_rule_class=0.5,
                           verbose=FALSE,
                           serialize=FALSE)
sink()

var_select <- var_selection_by_permute(bart_init)

sign_feats$sign <- labels
pdf("img/regsite-sign-boxplots.pdf")
boxplot(log10_hotspot_pval_min ~ sign, data=sign_feats, ylab="Hotspot log10 pval")
boxplot(disopred_score ~ sign, data=sign_feats, ylab="DISOPRED score")
boxplot(pos_in_domain ~ sign, data=sign_feats, ylab="Position in domain")
boxplot(prot_length ~ sign, data=sign_feats, ylab="Protein length")
boxplot(sample(phosfun_nonreg$phosfun_score, 298),
        sign_feats[sign_feats$sign=="down", "phosfun_score"],
        sign_feats[sign_feats$sign=="up", "phosfun_score"],
        names=c("non-regulatory", "inhibiting", "activating"),
        ylab="Functional Score")
dev.off()

## final_feats <- c("log10_hotspot_pval_min", "exp3d_max_ddG_effect",
##                  "exp3d_ala_ddG_effect", "disopred_score",
##                  "pos_in_domain", "pfam")
final_feats <- c("log10_hotspot_pval_min", "disopred_score",
                 "pos_in_domain", "pos_perc", "pfam", "residue",
                 "is_tyr_kin")
## final_feats <- c("log10_hotspot_pval_min", "pos_in_domain", "pfam",
##                  "disopred_score")
## preds <- sign_feats[final_feats]
## rownames(preds) <- rownames(sign_feats)
## preds$region_uniprot <- as.factor(preds$region_uniprot)

k <- 5
n <- 5
all_probs <- NULL
all_labels <- NULL
for (j in 1:n){
    eval_rows <- c(sample(which(labels=="down"), 50),
                   sample(which(labels=="up"), 50))
    eval_rows <- sample(eval_rows)
    labels_sub <- labels[eval_rows]
    preds <- sign_feats[eval_rows, final_feats]
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
        val_labels <- labels_sub[val_rows]
        train_set <- preds[train_rows,]
        train_labels <- labels_sub[train_rows]
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
        val_preds <- 1.0-bart_machine_get_posterior(bart, val_set)$y_hat
        sink()
        if (is.null(all_probs)){
            all_probs <- val_preds
            all_labels <- val_labels
        }else{
            if (length(val_preds) < val_n){
                num.na <- val_n - length(val_preds)
                val_preds <- c(val_preds, rep(NA, num.na))
            }
            ## BART seems to be assigning P>cutoff to the first level, which is
            ## "activates", whereas ROCR expects P<cutoff to be the first level
            all_probs <- cbind(all_probs, val_preds)
            all_labels <- cbind(all_labels, val_labels)
        }
    }
}

rocr_pred <- prediction(all_probs, all_labels)
rocr_roc <- performance(rocr_pred, measure="tpr", x.measure="fpr")
rocr_prec <- performance(rocr_pred, measure="prec", x.measure="rec")
rocr_auc <- performance(rocr_pred, "auc")
rocr_mcc <- performance(rocr_pred, "mat")
mean_auc <- mean(unlist(rocr_auc@y.values), na.rm=TRUE)
se_auc <- sd(unlist(rocr_auc@y.values), na.rm=TRUE)/sqrt(k)
fprs <- rocr_roc@x.values
fpr0.05s <- unlist(lapply(fprs,
                        function(f){
                            fpr0.05 <-max(which(f<0.05))
                            rocr_roc@alpha.values[[1]][fpr0.05]
                        }))
cutoff <- mean(fpr0.05s)

rocr_balance <- performance(rocr_pred, measure="tpr", x.measure="tnr")

avg_mccs <- perf_vert_avg(rocr_mcc)
max_mcc_cutoff <- avg_mccs@x.values[[1]][which.max(avg_mccs@y.values[[1]])]
avg_balance <- perf_thresh_avg(rocr_balance)
approx_max_mcc_cutoff_i <- which.min(abs(avg_balance@alpha.values[[1]]-max_mcc_cutoff))
max_mcc_tnr <- avg_balance@x.values[[1]][approx_max_mcc_cutoff_i]
max_mcc_tpr <- avg_balance@y.values[[1]][approx_max_mcc_cutoff_i]

message(paste("FPR 0.05 cutoff: ", cutoff))
pdf("img/regsite-sign-bart-eval.pdf")
## plot(rocr_roc, spread.estimate="boxplot", avg="vertical",
##      main=paste0("BART\n(mean AUC=", format(mean_auc, digits=2),
##                  ", S.E.M=", format(se_auc, digits=2), ")"))
## abline(v=0.05, lty=2, col="blue")
## legend("bottomright",
##        legend=paste0("FPR=0.05 -> ", format(cutoff, digits=2)),
##        lty=2, col="blue")
## plot(rocr_prec, spread.estimate="boxplot", avg="vertical",
##      main=paste0("BART"))
plot(rocr_mcc, spread.estimate="stderror", avg="vertical")
abline(v=max_mcc_cutoff, col="blue", lty=2)
plot(rocr_balance, spread.estimate="stderror", avg="threshold",
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

sink("/dev/null")
bart <- bartMachineCV(sign_feats[final_feats], labels,
                      use_missing_data=TRUE,
                      use_missing_data_dummies_as_covars=TRUE,
                      replace_missing_data_with_x_j_bar=FALSE,
                      impute_missingness_with_x_j_bar_for_lm=FALSE,
                      prob_rule_class=(1-max_mcc_cutoff),
                      verbose=FALSE,
                      serialize=FALSE)
sink()

phospho_feats <- read.delim("data/phosfun-features.tsv")

new_preds <- predict(bart, phospho_feats[final_feats],
                     type="class", prob_rule_class=(1-max_mcc_cutoff))
new_probs <- 1.0-bart_machine_get_posterior(bart, phospho_feats[final_feats])$y_hat

out_tbl <- data.frame(kinase=phospho_feats$kinase,
                      position=phospho_feats$position,
                      prob.act=(new_probs-max_mcc_cutoff))
out_tbl <- out_tbl[order(out_tbl$kinase, out_tbl$position),]

write.table(out_tbl, "out/reg-site-bart-sign-preds.tsv", quote=FALSE,
            row.names=FALSE, col.names=TRUE, sep="\t")
