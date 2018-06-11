library(dplyr)
argv <- commandArgs(TRUE)
if (length(argv) != 6){
    stop("USAGE: <script> ASSOC_PREDS DIRECT_PREDS SIGN_PREDS ASSOC_VAL_SET DIRECT_VAL_SET SIGN_VAL_SET")
}

assoc_preds_file <- argv[1]
direct_preds_file <- argv[2]
sign_preds_file <- argv[3]
assoc_val_set_file <- argv[4]
direct_val_set_file <- argv[5]
sign_val_set_file <- argv[6]

assoc_preds <- read.delim(assoc_preds_file, as.is=TRUE)
direct_preds <- read.delim(direct_preds_file, as.is=TRUE)
sign_preds <- read.delim(sign_preds_file, as.is=TRUE)
assoc_val_set <- read.table(assoc_val_set_file, as.is=TRUE)
direct_val_set <- read.table(direct_val_set_file, as.is=TRUE)
sign_val_set <- read.table(sign_val_set_file, as.is=TRUE)

assoc_preds <- assoc_preds[,c("prot1", "prot2", "bart.pred.mean", "bart.pred.sd")]
names(assoc_preds) <- c("prot1", "prot2", "assoc.pred.mean", "assoc.pred.sd")
direct_preds <- direct_preds[,c("prot1", "prot2", "bart.pred.mean", "bart.pred.sd")]
names(direct_preds) <- c("prot1", "prot2", "direct.pred.mean", "direct.pred.sd")

final_tbl <- merge(assoc_preds, direct_preds,
                   by.x=c("prot1", "prot2"), all.x=TRUE,
                   by.y=c("prot1", "prot2"), all.y=FALSE)
final_tbl <- merge(final_tbl, sign_preds,
                   by.x=c("prot1", "prot2"), all.x=TRUE,
                   by.y=c("prot1", "prot2"), all.y=FALSE)

rownames(final_tbl) <- paste(final_tbl$prot1, final_tbl$prot2, sep="-")
rownames(assoc_val_set) <- paste(assoc_val_set[,1], assoc_val_set[,2], sep="-")
rownames(direct_val_set) <- paste(direct_val_set[,1], direct_val_set[,2], sep="-")
rownames(sign_val_set) <- paste(sign_val_set[,1], sign_val_set[,2], sep="-")

final_tbl$assoc.is.true <- rownames(final_tbl) %in% rownames(assoc_val_set)
final_tbl$direct.is.true <- rownames(final_tbl) %in% rownames(direct_val_set)
final_tbl$direct.is.true[is.na(final_tbl$direct.pred.mean)] <- NA

final_tbl$pred_sign <- sapply(final_tbl$prob.act.mean,
                              function(score){
                                  if (is.na(score))
                                      return (NA)
                                  if (score >= 0.5)
                                      return ("up")
                                  if (score < 0.5)
                                      return ("down")
                                  return ("ambig")
                              })

final_tbl$sign.is.correct <- sapply(rownames(final_tbl),
                                 function(pair){
                                     if (!(pair %in% rownames(sign_val_set)))
                                         return (NA)
                                     if (is.na(final_tbl[pair, "prob.act.mean"]) |
                                         is.na(final_tbl[pair, "prob.act.mean"]))
                                         return (NA)
                                     ## if (final_tbl[pair, "prob.act.mean"] < 0.7 &
                                     ##     final_tbl[pair, "prob.act.mean"] > 0.3)
                                     ##     return (NA)
                                     return ((final_tbl[pair, "prob.act.mean"] >= 0.5 &&
                                              sign_val_set[pair, 3] == "activates") ||
                                             (final_tbl[pair, "prob.act.mean"] < 0.5 &&
                                              sign_val_set[pair, 3] == "inhibits"))
                                 })

final_tbl <- final_tbl[order(rownames(final_tbl)),]

erk_akt1_kins <- c("BRAF", "MAP2K1", "MAPK1", "MTOR", "PDPK1", "AKT1", "GSK3B",
                   "RPS6KB1")
erk_akt1_tbl <- subset(final_tbl, prot1 %in% erk_akt1_kins &
                                  prot2 %in% erk_akt1_kins &
                                  assoc.pred.mean > 0.5) &
                                  (is.na(direct.pred.mean) | direct.pred.mean > 0.7))

summ_tbl_a <- final_tbl %>% group_by(prot1) %>%
    summarise(min.score.out=min(assoc.pred.mean, na.rm=TRUE),
              mean.score.out=mean(assoc.pred.mean, na.rm=TRUE),
              med.score.out=median(assoc.pred.mean, na.rm=TRUE),
              max.score.out=max(assoc.pred.mean, na.rm=TRUE))
summ_tbl_b <- final_tbl %>% group_by(prot2) %>%
    summarise(min.score.in=min(assoc.pred.mean, na.rm=TRUE),
              mean.score.in=mean(assoc.pred.mean, na.rm=TRUE),
              med.score.in=median(assoc.pred.mean, na.rm=TRUE),
              max.score.in=max(assoc.pred.mean, na.rm=TRUE))
summ_tbl <- merge(summ_tbl_a, summ_tbl_b, by.x="prot1", by.y="prot2")

## file_base <- strsplit(basename(assoc_preds_file), split="\\.")[[1]][1]
## pred_ver <- sub("kinase-merged-pred-", "", sub("-bart-preds", "", file_base))
## if (pred_ver != "kinase-merged-pred"){
##     fn_mod <- paste0("-", pred_ver)
## }else{
##     fn_mod <- ""
## }

write.table(final_tbl,
            ## paste0("out/full-annotated-predictor", fn_mod, ".tsv"),
            "out/full-annotated-predictor.tsv",
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(subset(final_tbl, assoc.pred.mean>0.92),
            ## paste0("out/full-annotated-predictor", fn_mod, "-hiconf.tsv"),
            "out/full-annotated-predictor-hiconf.tsv",
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(erk_akt1_tbl,
            ## paste0("out/full-annotated-predictor", fn_mod, "-erk-akt1.tsv"),
            "out/full-annotated-predictor-erk-akt1.tsv",
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(summ_tbl,
            ## paste0("out/full-annotated-predictor", fn_mod, "-summary.tsv"),
            "out/full-annotated-predictor-summary.tsv",
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
