library(dplyr)
argv <- commandArgs(TRUE)
if (length(argv) != 4){
    stop("USAGE: <script> REG_PREDS SIGN_PREDS VAL_SET SIGN_VAL_SET")
}

reg_preds_file <- argv[1]
sign_preds_file <- argv[2]
reg_val_set_file <- argv[3]
sign_val_set_file <- argv[4]

reg_preds <- read.delim(reg_preds_file, as.is=TRUE)
sign_preds <- read.delim(sign_preds_file, as.is=TRUE)
reg_val_set <- read.table(reg_val_set_file, as.is=TRUE)
sign_val_set <- read.table(sign_val_set_file, as.is=TRUE)

final_tbl <- merge(reg_preds, sign_preds,
                   by.x=c("prot1", "prot2"), all.x=TRUE,
                   by.y=c("prot1", "prot2"), all.y=FALSE)

rownames(final_tbl) <- paste(final_tbl$prot1, final_tbl$prot2, sep="-")
rownames(reg_val_set) <- paste(reg_val_set[,1], reg_val_set[,2], sep="-")
rownames(sign_val_set) <- paste(sign_val_set[,1], sign_val_set[,2], sep="-")

final_tbl$reg.is.true <- rownames(final_tbl) %in% rownames(reg_val_set)

final_tbl$pred_sign <- sapply(final_tbl$prob.act,
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
                                     if (is.na(final_tbl[pair, "prob.act"]) |
                                         is.na(final_tbl[pair, "prob.act"]))
                                         return (NA)
                                     ## if (final_tbl[pair, "prob.act"] < 0.7 &
                                     ##     final_tbl[pair, "prob.act"] > 0.3)
                                     ##     return (NA)
                                     return ((final_tbl[pair, "prob.act"] >= 0.5 &&
                                              sign_val_set[pair, 3] == "activates") ||
                                             (final_tbl[pair, "prob.act"] < 0.5 &&
                                              sign_val_set[pair, 3] == "inhibits"))
                                 })

final_tbl <- final_tbl[order(rownames(final_tbl)),]

erk_akt1_kins <- c("BRAF", "MAP2K1", "MAPK1", "MTOR", "PDPK1", "AKT1", "GSK3B",
                   "RPS6KB1")
erk_akt1_tbl <- subset(final_tbl, prot1 %in% erk_akt1_kins &
                                  prot2 %in% erk_akt1_kins &
                                  bart.pred > 0.8)

summ_tbl_a <- final_tbl %>% group_by(prot1) %>%
    summarise(min.score.out=min(bart.pred, na.rm=TRUE),
              mean.score.out=mean(bart.pred, na.rm=TRUE),
              med.score.out=median(bart.pred, na.rm=TRUE),
              max.score.out=max(bart.pred, na.rm=TRUE))
summ_tbl_b <- final_tbl %>% group_by(prot2) %>%
    summarise(min.score.in=min(bart.pred, na.rm=TRUE),
              mean.score.in=mean(bart.pred, na.rm=TRUE),
              med.score.in=median(bart.pred, na.rm=TRUE),
              max.score.in=max(bart.pred, na.rm=TRUE))
summ_tbl <- merge(summ_tbl_a, summ_tbl_b, by.x="prot1", by.y="prot2")

if (grepl("no-dcg", reg_preds_file)){
    dcg_fn_mod <- "-no-dcg"
}else{
    dcg_fn_mod <- ""
}

write.table(final_tbl,
            paste0("out/full-annotated-predictor", dcg_fn_mod, ".tsv"),
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(subset(final_tbl, bart.pred>0.9),
            paste0("out/full-annotated-predictor", dcg_fn_mod, "-hiconf.tsv"),
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(erk_akt1_tbl,
            paste0("out/full-annotated-predictor", dcg_fn_mod, "-erk-akt1.tsv"),
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(summ_tbl,
            paste0("out/full-annotated-predictor", dcg_fn_mod, "-summary.tsv"),
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
