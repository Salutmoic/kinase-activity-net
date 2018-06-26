library(dplyr)
argv <- commandArgs(TRUE)
if (length(argv) != 3){
    stop("USAGE: <script> ASSOC_PREDS DIRECT_PREDS SIGN_PREDS")
}

assoc_preds_file <- argv[1]
direct_preds_file <- argv[2]
sign_preds_file <- argv[3]

assoc_preds <- read.delim(assoc_preds_file, as.is=TRUE)
direct_preds <- read.delim(direct_preds_file, as.is=TRUE)
sign_preds <- read.delim(sign_preds_file, as.is=TRUE)

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

final_tbl <- final_tbl[order(rownames(final_tbl)),]

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

write.table(final_tbl,
            "out/final-predictors/full-predictor.tsv",
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(summ_tbl,
            "out/final-predictors/full-predictor-summary.tsv",
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
