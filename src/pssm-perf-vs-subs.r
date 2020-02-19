library(dplyr)

merged_pred_file  <- "out/kinase-merged-pred.tsv"
bart_file <- "out/kinase-merged-pred-no-mfs-rand-negs-bart-preds.tsv"
kin_sub_tbl_file <- "data/psiteplus-kinase-substrates.tsv"
val_set_file <- "data/validation-set-omnipath.tsv"

merged_pred <- read.delim(merged_pred_file, as.is=TRUE)
names(merged_pred)[1:2] <- c("prot1", "prot2")
bart <- read.delim(bart_file, as.is=TRUE)
kin_sub_tbl <- read.delim(kin_sub_tbl_file)
val_set <- read.table(val_set_file, as.is=TRUE)
names(val_set) <- c("prot1", "prot2")

merged_pred_summ <- merged_pred %>% group_by(prot1) %>%
    summarise(min.pssm.out=min(max.pssm.score, na.rm=TRUE),
              max.pssm.out=max(max.pssm.score, na.rm=TRUE),
              med.pssm.out=median(max.pssm.score, na.rm=TRUE),
              mean.pssm.out=mean(max.pssm.score, na.rm=TRUE),
              mean.func.score=mean(max.func.score, na.rm=TRUE),
              kin.type=unique(kinase.type))
merged_pred_summ2 <- merged_pred %>% group_by(prot2) %>%
    summarise(min.pssm.in=min(max.pssm.score, na.rm=TRUE),
              max.pssm.in=max(max.pssm.score, na.rm=TRUE),
              med.pssm.in=median(max.pssm.score, na.rm=TRUE),
              mean.pssm.in=mean(max.pssm.score, na.rm=TRUE),
              mean.func.score=mean(max.func.score, na.rm=TRUE),
              kin.type=unique(sub.kinase.type))

bart_summ <- bart %>% group_by(prot1) %>%
    summarise(min.bart.out=min(bart.pred, na.rm=TRUE),
              max.bart.out=max(bart.pred, na.rm=TRUE),
              med.bart.out=median(bart.pred, na.rm=TRUE),
              mean.bart.out=mean(bart.pred, na.rm=TRUE))
bart_summ2 <- bart %>% group_by(prot2) %>%
    summarise(min.bart.in=min(bart.pred, na.rm=TRUE),
              max.bart.in=max(bart.pred, na.rm=TRUE),
              med.bart.in=median(bart.pred, na.rm=TRUE),
              mean.bart.in=mean(bart.pred, na.rm=TRUE))


merged_pred_summ$num.subs <- log10(sapply(merged_pred_summ$prot1,
                             function(kin){
                                 nrow(subset(kin_sub_tbl, GENE==kin))
                             }))

merged_pred_summ$num.known.out <- sapply(merged_pred_summ$prot1,
                                  function(kin){
                                      nrow(subset(val_set, prot1==kin))
                                  })
merged_pred_summ2$num.kins <- log10(sapply(merged_pred_summ2$prot2,
                                    function(kin){
                                        nrow(subset(kin_sub_tbl, SUB_GENE==kin))
                                    }))
merged_pred_summ2$num.known.in <- sapply(merged_pred_summ2$prot2,
                                  function(kin){
                                      nrow(subset(val_set, prot2==kin))
                                  })

full_summ <- merge(merged_pred_summ, bart_summ)
full_summ2 <- merge(merged_pred_summ2, bart_summ2)

full_summ_st <- subset(full_summ, kin.type=="ST")
full_summ_y <- subset(full_summ, kin.type=="Y")

full_summ2_st <- subset(full_summ2, kin.type=="ST")
full_summ2_y <- subset(full_summ2, kin.type=="Y")


pdf("img/merged_pred-perf-vs-subs.pdf")
## number of substrates vs max MERGED_PRED score (regulator)
par(mfrow=c(2,2))
plot(full_summ_st$num.subs, full_summ_st$min.pssm.out,
     pch=19, xlab="log10(# substrates)", ylab="min(max. PSSM score)",
     main="ST Kinases", ylim=c(0, 1))
plot(full_summ_st$num.subs, full_summ_st$max.pssm.out,
     pch=19, xlab="log10(# substrates)", ylab="max(max. PSSM score)",
     main="ST Kinases", ylim=c(0, 1))
plot(full_summ_st$num.subs, full_summ_st$med.pssm.out,
     pch=19, xlab="log10(# substrates)", ylab="median(max. PSSM score)",
     main="ST Kinases", ylim=c(0, 1))
plot(full_summ_st$num.subs, full_summ_st$mean.pssm.out,
     pch=19, xlab="log10(# substrates)", ylab="mean(max. PSSM score)",
     main="ST Kinases", ylim=c(0, 1))
par(mfrow=c(2,2))
plot(full_summ_y$num.subs, full_summ_y$min.pssm.out,
     pch=19, xlab="log10(# substrates)", ylab="min(max. PSSM score)",
     main="Y Kinases", ylim=c(0, 1), xlim=c(0, 350))
plot(full_summ_y$num.subs, full_summ_y$max.pssm.out,
     pch=19, xlab="log10(# substrates)", ylab="max(max. PSSM score)",
     main="Y Kinases", ylim=c(0, 1), xlim=c(0, 350))
plot(full_summ_y$num.subs, full_summ_y$med.pssm.out,
     pch=19, xlab="log10(# substrates)", ylab="median(max. PSSM score)",
     main="Y Kinases", ylim=c(0, 1), xlim=c(0, 350))
plot(full_summ_y$num.subs, full_summ_y$mean.pssm.out,
     pch=19, xlab="log10(# substrates)", ylab="mean(max. PSSM score)",
     main="Y Kinases", ylim=c(0, 1), xlim=c(0, 350))


## number of known interactions vs max PSSM score (regulator)
par(mfrow=c(2,2))
plot(full_summ_st$num.known.out, full_summ_st$min.pssm.out,
     pch=19, xlab="Num. Known Interactions (out)", ylab="min(max. PSSM score)",
     main="ST Kinases", ylim=c(0, 1))
plot(full_summ_st$num.known.out, full_summ_st$max.pssm.out,
     pch=19, xlab="Num. Known Interactions (out)", ylab="max(max. PSSM score)",
     main="ST Kinases", ylim=c(0, 1))
plot(full_summ_st$num.known.out, full_summ_st$med.pssm.out,
     pch=19, xlab="Num. Known Interactions (out)", ylab="median(max. PSSM score)",
     main="ST Kinases", ylim=c(0, 1))
plot(full_summ_st$num.known.out, full_summ_st$mean.pssm.out,
     pch=19, xlab="Num. Known Interactions (out)", ylab="mean(max. PSSM score)",
     main="ST Kinases", ylim=c(0, 1))
par(mfrow=c(2,2))
plot(full_summ_y$num.known.out, full_summ_y$min.pssm.out,
     pch=19, xlab="Num. Known Interactions (out)", ylab="min(max. PSSM score)",
     main="Y Kinases", ylim=c(0, 1))
plot(full_summ_y$num.known.out, full_summ_y$max.pssm.out,
     pch=19, xlab="Num. Known Interactions (out)", ylab="max(max. PSSM score)",
     main="Y Kinases", ylim=c(0, 1))
plot(full_summ_y$num.known.out, full_summ_y$med.pssm.out,
     pch=19, xlab="Num. Known Interactions (out)", ylab="median(max. PSSM score)",
     main="Y Kinases", ylim=c(0, 1))
plot(full_summ_y$num.known.out, full_summ_y$mean.pssm.out,
     pch=19, xlab="Num. Known Interactions (out)", ylab="mean(max. PSSM score)",
     main="Y Kinases", ylim=c(0, 1))

## number of substrates vs BART prob (regulator)
par(mfrow=c(2, 2))
plot(full_summ_st$num.subs, full_summ_st$min.bart.out,
     pch=19, xlab="log10(# substrates)", ylab="min(BART prob.)",
     main="ST Kinases", ylim=c(0, 1))
plot(full_summ_st$num.subs, full_summ_st$max.bart.out,
     pch=19, xlab="log10(# substrates)", ylab="max(BART prob.)",
     main="ST Kinases", ylim=c(0, 1))
plot(full_summ_st$num.subs, full_summ_st$med.bart.out,
     pch=19, xlab="log10(# substrates)", ylab="median(BART prob.)",
     main="ST Kinases", ylim=c(0, 1))
plot(full_summ_st$num.subs, full_summ_st$mean.bart.out,
     pch=19, xlab="log10(# substrates)", ylab="mean(BART prob.)",
     main="ST Kinases", ylim=c(0, 1))
par(mfrow=c(2,2))
plot(full_summ_y$num.subs, full_summ_y$min.bart.out,
     pch=19, xlab="log10(# substrates)", ylab="min(BART prob.)",
     main="Y Kinases", ylim=c(0, 1), xlim=c(0, 350))
plot(full_summ_y$num.subs, full_summ_y$max.bart.out,
     pch=19, xlab="log10(# substrates)", ylab="max(BART prob.)",
     main="Y Kinases", ylim=c(0, 1), xlim=c(0, 350))
plot(full_summ_y$num.subs, full_summ_y$med.bart.out,
     pch=19, xlab="log10(# substrates)", ylab="median(BART prob.)",
     main="Y Kinases", ylim=c(0, 1), xlim=c(0, 350))
plot(full_summ_y$num.subs, full_summ_y$mean.bart.out,
     pch=19, xlab="log10(# substrates)", ylab="mean(BART prob.)",
     main="Y Kinases", ylim=c(0, 1), xlim=c(0, 350))

## number of known interactions vs BART prob (regulator)
par(mfrow=c(2,2))
plot(full_summ_st$num.known.out, full_summ_st$min.bart.out,
     pch=19, xlab="Num. Known Interactions (out)", ylab="min(BART prob.)",
     main="ST Kinases", ylim=c(0, 1))
plot(full_summ_st$num.known.out, full_summ_st$max.bart.out,
     pch=19, xlab="Num. Known Interactions (out)", ylab="max(BART prob.)",
     main="ST Kinases", ylim=c(0, 1))
plot(full_summ_st$num.known.out, full_summ_st$med.bart.out,
     pch=19, xlab="Num. Known Interactions (out)", ylab="median(BART prob.)",
     main="ST Kinases", ylim=c(0, 1))
plot(full_summ_st$num.known.out, full_summ_st$mean.bart.out,
     pch=19, xlab="Num. Known Interactions (out)", ylab="mean(BART prob.)",
     main="ST Kinases", ylim=c(0, 1))
par(mfrow=c(2,2))
plot(full_summ_y$num.known.out, full_summ_y$min.bart.out,
     pch=19, xlab="Num. Known Interactions (out)", ylab="min(BART prob.)",
     main="Y Kinases", ylim=c(0, 1))
plot(full_summ_y$num.known.out, full_summ_y$max.bart.out,
     pch=19, xlab="Num. Known Interactions (out)", ylab="max(BART prob.)",
     main="Y Kinases", ylim=c(0, 1))
plot(full_summ_y$num.known.out, full_summ_y$med.bart.out,
     pch=19, xlab="Num. Known Interactions (out)", ylab="median(BART prob.)",
     main="Y Kinases", ylim=c(0, 1))
plot(full_summ_y$num.known.out, full_summ_y$mean.bart.out,
     pch=19, xlab="Num. Known Interactions (out)", ylab="mean(BART prob.)",
     main="Y Kinases", ylim=c(0, 1))

## number of known interactions vs BART prob (substrate)
par(mfrow=c(2,2))
plot(full_summ2_st$num.known.in, full_summ2_st$min.bart.in,
     pch=19, xlab="Num. Known Interactions (in)", ylab="min(BART prob)",
     main="ST Kinases", ylim=c(0, 1))
plot(full_summ2_st$num.known.in, full_summ2_st$max.bart.in,
     pch=19, xlab="Num. Known Interactions (in)", ylab="max(BART prob)",
     main="ST Kinases", ylim=c(0, 1))
plot(full_summ2_st$num.known.in, full_summ2_st$med.bart.in,
     pch=19, xlab="Num. Known Interactions (in)", ylab="median(BART prob)",
     main="ST Kinases", ylim=c(0, 1))
plot(full_summ2_st$num.known.in, full_summ2_st$mean.bart.in,
     pch=19, xlab="Num. Known Interactions (in)", ylab="mean(BART prob)",
     main="ST Kinases", ylim=c(0, 1))
par(mfrow=c(2,2))
plot(full_summ2_y$num.known.in, full_summ2_y$min.bart.in,
     pch=19, xlab="Num. Known Interactions (in)", ylab="min(BART prob)",
     main="Y Kinases", ylim=c(0, 1))
plot(full_summ2_y$num.known.in, full_summ2_y$max.bart.in,
     pch=19, xlab="Num. Known Interactions (in)", ylab="max(BART prob)",
     main="Y Kinases", ylim=c(0, 1))
plot(full_summ2_y$num.known.in, full_summ2_y$med.bart.in,
     pch=19, xlab="Num. Known Interactions (in)", ylab="median(BART prob)",
     main="Y Kinases", ylim=c(0, 1))
plot(full_summ2_y$num.known.in, full_summ2_y$mean.bart.in,
     pch=19, xlab="Num. Known Interactions (in)", ylab="mean(BART prob)",
     main="Y Kinases", ylim=c(0, 1))

## max func. score vs BART prob (regulator)
par(mfrow=c(2,2))
plot(full_summ_st$mean.func.score, full_summ_st$min.bart.out,
     pch=19, xlab="max functional score", ylab="min(BART prob)",
     main="ST Kinases", ylim=c(0, 1))
plot(full_summ_st$mean.func.score, full_summ_st$max.bart.out,
     pch=19, xlab="max functional score", ylab="max(BART prob)",
     main="ST Kinases", ylim=c(0, 1))
plot(full_summ_st$mean.func.score, full_summ_st$med.bart.out,
     pch=19, xlab="max functional score", ylab="median(BART prob)",
     main="ST Kinases", ylim=c(0, 1))
plot(full_summ_st$mean.func.score, full_summ_st$mean.bart.out,
     pch=19, xlab="max functional score", ylab="mean(BART prob)",
     main="ST Kinases", ylim=c(0, 1))
par(mfrow=c(2,2))
plot(full_summ_y$mean.func.score, full_summ_y$min.bart.out,
     pch=19, xlab="max functional score", ylab="min(BART prob)",
     main="Y Kinases", ylim=c(0, 1))
plot(full_summ_y$mean.func.score, full_summ_y$max.bart.out,
     pch=19, xlab="max functional score", ylab="max(BART prob)",
     main="Y Kinases", ylim=c(0, 1))
plot(full_summ_y$mean.func.score, full_summ_y$med.bart.out,
     pch=19, xlab="max functional score", ylab="median(BART prob)",
     main="Y Kinases", ylim=c(0, 1))
plot(full_summ_y$mean.func.score, full_summ_y$mean.bart.out,
     pch=19, xlab="max functional score", ylab="mean(BART prob)",
     main="Y Kinases", ylim=c(0, 1))

## max func. score vs BART prob (substrate)
par(mfrow=c(2,2))
plot(full_summ2_st$mean.func.score, full_summ2_st$min.bart.in,
     pch=19, xlab="max functional score", ylab="min(BART prob)",
     main="ST Kinases", ylim=c(0, 1))
plot(full_summ2_st$mean.func.score, full_summ2_st$max.bart.in,
     pch=19, xlab="max functional score", ylab="max(BART prob)",
     main="ST Kinases", ylim=c(0, 1))
plot(full_summ2_st$mean.func.score, full_summ2_st$med.bart.in,
     pch=19, xlab="max functional score", ylab="median(BART prob)",
     main="ST Kinases", ylim=c(0, 1))
plot(full_summ2_st$mean.func.score, full_summ2_st$mean.bart.in,
     pch=19, xlab="max functional score", ylab="mean(BART prob)",
     main="ST Kinases", ylim=c(0, 1))
par(mfrow=c(2,2))
plot(full_summ2_y$mean.func.score, full_summ2_y$min.bart.in,
     pch=19, xlab="max functional score", ylab="min(BART prob)",
     main="Y Kinases", ylim=c(0, 1))
plot(full_summ2_y$mean.func.score, full_summ2_y$max.bart.in,
     pch=19, xlab="max functional score", ylab="max(BART prob)",
     main="Y Kinases", ylim=c(0, 1))
plot(full_summ2_y$mean.func.score, full_summ2_y$med.bart.in,
     pch=19, xlab="max functional score", ylab="median(BART prob)",
     main="Y Kinases", ylim=c(0, 1))
plot(full_summ2_y$mean.func.score, full_summ2_y$mean.bart.in,
     pch=19, xlab="max functional score", ylab="mean(BART prob)",
     main="Y Kinases", ylim=c(0, 1))

dev.off()

pdf("tmp/kinase-bart-prob-dists.pdf")
akt1 <- subset(bart, prot1=="AKT1")
plot(density(akt1$bart.pred))

mapk1 <- subset(bart, prot1=="MAPK1")
plot(density(mapk1$bart.pred))

s6k <- subset(bart, prot1=="RPS6KB1")
plot(density(s6k$bart.pred))

jak2 <- subset(bart, prot1=="JAK2")
plot(density(jak2$bart.pred))
dev.off()

