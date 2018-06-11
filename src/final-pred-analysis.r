library(ggplot2)
library(tidyr)
library(dplyr)
library(classInt)
library(scales)

argv <- commandArgs(TRUE)
if (length(argv) != 4){
    stop("USAGE: <script> PRED_FILE FEAT_FILE VAL_SET MAX_F_CUTOFF")
}

pred_file <- argv[1]
feat_file <- argv[2]
val_set_file <- argv[3]
max_f_cutoff <- as.numeric(argv[4])

file_base <- strsplit(basename(pred_file), split="\\.")[[1]][1]

preds <- read.delim(pred_file, as.is=TRUE)
feats <- read.delim(feat_file, as.is=TRUE)

val_set <- read.table(val_set_file, as.is=TRUE)
names(val_set) <- c("prot1", "prot2")

out_deg_val_df <- as.data.frame(table(val_set$prot1))
names(out_deg_val_df) <- c("prot", "out.degree.val")
in_deg_val_df <- as.data.frame(table(val_set$prot2))
names(in_deg_val_df) <- c("prot", "in.degree.val")
deg_val_df <- merge(out_deg_val_df, in_deg_val_df, by=c("prot"), all=TRUE)
deg_val_df[is.na(deg_val_df)] <- 0
deg_val_df <- deg_val_df %>%
    gather(direction, val.degree, out.degree.val:in.degree.val) %>%
    mutate(direction=gsub("\\.degree\\.val", "", direction)) %>%
    arrange(prot, direction)

cutoffs <- c(max_f_cutoff, seq(round(max_f_cutoff, digits=1), 0.0, -0.1))

pdf(paste0("img/", file_base, "-val-deg-vs-pred-deg.pdf"), height=3.5, width=7)
for (cutoff in cutoffs){
    pred_net <- subset(preds, bart.pred.mean > cutoff)
    out_deg_pred_df <- as.data.frame(table(pred_net$prot1))
    names(out_deg_pred_df) <- c("prot", "out.degree.pred")
    in_deg_pred_df <- as.data.frame(table(pred_net$prot2))
    names(in_deg_pred_df) <- c("prot", "in.degree.pred")
    deg_pred_df <- merge(out_deg_pred_df, in_deg_pred_df, by=c("prot"), all=TRUE)
    deg_pred_df[is.na(deg_pred_df)] <- 0
    deg_pred_df <- deg_pred_df %>%
        gather(direction, pred.degree, out.degree.pred:in.degree.pred) %>%
        mutate(direction=gsub("\\.degree\\.pred", "", direction)) %>%
        arrange(prot, direction)
    deg_df <- merge(deg_val_df, deg_pred_df, by=c("prot", "direction"), all=TRUE)
    deg_df[is.na(deg_df)] <- 0
    print(dim(deg_df))
    print(ggplot(deg_df, aes(val.degree, pred.degree)) +
          geom_count() +
          geom_abline(slope=1, intercept=0) +
          facet_grid(. ~ direction) +
          theme_bw() +
          labs(title=paste("Edge cutoff =", format(cutoff, digits=2)),
               x="degree in validation set",
               y="degree in predicted network"))
}
dev.off()

kin.citations <- read.table("data/kinome-entrez-data.tsv", as.is=TRUE)
names(kin.citations) <- c("kinase", "citations", "citations.filt")

num.cite.ints <- classIntervals(kin.citations$citations.filt, 10,
                                style="quantile")
num.cite.ints$brks <- round(num.cite.ints$brks)
kin.citations$cite.class <- cut(kin.citations$citations.filt,
                                num.cite.ints$brks,
                                include.lowest=TRUE, dig.lab=4)
rownames(kin.citations) <- kin.citations$kinase

preds$cite.class1 <- sapply(preds$prot1,
                            function(kin){
                                kin.citations[kin, "cite.class"]
                            })
preds$cite.class2 <- sapply(preds$prot2,
                            function(kin){
                                kin.citations[kin, "cite.class"]
                            })

pdf(paste0("img/", file_base, "-val-cites-vs-pred-cites.pdf"))
for (cutoff in cutoffs){
    pred_net <- subset(preds, bart.pred.mean > cutoff)
    print(ggplot(pred_net, aes(x=cite.class1, y=cite.class2)) +
          geom_bin2d() +
          scale_fill_distiller(direction=1, palette="Blues", na.value="white",
                               guide_colourbar(title="number of\ninteracting\npairs"),
                               limits=c(1,3500)) +
          scale_x_discrete(drop=FALSE) +
          scale_y_discrete(drop=FALSE) +
          theme_bw() +
          theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank()) +
          labs(title=paste("Edge cutoff =", format(cutoff, digits=2)),
               x="number of publications (kinase)",
               y="number of publications (substrate kinase)"))
}
dev.off()

relevant_feats <- setdiff(colnames(feats),
                          c("node1", "node2", "kinase.type", "sub.kinase.type",
                            "celline.sel2", "tissue.sel2"))

merged_dat <- merge(preds, feats, by.x=c("prot1", "prot2"),
                    by.y=c("node1", "node2"))

pred_vs_feats <- merged_dat %>%
    gather(feature, feature.value, relevant_feats) %>%
    select(prot1, prot2, feature, bart.pred.mean, feature.value) %>%
    filter(prot1 != prot2 & !is.na(feature.value) & !is.na(bart.pred.mean))
pred_vs_feats$cite.class1 <- sapply(pred_vs_feats$prot1,
                               function(kin){
                                   kin.citations[kin, "cite.class"]
                               })
pred_vs_feats$cite.class2 <- sapply(pred_vs_feats$prot2,
                               function(kin){
                                   kin.citations[kin, "cite.class"]
                               })
pred_vs_feats$cites1 <- sapply(pred_vs_feats$prot1,
                               function(kin){
                                   kin.citations[kin, "citations.filt"]
                               })
pred_vs_feats$cites2 <- sapply(pred_vs_feats$prot2,
                               function(kin){
                                   kin.citations[kin, "citations.filt"]
                               })

pdf(paste0("img/", file_base, "-citations-vs-pred.pdf"))
ggplot(pred_vs_feats, aes(cite.class1, bart.pred.mean)) +
    geom_violin() +
    theme_bw() +
    theme(axis.text.x=element_text(angle=45, hjust=1))
dev.off()

pred_vs_feats_sum <- pred_vs_feats %>%
    group_by(prot1, cites1) %>%
    summarise(mean.pred=mean(bart.pred.mean, na.rm=TRUE),
              sd.pred=sd(bart.pred.mean, na.rm=TRUE)) %>%
    arrange(desc(mean.pred))
pred_vs_feats_sum$prot1 <- factor(pred_vs_feats_sum$prot1,
                                  levels=unique(pred_vs_feats_sum$prot1))
pdf(paste0("img/", file_base, "-ranked-pred-scores.pdf"), width=84)
ggplot(pred_vs_feats_sum, aes(prot1, mean.pred, colour=log10(cites1))) +
    geom_pointrange(aes(ymin=(mean.pred-sd.pred), ymax=(mean.pred+sd.pred))) +
    scale_colour_gradient(low="#C6DBEF", high="#08306B") +
    theme_bw() +
    theme(axis.text.x=element_text(angle=45, hjust=1))
dev.off()


pred_vs_feats_sum <- pred_vs_feats %>%
    group_by(prot1, cite.class1) %>%
    summarise(mean.pred=mean(bart.pred.mean, na.rm=TRUE))
pdf(paste0("img/", file_base, "-citations-vs-mean-pred.pdf"))
ggplot(pred_vs_feats_sum, aes(cite.class1, mean.pred)) +
    geom_violin() +
    theme_bw() +
    theme(axis.text.x=element_text(angle=45, hjust=1))
dev.off()

pdf(paste0("img/", file_base, "-citations-vs-feats.pdf"))
ggplot(pred_vs_feats, aes(cite.class1, feature.value)) +
    geom_violin() +
    facet_wrap(~feature, ncol=3) +
    theme_bw() +
    theme(axis.text.x=element_text(angle=45, hjust=1))
dev.off()
