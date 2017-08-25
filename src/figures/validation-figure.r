library(dplyr)
library(ggplot2)
library(readr)
library(tibble)
library(tidyr)
library(scales)
library(cowplot)
library(viridis)
library(classInt)
library(grid)
library(RColorBrewer)
source("src/rocr-ggplot2.r")

theme_set(theme_cowplot(font_size=10))

trim_tails <- function(range = c(-Inf, Inf), n){
    trans_new("trim_tails",
              transform = function(x) {
                  force(range)
                  desired_breaks <- extended_breaks(n = n)(x[x >= range[1] & x <= range[2]])
                  break_increment <- diff(desired_breaks)[1]
                  x[x < range[1]] <- range[1] ## - break_increment
                  x[x > range[2]] <- range[2] ## + break_increment
                  return (x)
              },
              inverse = function(x) return (x),

              breaks = function(x) {
                  force(range)
                  return (extended_breaks(n = n)(x))
              },
              format = function(x) {
                  force(range)
                  ## x[1] <- paste("<", range[1])
                  x_tmp <- sapply(x, paste0, "%")
                  x_tmp[which.max(x)] <- paste0(range[2], "%+")
                  x_tmp[length(x_tmp)] <- ""
                  return (x_tmp)
              })
}

file_base1 <- "img/rds/kinase-merged-feats-clean-val-rand-negs"
file_base2 <- "img/rds/kinase-merged-feats-clean-hinegs-bart-preds-final-val-rand-negs"
file_base3 <- "img/rds/kinase-merged-feats-clean-hinegs-bart-rand-negs-cv"

feat_roc_data <- readRDS(paste0(file_base1, "-roc-data.rds"))
roc_data <- readRDS(paste0(file_base3, "-roc-data.rds"))

full_roc_data <- rbind(feat_roc_data, roc_data)

full_roc_data$feature <- factor(full_roc_data$feature,
                                levels=levels(full_roc_data$feature),
                                labels=c("coregulation (CPTAC)",
                                         "coregulation (Wilkes et al.)",
                                         "coexpression (Prot. Atlas)",
                                         "coexpression (GTEx)",
                                         "max. PSSM score",
                                         "DCG score",
                                         "max. functional score",
                                         "regulator tissue spec.",
                                         "substrate tissue spec.",
                                         "final predictor"))

roc_plot <- rocr2d.ggplot2(full_roc_data, "false positive rate",
                           "true positive rate", 1, c(0, 1)) +
    coord_fixed() +
    ## scale_colour_brewer(palette="Set3") +
    scale_colour_manual(values=c(brewer.pal(name="Set3", n=9), "#666666")) +
    scale_x_continuous(expand=c(0.001, 0, 0.001, 0), limits=c(0, 1),
                       breaks=c(0, 0.5, 1.0)) +
    scale_y_continuous(expand=c(0.001, 0, 0.001, 0), limits=c(0, 1),
                       breaks=c(0, 0.5, 1.0)) +
    ## guides(colour=guide_legend(keywidth=0.2, keyheight=0.15, default.unit="inch")) +
    guides(colour=FALSE) +
    theme(legend.title=element_blank(), legend.text=element_text(size=7),
           plot.margin=margin(0, 0, 0, 0))

save_plot("img/figures/panels/validation-roc.pdf", roc_plot)

feat_prec_data <- readRDS(paste0(file_base1, "-prec-data.rds"))
prec_data <- readRDS(paste0(file_base3, "-prec-data.rds"))

full_prec_data <- rbind(feat_prec_data, prec_data)

full_prec_data$feature <- factor(full_prec_data$feature,
                                levels=levels(full_prec_data$feature),
                                labels=c("coregulation (CPTAC)",
                                         "coregulation (Wilkes et al.)",
                                         "coexpression (Prot. Atlas)",
                                         "coexpression (GTEx)",
                                         "max. PSSM score",
                                         "DCG score",
                                         "max. functional score",
                                         "regulator tissue spec.",
                                         "substrate tissue spec.",
                                         "final predictor"))

prec_plot <- rocr2d.ggplot2(full_prec_data, "recall", "precision", 1, c(1/9, 0)) +
    coord_fixed() +
    scale_colour_manual(values=c(brewer.pal(name="Set3", n=9), "#666666")) +
    scale_x_continuous(expand=c(0.001, 0, 0.001, 0), limits=c(0, 1),
                       breaks=c(0, 0.5, 1.0)) +
    scale_y_continuous(expand=c(0.001, 0, 0.001, 0), limits=c(0, 1),
                       breaks=c(0, 0.5, 1.0)) +
    guides(colour=guide_legend(keywidth=0.2, keyheight=0.15, default.unit="inch")) +
    theme(legend.title=element_blank(), legend.text=element_text(size=7),
           plot.margin=margin(0, 0, 0, 0))

save_plot("img/figures/panels/validation-prec.pdf", prec_plot)

feat_auc_plot <- readRDS(paste0(file_base1, "-auc.rds"))

feat_auc_data <- readRDS(paste0(file_base1, "-auc-data.rds"))
pred_auc_data <- readRDS(paste0(file_base3, "-auc-data.rds"))

table(feat_auc_data$feature)
table(pred_auc_data$feature)

auc_data <- rbind(feat_auc_data, pred_auc_data)
auc_data$feature <- factor(auc_data$feature,
                           levels=levels(auc_data$feature),
                           labels=c("coregulation (CPTAC)",
                                    "coregulation (Wilkes et al.)",
                                    "coexpression (Prot. Atlas)",
                                    "coexpression (GTEx)",
                                    "max. PSSM score",
                                    "DCG score",
                                    "max. functional score",
                                    "regulator tissue spec.",
                                    "substrate tissue spec.",
                                    "final predictor"))
auc_plot <- rocr.sum.ggplot2(auc_data, "AUC", 100, c(0.5, 0))  +
    scale_colour_manual(values=c(brewer.pal(name="Set3", n=9), "#666666")) +
    ylim(0.5, 1.0) +
    guides(colour=FALSE, alpha=FALSE, fill=FALSE)

save_plot("img/figures/panels/validation-auc.pdf", auc_plot)

net_dens_plot <- readRDS("img/rds/merged-predictor-net-density-by-cites.rds") +
    coord_fixed() +
    scale_colour_brewer("min. #\ncitations", palette="Set3") +
    guides(colour=guide_legend(keywidth=0.2, keyheight=0.2, default.unit="inch"),
           legend.title=element_blank(),
           legend.text=element_text(size=7))

save_plot("img/figures/panels/validation-net-density.pdf", net_dens_plot)

site_sign_data <- read_tsv("img/rds/site-sign-bart-pred-mcc-data.tsv")
reg_sign_data <- read_tsv("img/rds/sign-bart-pred-mcc-data.tsv")

sign_data <- rbind(site_sign_data, reg_sign_data)

sign_data$feature <- factor(sign_data$feature,
                            levels=c("site sign prob.",
                                     "regulation sign prob."),
                            labels=c("phosphosite sign",
                                     "kinase-kinase\nregulation sign"))

## 3*20 --> 3-fold cross validation, 20 reps
sign_plot <- rocr2d.ggplot2(sign_data, "\"activation\" probability cutoff",
                            "Matthews correlation coefficient", 3*20, NULL) +
    facet_wrap(~feature, nrow=1) +
    scale_colour_manual(values=c("#000000", "#000000")) +
    guides(colour=FALSE) +
    ylim(0, 0.5) +
    scale_x_continuous(limits=c(0, 1), breaks=c(0, 0.5, 1.0)) +
    theme(strip.background=element_blank(), strip.placement="inside",
          panel.spacing=unit(0.5, "cm"),
          strip.text=element_text(vjust=0))

save_plot("img/figures/panels/validation-sign.pdf", sign_plot)

lo_kin_dir <- "out/kinase-merged-feats-clean-lo-kin/"
lo_out_data <- read_tsv(paste0(lo_kin_dir,
                               "kinase-merged-feats-clean-lo-kin-out.tsv"),
                        col_names=FALSE)
names(lo_out_data) <- c("node1", "node2", "lo.kin", "bart.pred.mean")
lo_out_data <- lo_out_data %>% mutate(id=paste(node1, node2, sep="-"))
lo_in_data <- read_tsv(paste0(lo_kin_dir,
                              "kinase-merged-feats-clean-lo-kin-in.tsv"),
                       col_names=FALSE)
names(lo_in_data) <- c("node1", "node2", "lo.kin", "bart.pred.mean")
lo_in_data <- lo_in_data %>% mutate(id=paste(node1, node2, sep="-"))

true_pos <- read_tsv("data/direct-validation-set-omnipath.tsv", col_names=FALSE)
names(true_pos) <- c("node1", "node2")
true_pos <- true_pos %>% mutate(id=paste(node1, node2, sep="-"))
citations <- read.delim("out/kinase-citations.tsv", as.is=TRUE)
citations <- subset(citations, count.type=="num.subs")

lo_out_ranks <- lo_out_data %>%
    group_by(lo.kin) %>%
    arrange(desc(bart.pred.mean, .by_group=TRUE)) %>%
    mutate(rank=rank(-bart.pred.mean, ties.method="first")) %>%
    filter(id %in% true_pos$id)  %>%
    group_by(node1) %>%
    summarise(best.rank=min(rank))
lo_in_ranks <- lo_in_data %>%
    group_by(lo.kin) %>%
    arrange(desc(bart.pred.mean, .by_group=TRUE)) %>%
    mutate(rank=rank(-bart.pred.mean, ties.method="first")) %>%
    filter(id %in% true_pos$id)  %>%
    group_by(node2) %>%
    summarise(best.rank=min(rank))
lo_out_ranks$direction <- rep("regulator", nrow(lo_out_ranks))
lo_in_ranks$direction <- rep("substrate", nrow(lo_in_ranks))
names(lo_out_ranks) <- c("kinase", "best.rank", "direction")
names(lo_in_ranks) <- c("kinase", "best.rank", "direction")
lo_ranks <- rbind(lo_out_ranks, lo_in_ranks)
lo_ranks$direction <- factor(lo_ranks$direction)
lo_ranks <- lo_ranks %>%
    inner_join(citations, by="kinase")
lo_ranks$is.random <- rep(FALSE, nrow(lo_ranks))
lo_ranks$class <- rep("predictions", nrow(lo_ranks))

lo_out_ranks_rand <- lo_out_data %>%
    group_by(lo.kin) %>%
    mutate(bart.pred.mean.rand = sample(bart.pred.mean)) %>%
    arrange(desc(bart.pred.mean.rand, .by_group=TRUE)) %>%
    mutate(rank=rank(-bart.pred.mean.rand, ties.method="first")) %>%
    filter(id %in% true_pos$id) %>%
    group_by(node1) %>%
    summarise(best.rank=min(rank))
lo_in_ranks_rand <- lo_in_data %>%
    group_by(lo.kin) %>%
    mutate(bart.pred.mean.rand = sample(bart.pred.mean)) %>%
    arrange(desc(bart.pred.mean.rand, .by_group=TRUE)) %>%
    mutate(rank=rank(-bart.pred.mean.rand, ties.method="first")) %>%
    filter(id %in% true_pos$id) %>%
    group_by(node2) %>%
    summarise(best.rank=min(rank))
lo_out_ranks_rand$direction <- rep("regulator", nrow(lo_out_ranks_rand))
lo_in_ranks_rand$direction <- rep("substrate", nrow(lo_in_ranks_rand))
names(lo_out_ranks_rand) <- c("kinase", "best.rank", "direction")
names(lo_in_ranks_rand) <- c("kinase", "best.rank", "direction")
lo_ranks_rand <- rbind(lo_out_ranks_rand, lo_in_ranks_rand)
lo_ranks_rand$direction <- factor(lo_ranks_rand$direction)
lo_ranks_rand <- lo_ranks_rand %>%
    inner_join(citations, by="kinase")
lo_ranks_rand$is.random <- rep(TRUE, nrow(lo_ranks_rand))
lo_ranks_rand$class <- rep("random", nrow(lo_ranks_rand))

lo_ranks <- rbind(lo_ranks, lo_ranks_rand)
lo_ranks$class <- as.factor(lo_ranks$class)

rank_plot <- lo_ranks %>%
    arrange(num.cites.filt) %>%
    ggplot(aes(x=direction, y=best.rank, fill=class)) +
    geom_violin(draw_quantiles=c(0.25, 0.5, 0.75), adjust=0.25) +
    ## geom_boxplot() +
    ## geom_jitter() +
    ## facet_wrap(~direction) +
    scale_y_log10("top rank of known\n relationship",
                  expand=c(0.05, 0, 0.05, 0)) +
    ## annotation_logticks(sides="b") +
    scale_x_discrete("role of kinase") +
    scale_fill_manual(values=c("#8f2150", "#b3b3b3")) +
    ## coord_flip() +
    guides(fill=guide_legend(keywidth=0.1, keyheight=0.1, default.unit="inch")) +
    theme(legend.title=element_blank(), legend.text=element_text(size=7),
          legend.position="top", legend.direction="horizontal",
          legend.margin=margin(0, 0, 1, 0),
          axis.text.y=element_text(angle=90, hjust=0.5))

save_plot("img/figures/panels/validation-rank-violin.pdf", rank_plot)

lo_ranks %>% filter(!is.random)

rank_plot2 <- lo_ranks %>%
    filter(!is.random) %>%
    ggplot(aes(x=num.cites.filt, y=best.rank)) +
    geom_point() +
    geom_smooth(method=lm) +
    scale_x_log10("# citations") +
    scale_y_log10("top rank of known\n relationship") +
    facet_wrap("direction") +
    theme(strip.background=element_blank(), strip.placement="outside",
          panel.spacing=unit(0.5, "cm"),
          strip.text=element_text(vjust=1))

save_plot("img/figures/panels/validation-rank-vs-citations.pdf", rank_plot2)

cutoff <- 0.5
final_preds <- read_tsv("out/final-predictors/merged-predictor.tsv") %>%
    filter(bart.pred.mean >= cutoff) %>%
    select(prot1, prot2)

final_preds$source <- rep("prediction", nrow(final_preds))
true_net <- true_pos
colnames(true_net)[1:2] <- c("prot1", "prot2")
true_net$source <- rep("literature", nrow(true_net))
true_net <- select(true_net, prot1, prot2, source)

net_dens_data <- rbind(final_preds, true_net)

num.cite.ints <- classIntervals(citations$num.cites.filt, 10, style="quantile")
num.cite.ints$brks <- round(num.cite.ints$brks)
citations$cite.class <- cut(citations$num.cites.filt, num.cite.ints$brks,
                            include.lowest=TRUE, dig.lab=4)

net_dens_data$cite.class1 <- sapply(net_dens_data$prot1,
                                   function(prot){
                                       kin_cites <- subset(citations, kinase==prot)
                                       return(kin_cites$cite.class)
                                   })
net_dens_data$cite.class2 <- sapply(net_dens_data$prot2,
                                   function(prot){
                                       kin_cites <- subset(citations, kinase==prot)
                                       return(kin_cites$cite.class)
                                   })

## print(as.data.frame(subset(net_dens_data, source=="prediction" &
##                                           cite.class1 %in% c("(39,57]", "(57,95]"))))

net_dens_data$int.id <- paste(net_dens_data$prot1, net_dens_data$prot2)
net_dens_data$cite.class.id <- paste(net_dens_data$cite.class1,
                                     net_dens_data$cite.class2)
final_preds_w_class <- filter(net_dens_data, source=="prediction", prot1 != prot2)
true_pos_w_class <- filter(net_dens_data, source=="literature", prot1 != prot2)

final_preds_contrib <- final_preds_w_class %>%
    mutate(is.known = final_preds_w_class$int.id %in% true_pos_w_class$int.id) %>%
    group_by(cite.class.id) %>%
    summarise(perc.known = length(which(is.known))/n()) %>%
    separate(cite.class.id, sep=" ", into=c("cite.class1", "cite.class2"))

final_preds_counts <- as_tibble(table(final_preds_w_class$cite.class1,
                                      final_preds_w_class$cite.class2))
names(final_preds_counts) <- c("cite.class1", "cite.class2", "n")

final_preds_data <- final_preds_contrib %>%
    inner_join(final_preds_counts)

blank_spaces <- sapply(1:10, function(n) paste0(rep(" ", n), collapse=""))
final_preds_data$cite.class1 <- factor(final_preds_data$cite.class1,
                                       levels=levels(citations$cite.class),
                                       labels=blank_spaces)
final_preds_data$cite.class2 <- factor(final_preds_data$cite.class2,
                                       levels=levels(citations$cite.class),
                                       labels=blank_spaces)

cite_plot <- ggplot(final_preds_data,
                    aes(x=cite.class1, y=cite.class2, colour=perc.known*100,
                        size=n)) +
    geom_point() +
    scale_colour_viridis("known\nrelationships",
                         trans=trim_tails(range=c(0, 20), n=5),
                         breaks=c(0, 5, 10, 15, 20, 25),
                         option="plasma") +
    scale_size("# predicted\nrelationships", range=c(0, 5)) +
    scale_x_discrete(drop=FALSE) +
    scale_y_discrete(drop=FALSE) +
    coord_fixed() +
    theme(axis.title=element_text(size=8),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          legend.box = "vertical", legend.title=element_text(size=8),
          legend.text=element_text(size=7), axis.text=element_text(size=7),
          legend.key.height=unit(0.35, "cm"),
          legend.key.width=unit(0.3, "cm"),
          legend.spacing=unit(0.2, "cm"),
          axis.ticks=element_blank()) +
    xlab("# publications\n(\"upstream\" protein kinase)") +
    ylab("# publications\n(\"downstream\" substrate kinase)") +
    annotation_custom(grob=polygonGrob(x=c(0.5, 11.5, 11.5, 0.5),
                                       y=c(2, 2, 0, 2),
                                       gp=gpar(fill="grey75", col="grey75")),
                      xmin=0.1, xmax=1, ymin=-1, ymax=-0.5) +
    annotation_custom(grob=polygonGrob(y=c(0.5, 11.5, 11.5, 0.5),
                                       x=c(2, 2, 0, 2),
                                       gp=gpar(fill="grey75", col="grey75")),
                      ymin=0.1, ymax=1, xmin=-1, xmax=-0.5)

cite_plot_gt <- ggplot_gtable(ggplot_build(cite_plot))
cite_plot_gt$layout$clip[cite_plot_gt$layout$name=="panel"] <- "off"

save_plot("img/figures/panels/validation-citations.pdf", cite_plot_gt)

true_pos_1src <- read_tsv("data/direct-validation-set-omnipath-1-src.tsv",
                        col_names=FALSE)
colnames(true_pos_1src) <- c("prot1", "prot2")
true_pos <- true_pos %>% select(node1, node2, id)
colnames(true_pos) <- c("prot1", "prot2", "id")
true_pos_1src <- true_pos_1src %>% mutate(id=paste(prot1, prot2, sep="-"))

true_pos_1src <- subset(true_pos_1src, !(id %in% true_pos$id))

pred_net <- read_tsv("out/final-predictors/merged-predictor.tsv") %>%
    select(prot1, prot2, bart.pred.mean) %>%
    mutate(id = paste(prot1, prot2, sep="-"))

pred_net$val_class <- rep(NA, nrow(pred_net))
pred_net$val_class[pred_net$id %in% true_pos$id] <- "2+"
pred_net$val_class[pred_net$id %in% true_pos_1src$id] <- "1"
pred_net$val_class[!(pred_net$id %in% true_pos$id) &
                   !(pred_net$id %in% true_pos_1src$id)] <- "0"
pred_net$val_class <- factor(pred_net$val_class)

pred_net <- pred_net %>%
    group_by(val_class) %>%
    mutate(outlier = bart.pred.mean > (quantile(bart.pred.mean, 0.75, na.rm=TRUE) + IQR(bart.pred.mean, na.rm=TRUE)*1.5)) %>%
    ungroup()

val_class_plot <- pred_net %>%
    ggplot(aes(x=val_class, y=bart.pred.mean, fill=val_class)) +
    ## geom_violin(draw_quantiles=c(0.25, 0.5, 0.75), scale="width") +
    geom_boxplot(outlier.shape=NA) +
    geom_jitter(data = function(x) dplyr::filter_(x, ~ outlier), alpha=0.1,
                size=0.25, width=0.15) +
    scale_fill_manual(values=c("#b3b3b3", "#21508f", "#8f2150")) +
    guides(fill=FALSE) +
    ## coord_flip() +
    scale_x_discrete("# sources in Omnipath reporting\nregulatory relationship") +
    scale_y_continuous("predicted posterior probability\nof regulatory relationship") +
    theme(legend.title=element_blank(), legend.text=element_text(size=7),
          legend.position="top", legend.direction="horizontal",
          axis.text.y=element_text(angle=90, hjust=0.5))

save_plot("img/figures/panels/validation-prob-by-sources.pdf", val_class_plot)

## clust_enrich_plot <- readRDS("img/rds/cluster-enrichment-ggplot.rds")
pathways = read_tsv("out/final-predictors/cluster-analysis/pathway-enrichment-plot-data.tsv") %>%
    filter(file == "scaled-Louvain-0.5.tsv") %>%
    arrange(cluster, p.values)

tbl_sub = pathways %>%
    group_by(cluster) %>%
    top_n(n = -5, wt = p.values) %>%
    group_by(cluster) %>%
    top_n(n = -5, wt = pathways) %>%
    mutate(perc_annotated = 100*no_genes_in_top_enriched/no_genes_in_cluster)

tbl_sub$pathways <- reorder(factor(tbl_sub$pathways), tbl_sub$cluster)
tbl_sub$cluster <- factor(tbl_sub$cluster,
                          levels = sort(unique(tbl_sub$cluster)),
                          labels = paste("cluster", sort(unique(tbl_sub$cluster))))

cluster_ns_df <- unique(tbl_sub[,c("cluster", "no_genes_in_cluster")])
cluster_ns <- list()
for (i in 1:nrow(cluster_ns_df)){
    index <- as.character(cluster_ns_df$cluster[i])
    cluster_ns[[index]] <- cluster_ns_df$no_genes_in_cluster[i]
}

cluster_labeller <- function(cluster_names){
    return (sapply(cluster_names,
                   function(n){
                       paste0(n, "\n(n=", cluster_ns[[n]], ")")
                   }))
}
cluster_labeller <- as_labeller(cluster_labeller)

clust_enrich_plot <- tbl_sub %>%
    group_by(cluster) %>%
    mutate(num_path = length(unique(pathways))) %>%
    ungroup() %>%
    mutate(width = num_path/max(num_path)) %>%
    ggplot(aes(x = pathways, y = perc_annotated, fill=-log10(p.values))) +
    geom_bar(aes(width= 0.75*width), stat="identity") +
    scale_y_continuous("% of cluster annotated",
                       expand=c(0.01, 0.0, 0.01, 0.0)) +
    scale_fill_viridis(expression(atop("enrichment",
                                       paste(-log[10], "(", italic(p), ")")))) +
    ## xlab("annotated pathway") +
    coord_flip() +
    facet_wrap(~cluster, nrow=4, strip.position="right", scales="free_y",
               labeller=cluster_labeller) +
    theme(axis.title.y=element_blank(),
          legend.title=element_text(size=8),
          legend.text=element_text(size=7), axis.text=element_text(size=7),
          legend.key.height=unit(0.35, "cm"),
          legend.key.width=unit(0.3, "cm"),
          strip.background=element_blank(), strip.placement="outside",
          strip.text=element_text(size=8),
          panel.spacing=unit(0.1, "cm"))

save_plot("img/figures/panels/validation-clust-enrich.pdf", clust_enrich_plot)

top_grid <- plot_grid(roc_plot, prec_plot, rank_plot,
                      nrow=1, scale=0.9, align="h", axis="b",
                      rel_widths=c(0.54, 1.0, 0.5), labels=c("", "", "b"))
mid_grid <- plot_grid(val_class_plot, clust_enrich_plot,
                      nrow=1, scale=0.9, rel_widths=c(0.4, 1.0),
                      labels=c("", "d"), align="h", axis="b")
bot_grid <- plot_grid(cite_plot_gt, sign_plot, nrow=1,
                      scale=0.9, rel_widths=c(1.0, 1.0),
                      labels=c("", "f"))
figure <- plot_grid(top_grid, mid_grid, bot_grid, nrow=3,
                    rel_heights=c(1.0, 1.3, 1.2),
                    labels=c("a", "c", "e"))
save_plot("img/figures/validation-figure.pdf", figure, ncol=3, nrow=3,
          base_height=2.5)

supp_figure <- plot_grid(auc_plot,
                         net_dens_plot,
                         rank_plot2,
                         nrow=3, ncol=1, labels=c("a", "b", "c"), scale=0.9)
save_plot("img/figures/validation-supp-figure.pdf", supp_figure, ncol=2, nrow=4,
          base_height=2)
