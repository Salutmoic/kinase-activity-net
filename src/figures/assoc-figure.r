library(dplyr)
library(ggplot2)
library(readr)
library(tibble)
library(tidyr)
library(scales)
library(viridis)
library(cowplot)
library(gridExtra)
library(e1071)

theme_set(theme_cowplot(font_size=10))

extract_data <- function(kin1, kin2, cor_data, psites){
    excor_i <- which(cor_data$node1 == kin1 &
                         cor_data$node2 == kin2)
    excor_sites1_str <- cor_data[excor_i, ]$kin1.pos
    excor_sites2_str <- cor_data[excor_i, ]$kin2.pos
    excor_sites1 <- as.integer(strsplit(excor_sites1_str, split=",")[[1]])
    excor_sites2 <- as.integer(strsplit(excor_sites2_str, split=",")[[1]])
    excor_scores1_str <- cor_data[excor_i, ]$all.weights1
    excor_scores2_str <- cor_data[excor_i, ]$all.weights2
    excor_scores1 <- as.numeric(strsplit(excor_scores1_str, split=",")[[1]])
    excor_scores2 <- as.numeric(strsplit(excor_scores2_str, split=",")[[1]])
    excor_pvals_str <- cor_data[excor_i, ]$all.p.vals
    excor_pvals <- as.numeric(strsplit(excor_pvals_str, split=",")[[1]])
    excor_cor_data <- as_tibble(data.frame(kin1=rep(kin1, length(excor_sites1)),
                                           kin2=rep(kin2, length(excor_sites2)),
                                           site1=excor_sites1,
                                           site2=excor_sites2,
                                           score1=excor_scores1,
                                           score2=excor_scores2,
                                           pval=excor_pvals))
    excor_cor_data$site1 <- reorder(excor_cor_data$site1,
                                    -excor_cor_data$score1, order=TRUE)
    excor_cor_data$site2 <- reorder(excor_cor_data$site2,
                                    excor_cor_data$score2, order=TRUE)
    excor_cor_data <- excor_cor_data %>%
        mutate(coreg.score=pval*score1*score2)
    excor_cor_data_tmp <- merge(excor_cor_data, psites, by.x=c("kin1", "site1"),
                                by.y=c("kin", "site"))
    excor_cor_data <- merge(excor_cor_data_tmp, psites, by.x=c("kin2", "site2"),
                            by.y=c("kin", "site"), suffixes=c("1", "2"))
    excor_cor_data$pos1 <- paste0(excor_cor_data$residue1, excor_cor_data$site1)
    excor_cor_data$pos2 <- paste0(excor_cor_data$residue2, excor_cor_data$site2)
    return(excor_cor_data)
}

build_func_data <- function(excor_cor_data){
    excor_func_data1 <- unique(excor_cor_data[, c("kin1", "site1", "score1")])
    names(excor_func_data1) <- c("kin", "site", "score")
    excor_func_data2 <- unique(excor_cor_data[, c("kin2", "site2", "score2")])
    names(excor_func_data2) <- c("kin", "site", "score")
    excor_func_data <- rbind(excor_func_data1, excor_func_data2)
    return(excor_func_data)
}

extract_func_data <- function(exkin, excor_func_data, order_sign, psites){
    kin_excor_func_data <- excor_func_data %>%
        filter(kin==exkin)
    kin_excor_func_data <- merge(kin_excor_func_data, psites)
    kin_excor_func_data$pos <- as.factor(paste0(kin_excor_func_data$residue,
                                                kin_excor_func_data$site))
    kin_excor_func_data$pos <- reorder(kin_excor_func_data$pos,
                                        order_sign*kin_excor_func_data$score,
                                        order=TRUE)
    return(kin_excor_func_data)
}

do_func_cor_plot <- function(kin1, kin2, kin3, excor_cor_data1, excor_cor_data2,
                              excor_func_data1, excor_func_data2, excor_func_data3){
    empty_df <- data.frame()
    empty_plot <- ggplot(empty_df) + geom_blank()

    heat_palette <- "viridis"
    max_val <- max(excor_cor_data1$pval, excor_cor_data2$pval, na.rm=TRUE)
    heat_plot1 <- excor_cor_data1 %>%
        ggplot(aes(x=site1, y=site2, fill=pval)) +
        geom_tile(colour="white", size=1.1) +
        scale_x_discrete(expand=c(0.01, 0.0, 0.01, 0.0)) +
        scale_y_discrete(expand=c(0.01, 0.0, 0.01, 0.0)) +
        scale_fill_viridis(limits=c(0, max_val), option=heat_palette) +
        ## scale_fill_continuous(limits=c(0, 1)) +
        theme_void() +
        theme(plot.margin=unit(c(0,0,0,0), "cm")) +
        labs(x=element_blank(), y=element_blank()) +
        guides(fill=FALSE)
        ## coord_fixed() +
    heat_plot2 <- excor_cor_data2 %>%
        ggplot(aes(x=site1, y=site2, fill=pval)) +
        geom_tile(colour="white", size=1.1) +
        scale_x_discrete(expand=c(0.01, 0.0, 0.01, 0.0)) +
        scale_y_discrete(expand=c(0.01, 0.0, 0.01, 0.0)) +
        scale_fill_viridis(limits=c(0, max_val), 
                           option=heat_palette) +
        ## scale_fill_continuous(limits=c(0, 1), breaks=c(0, 0.5, 1.0)) +
        theme_void() +
        theme(plot.margin=unit(c(0,0,0,0), "cm")) +
        labs(x=element_blank(), y=element_blank()) +
        ## coord_fixed() +
        guides(fill=guide_colourbar(expression(atop("cophospho.",
                                                    paste(-log[10], "(", italic(p), ")"))),
                                    title.position="top"))
    bar_fill <- "#000000"
    excor_func_plot1 <- excor_func_data1 %>%
        ggplot(aes(x=pos, y=score)) +
        geom_bar(stat="identity", fill=bar_fill, colour="white", size=1.1,
                 width=1.0) +
        theme(axis.line.x=element_blank(), plot.margin=unit(c(0, 0, 0, 0), "cm"),
              axis.ticks.x=element_blank(), axis.line.y=element_blank(),
              axis.ticks.y=element_blank(), axis.text.y=element_blank(),
              axis.title.y=element_blank()) +
        scale_x_discrete(kin1, position="top",
                         expand=c(0.01, 0.0, 0.01, 0.0)) +
        scale_y_continuous(expand=c(0.01, 0.0, 0.01, 0.0)) +
        ## coord_fixed(ratio=1/1.5) +
        labs(x=NULL, y="functional\nscore")
    excor_func_plot2 <- excor_func_data2 %>%
        ggplot(aes(x=pos, y=score)) +
        geom_bar(stat="identity", fill=bar_fill, colour="white", size=1.1,
                 width=1.0) +
        theme(axis.line.x=element_blank(), plot.margin=unit(c(0, 0, 0, 0), "cm"),
              axis.ticks.x=element_blank()) +
        scale_x_discrete(kin2, position="top",
                         expand=c(0.01, 0.0, 0.01, 0.0)) +
        scale_y_continuous(position="right", limits=c(0, 1), breaks=c(1.0, 0.0),
                           expand=c(0.01, 0.0, 0.01, 0.0)) +
        ## coord_fixed(ratio=1/1.5) +
        labs(x=NULL, y="functional\nscore")
    excor_func_plot3 <- excor_func_data3 %>%
        ggplot(aes(x=pos, y=score)) +
        geom_bar(stat="identity", fill=bar_fill, colour="white", size=1.1,
                 width=1.0) +
        scale_x_discrete(kin3, expand=c(0.01, 0.0, 0.01, 0.0)) +
        scale_y_reverse(position="left", limits=c(1, 0), breaks=c(1.0, 0.0),
                        expand=c(0.01, 0.0, 0.01, 0.0)) +
        labs(y="functional\nscore") +
        coord_flip() +
        theme(axis.line.y=element_blank(), plot.margin=unit(c(0, 0, 0, 0), "cm"),
              axis.ticks.y=element_blank(),
              axis.text.y=element_text(angle = 90, hjust = 0.5)) ## , aspect.ratio=4)
    plot_alignment1 <- align_plots(excor_func_plot1, heat_plot1, align="v",
                                   axis="lr")
    plot_alignment2 <- align_plots(excor_func_plot2, heat_plot2, align="v",
                                   axis="lr")
    plot_alignment3 <- align_plots(excor_func_plot3, plot_alignment1[[2]],
                                   align="h")
    bottom_row <- plot_grid(excor_func_plot3, plot_alignment1[[2]],
                            plot_alignment2[[2]], align="h",
                            ncol=3, rel_widths=c(0.45, 0.52, 1.0))
    reg_site_cor_plot <- grid.arrange(empty_plot, plot_alignment1[[1]],
                                      plot_alignment2[[1]], bottom_row,
                                      layout_matrix=matrix(c(1, 2, 3, 4, 4, 4),
                                                           nrow=2, byrow=TRUE),
                                      widths=c(0.45, 0.52, 1.0),
                                      heights=c(0.25, 1.0))
    return(reg_site_cor_plot)
}

reg_cor_dir <- "out/reg-site-cor-by-exp/"
reg_cor_272_file <- paste0(reg_cor_dir, "reg-site-cor-exp-272-annot.tsv")
reg_cor_293_file <- paste0(reg_cor_dir, "reg-site-cor-exp-293-annot.tsv")
reg_cor_272 <- read_tsv(reg_cor_272_file,
                        col_types=cols(
                            node1 = col_character(),
                            node2 = col_character(),
                            reg.site.cor.p = col_double(),
                            all.p.vals = col_character(),
                            all.weights1 = col_character(),
                            all.weights2 = col_character(),
                            kin1.pos = col_character(),
                            kin2.pos = col_character(),
                            cor.p.272 = col_double()
                        ))
reg_cor_293 <- read_tsv(reg_cor_293_file,
                        col_types=cols(
                            node1 = col_character(),
                            node2 = col_character(),
                            reg.site.cor.p = col_double(),
                            all.p.vals = col_character(),
                            all.weights1 = col_character(),
                            all.weights2 = col_character(),
                            kin1.pos = col_character(),
                            kin2.pos = col_character(),
                            cor.p.293 = col_double()
                        ))
true_pos <- read_tsv("data/direct-validation-set-omnipath.tsv", col_names=FALSE)
names(true_pos) <- c("node1", "node2")
true_pos <- true_pos %>% mutate(id=paste(node1, node2, sep="-"))
pride_sites <- read_tsv("data/pride-phosphosites.tsv", col_names=FALSE)
names(pride_sites) <- c("kin", "site", "residue", "motif")

atlas_coexp <- read_tsv("out/rna-tissue-cor.tsv")
gtex_coexp <- read_tsv("out/gtex-tissue-coexpr.tsv")

## Coregulation heatmaps
exkin1 <- "MAPK3"
exkin2 <- "MAPK4"
exkin3 <- "RPS6KA3"
excor_cor_data1 <- extract_data(exkin1, exkin3, reg_cor_293, pride_sites)
excor_cor_data2 <- extract_data(exkin2, exkin3, reg_cor_293, pride_sites)
excor_func_data1 <- build_func_data(excor_cor_data1)
excor_func_data2 <- build_func_data(excor_cor_data2)
kin1_excor_func_data <- extract_func_data(exkin1, excor_func_data1, -1,
                                          pride_sites)
kin2_excor_func_data <- extract_func_data(exkin2, excor_func_data2, -1,
                                           pride_sites)
kin3_excor_func_data <- extract_func_data(exkin3, excor_func_data1, 1,
                                           pride_sites)
reg_site_cor_plot <- do_func_cor_plot(exkin1, exkin2, exkin3, excor_cor_data1,
                                       excor_cor_data2, kin1_excor_func_data,
                                       kin2_excor_func_data,kin3_excor_func_data)

save_plot("img/figures/panels/assoc-reg-site-heatmap.pdf", reg_site_cor_plot)

## Coregulation scores
## excor_cor_data <- rbind(excor_cor_data1, excor_cor_data2) %>%
##     mutate(pair=paste(kin1, kin2, sep="-"))
## excor_cor_data$pair <- factor(excor_cor_data$pair,
##                               levels=c(paste(exkin1, exkin3, sep="-"),
##                                        paste(exkin2, exkin3, sep="-")))
excor_cor_data <- rbind(excor_cor_data1, excor_cor_data2) %>%
    group_by(kin1) %>%
    mutate(pair=paste0(kin1, "-", pos1, " / ", kin2, "-", pos2),
           is.max=(coreg.score == max(coreg.score, na.rm=TRUE)))
excor_cor_data$labelsize <- rep(NA, nrow(excor_cor_data))
excor_cor_data[which(excor_cor_data$is.max), "labelsize"] <- 2.5
excor_cor_data[which(!excor_cor_data$is.max), "labelsize"] <- 0

coreg_score_plot <- excor_cor_data %>%
    ggplot(aes(x=kin1, y=coreg.score, colour=kin1, label=pair)) +
    geom_text(aes(size=labelsize, hjust=(1+sign(coreg.score-1))/2,
                  vjust=(1+sign(1-coreg.score))/2, nudge_x=2, nudge_y=2),
              position=position_jitter()) +
    geom_jitter(aes(width=0.4, height=0.0, size=1 + as.integer(is.max))) +
    scale_colour_manual(values=c("#21908CFF", "#440154FF")) +
    scale_x_discrete(limits=rev(levels(excor_cor_data$kin1))) +
    scale_size_identity() +
    labs(x="kinase", y=paste("coregulation score w/", exkin3)) +
    coord_flip() +
    guides(colour=FALSE) +
    theme(axis.title.y=element_blank())

save_plot("img/figures/panels/assoc-coreg-score.pdf", coreg_score_plot)

## Coexpression
exkin1 <- "LYN"
exkin2 <- "YES1"
exkin3 <- "BTK"

atlas_exp <- read_tsv("data/rna_tissue.tsv")

kin_cor <- data.frame(kinase = c(exkin1, exkin2))
for (i in 1:2){
    exkin_values <- atlas_exp %>% filter(kinase == c(exkin1, exkin2)[i]) %>% select(value)
    exkin3_values <- atlas_exp %>% filter(kinase == exkin3) %>% select(value)
    cor_res <- cor(exkin_values$value, exkin3_values$value, method="spearman")
    kin_cor$text[i] <- paste("r =", round(cor_res, digits=2))
}

coexp_plot <- atlas_exp %>%
    filter(kinase %in% c(exkin1, exkin2)) %>%
    mutate(exkin3.value=rep(subset(atlas_exp, kinase==exkin3)$value, 2)) %>%
    ggplot(aes(x=value, y=exkin3.value)) +
    geom_smooth(method="lm", formula=y~x, colour="#8f2150") +
    geom_point() +
    facet_wrap("kinase", scales="free_x") +
    labs(x="RNA expression (TPM)", y=paste(exkin3, "RNA\nexpression (TPM)")) +
    theme(strip.background=element_blank(), strip.placement="outside",
          panel.spacing=unit(0.5, "cm")) +
    geom_text(data=kin_cor, aes(0, -15, label=text, hjust=0))

save_plot("img/figures/panels/assoc-coexp.pdf", coexp_plot)

spec_plot1 <- atlas_exp %>%
    filter(kinase %in% c(exkin1, exkin2, exkin3)) %>%
    ggplot(aes(x=value, fill=kinase, colour=kinase, alpha=0.7)) +
    geom_density() +
    xlab("RNA expression (TPM)") +
    guides(alpha=FALSE) +
    scale_fill_viridis(discrete=TRUE, direction=-1) +
    scale_colour_viridis(discrete=TRUE, direction=-1) +
    scale_x_continuous(expand=c(0.01, 0, 0.01, 0)) +
    scale_y_continuous(expand=c(0.01, 0, 0.01, 0)) +
    theme(legend.position=c(1, 1), legend.justification=c(1, 1))

save_plot("img/figures/panels/assoc-expr-spec-density.pdf", spec_plot1)

spec_dat_raw <- atlas_exp %>%
    group_by(kinase) %>%
    summarise(specificity=skewness(value))

spec_dat <- read_tsv("out/rna-tissue-selectivity.tsv")

spec_dat_sub <- spec_dat_raw %>%
    filter(kinase %in% c(exkin1, exkin2, exkin3))
spec_dat_sub$kinase <- factor(spec_dat_sub$kinase)
spec_plot2 <- spec_dat_sub %>%
    ggplot(aes(x=kinase, y=specificity, fill=kinase, colour=kinase)) +
    geom_bar(stat="identity") +
    scale_y_continuous("tissue specificity", expand=c(0, 0, 0, 0)) +
    scale_x_discrete(limits=rev(levels(spec_dat_sub$kinase))) +
    scale_fill_viridis(discrete=TRUE, direction=-1) +
    scale_colour_viridis(discrete=TRUE, direction=-1) +
    guides(fill=FALSE, colour=FALSE)  +
    theme(axis.title.y=element_blank()) +
    coord_flip()

save_plot("img/figures/panels/assoc-expr-spec-score.pdf", spec_plot2)

## Performance plot
coreg_293_perf_data <- reg_cor_293 %>%
    mutate(id=paste(node1, node2, sep="-"),
           is.true=id %in% true_pos$id,
           score=rep("coreg.293.score", length(node1))) %>%
    select(node1, node2, is.true, score, cor.p.293) %>%
    filter(!is.na(cor.p.293))
names(coreg_293_perf_data)[5] <- "value"
coreg_272_perf_data <- reg_cor_272 %>%
    mutate(id=paste(node1, node2, sep="-"),
           is.true=id %in% true_pos$id,
           score=rep("coreg.272.score", length(node1))) %>%
    select(node1, node2, is.true, score, cor.p.272) %>%
    filter(!is.na(cor.p.272))
names(coreg_272_perf_data)[5] <- "value"

min.coex <- min(atlas_coexp$prot.atlas.tissue, na.rm=TRUE)
max.coex <- max(atlas_coexp$prot.atlas.tissue, na.rm=TRUE)
atlas_coexp$prot.atlas.tissue.norm <- (atlas_coexp$prot.atlas.tissue-min.coex)/(max.coex-min.coex)

min.coex <- min(gtex_coexp$tissue.coexpr, na.rm=TRUE)
max.coex <- max(gtex_coexp$tissue.coexpr, na.rm=TRUE)
gtex_coexp$tissue.coexpr.norm <- (gtex_coexp$tissue.coexpr-min.coex)/(max.coex-min.coex)

tis_coexp_perf_data1 <- atlas_coexp %>%
    mutate(id=paste(node1, node2, sep="-"),
           is.true=id %in% true_pos$id,
           score=rep("prot.atlas", length(node1))) %>%
    select(node1, node2, is.true, score, prot.atlas.tissue) %>%
    filter(!is.na(prot.atlas.tissue))
names(tis_coexp_perf_data1)[5] <- "value"
tis_coexp_perf_data2 <- gtex_coexp %>%
    mutate(id=paste(node1, node2, sep="-"),
           is.true=id %in% true_pos$id,
           score=rep("gtex", length(node1))) %>%
    select(node1, node2, is.true, score, tissue.coexpr) %>%
    filter(!is.na(tissue.coexpr))
names(tis_coexp_perf_data2)[5] <- "value"
tis_spec_perf_data <- spec_dat %>%
    mutate(id=paste(node1, node2, sep="-"),
           is.true=id %in% true_pos$id,
           specificity=abs(tissue.sel1-tissue.sel2),
           score=rep("spec", length(node1))) %>%
    select(node1, node2, is.true, score, specificity) %>%
    filter(!is.na(specificity))
names(tis_spec_perf_data)[5] <- "value"
perf_data <- rbind(tis_coexp_perf_data1, tis_coexp_perf_data2,
                   tis_spec_perf_data, coreg_293_perf_data,
                   coreg_272_perf_data)
perf_data$is.true <- factor(perf_data$is.true,
                            levels=c(FALSE, TRUE),
                            labels=c("no known relationship",
                                     "known relationship"))
perf_data$score <- factor(perf_data$score,
                          levels=c("coreg.293.score", "coreg.272.score",
                                   "prot.atlas", "gtex", "spec"),
                          labels=c("max. coregulation score\n(CPTAC breast cancer)",
                                   "max. coregulation score\n(Wilkes et al. 2015)",
                                   "Protein Atlas tissue\nRNA coexpression",
                                   "GTEx RNA coexpression",
                                   "RNA expression tissue\nspecificity (absolute\ndifference)"))
perf_plot <- perf_data %>%
    ggplot(aes(x=1, y=value, fill=is.true)) +
    geom_boxplot(outlier.size=0.25) +
    facet_wrap("score", scales="free_x", nrow=1, strip.position="bottom") +
    scale_fill_manual("", values=c("#b3b3b3", "#8f2150"),
                      guide=guide_legend(direction="horizontal")) +
    theme(legend.position="bottom", legend.title=element_blank(),
          axis.ticks.y=element_blank(), axis.text.y=element_blank(),
          axis.line.y=element_blank(), axis.title=element_blank(),
          strip.background=element_blank(), strip.placement="outside",
          strip.text=element_text(vjust=1),
          panel.spacing=unit(0.5, "cm"),
          legend.margin=margin(0,0,0,0)) +
    coord_flip()

save_plot("img/figures/panels/assoc-perf.pdf", perf_plot, ncol=3)

top_left_grid <- plot_grid(reg_site_cor_plot, coreg_score_plot,
                           ncol=1, nrow=2, rel_heights=c(1.0, 0.33),
                           scale=0.9)
exp_align <- align_plots(coexp_plot, spec_plot1, align="v",
                         axis="lr")
spec_grid <- plot_grid(exp_align[[2]], spec_plot2, nrow=2, ncol=1, align="v",
                       scale=0.9)
top_right_grid <- plot_grid(exp_align[[1]], spec_grid, nrow=2,
                            ncol=1, labels=c("b", "c"), rel_heights=c(0.8, 1.0),
                            scale=0.9)
final_plot_top <- plot_grid(top_left_grid, top_right_grid,
                            nrow=1, ncol=2, labels=c("", ""))
final_plot <- plot_grid(final_plot_top, perf_plot, nrow=2, ncol=1,
                        labels=c("a", "d"), rel_heights=c(1.0, 0.29),
                        scale=c(1.0, 0.9))


save_plot("img/figures/assoc-figure.pdf", final_plot, nrow=2, ncol=2)
