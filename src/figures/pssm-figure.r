library(dplyr)
library(ggplot2)
library(readr)
library(tibble)
library(tidyr)
library(scales)
library(viridis)
library(cowplot)

theme_set(theme_cowplot(font_size=10))

calc_dcg <- function(func_scores, pssm_scores){
    n <- length(func_scores)
    func_scores_sort <- func_scores[order(pssm_scores, decreasing=TRUE)]
    raw_dcg <- cumsum(func_scores_sort/(log2((1:n)+1)))
    raw_dcg <- raw_dcg[order(order(pssm_scores, decreasing=TRUE))]
    return(raw_dcg)
}

pssms <- as_tibble(read.delim("out/kinase-pssm-annotated.tsv", as.is=TRUE))
pssms$pssm.source <- factor(pssms$pssm.source,
                            levels=c("sdr", "family", "own"),
                            labels=c("SDR similarity",
                                     "kinase family",
                                     "own"))
pssms$kinase.type <- factor(pssms$kinase.type,
                            levels=c("ST", "Y"),
                            labels=c("S/T", "Y"))

## Family/SDR assignments from David

sdr_sim_list <- readRDS("data/external/pssm-assignment/Brandon1.rds")
sdr_sim_lens <- unlist(lapply(sdr_sim_list, length))
sdr_sim_class_labs <- c('>-0.4','>-0.3','>-0.2','>-0.1','>0.0','>0.1','>0.2',
                        '>0.3','>0.4','>0.5','>0.6','>0.7','>0.8','>0.9',
                        'random','duplicate')
sdr_sim_classes <- c()
for (i in 1:length(sdr_sim_list)){
     print(i)
     lab <- sdr_sim_class_labs[i]
     len <- sdr_sim_lens[i]
     sdr_sim_classes <- c(sdr_sim_classes, rep(lab, len))
}
sdr_sim <- new_tibble(list(sdr.sim.class=as.factor(sdr_sim_classes),
                           pssm.distance=unlist(sdr_sim_list))) %>%
    filter(sdr.sim.class %in% c(">-0.4", ">-0.2", ">0.0", ">0.2",
                                ">0.4", ">0.6", ">0.8", "duplicate", "random"))

sdr_sim$sdr.sim.class <- factor(sdr_sim$sdr.sim.class,
                                levels=c(">-0.4", ">-0.2", ">0.0", ">0.2",
                                ">0.4", ">0.6", ">0.8", "duplicate", "random"))

sdr_sim_boxplot <- sdr_sim %>%
    ggplot(aes(x=sdr.sim.class, y=pssm.distance)) +
    geom_hline(yintercept=1.112, colour="#8F2150", alpha=1.0) +
    geom_boxplot(outlier.size=0.25, alpha=0.75) +
    labs(x="SDR similarity", y="PSSM distance") +
    theme(axis.text.x=element_text(angle=45, hjust=1))

save_plot("img/figures/panels/pssms-sdr-sim.pdf", sdr_sim_boxplot)

sdr_vs_fam_dat <- readRDS("data/external/pssm-assignment/Brandon3.rds")
sdr_vs_fam <- new_tibble(list(sdr.dist=sdr_vs_fam_dat[[1]],
                              fam.dist=sdr_vs_fam_dat[[2]]))

sdr_vs_fam_plot <- sdr_vs_fam %>%
    ggplot(aes(x=sdr.dist, y=fam.dist)) +
    geom_abline(slope=1, intercept=0, colour="#8F2150") +
    geom_point(alpha=0.75) +
    ylim(0, 1.6) +
    xlim(0, 1.6) +
    labs(x="PSSM distance\n(SDR)",
         y="PSSM distance\n(family)")

save_plot("img/figures/panels/pssms-sdr-vs-fam.pdf", sdr_vs_fam_plot)

## PSSM sources: own/computed, assigned by family, assigned by SDR
sources_plot <- pssms %>%
    select(node1, kinase.type, pssm.source) %>%
    filter(!is.na(kinase.type), !is.na(pssm.source)) %>%
    unique %>%
    ggplot(aes(x=kinase.type, fill=pssm.source)) +
    geom_bar(na.rm=TRUE) +
    ## scale_fill_manual("PSSM source", values=c("#e5c494", "#ffd92f", "#a6d854")) +
    scale_fill_viridis("PSSM source:", discrete=TRUE, direction=-1) +
    scale_y_continuous(label=comma, expand=c(0.0, 0.0, 0.01, 0.0), limits=c(0, 300)) +
    labs(x="kinase type") +
    theme(legend.position=c(1, 1), legend.justification=c(1, 1),
          legend.direction="horizontal",
          legend.margin=margin(0,0,0,0)) +
    coord_flip()

save_plot("img/figures/panels/pssms-sources.pdf", sources_plot)

## exsub_id <- "MAPKAPK2"
## exreg1_id <- "MAPK14"
## exreg2_id <- "CAMK1"
exsub_id <- "AKT1"
exreg1_id <- "PDPK1"
exreg2_id <- "PRKCG" #"GRK2"
exreg1_exsub_i <- which(pssms$node1 == exreg1_id &
                        pssms$node2 == exsub_id)
exsub_sites_str <- pssms[exreg1_exsub_i, ]$all.sites
exsub_sites <- strsplit(exsub_sites_str, split=",")[[1]]
exsub_pos <- as.numeric(gsub("[TS]", "", exsub_sites))
exsub_func_scores_str <- pssms[exreg1_exsub_i, ]$all.func.scores
exsub_func_scores <- as.numeric(strsplit(exsub_func_scores_str, split=",")[[1]])
exreg1_exsub_scores_str <- pssms[exreg1_exsub_i, ]$all.pssm.scores
exreg1_exsub_scores <- as.numeric(strsplit(exreg1_exsub_scores_str, split=",")[[1]])
exreg2_exsub_i <- which(pssms$node1 == exreg2_id &
                      pssms$node2 == exsub_id)
exreg2_exsub_scores_str <- pssms[exreg2_exsub_i, ]$all.pssm.scores
exreg2_exsub_scores <- as.numeric(strsplit(exreg2_exsub_scores_str, split=",")[[1]])
exsub_data_full <- as_tibble(data.frame(kinase=c(rep(exreg1_id, length(exsub_sites)),
                                                 rep(exreg2_id, length(exsub_sites))),
                                        site=c(exsub_sites, exsub_sites),
                                        score=c(exreg1_exsub_scores,
                                                exreg2_exsub_scores)))
exsub_func_data <- as_tibble(data.frame(site=exsub_sites,
                                        func.score=exsub_func_scores)) %>%
    filter(!is.na(func.score))
exsub_func_data$site <- as.factor(exsub_func_data$site)
exsub_func_data$site <- reorder(exsub_func_data$site,
                                -exsub_func_data$func.score)
exsub_data <- exsub_data_full %>%
    filter(site %in% exsub_func_data$site)
exsub_data$site <- factor(exsub_data$site,
                          levels=levels(exsub_func_data$site))
## levels(exsub_data$site) <- levels(exsub_func_data$site)

cols <- c(" "="#b3b3b3")
exsub_site_plot <- exsub_func_data %>%
    ggplot(aes(x=site, y=func.score)) +
    geom_bar(data=exsub_data,
             aes(x=site, y=score, fill=kinase),
             stat="identity",
             position=position_dodge()) +
             ## width=1.5) +
    ## geom_segment(aes(x=site, y=0, xend=site, yend=func.score),
    ##              colour="#b3b3b3", size=1.25) +
    geom_segment(aes(x=site, y=0, xend=site, yend=func.score,
                     colour=" "), size=1.25) +
    ## geom_point(colour="#b3b3b3", size=2.5) +
    geom_point(aes(colour=" "), size=2.5) +
    scale_colour_manual("functional score:", values=cols) +
    scale_y_continuous("functional score",
                       sec.axis=sec_axis(~., name="PSSM score"),
                       limits=c(0, 1.1),
                       expand=c(0.0, 0.0, 0.01, 0.0)) +
    scale_fill_manual("PSSM score:", values=c("#21908CFF", "#440154FF")) +
    guides(fill=guide_legend(title.position="left"),
           colour=guide_legend(title.position="left")) +
    xlab(paste(exsub_id, "phosphosite")) +
    theme(legend.position=c(1, 1), legend.justification=c(1, 1),
          ## legend.title=element_blank(),
          legend.spacing=unit(0.01, "mm"),
          legend.margin=margin(t=0, r=10, b=0, l=0),
          legend.direction="horizontal")

save_plot("img/figures/panels/pssms-example-sites.pdf", exsub_site_plot, ncol=3)

exsub_final_data <- pssms %>%
    filter(node1 %in% c(exreg1_id, exreg2_id),
           node2 == exsub_id) %>%
    select(node1, max.pssm.score, dcg) %>%
    gather(score, value, -node1, convert=TRUE)
exsub_final_data$score <- factor(exsub_final_data$score,
                                 levels=c("max.pssm.score", "dcg"),
                                 labels=c("max. PSSM score",
                                          "DCG score"))

exsub_final_plot <- exsub_final_data %>%
    ggplot(aes(x=score, y=value, fill=node1)) +
    geom_bar(stat="identity", position=position_dodge()) +
    ## scale_fill_brewer(palette="Set2") +
    scale_fill_manual(values=c("#21908CFF", "#440154FF")) +
    scale_y_continuous("score value", limits=c(0, 1),
                       expand=c(0.0, 0.0, 0.0, 0.0)) +
    xlab(paste("regulation of", exsub_id)) +
    guides(fill=guide_legend("kinase")) +
    theme(legend.title=element_blank(),
          legend.position=c(1, 1), legend.justification=c(1, 1),
          legend.direction="horizontal",
          legend.margin=margin(0,0,0,0))

save_plot("img/figures/panels/pssms-example-scores.pdf", exsub_final_plot)

dummy_dcg_data <- as_tibble(data.frame(kinase=c(rep("best case",
                                                    nrow(exsub_func_data)),
                                                rep("worst case",
                                                    nrow(exsub_func_data))),
                                       site=exsub_data$site,
                                       score=c(rank(exsub_func_data$func.score),
                                               rank(-exsub_func_data$func.score))))

dcg_data <- rbind(exsub_data, dummy_dcg_data) %>%
    group_by(kinase) %>%
    mutate(dcg=calc_dcg(exsub_func_data$func.score, score),
           is.dummy=!(kinase %in% c(exreg1_id, exreg2_id))) %>%
    arrange(kinase, -score) %>%
    group_by(kinase)

dcg_data[which(!dcg_data$is.dummy), "linetype"] <- "solid"
dcg_data[which(dcg_data$is.dummy), "linetype"] <- "22"

dcg_data$kinase <- factor(dcg_data$kinase,
                          levels=c(exreg1_id, exreg2_id, "best case", "worst case"))

dcg_plot <- dcg_data %>%
    ggplot(aes(x=rep(1:nrow(exsub_func_data), 4), y=dcg, colour=kinase,
               group=kinase)) +
    geom_point(aes(size=rep(c(rep(0, 12), 2.5), 4))) +
    ## geom_step(aes(alpha=rep(1:nrow(exsub_func_data), 4))) +
    geom_step(aes(linetype=linetype)) +
    scale_shape_identity() +
    scale_size_identity() +
    scale_linetype_identity() +
    guides(linetype=FALSE, colour=guide_legend(keyheight=0.7)) +
    ## scale_alpha_continuous(guide=FALSE) +
    ## scale_colour_brewer(palette="Set2") +
    scale_colour_manual(values=c("#21908CFF", "#440154FF",
                                 "#8F2125FF", "#5B8F21FF")) +
    labs(x=paste("PSSM-based rank of", exsub_id, "phosphosite"),
         y="DCG of phosphosite\nfunctional score") +
    theme(legend.position=c(1, 0), legend.justification=c(1, 0),
          legend.title=element_blank())

save_plot("img/figures/panels/pssms-dcg.pdf", dcg_plot)

true_pos <- read_tsv("data/direct-validation-set-omnipath.tsv", col_names=FALSE)
names(true_pos) <- c("node1", "node2")
true_pos <- true_pos %>% mutate(id=paste(node1, node2, sep="-"))

perf_data <- pssms %>%
    mutate(id=paste(node1, node2, sep="-"),
           is.true=id %in% true_pos$id) %>%
    select(node1, node2, max.pssm.score, dcg, max.func.score, is.true) %>%
    mutate(interaction=max.pssm.score*dcg*max.func.score) %>%
    filter(!is.na(max.pssm.score)) %>%
    gather(score, value, -node1, -node2, -is.true)
## pvals <- c()
## for (s in sort(unique(perf_data$score))){
##     perf_data_sub <- subset(perf_data, score==s)
##     k <- kruskal.test(value~is.true, data=perf_data_sub)
##     pvals <- c(pvals, k$p.value)
## }

perf_data$score <- factor(perf_data$score,
                          levels=c("max.pssm.score", "max.func.score", "dcg",
                                   "interaction"),
                          labels=c("max. PSSM score",
                                   "max. functional score",
                                   "DCG score",
                                   "three-way product"))
perf_data$is.true <- factor(perf_data$is.true,
                            levels=c(FALSE, TRUE),
                            labels=c("no known relationship",
                                     "known relationship"))


perf_plot <- perf_data %>%
    filter(score != "three-way product") %>%
    ggplot(aes(x=1, y=value, fill=is.true)) +
    geom_boxplot(outlier.size=0.25) +
    facet_wrap("score", scales="free_x", nrow=1, strip.position="bottom") +
    scale_fill_manual("", values=c("#b3b3b3", "#8f2150"),
                      guide=guide_legend(direction="horizontal")) +
    theme(legend.position="bottom", legend.title=element_blank(),
          axis.ticks.y=element_blank(), axis.text.y=element_blank(),
          axis.line.y=element_blank(), axis.title=element_blank(),
          strip.background=element_blank(), strip.placement="outside",
          panel.spacing=unit(0.5, "cm"),
          strip.text=element_text(vjust=1),
          legend.margin=margin(0,0,0,0)) +
    coord_flip()

save_plot("img/figures/panels/pssms-perf.pdf", perf_plot, ncol=3)

final_plot <- plot_grid(plot_grid(sdr_sim_boxplot, sdr_vs_fam_plot,
                                  ncol=2, labels=c("", "b"), scale=0.9),
                        sources_plot,
                        exsub_site_plot,
                        plot_grid(dcg_plot, exsub_final_plot, ncol=2,
                                  labels=c("", "f"), align="h",
                                  scale=0.9),
                        perf_plot,
                        ncol=1, nrow=5, labels=c("a", "c", "d", "e", "g"),
                        rel_heights=c(1.0, 0.5, 1.0,  1.0, 0.65),
                        scale=c(1.0, 0.9, 0.9, 1.0, 0.9))
save_plot("img/figures/pssms-figure.pdf", final_plot, ncol=3, nrow=5,
          base_height=2.0)
