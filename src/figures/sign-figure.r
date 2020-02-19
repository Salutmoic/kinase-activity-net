library(dplyr)
library(ggplot2)
library(readr)
library(tibble)
library(tidyr)
library(scales)
library(viridis)
library(cowplot)
library(grid)
library(gridExtra)

theme_set(theme_cowplot(font_size=10))

calc_dcg <- function(func_scores, pssm_scores){
    n <- length(func_scores)
    func_scores_sort <- func_scores[order(pssm_scores, decreasing=TRUE)]
    raw_dcg <- cumsum(func_scores_sort/(log2((1:n)+1)))
    ## raw_dcg <- cumsum((2^func_scores_sort - 1)/(log2((1:n)+1)))
    raw_dcg <- raw_dcg[order(order(pssm_scores, decreasing=TRUE))]
    return(raw_dcg)
}
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
    excor_signs1_str <- cor_data[excor_i, ]$all.signs1
    excor_signs2_str <- cor_data[excor_i, ]$all.signs2
    excor_signs1 <- as.numeric(strsplit(excor_signs1_str, split=",")[[1]])
    excor_signs2 <- as.numeric(strsplit(excor_signs2_str, split=",")[[1]])
    excor_ests_str <- cor_data[excor_i, ]$all.cor.ests
    excor_ests <- as.numeric(strsplit(excor_ests_str, split=",")[[1]])
    excor_cor_data <- as_tibble(data.frame(kin1=rep(kin1, length(excor_sites1)),
                                           kin2=rep(kin2, length(excor_sites2)),
                                           site1=excor_sites1,
                                           site2=excor_sites2,
                                           score1=excor_scores1,
                                           score2=excor_scores2,
                                           sign1=excor_signs1,
                                           sign2=excor_signs2,
                                           cor.est=excor_ests))
    excor_cor_data$site1 <- reorder(excor_cor_data$site1,
                                    -excor_cor_data$score1, order=TRUE)
    excor_cor_data$site2 <- reorder(excor_cor_data$site2,
                                    excor_cor_data$score2, order=TRUE)
    excor_cor_data <- excor_cor_data %>%
        mutate(coreg.score=cor.est*score1*score2*sign(sign1)*sign(sign2))
    excor_cor_data_tmp <- merge(excor_cor_data, psites, by.x=c("kin1", "site1"),
                                by.y=c("kin", "site"))
    excor_cor_data <- merge(excor_cor_data_tmp, psites, by.x=c("kin2", "site2"),
                            by.y=c("kin", "site"), suffixes=c("1", "2"))
    excor_cor_data$pos1 <- paste0(excor_cor_data$residue1, excor_cor_data$site1)
    excor_cor_data$pos2 <- paste0(excor_cor_data$residue2, excor_cor_data$site2)
    return(excor_cor_data)
}

build_func_data <- function(excor_cor_data){
    excor_func_data1 <- unique(excor_cor_data[, c("kin1", "site1", "score1",
                                                  "sign1")])
    names(excor_func_data1) <- c("kin", "site", "score", "sign")
    excor_func_data2 <- unique(excor_cor_data[, c("kin2", "site2", "score2",
                                                  "sign2")])
    names(excor_func_data2) <- c("kin", "site", "score", "sign")
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

do_func_cor_plot_vert <- function(reg_kin, sub_kin1, sub_kin2, excor_cor_data1,
                                  excor_cor_data2, excor_func_data_reg,
                                  excor_func_data_sub1, excor_func_data_sub2){
    empty_df <- data.frame()
    empty_plot <- ggplot(empty_df) + geom_blank()

    min_cor <- min(min(excor_cor_data1$cor.est, na.rm=TRUE),
                   min(excor_cor_data2$cor.est, na.rm=TRUE))
    max_cor <- max(max(excor_cor_data1$cor.est, na.rm=TRUE),
                   max(excor_cor_data2$cor.est, na.rm=TRUE))
    cor_lim <- ceiling(max(abs(min_cor), max_cor))
    heat_palette <- "PRGn"
    heat_plot1 <- excor_cor_data1 %>%
        ggplot(aes(x=site1, y=site2, fill=cor.est)) +
        geom_tile(colour="white", size=1.1) +
        scale_x_discrete(expand=c(0.01, 0.0, 0.01, 0.0)) +
        scale_y_discrete(expand=c(0.01, 0.0, 0.01, 0.0)) +
        scale_fill_distiller(limits=c(-cor_lim, cor_lim),
                             palette=heat_palette,
                             breaks=c(-cor_lim, 0, cor_lim),
                             na.value="darkgrey") +
        ## scale_fill_continuous(limits=c(0, 1)) +
        theme_void() +
        theme(plot.margin=unit(c(0,0,0,0), "cm")) +
        labs(x=element_blank(), y=element_blank()) +
        guides(fill=FALSE)
        ## coord_fixed() +
    heat_plot2 <- excor_cor_data2 %>%
        ggplot(aes(x=site1, y=site2, fill=cor.est)) +
        geom_tile(colour="white", size=1.1) +
        scale_x_discrete(expand=c(0.01, 0.0, 0.01, 0.0)) +
        scale_y_discrete(expand=c(0.01, 0.0, 0.01, 0.0)) +
        scale_fill_distiller(limits=c(-cor_lim, cor_lim),
                             palette=heat_palette,
                             breaks=c(-cor_lim, 0, cor_lim),
                             na.value="darkgrey") +
        ## scale_fill_continuous(limits=c(0, 1), breaks=c(0, 0.5, 1.0)) +
        theme_void() +
        theme(plot.margin=unit(c(0,0,0,0), "cm"),
              legend.title=element_text(size=10),
              legend.text=element_text(size=8),
              legend.direction="vertical",
              legend.margin=unit(c(0,0,0,0), "cm"),
              legend.position="right",
              ## legend.position=c(1, 0.325),
              ## legend.justification=c(0, 0.5),
              legend.key.height=unit(0.35, "cm"),
              legend.key.width=unit(0.3, "cm")) +
        labs(x=element_blank(), y=element_blank()) +
        ## coord_fixed() +
        guides(fill=guide_colourbar("coregulation\ncoefficient",
                                    title.position="top"))
    sign_palette <- "RdYlBu"
    excor_func_plot_reg <- excor_func_data_reg %>%
        ggplot(aes(x=pos, y=score, fill=sign)) +
        geom_bar(stat="identity", colour="white", size=1.1,
                 width=1.0) +
        theme(axis.line.x=element_blank(), plot.margin=unit(c(0, 0, 0, 0), "cm"),
              axis.ticks.x=element_blank(), ## axis.title=element_text(size=10),
              ## axis.text=element_text(size=8),
              axis.text.x=element_text(angle=270, hjust=1.0),
              axis.text.x.top = element_text(vjust = 0.5),
              ## legend.title=element_text(size=10),
              ## legend.text=element_text(size=8),
              legend.key.height=unit(0.35, "cm"),
              legend.key.width=unit(0.3, "cm"),
              legend.margin=unit(c(0,0,0,0), "cm"),
              legend.position=c(1.125,-0.4),
              legend.justification=c(0, 1)) +
        scale_x_discrete(reg_kin, position="top",
                         expand=c(0.01, 0.0, 0.01, 0.0)) +
        scale_y_continuous(position="right", limits=c(0, 1), breaks=c(1.0, 0.0),
                           expand=c(0.01, 0.0, 0.01, 0.0)) +
        scale_fill_distiller("site sign\nprediction", palette=sign_palette) +
        ## coord_fixed(ratio=1/1.5) +
        labs(x=NULL, y="functional\nscore")
    excor_func_plot_sub1 <- excor_func_data_sub1 %>%
        ggplot(aes(x=pos, y=score, fill=sign)) +
        geom_bar(stat="identity", colour="white", size=1.1,
                 width=1.0) +
        scale_x_discrete(sub_kin1, expand=c(0.01, 0.0, 0.01, 0.0)) +
        scale_y_reverse(position="left", limits=c(1, 0), breaks=c(1.0, 0.0),
                        expand=c(0.01, 0.0, 0.01, 0.0)) +
        scale_fill_distiller(guide=FALSE, palette=sign_palette) +
        coord_flip() +
        theme(axis.line=element_blank(), plot.margin=unit(c(0, 0, 0, 0), "cm"),
              axis.ticks=element_blank(),
              ## axis.text.y=element_text(size = 8, angle = 0),
              axis.text.y=element_text(angle = 0),
              axis.text.x=element_blank(),
              ## axis.title.y=element_text(size=10),
              axis.title.x=element_blank()) ## , aspect.ratio=4)
    excor_func_plot_sub2 <- excor_func_data_sub2 %>%
        ggplot(aes(x=pos, y=score, fill=sign)) +
        geom_bar(stat="identity", colour="white", size=1.1,
                 width=1.0) +
        scale_x_discrete(sub_kin2, expand=c(0.01, 0.0, 0.01, 0.0)) +
        scale_y_reverse(position="left", limits=c(1, 0), breaks=c(1.0, 0.0),
                        expand=c(0.01, 0.0, 0.01, 0.0)) +
        scale_fill_distiller(guide=FALSE, palette=sign_palette) +
        labs(y="functional\nscore") +
        coord_flip() +
        theme(axis.line.y=element_blank(), plot.margin=unit(c(0, 0, 0, 0), "cm"),
              axis.ticks.y=element_blank(),
              axis.text.y=element_text(angle = 0))## ,
              ## axis.title=element_text(size=10),
              ## axis.text=element_text(size=8)) ## , aspect.ratio=4)
    plot_alignment1 <- align_plots(excor_func_plot_reg, heat_plot1,
                                   heat_plot2,
                                   align="v", axis="lr")
    plot_alignment2 <- align_plots(excor_func_plot_sub1, plot_alignment1[[2]],
                                   align="h")
    plot_alignment3 <- align_plots(excor_func_plot_sub2, plot_alignment1[[3]],
                                   align="h")
    ## bottom_row <- plot_grid(excor_func_plot3, plot_alignment1[[2]],
    ##                         plot_alignment2[[2]], align="h",
    ##                         ncol=3, rel_widths=c(0.45, 0.52, 1.0))
    reg_site_cor_plot <- grid.arrange(empty_plot, plot_alignment1[[1]],
                                      plot_alignment2[[1]], plot_alignment2[[2]],
                                      plot_alignment3[[1]], plot_alignment3[[2]],
                                      layout_matrix=matrix(c(1, 2, 3, 4, 5, 6),
                                                           nrow=3, byrow=TRUE),
                                      widths=c(0.5, 1.0),
                                      heights=c(0.75, 1.0, 1.2))
    return(reg_site_cor_plot)
}

pssms <- as_tibble(read.delim("out/kinase-pssm-signed-annotated.tsv", as.is=TRUE))
pssms$kinase.type <- factor(pssms$kinase.type,
                            levels=c("ST", "Y"),
                            labels=c("S/T", "Y"))

exsub_id <- "MAP2K1"
exreg1_id <- "BRAF"
exreg2_id <- "MAPK3"
## exsub_id <- "ABL1"
## exreg1_id <- "PAK2"
## exreg2_id <- "STK4"
exreg1_exsub_i <- which(pssms$node1 == exreg1_id &
                        pssms$node2 == exsub_id)
exsub_sites_str <- pssms[exreg1_exsub_i, ]$all.sites
exsub_sites <- strsplit(exsub_sites_str, split=",")[[1]]
exsub_pos <- as.numeric(gsub("[TS]", "", exsub_sites))
exsub_func_scores_str <- pssms[exreg1_exsub_i, ]$all.func.scores
exsub_func_scores <- as.numeric(strsplit(exsub_func_scores_str, split=",")[[1]])
exsub_signs_str <- pssms[exreg1_exsub_i, ]$all.signs
exsub_signs <- as.numeric(strsplit(exsub_signs_str, split=",")[[1]])
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
                                        func.score=exsub_func_scores,
                                        sign=exsub_signs)) %>%
    mutate(signed.func.score=sign(sign)*func.score) %>%
    filter(!is.na(signed.func.score), !is.na(sign))
exsub_func_data$site <- as.factor(exsub_func_data$site)
exsub_func_data$site <- reorder(exsub_func_data$site,
                                -exsub_func_data$signed.func.score)
exsub_data <- exsub_data_full %>%
    filter(site %in% exsub_func_data$site)
exsub_data$site <- factor(exsub_data$site,
                          levels=levels(exsub_func_data$site))

cols <- c(" "="#b3b3b3")
exsub_site_plot <- exsub_func_data %>%
    ggplot(aes(x=site, y=signed.func.score)) +
    geom_bar(data=exsub_data,
             aes(x=site, y=score, fill=kinase),
             stat="identity",
             position=position_dodge()) +
    geom_segment(aes(x=site, y=0, xend=site, yend=signed.func.score,
                     colour=" "), size=1.25) +
    geom_point(aes(colour=" "), size=2.5) +
    scale_colour_manual("signed functional score:", values=cols) +
    scale_fill_manual("PSSM score:", values=c("#440154FF", "#21908CFF")) +
    scale_y_continuous("signed functional score",
                       sec.axis=sec_axis(~., name="PSSM score")) +
    xlab(paste(exsub_id, "phosphosite")) +
    guides(fill=guide_legend(title.position="left"),
           colour=guide_legend(title.position="left")) +
    theme(legend.position=c(0, 0), legend.justification=c(0, 0),
          axis.text.x=element_text(angle=45, hjust=1.0),
          legend.spacing=unit(0.01, "mm"),
          legend.margin=margin(t=0, r=0, b=5, l=10),
          legend.direction="horizontal")

save_plot("img/figures/panels/signed-pssms-example-sites.pdf", exsub_site_plot, ncol=3)

exsub_final_data <- pssms %>%
    filter(node1 %in% c(exreg1_id, exreg2_id),
           node2 == exsub_id) %>%
    select(node1, signed.dcg)
##     gather(score, value, -node1, convert=TRUE)
## exsub_final_data$score <- factor(exsub_final_data$score,
##                                  levels=c("signed.dcg"),
##                                  labels=c("signed DCG score"))

exsub_final_plot <- exsub_final_data %>%
    ggplot(aes(x=node1, y=signed.dcg, fill=node1)) +
    geom_hline(yintercept=0.0, linetype=2, colour="grey") +
    geom_bar(stat="identity", position=position_dodge()) +
    scale_fill_manual(values=c("#440154FF", "#21908CFF")) +
    scale_y_continuous(paste("signed DCG score for", exsub_id),
                       limits=c(-1, 1), expand=c(0.01, 0.0, 0.01, 0.0)) +
    scale_x_discrete("", expand=c(0.01, 0.0, 0.01, 0.0)) +
    guides(fill=FALSE) +
    ## coord_flip() +
    theme(legend.title=element_blank(),
          legend.position=c(1, 1), legend.justification=c(1, 1),
          legend.direction="horizontal",
          axis.text.x=element_text(angle=45, hjust=1.0),
          plot.margin=margin(0,0,0,0))

save_plot("img/figures/panels/sign-pssms-example-scores.pdf", exsub_final_plot)

dummy_dcg_data <- as_tibble(data.frame(kinase=c(rep("best positive",
                                                    nrow(exsub_func_data)),
                                                rep("best negative",
                                                    nrow(exsub_func_data))),
                                       site=exsub_data$site,
                                       score=c(rank(exsub_func_data$signed.func.score),
                                               rank(-exsub_func_data$signed.func.score))))

dcg_data <- rbind(exsub_data, dummy_dcg_data) %>%
    group_by(kinase) %>%
    mutate(dcg=calc_dcg(exsub_func_data$signed.func.score, score),
           is.dummy=!(kinase %in% c(exreg1_id, exreg2_id)),
           is.max=abs(dcg)==max(abs(dcg))) %>%
    arrange(kinase, -score) %>%
    group_by(kinase)

dcg_data$used <- rep(NA, nrow(dcg_data))
dcg_data$x <- rep(NA, nrow(dcg_data))
for (kin in dcg_data$kinase){
    max_row <- which(dcg_data$kinase == kin & dcg_data$is.max)
    j <- 1
    for (i in which(dcg_data$kinase == kin)){
        dcg_data$used[i] <- i <= max_row
        dcg_data$x[i] <- j
        j <- j+1
    }
}
dcg_data_dup <- dcg_data[which(dcg_data$is.max),]
dcg_data_dup$used <- rep(FALSE, nrow(dcg_data_dup))
dcg_data_dup$is.max <- rep(FALSE, nrow(dcg_data_dup))
dcg_data <- rbind(dcg_data, dcg_data_dup) %>% arrange(kinase, -score)

dcg_data$linetype <- rep(NA, nrow(dcg_data))
dcg_data[which(dcg_data$used & !dcg_data$is.dummy), "linetype"] <- "solid"
dcg_data[which(dcg_data$used & dcg_data$is.dummy), "linetype"] <- "32"
dcg_data[which(!dcg_data$used), "linetype"] <- "12"
dcg_data$pointsize <- rep(NA, nrow(dcg_data))
dcg_data[which(dcg_data$is.max), "pointsize"] <- 2.5
dcg_data[which(!dcg_data$is.max), "pointsize"] <- 0

dcg_data$kinase <- factor(dcg_data$kinase,
                          levels=c(exreg2_id, exreg1_id, "best negative",
                                   "best positive"))
dcg_data$cases <- paste(dcg_data$kinase, dcg_data$used)

dcg_plot <- dcg_data %>%
    ggplot(aes(x=x, y=dcg, colour=kinase,
               group=cases)) +
    geom_hline(yintercept=0.0, linetype=2, colour="grey") +
    geom_step(aes(linetype=linetype)) +
    geom_point(aes(size=pointsize)) +
    scale_size_identity() +
    scale_linetype_identity() +
    guides(colour=guide_legend(keyheight=0.7)) + ## , nrow=2, ncol=2)) +
    scale_colour_manual(values=c("#21908CFF", "#440154FF",
                                 "#5B8F21FF", "#8F2125FF")) +
    labs(x=paste0("PSSM-based rank of ", exsub_id, " phosphosite"),
         y="DCG of phosphosite\nsigned functional score") +
    ## theme(legend.position=c(0.55, 0.6), legend.justification=c(0.5, 1),
    ##       legend.title=element_blank(), plot.margin=margin(0,0,0,0)) + ## ,
    theme(legend.position=c(1, 1), legend.justification=c(0.5, 1),
          legend.title=element_blank(), plot.margin=margin(0,0,0,0)) + ## ,
    scale_x_continuous(limits=c(0, 15), expand=c(0.01, 0.0, 0.01, 0.0))

save_plot("img/figures/panels/sign-dcg.pdf", dcg_plot)

reg_cor_dir <- "out/reg-site-cor-by-exp/"
reg_cor_293_file <- paste0(reg_cor_dir, "reg-site-cor-sign-exp-293-annot.tsv")
reg_cor_293 <- read_tsv(reg_cor_293_file,
                        col_types=cols(
                            node1 = col_character(),
                            node2 = col_character(),
                            cor.est.293 = col_double(),
                            all.cor.ests = col_character(),
                            all.weights1 = col_character(),
                            all.weights2 = col_character(),
                            all.signs1 = col_character(),
                            all.signs2 = col_character(),
                            kin1.pos = col_character(),
                            kin2.pos = col_character()
                        ))
reg_cor_272_file <- paste0(reg_cor_dir, "reg-site-cor-sign-exp-272-annot.tsv")
reg_cor_272 <- read_tsv(reg_cor_272_file,
                        col_types=cols(
                            node1 = col_character(),
                            node2 = col_character(),
                            cor.est.272 = col_double(),
                            all.cor.ests = col_character(),
                            all.weights1 = col_character(),
                            all.weights2 = col_character(),
                            all.signs1 = col_character(),
                            all.signs2 = col_character(),
                            kin1.pos = col_character(),
                            kin2.pos = col_character()
                        ))
pride_sites <- read_tsv("data/pride-phosphosites.tsv", col_names=FALSE)
names(pride_sites) <- c("kin", "site", "residue", "motif")

exreg_id <- "CDK1"
exinhib_id <- "MAP2K2"
exact_id <- "MAPK6"
## exreg_id <- "PLK1"
## exact_id <- "BUB1B"
## exinhib_id <- "PKMYT1"
excor_cor_data1 <- extract_data(exreg_id, exinhib_id, reg_cor_293, pride_sites)
excor_cor_data2 <- extract_data(exreg_id, exact_id, reg_cor_293, pride_sites)
excor_func_data1 <- build_func_data(excor_cor_data1)
excor_func_data2 <- build_func_data(excor_cor_data2)
reg_excor_func_data <- extract_func_data(exreg_id, excor_func_data1, -1,
                                           pride_sites)
subkin1_excor_func_data <- extract_func_data(exinhib_id, excor_func_data1, 1,
                                             pride_sites)
subkin2_excor_func_data <- extract_func_data(exact_id, excor_func_data2, 1,
                                             pride_sites)
reg_site_cor_plot <- do_func_cor_plot_vert(exreg_id, exinhib_id, exact_id,
                                           excor_cor_data1, excor_cor_data2,
                                           reg_excor_func_data,
                                           subkin1_excor_func_data,
                                           subkin2_excor_func_data)

save_plot("img/figures/panels/sign-assoc-reg-site-cor.pdf", reg_site_cor_plot)

## Coregulation scores
excor_cor_data <- rbind(excor_cor_data1, excor_cor_data2) %>%
    group_by(kin1) %>%
    mutate(pair=paste0(pos2, " (", kin1, " ", pos1, ")"),
           is.max=(coreg.score == max(coreg.score, na.rm=TRUE) |
                   coreg.score == min(coreg.score, na.rm=TRUE)))
excor_cor_data$labelsize <- rep(NA, nrow(excor_cor_data))
excor_cor_data[which(excor_cor_data$is.max), "labelsize"] <- 2.5
excor_cor_data[which(!excor_cor_data$is.max), "labelsize"] <- 0

coreg_score_plot <- excor_cor_data %>%
    ggplot(aes(x=kin2, y=coreg.score, colour=kin2, label=pair)) +
    geom_hline(yintercept=0, linetype=2, colour="grey") +
    geom_text(aes(size=labelsize, hjust=(1+sign(coreg.score))/2,
                  vjust=(1+sign(coreg.score))/2, nudge_x=1, nudge_y=2),
              position=position_jitter()) +
    geom_jitter(aes(width=0.4, height=0.0, size=1 + as.integer(is.max))) +
    scale_size_identity() +
    scale_colour_manual(values=c("#21908CFF", "#440154FF")) +
    labs(x="", y=paste(exreg_id, "coregulation sign score")) +
    coord_flip() +
    guides(colour=FALSE) +
    theme(axis.title.y=element_blank(),
          plot.margin=margin(0,0,0,0))

save_plot("img/figures/panels/sign-coreg-score.pdf", coreg_score_plot)

val_set <- read_tsv("data/sign-validation-set-omnipath.tsv", col_names=FALSE)
names(val_set) <- c("node1", "node2", "sign")
val_set <- val_set %>% mutate(id=paste(node1, node2, sep="-"))

pssms_perf_data <- pssms %>%
    mutate(id=paste(node1, node2, sep="-")) %>%
    inner_join(val_set, by=c("node1", "node2")) %>%
    filter(!is.na(signed.dcg)) %>%
    select(node1, node2, signed.dcg, signed.func.score, sign) %>%
    gather(score, value, -node1, -node2, -sign)

reg_cor_293_perf_data <- reg_cor_293 %>%
    mutate(id=paste(node1, node2, sep="-")) %>%
    inner_join(val_set, by=c("node1", "node2")) %>%
    filter(!is.na(cor.est.293), abs(cor.est.293)<10) %>%
    select(node1, node2, cor.est.293, sign) %>%
    gather(score, value, -node1, -node2, -sign)

reg_cor_272_perf_data <- reg_cor_272 %>%
    mutate(id=paste(node1, node2, sep="-")) %>%
    inner_join(val_set, by=c("node1", "node2")) %>%
    filter(!is.na(cor.est.272)) %>%
    select(node1, node2, cor.est.272, sign) %>%
    gather(score, value, -node1, -node2, -sign)

perf_data <- rbind(pssms_perf_data, reg_cor_293_perf_data, reg_cor_272_perf_data)
perf_data$score <- factor(perf_data$score,
                          levels=c("cor.est.293", "cor.est.272",
                                   "signed.func.score", "signed.dcg"),
                          labels=c("signed coregulation score\n(CPTAC breast cancer)",
                                   "signed coregulation score\n(Wilkes et al. 2015)",
                                   "strongest signed functional score",
                                   "signed DCG"))
perf_data$sign <- factor(perf_data$sign,
                         levels=c("activates", "inhibits"),
                         labels=c("activating relationship",
                                  "inhibitory relationship"))

perf_plot <- perf_data %>%
    ggplot(aes(x=1, y=value, fill=sign)) +
    geom_boxplot(outlier.size=0.25) +
    facet_wrap("score", scales="free_x", nrow=1, strip.position="bottom") +
    scale_fill_manual("", values=c("#7b4c87", "#82cfcc"),
                      guide=guide_legend(direction="horizontal")) +
    theme(legend.position="top", legend.title=element_blank(),
          axis.ticks.y=element_blank(), axis.text.y=element_blank(),
          axis.line.y=element_blank(), axis.title=element_blank(),
          strip.background=element_blank(), strip.placement="outside",
          panel.spacing=unit(0.5, "cm"),
          strip.text=element_text(vjust=1),
          legend.margin=margin(0,0,0,0)) +
    coord_flip()

save_plot("img/figures/panels/sign-perf.pdf", perf_plot, ncol=3)

reg_site_val_set <- read_tsv("data/psiteplus-reg-sites.tsv")
reg_site_val_set <- reg_site_val_set %>%
    filter(!(grepl("activity, induced", reg_site_val_set$ON_FUNCTION) &
             grepl("activity, inhibited", reg_site_val_set$ON_FUNCTION)) &
           (grepl("activity, induced", reg_site_val_set$ON_FUNCTION) |
            grepl("activity, inhibited", reg_site_val_set$ON_FUNCTION))) %>%
    select(GENE, residue, position, ON_FUNCTION)
reg_site_val_set$sign <- factor(sapply(1:nrow(reg_site_val_set),
                                       function(n){
                                           act <- reg_site_val_set$ON_FUNCTION[n]
                                           if (grepl("activity, induced", act)){
                                               return ("up")
                                           }else if (grepl("activity, inhibited", act)){
                                               return ("down")
                                           }else{
                                               return (NA)
                                           }
                                       }),
                                levels=c("up", "down"),
                                labels=c("activating site", "inhibitory site"))
names(reg_site_val_set)[1] <- "kinase"

reg_site_data <- read_tsv("out/reg-site-sign-feats.tsv") %>%
    inner_join(reg_site_val_set, by=c("kinase", "position",  "residue"))  %>%
    select(kinase, position, residue, log10_hotspot_pval_min, disopred_score,
           pos_in_domain, pos_perc, sign)  %>%
    gather(score, value, -kinase, -position, -residue, -sign)

reg_site_data$score <- factor(reg_site_data$score,
                              levels=c("log10_hotspot_pval_min",
                                       "disopred_score",
                                       "pos_in_domain",
                                       "pos_perc"),
                              labels=c("phosphorylation hotspot\n(-log10 p-value)",
                                       "DISOPRED score",
                                       "% position in kinase domain",
                                       "% position in protein"))

reg_site_perf_plot <- reg_site_data %>%
    ggplot(aes(x=1, y=value, fill=sign)) +
    geom_boxplot(outlier.size=0.25) +
    facet_wrap("score", scales="free_x", nrow=1, strip.position="bottom") +
    scale_fill_manual("", values=c("#7b4c87", "#82cfcc"),
                      guide=guide_legend(direction="horizontal")) +
    theme(legend.position="top", legend.title=element_blank(),
          axis.ticks.y=element_blank(), axis.text.y=element_blank(),
          axis.line.y=element_blank(), axis.title=element_blank(),
          strip.background=element_blank(), strip.placement="outside",
          panel.spacing=unit(0.5, "cm"),
          strip.text=element_text(vjust=1),
          legend.margin=margin(0,0,0,0)) +
    coord_flip()

save_plot("img/figures/panels/sign-regsite-perf.pdf", reg_site_perf_plot, ncol=3)

mid_align <- align_plots(coreg_score_plot, dcg_plot, exsub_final_plot, align="h")
reg_cor_grid <- plot_grid(reg_site_cor_plot, mid_align[[1]], ncol=1,
                          scale=0.9, rel_heights=c(1.0, 0.3))
dcg_align <- align_plots(exsub_site_plot, mid_align[[2]], align="v", axis="l")
dcg_grid <- plot_grid(dcg_align[[2]], mid_align[[3]], ncol=2,
                      rel_widths=c(1.0, 0.35), scale=0.9)
pssms_grid <- plot_grid(dcg_align[[1]], dcg_grid, ncol=1, scale=0.9,
                        labels=c("c", "d"))
ex_grid <- plot_grid(reg_cor_grid, pssms_grid, ncol=2, rel_widths=c(0.7, 1.0))

final_plot <- plot_grid(reg_site_perf_plot,
                        ex_grid,
                        perf_plot,
                        labels=c("a", "b", "e"), ncol=1,
                        rel_heights=c(0.35, 1.5, 0.35),
                        scale=c(0.95, 1.0, 0.95))

save_plot("img/figures/sign-figure.pdf", final_plot, ncol=4, nrow=4,
          base_height=2.0)
