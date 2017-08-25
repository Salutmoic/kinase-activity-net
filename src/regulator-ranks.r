suppressPackageStartupMessages(library(data.table))

suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))

argv <- commandArgs(TRUE)
if (length(argv) != 2){
    stop("USAGE: <script> MERGED_NETWORK DIRECT_VAL_SET")
}

net_file <- argv[1]
true_intxns_file <- argv[2]

network <- read.delim(net_file, as.is=TRUE)
names(network) <- c("prot1", "prot2", "lo.kin", "bart.pred.mean")
true_intxns_tbl <- read.table(true_intxns_file, as.is=TRUE)

network$id <- paste(network$prot1, network$prot2, sep="-")

possible_true_intxns <- paste(true_intxns_tbl[, 1], true_intxns_tbl[, 2],
                              sep = "-")
possible_true_intxns <- intersect(possible_true_intxns, network$id)

citations <- read.delim("out/kinase-citations.tsv", as.is=TRUE)
citations <- subset(citations, count.type=="num.subs")


ranks_by_kin <- network %>%
    group_by(lo.kin) %>%
    arrange(desc(bart.pred.mean), .by_group=TRUE) %>%
    mutate(rank=rank(-bart.pred.mean, ties.method="first"),
           is.known=(id %in% possible_true_intxns))

regulator_ranks <- ranks_by_kin %>%
    filter(is.known)

min_rank <- regulator_ranks %>%
    summarise(min.rank=min(rank))
min_rank <- merge(min_rank, citations, by.x=c("lo.kin"), by.y=c("kinase"))

wkins <- c()
pvals <- c()
for (kin in unique(regulator_ranks$lo.kin)){
    kin.known <- subset(ranks_by_kin, lo.kin==kin & is.known)
    if (nrow(kin.known) < 5)
        next
    kin.unknown <- subset(ranks_by_kin, lo.kin==kin & !is.known)
    w <- wilcox.test(kin.known$bart.pred.mean, kin.unknown$bart.pred.mean,
                     alternative="greater")
    wkins <- c(wkins, kin)
    pvals <- c(pvals, w$p.value)
}
pred_wilcox <- data.frame(kinase=wkins, p=pvals, p.adj=p.adjust(pvals, method="fdr"))
pred_wilcox <- merge(pred_wilcox, citations, by.x=c("kinase"), by.y=c("kinase"))

file_base <- strsplit(net_file, split="\\.")[[1]][1]
img_file_base <- basename(file_base)
img_file <- paste0("img/", img_file_base, ".pdf")

pdf(img_file, width=28, height=7)
ggplot(regulator_ranks, aes(reorder(lo.kin, rank, FUN=median), rank)) +
    geom_boxplot() +
    geom_hline(yintercept=10, linetype=2, colour="grey") +
    scale_y_log10() +
    theme_bw(base_size=12) +
    labs(x="kinase", y="regulator prediction rank") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggplot(min_rank, aes(reorder(lo.kin, min.rank), min.rank)) +
    geom_point(aes(colour=log10(num.cites.filt)), size=3) +
    scale_colour_distiller("log10(citations)", palette="Spectral") +
    geom_hline(yintercept=10, linetype=2, colour="grey") +
    scale_y_log10() +
    theme_bw(base_size=20) +
    labs(x="kinase", y="best regulator prediction rank") +
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=9))
ggplot(pred_wilcox, aes(reorder(kinase, p.adj), p.adj)) +
    geom_point(aes(colour=log10(num.cites.filt)), size=3) +
    scale_colour_distiller("log10(citations)", palette="Spectral") +
    geom_hline(yintercept=0.05, linetype=2, colour="grey") +
    scale_y_log10() +
    theme_bw(base_size=20) +
    labs(x="kinase", y="Wilcox test p-value") +
    theme(axis.text.x = element_text(angle=45, hjust=1, size=15))
dev.off()
