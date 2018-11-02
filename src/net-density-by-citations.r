suppressMessages(library(ROCR))
suppressMessages(library(dplyr))
suppressMessages(library(readr))
suppressMessages(library(ggplot2))
source("src/rocr-ggplot2.r")
suppressMessages(library(cowplot))

argv <- commandArgs(TRUE)
if (length(argv) != 2){
    stop("USAGE: <script> NET_FILE VAL_SET")
}

net_file <- argv[1]
val_set_file <- argv[2]
net <- read_tsv(net_file) %>%
    mutate(id = paste(prot1, prot2))
val_set <- read_tsv(val_set_file, col_names=FALSE) %>%
    mutate(id = paste(X1, X2))

citations <- read_tsv("out/kinase-citations.tsv") %>%
    filter(count.type == "num.subs")

prs_full <- NULL
dens_full <- NULL
for (n in c(0, 10*2^c(0, 2, 4, 6, 7))){
    print(n)
    kins <- filter(citations, num.cites.filt >= n)$kinase
    net_sub <- filter(net, prot1 %in% kins, prot2 %in% kins,
                      prot1 != prot2)
    ## net_sub <- net_sub[complete.cases(net_sub),]
    net_n <- nrow(net_sub)
    true_intxns <- which(net_sub$id %in% val_set$id)
    unknown_intxns <- which(!(net_sub$id %in% val_set$id))
    intxns <- c(true_intxns, unknown_intxns)
    bart_score <- net_sub[intxns, "bart.pred.mean"]
    labels <- c(rep(TRUE, length(true_intxns)),
                rep(FALSE, length(unknown_intxns)))
    pred <- prediction(bart_score, labels)
    perf_pr <- performance(pred, measure = "prec", x.measure = "rec")
    if (is.null(prs_full)){
        prs_full <- data.frame(x = perf_pr@x.values[[1]],
                               y = perf_pr@y.values[[1]],
                               feature = rep(n, length(perf_pr@x.values[[1]])))
        dens_full <- data.frame(alpha = perf_pr@alpha.values[[1]],
                                y = perf_pr@y.values[[1]],
                                feature = rep(n, length(perf_pr@y.values[[1]])))
        dens_full$x <- sapply(dens_full$alpha,
                                    function(alpha){
                                        return(length(which(net_sub$bart.pred.mean >= alpha))/net_n)
                                    })

    }else{
        new_pr <- data.frame(x = perf_pr@x.values[[1]],
                             y = perf_pr@y.values[[1]],
                             feature = rep(n, length(perf_pr@x.values[[1]])))
        new_dens <- data.frame(alpha = perf_pr@alpha.values[[1]],
                               y = perf_pr@y.values[[1]],
                               feature = rep(n, length(perf_pr@y.values[[1]])))
        new_dens$x <- sapply(new_dens$alpha,
                                    function(alpha){
                                        return(length(which(net_sub$bart.pred.mean >= alpha))/net_n)
                                    })
        prs_full <- rbind(prs_full, new_pr)
        dens_full <- rbind(dens_full, new_dens)
    }
}

prs_full$feature <- as.factor(prs_full$feature)
prs_full$y[which(is.nan(prs_full$y))] <- NA
prs_full$y[which(is.infinite(prs_full$y))] <- NA
prs_full$x[which(is.nan(prs_full$x))] <- NA
prs_full$x[which(is.infinite(prs_full$x))] <- NA

dens_full$feature <- as.factor(dens_full$feature)
dens_full$y[which(is.nan(dens_full$y))] <- NA
dens_full$y[which(is.infinite(dens_full$y))] <- NA
dens_full <- dens_full[,c("x", "y", "feature")]

pdf("img/net-density-by-citations.pdf")
prec_full_plot <- rocr2d.ggplot2(prs_full, "recall", "precision", 1, NULL) +
    scale_colour_discrete("min. # citations")
print(prec_full_plot)
density_plot <- rocr2d.ggplot2(dens_full, "network density",
                               "fraction of\npredicted edges known", 1, NULL) +
    scale_colour_discrete("min. # citations")
print(density_plot)
dev.off()

data_basename <- strsplit(basename(net_file), split = "\\.")[[1]][1]
rds_basename <- paste0("img/rds/", data_basename)

saveRDS(prec_full_plot, file=paste0(data_basename, "-prec-full-by-cites.rds"))
saveRDS(density_plot, file=paste0(rds_basename, "-net-density-by-cites.rds"))
