library(dplyr)
library(readr)
library(tidyr)
library(VGAM)
library(igraph)
library(Bivariate.Pareto)
library(ggplot2)
library(copula)

tweak_deg_dist <- function(deg_dist){
    j <- which.max(c(sum(deg_dist$degree.out), sum(deg_dist$degree.in)))
    k <- which.min(c(sum(deg_dist$degree.out), sum(deg_dist$degree.in)))
    if (sum(deg_dist[,j]) %% 2){
        i <- sample.int(nrow(deg_dist), size=1)
        deg_dist[i, j] <- deg_dist[i, j] + 1
    }
    while (sum(deg_dist$degree.out) != sum(deg_dist$degree.in)){
        i <- sample.int(nrow(deg_dist), size=1)
        deg_dist[i, k] <- deg_dist[i, k] + 1
    }
    return (deg_dist)
}

true_pos <- read_tsv("data/direct-validation-set-omnipath.tsv",
                     col_names=FALSE)
names(true_pos) <- c("prot1", "prot2")

network <- read_tsv("out/final-predictors/merged-predictor.tsv") %>%
    filter(prot1 != prot2, !is.na(bart.pred.mean))

kinome <- as.character(read_tsv("data/human-kinome.txt", col_names=FALSE)$X1)
num_kin <- length(kinome)

true_pos$id <- paste(true_pos$prot1, true_pos$prot2,  sep="\t")
network$id <- paste(network$prot1, network$prot2, sep="\t")

true_pos_net <- filter(network, id %in% true_pos$id)
true_pos_ecdf <- ecdf(true_pos_net$bart.pred.mean)

net_unknown <- filter(network, !(id %in% true_pos$id))
net_unknown$prob.true <- true_pos_ecdf(net_unknown$bart.pred.mean)

n <- nrow(true_pos)
N_trials <- floor(nrow(network)/2^c(2, 4, 6, 8))
for (N in c(N_trials, nrow(network)/9){
    print(N)
    true_pos_sim <- c(sample(net_unknown$id, size=(N-n), prob=net_unknown$prob.true),
                      true_pos$id)
    true_pos_sim_df <- as_tibble(as.data.frame(true_pos_sim)) %>%
        separate(true_pos_sim, into=c("prot1", "prot2"))
    write_tsv(true_pos_sim_df, paste0("out/sim-valsets/direct-valset-", N, ".tsv"),
              col_names=FALSE)
}

true_pos_kins <- unique(c(true_pos$prot1, true_pos$prot2))

out_deg_dist <- as.data.frame(table(true_pos$prot1))
names(out_deg_dist) <- c("kinase", "degree")
out_deg_dist$kinase <- as.character(out_deg_dist$kinase)
## extra_kins <- setdiff(true_pos_kins, out_deg_dist$kinase)
## extra_kins_deg <- data.frame(kinase=extra_kins, degree=rep(0, length(extra_kins)))
## out_deg_dist <- rbind(out_deg_dist, extra_kins_deg)

in_deg_dist <- as.data.frame(table(true_pos$prot2))
names(in_deg_dist) <- c("kinase", "degree")
in_deg_dist$kinase <- as.character(in_deg_dist$kinase)
## extra_kins <- setdiff(true_pos_kins, in_deg_dist$kinase)
## extra_kins_deg <- data.frame(kinase=extra_kins, degree=rep(0, length(extra_kins)))
## in_deg_dist <- rbind(in_deg_dist, extra_kins_deg)

deg_dist <- merge(out_deg_dist, in_deg_dist, by="kinase",
                  suffixes=c(".out", ".in"), all=TRUE)
deg_dist <- subset(deg_dist, degree.out > 0 & degree.in > 0)

cor.test(deg_dist$degree.in, deg_dist$degree.out, method="spearman")

## Fit the degree distribution of the known true-positive network with
## a Pareto distribution.
out_deg_fit <- vglm(degree ~ 1, truncpareto(lower=0.1, upper=(num_kin-1)),
                                            data=out_deg_dist)
out_deg_shape <- Coef(out_deg_fit)[1]

in_deg_fit <- vglm(degree ~ 1, truncpareto(lower=0.1, upper=(num_kin-1)),
                                            data=in_deg_dist)
in_deg_shape <- Coef(in_deg_fit)[1]

deg_fit <- vglm(c(degree.out, degree.in) ~ 1, truncpareto(lower=0.1, upper=(num_kin-1)),
                data=deg_dist)
deg_shape <- Coef(deg_fit)[1]

for (i in 1:5){
    foo <- round(rtruncpareto(num_kin, lower=0.1, upper=(num_kin-1),
                              shape=out_deg_shape))
    print(c(min(foo), mean(foo), max(foo), sum(foo)))
    foo <- round(rtruncpareto(num_kin, lower=0.1, upper=(num_kin-1),
                              shape=in_deg_shape))
    print(c(min(foo), mean(foo), max(foo), sum(foo)))
    foo <- round(rtruncpareto(num_kin, lower=0.1, upper=(num_kin-1),
                              shape=deg_shape))
    print(c(min(foo), mean(foo), max(foo), sum(foo)))
    message()
}

## frank_pareto_theta <- 6
## bivar_deg_fit <- MLE.Frank.Pareto(c(deg_dist$degree.in, deg_dist$degree.out),
##                                   c(rep(1, length(deg_dist$degree.in)),
##                                     rep(0, length(deg_dist$degree.out))),
##                                   c(rep(0, length(deg_dist$degree.in)),
##                                     rep(1, length(deg_dist$degree.out))),
##                                   Theta = frank_pareto_theta)

## ## Do multiple iterations of this
## sim_deg_dist <- as.data.frame(round(Frank.Pareto(num_kin, Theta = frank_pareto_theta,
##                                                  Alpha1 = bivar_deg_fit$Alpha1[1],
##                                                  Alpha2 = bivar_deg_fit$Alpha2[1],
##                                                  Gamma1 = bivar_deg_fit$Gamma1[1],
##                                                  Gamma2 = bivar_deg_fit$Gamma2[1])))
## names(sim_deg_dist) <- c("degree.in", "degree.out")

dtruncpareto2 <- function(x, shape, log = FALSE){
    dtruncpareto(x=x, lower=1, upper=(num_kin-1), shape=shape, log=log)
}
ptruncpareto2 <- function(q, shape, lower.tail = TRUE, log.p = FALSE){
    ptruncpareto(q=q, lower=1, upper=(num_kin-1), shape=shape,
                 lower.tail=lower.tail, log.p=log.p)
}
qtruncpareto2 <- function(p, shape){
    qtruncpareto(p=p, lower=1, upper=(num_kin-1), shape=shape)
}
rtruncpareto2 <- function(n, shape){
    rtruncpareto(n=n, lower=1, upper=(num_kin-1), shape=shape)
}

bipareto_mvdc <- mvdc(gumbelCopula(dim=2), c("truncpareto2", "truncpareto2"),
                      paramMargins = list(list(shape=NA), list(shape=NA)),
                      marginsIdentical = FALSE)
bipareto_fit <- fitMvdc(as.matrix(deg_dist[2:3]+1e-15), bipareto_mvdc, start=c(1, 1, 10),
                        method="BFGS", optim.control=list(trace=TRUE, REPORT=2))
bipareto_fit_est <- bipareto_fit@estimate
deg_bipareto <- mvdc(gumbelCopula(param=bipareto_fit_est[3], dim=2),
                     c("truncpareto2", "truncpareto2"),
                     paramMargins = list(list(shape=bipareto_fit_est[1]),
                                         list(shape=bipareto_fit_est[2])),
                     marginsIdentical = FALSE)
## deg_bipareto <- mvdc(gumbelCopula(param=bipareto_fit_est[3], dim=2),
##                      c("truncpareto2", "truncpareto2"),
##                      paramMargins = list(list(shape=0.25),
##                                          list(shape=0.25)),
##                      marginsIdentical = FALSE)
sim_deg_dist <- as.data.frame(round(rMvdc(num_kin, deg_bipareto)))
names(sim_deg_dist) <- names(deg_dist)[2:3]

gamma_mvdc <- mvdc(gumbelCopula(dim=2), c("gamma", "gamma"),
                   paramMargins = list(list(shape=NA, scale=NA),
                                       list(shape=NA, scale=NA)),
                      marginsIdentical = FALSE)
gamma_fit <- fitMvdc(as.matrix(deg_dist[2:3]+1e-15), gamma_mvdc,
                        start=c(1, 1, 1, 1, 2),
                        method="BFGS", optim.control=list(trace=TRUE, REPORT=2))
gamma_fit_est <- gamma_fit@estimate
## deg_gamma <- mvdc(gumbelCopula(param=gamma_fit_est[5], dim=2),
##                      c("gamma", "gamma"),
##                   paramMargins = list(list(shape=gamma_fit_est[1],
##                                            scale=gamma_fit_est[2]),
##                                       list(shape=gamma_fit_est[3],
##                                            scale=gamma_fit_est[4])),
##                      marginsIdentical = FALSE)
deg_gamma <- mvdc(gumbelCopula(param=gamma_fit_est[5], dim=2),
                     c("gamma", "gamma"),
                  paramMargins = list(list(shape=15,
                                           scale=5),
                                      list(shape=15,
                                           scale=5)),
                     marginsIdentical = FALSE)
sim_deg_dist <- as.data.frame(round(rMvdc(num_kin, deg_gamma)))
names(sim_deg_dist) <- names(deg_dist)[2:3]

sim_deg_dist <- tweak_deg_dist(sim_deg_dist)
cor.test(sim_deg_dist$degree.in, sim_deg_dist$degree.out, method="spearman")
pdf("tmp/deg-dist.pdf")
ggplot(deg_dist, aes(x=degree.in, y=degree.out)) + geom_bin2d()
ggplot(sim_deg_dist, aes(x=degree.in, y=degree.out)) + geom_bin2d()
dev.off()

sim_deg_assign <- c()
sim_deg_kinase <- c()
kinome_sorted <- deg_dist$kinase[order(deg_dist$degree.out, decreasing=TRUE)]
kinome_sorted <- c(kinome_sorted, sample(setdiff(kinome, kinome_sorted)))
for (kin in kinome_sorted){
    if (kin %in% deg_dist$kinase){
        out_deg <- subset(deg_dist, kinase==kin)$degree.out
        in_deg <- subset(deg_dist, kinase==kin)$degree.in
    }else{
        out_deg <- 0
        in_deg <- 0
    }
    valid_is <- setdiff(1:nrow(sim_deg_dist), sim_deg_assign)
    sim_deg_is <- intersect(which(sim_deg_dist$degree.out >= out_deg &
                                  sim_deg_dist$degree.in >= in_deg),
                            valid_is)
    if (is.null(sim_deg_is) || length(sim_deg_is) == 0){
        sim_deg_i <- sample(valid_is, size=1)
    }else{
        sim_deg_i <- sample(sim_deg_is, size=1)
    }
    sim_deg_kinase <- c(sim_deg_kinase, kin)
    sim_deg_assign <- c(sim_deg_assign, sim_deg_i)
}

sim_deg_df <- data.frame(kinase=sim_deg_kinase,
                         degree.in=sim_deg_dist$degree.in[sim_deg_assign],
                         degree.out=sim_deg_dist$degree.out[sim_deg_assign])

true_pos_sim <- true_pos$id
for (kin in kinome_sorted){
    kin_out_deg <- subset(sim_deg_df, kinase==kin)$degree.out
    kin_out_true <- subset(network, prot1==kin & id %in% true_pos_sim)$id
    kin_out_unknown <- subset(network, prot1==kin & !(id %in% true_pos_sim))
    kin_out_unknown$prob.true <- true_pos_ecdf(kin_out_unknown$bart.pred.mean)
    i <- kin_out_deg - length(kin_out_true)
    if (i<=0){
        new_out_true <- NULL
    }else{
        num_pos <- length(which(kin_out_unknown$prob.true > 0))
        if (i > num_pos)
            i <- num_pos
        new_out_true <- sample(kin_out_unknown$id, size=i,
                               prob=kin_out_unknown$prob.true)
    }
    kin_in_deg <- subset(sim_deg_df, kinase==kin)$degree.in
    kin_in_true <- subset(network, prot2==kin & id %in% true_pos_sim)$id
    kin_in_unknown <- subset(network, prot2==kin & !(id %in% true_pos_sim))
    kin_in_unknown$prob.true <- true_pos_ecdf(kin_in_unknown$bart.pred.mean)
    i <- kin_in_deg - length(kin_in_true)
    if (i<=0){
        new_in_true <- NULL
    }else{
        num_pos <- length(which(kin_in_unknown$prob.true > 0))
        if (i > num_pos)
            i <- num_pos
        new_in_true <- sample(kin_in_unknown$id, size=i,
                              prob=kin_in_unknown$prob.true)
    }
    true_pos_sim <- c(true_pos_sim, new_out_true, new_in_true)
}

true_pos_sim_df <- as_tibble(as.data.frame(true_pos_sim)) %>%
    separate(true_pos_sim, into=c("prot1", "prot2"))

write_tsv(true_pos_sim_df,
          paste0("out/sim-valsets/direct-valset-estimated-gamma.tsv"),
          col_names=FALSE)

## Repeat for subset of true network with citations >100
