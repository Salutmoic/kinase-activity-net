library(ggplot2)

get.net.perf <- function(net, cutoffs, alpha){
    assoc.cos <- c()
    direct.cos <- c()
    f.vals <- c()
    precs <- c()
    recs <- c()
    p <- length(which(net$assoc.is.true))
    for (assoc.co in cutoffs){
        for (direct.co in cutoffs){
            net.sub <- subset(net, assoc.pred.mean > assoc.co &
                                   direct.pred.mean > direct.co)
            tp <- length(which(net.sub$assoc.is.true))
            fp <- length(which(!net.sub$assoc.is.true))
            fn <- p - tp
            prec <- tp/(tp+fp)
            rec <- tp/(tp+fn)
            f <- 1.0/(alpha*(1.0/prec)+(1.0-alpha)*(1.0/rec))
            assoc.cos <- c(assoc.cos, assoc.co)
            direct.cos <- c(direct.cos, direct.co)
            f.vals <- c(f.vals, f)
            precs <- c(precs, prec)
            recs <- c(recs, rec)
        }
    }
    results <- data.frame(assoc.co=assoc.cos, direct.co=direct.cos, f.val=f.vals,
                          precision=precs, recall=recs)
    return(results)
}

net.full <- read.delim("out/full-annotated-predictor.tsv", as.is=TRUE)
net.pruned <- read.delim("out/full-annotated-predictor-pruned.tsv", as.is=TRUE)

cutoffs <- seq(0.0, 1.0, 0.025)

net.full.perf.prec <- get.net.perf(net.full, cutoffs, 0.75)
net.full.perf.prec$balance <- rep("favor precision (0.75)", nrow(net.full.perf.prec))
net.full.perf.prec2 <- get.net.perf(net.full, cutoffs, 0.9)
net.full.perf.prec2$balance <- rep("favor precision (0.9)", nrow(net.full.perf.prec2))
net.full.perf.bal <- get.net.perf(net.full, cutoffs, 0.5)
net.full.perf.bal$balance <- rep("balanced (0.5)", nrow(net.full.perf.bal))
net.full.perf <- rbind(net.full.perf.prec, net.full.perf.prec2)
net.full.perf <- rbind(net.full.perf, net.full.perf.bal)
net.full.perf$balance <- factor(net.full.perf$balance,
                                levels=c("balanced (0.5)",
                                         "favor precision (0.75)",
                                         "favor precision (0.9)"))

message("Full network best case:")
max.f <- which.max(net.full.perf.prec2$f.val)
baseline <- net.full.perf.prec2[max.f, "assoc.co"]
direct.co.0 <- which(net.full.perf.prec2$assoc.co == baseline & net.full.perf.prec2$direct.co==0.0)
print(net.full.perf.prec2[direct.co.0,])
print(net.full.perf.prec2[max.f,])

net.pruned.perf.prec <- get.net.perf(net.pruned, cutoffs, 0.75)
net.pruned.perf.prec$balance <- rep("favor precision (0.75)", nrow(net.pruned.perf.prec))
net.pruned.perf.prec2 <- get.net.perf(net.pruned, cutoffs, 0.9)
net.pruned.perf.prec2$balance <- rep("favor precision (0.9)", nrow(net.pruned.perf.prec2))
net.pruned.perf.bal <- get.net.perf(net.pruned, cutoffs, 0.5)
net.pruned.perf.bal$balance <- rep("balanced (0.5)", nrow(net.pruned.perf.bal))
net.pruned.perf <- rbind(net.pruned.perf.prec, net.pruned.perf.prec2)
net.pruned.perf <- rbind(net.pruned.perf, net.pruned.perf.bal)
net.pruned.perf$balance <- factor(net.pruned.perf$balance,
                                levels=c("balanced (0.5)",
                                         "favor precision (0.75)",
                                         "favor precision (0.9)"))

message("Pruned network best case:")
max.f <- which.max(net.pruned.perf.prec2$f.val)
baseline <- net.pruned.perf.prec2[max.f, "assoc.co"]
direct.co.0 <- which(net.pruned.perf.prec2$assoc.co == baseline & net.pruned.perf.prec2$direct.co==0.0)
print(net.pruned.perf.prec2[direct.co.0,])
print(net.pruned.perf.prec2[max.f,])

pdf("img/network-performance.pdf", width=21)
ggplot(net.full.perf, aes(assoc.co, direct.co)) +
    geom_raster(aes(fill=f.val)) +
    facet_grid(~balance)
ggplot(net.pruned.perf, aes(assoc.co, direct.co)) +
    geom_raster(aes(fill=f.val)) +
    facet_grid(~balance)
dev.off()
