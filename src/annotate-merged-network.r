suppressMessages(library(dplyr))

argv <- commandArgs(TRUE)
if (!(length(argv) %in% c(4, 6, 7))){
    stop("USAGE: <script> MERGED_NETWORK ASSOC_VAL_SET DIRECT_VAL_SET SIGN_VAL_SET [SCORE_CUTOFF KEEP_KNOWN CELLINE]")
}

net_file <- argv[1]
assoc_val_set_file <- argv[2]
direct_val_set_file <- argv[3]
sign_val_set_file <- argv[4]
if (length(argv) %in% c(6, 7)){
    score_co <- as.numeric(argv[5])
    keep_known <- as.logical(argv[6])
    if (length(argv) == 7){
        celline <- argv[7]
    }else{
        celline <- NULL
    }
}else{
    score_co <- NULL
    keep_known <- NULL
    celline <- NULL
}

network <- read.delim(net_file, as.is=TRUE)
assoc_val_set <- read.table(assoc_val_set_file, as.is=TRUE)
direct_val_set <- read.table(direct_val_set_file, as.is=TRUE)
sign_val_set <- read.table(sign_val_set_file, as.is=TRUE)

rownames(network) <- paste(network$prot1, network$prot2, sep="-")
rownames(assoc_val_set) <- paste(assoc_val_set[,1], assoc_val_set[,2], sep="-")
rownames(direct_val_set) <- paste(direct_val_set[,1], direct_val_set[,2], sep="-")
rownames(sign_val_set) <- paste(sign_val_set[,1], sign_val_set[,2], sep="-")

network$assoc.is.known <- rownames(network) %in% rownames(assoc_val_set)
rev.direct <- paste(direct_val_set[,2], direct_val_set[,1], sep="-")
rev.direct <- setdiff(rev.direct, rownames(direct_val_set))
network$direct.is.true <- sapply(rownames(network),
                                   function(pair){
                                       if (pair %in% rownames(direct_val_set)){
                                           return(TRUE)
                                       }else if (pair %in% rev.direct){
                                           return(FALSE)
                                       }else{
                                           return(NA)
                                       }
                                   })

network$pred_sign <- sapply(network$prob.act.mean,
                              function(score){
                                  if (is.na(score))
                                      return (NA)
                                  if (score >= 0.5)
                                      return ("up")
                                  if (score < 0.5)
                                      return ("down")
                                  return ("ambig")
                              })
network$true_sign <- sapply(rownames(network),
                              function(pair){
                                  if (pair %in% rownames(sign_val_set)){
                                      if (sign_val_set[pair, 3] == "activates"){
                                          return("up")
                                      }else if (sign_val_set[pair, 3] == "inhibits"){
                                          return("down")
                                      }else{
                                          return("ambig")
                                      }
                                  }else{
                                      return("ambig")
                                  }
                              })
network$sign.is.correct <- sapply(rownames(network),
                                 function(pair){
                                     if (!(pair %in% rownames(sign_val_set)))
                                         return (NA)
                                     if (is.na(network[pair, "prob.act.mean"]) |
                                         is.na(network[pair, "prob.act.mean"]))
                                         return (NA)
                                     return ((network[pair, "prob.act.mean"] >= 0.5 &&
                                              sign_val_set[pair, 3] == "activates") ||
                                             (network[pair, "prob.act.mean"] < 0.5 &&
                                              sign_val_set[pair, 3] == "inhibits"))
                                 })

file_base <- strsplit(net_file, split="\\.")[[1]][1]

if (!is.null(celline)){
    suppressMessages(library(mclust))
    celline.rna <- read.delim("data/rna_celline.tsv", as.is=TRUE)
    if (!(celline %in% celline.rna$celline))
        stop("Invalid celline")
    celline.sub <- subset(celline.rna, celline==celline)
    celline.clust <- Mclust(celline.sub$value, G=2)
    up.reg.clust <- which.max(celline.clust$parameters[["mean"]])
    expr.prots <- celline.sub$kinase[celline.clust$classification==up.reg.clust]
    network <- subset(network, prot1 %in% expr.prots & prot2 %in% expr.prots)
    if (nrow(network) == 0)
        stop("No expressed proteins")
    file_base <- paste0(file_base, "-", celline)
}
known <- (network$assoc.is.known &
          !is.na(network$direct.is.true) &
          network$direct.is.true)
if (!is.null(score_co)){
    ## ## Compute a distinct cutoff for each kinase
    ## kinase.quant <- network %>%
    ##     group_by(prot1) %>%
    ##     summarise(co.quantile=quantile(bart.pred.mean, probs=score_co, na.rm=TRUE))
    ## rownames(kinase.quant) <- kinase.quant$prot1
    ## network.tmp <- network
    ## network.tmp$known <- known
    ## network.filt <- NULL
    ## for (kinase in kinase.quant$prot1){
    ##     network.kin <- subset(network.tmp, prot1==kinase)
    ##     kin.co <- kinase.quant[kinase, ]$co.quantile
    ##     if (keep_known){
    ##         network.kin$missing <- network.kin$known & (is.na(network.kin$bart.pred.mean) |
    ##                                                     network.kin$bart.pred.mean < kin.co)
    ##         network.kin$known <- NULL
    ##         network.kin.filt <- network.kin[which(network.kin$bart.pred.mean >= kin.co | network.kin$missing),]
    ##     }else{
    ##         network.kin$known <- NULL
    ##         network.kin.filt <- network.kin[which(network.kin$bart.pred.mean >= kin.co),]
    ##     }
    ##     if (is.null(network.filt)){
    ##         network.filt <- network.kin.filt
    ##     }else{
    ##         network.filt <- rbind(network.filt, network.kin.filt)
    ##     }
    ## }
    ## network <- network.filt

    ## co.p <- unlist(apply(network[,c("bart.pred.mean", "bart.pred.sd")], 1,
    ##                      function(row){
    ##                          pred.mean <- row[1]
    ##                          pred.sd <- row[2]
    ##                          return(pnorm(score_co, mean=pred.mean, sd=pred.sd))
    ##                      }))
    p.cutoff <- 0.05
    network$missing <- known & (is.na(network$bart.pred.mean) |
                                network$bart.pred.mean < score_co)
    ## network$missing <- known & (is.na(network$bart.pred.mean) |
    ##                             co.p>p.cutoff)
    if (keep_known){
        network <- subset(network, bart.pred.mean > score_co | missing)
        ## network <- network[which(co.p<p.cutoff | network$missing),]
    }else{
        network <- subset(network, bart.pred.mean > score_co)
        ## network <- network[which(co.p<p.cutoff),]
    }
    file_base <- paste0(file_base, "-filt", 100*score_co)
}else{
    network$missing <- known & is.na(network$bart.pred.mean)
}
network$edge.status <- sapply(1:nrow(network),
                              function(n){
                                  if (network$missing[n]){
                                      return("missing")
                                  }else if (network$assoc.is.known[n] &&
                                            !is.na(network$direct.is.true[n]) &&
                                            network$direct.is.true[n] &&
                                            !network$missing[n]){
                                      return("correct")
                                  }else{
                                      return("predicted")
                                  }
                              })

network$draw.arrow <- !is.na(network$bart.pred.mean) | network$missing
network$arrow.type <- sapply(1:nrow(network),
                               function(n){
                                   if (!network$draw.arrow[n]){
                                       return("none")
                                   }else{
                                       ## if (network$missing[n])
                                       ##     return(network$true_sign[n])
                                       return(network$pred_sign[n])
                                   }
                               })
network$arrow.status <- sapply(1:nrow(network),
                               function(n){
                                   ## if (network$missing[n]){
                                   ##     return("missing")
                                   ## }else if (is.na(network$direct.is.true[n])){
                                   ##     return("unknown")
                                   ## }else if (network$direct.is.true[n]){
                                   ##     return("good")
                                   ## }else {
                                   ##     return("bad")
                                   ## }
                                   if (network$true_sign[n] == "ambig" ||
                                       is.na(network$sign.is.correct[n])){
                                       return("ambig")
                                   }else if (network$sign.is.correct[n]){
                                       return("good")
                                   }else{
                                       return("bad")
                                   }
                               })

write.table(network, paste0(file_base, "-annot.tsv"), quote=FALSE,
            sep="\t", row.names=FALSE, col.names=TRUE)
