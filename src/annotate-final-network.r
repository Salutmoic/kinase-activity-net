library(dplyr)

collapse_table <- function(network, keep_known){
    network <- network[order(network$direct.pred.mean),]
    keep <- !is.na(network$direct.pred.mean) | !duplicated(network[,c("assoc.pred.mean")])
    if (keep_known){
        keep <- keep | network$missing
    }
    return(network[keep,])
}

argv <- commandArgs(TRUE)
if (!(length(argv) %in% c(4, 8))){
    stop("USAGE: <script> NETWORK ASSOC_VAL_SET DIRECT_VAL_SET SIGN_VAL_SET [ASSOC_CUTOFF DIRECT_CUTOFF COLLAPSE KEEP_KNOWN]")
}

net_file <- argv[1]
assoc_val_set_file <- argv[2]
direct_val_set_file <- argv[3]
sign_val_set_file <- argv[4]
if (length(argv) == 8){
    assoc_co <- as.numeric(argv[5])
    direct_co <- as.numeric(argv[6])
    collapse <- as.logical(argv[7])
    keep_known <- as.logical(argv[8])
}else{
    assoc_co <- NULL
    direct_co <- NULL
    collapse <- NULL
    keep_known <- NULL
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

known <- (network$assoc.is.known &
          !is.na(network$direct.is.true) &
          network$direct.is.true)

file_base <- strsplit(net_file, split="\\.")[[1]][1]

if (!is.null(assoc_co)){
    network$missing <- known & ((is.na(network$assoc.pred.mean) |
                                network$assoc.pred.mean < assoc_co) |
                                (is.na(network$direct.pred.mean) |
                                network$direct.pred.mean <= direct_co))
    if (keep_known){
        network <- subset(network, (assoc.pred.mean > assoc_co &
                                    (is.na(direct.pred.mean) | direct.pred.mean > direct_co)) |
                                   missing)
    }else{
        network <- subset(network, assoc.pred.mean > assoc_co &
                                   (is.na(direct.pred.mean) | direct.pred.mean > direct_co))
    }
    file_base <- paste0(file_base, "-filt")
}else{
    network$missing <- known & is.na(network$direct.pred.mean)
}

network$edge.status <- sapply(1:nrow(network),
                              function(n){
                                  if (network$missing[n]){
                                      return("missing")
                                  }else if (network$assoc.is.known[n] &&
                                            !is.na(network$direct.is.true[n]) &&
                                            network$direct.is.true[n] &&
                                            !network$missing[n]){
                                      return("good")
                                  }else{
                                      return("novel")
                                  }
                              })

network$draw.arrow <- !is.na(network$direct.pred.mean) | network$missing
network$arrow.type <- sapply(1:nrow(network),
                               function(n){
                                   if (!network$draw.arrow[n]){
                                       return("none")
                                   }else{
                                       if (network$missing[n])
                                           return(network$true_sign[n])
                                       return(network$pred_sign[n])
                                   }
                               })
network$arrow.status <- sapply(1:nrow(network),
                               function(n){
                                   if (network$missing[n]){
                                       return("missing")
                                   }else if (is.na(network$direct.is.true[n])){
                                       return("unknown")
                                   }else if (network$direct.is.true[n]){
                                       return("good")
                                   }else {
                                       return("bad")
                                   }
                               })

if (!is.null(collapse) && collapse){
    network <- collapse_table(network, keep_known)
    network <- network[order(network$prot1, network$prot2),]
    file_base <- paste0(file_base, "-coll")
}

write.table(network, paste0(file_base, "-annot.tsv"), quote=FALSE,
            sep="\t", row.names=FALSE, col.names=TRUE)
