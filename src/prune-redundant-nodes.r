library(dplyr)

is_mutual_overlap <- function(overlap, max_overlap){
    if (any(is.nan(overlap))){
        return (FALSE)
    }
    return (all(overlap > max_overlap))
}

is_partial_overlap <- function(overlap, i, max_overlap){
    if (any(is.nan(overlap))){
        return (FALSE)
    }
    return (overlap[i] > max_overlap)
}

choose_redundant <- function(kin1, kin2, kin1_mean, kin2_mean, sub_overlap,
                             reg_overlap, max_overlap){
    if (((is_mutual_overlap(sub_overlap, max_overlap) &&
          is_mutual_overlap(reg_overlap, max_overlap)) ||
         (is_partial_overlap(sub_overlap, 1, max_overlap) &&
          is_partial_overlap(reg_overlap, 2, max_overlap)) ||
         (is_partial_overlap(sub_overlap, 2, max_overlap) &&
          is_partial_overlap(reg_overlap, 1, max_overlap)))){
        if (kin1_mean > kin2_mean){
            return (kin2)
        }else if (kin1_mean < kin2_mean){
            return (kin1)
        }
    }else if ((is_partial_overlap(sub_overlap, 1, max_overlap) &&
               is_partial_overlap(reg_overlap, 1, max_overlap))){
        return (kin1)
    }else if ((is_partial_overlap(sub_overlap, 2, max_overlap) &&
               is_partial_overlap(reg_overlap, 2, max_overlap))){
        return (kin2)
    }
    return (NULL)
}


argv <- commandArgs(TRUE)
if (length(argv) != 2){
    stop("USAGE: <script> NETWORK MAX_OVERLAP")
}

## Remove nodes from a (pre-filtered) network if they share too many
## edges with another node in the network.

net_file <- argv[1]
max_overlap <- as.numeric(argv[2])

network <- read.delim(net_file, as.is=TRUE)

to_prune <- rep(NA, 2)
net_sub <- network
while(length(to_prune)>0){
    all_kins <- sample(c(unique(net_sub$prot1), unique(net_sub$prot2)))
    to_prune <- c()
    for (kin1 in all_kins){
        if (kin1 %in% to_prune)
            next
        kin1_out_net <- subset(net_sub, prot1==kin1)
        kin1_subs <- unique(kin1_out_net$prot2)
        kin1_in_net <- subset(net_sub, prot2==kin1)
        kin1_regs <- unique(kin1_in_net$prot1)
        kin1_mean <- mean(kin1_out_net$bart.pred.mean, na.rm=TRUE)
        kin1_subs_n <- length(kin1_subs)
        kin1_regs_n <- length(kin1_regs)
        for (kin2 in all_kins){
            if (kin1==kin2 || kin2 %in% to_prune)
                next
            kin2_out_net <- subset(net_sub, prot1==kin2)
            kin2_subs <- unique(kin2_out_net$prot2)
            kin2_in_net <- subset(net_sub, prot2==kin2)
            kin2_regs <- unique(kin2_in_net$prot1)
            kin2_mean <- mean(kin2_out_net$bart.pred.mean, na.rm=TRUE)
            kin2_subs_n <- length(kin2_subs)
            kin2_regs_n <- length(kin2_regs)
            sub_intsct <- length(intersect(kin1_subs, kin2_subs))
            reg_intsct <- length(intersect(kin1_regs, kin2_regs))
            sub_overlap <- c(sub_intsct/kin1_subs_n, sub_intsct/kin2_subs_n)
            reg_overlap <- c(reg_intsct/kin1_regs_n, reg_intsct/kin2_regs_n)
            redundant <- choose_redundant(kin1, kin2, kin1_mean, kin2_mean,
                                          sub_overlap, reg_overlap, max_overlap)
            if (is.null(redundant))
                next
            if (!(redundant %in% to_prune))
                to_prune <- c(to_prune, redundant)
            if (redundant == kin1)
                break
       }
    }
    print(to_prune)
    net_sub <- subset(net_sub, !(prot1 %in% to_prune) &
                               !(prot2 %in% to_prune))
}

print(c(length(c(unique(network$prot1), unique(network$prot2))),
        length(c(unique(net_sub$prot1), unique(net_sub$prot2)))))

file_base <- strsplit(net_file, split="\\.")[[1]][1]

write.table(net_sub, paste0(file_base, "-pruned", (max_overlap*100), ".tsv"),
            quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)
