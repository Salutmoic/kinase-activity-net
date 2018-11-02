## "Borrowed" from the ROCR source code

prepare_perf <- function(perf){
    ## for infinite cutoff, assign maximal finite cutoff + mean difference
    ## between adjacent cutoff pairs
    if (length(perf@alpha.values)!=0) perf@alpha.values <-
                                          lapply(perf@alpha.values,
                                                 function(x) { isfin <- is.finite(x);
                                                     x[is.infinite(x)] <-
                                                         (max(x[isfin]) +
                                                          mean(abs(x[isfin][-1] -
                                                                   x[isfin][-length(x[isfin])])));
                                                     x } )
    ## remove samples with x or y not finite
    for (i in 1:length(perf@x.values)) {
        ind.bool <- (is.finite(perf@x.values[[i]]) &
                     is.finite(perf@y.values[[i]]))

        if (length(perf@alpha.values)>0)
            perf@alpha.values[[i]] <- perf@alpha.values[[i]][ind.bool]

        perf@x.values[[i]] <- perf@x.values[[i]][ind.bool]
        perf@y.values[[i]] <- perf@y.values[[i]][ind.bool]
    }
    return(perf)
}

perf_vert_avg <- function(perf){
    perf <- prepare_perf(perf)
    perf.avg <- perf
    x.values <- seq(min(unlist(perf@x.values)), max(unlist(perf@x.values)),
                    length=max( sapply(perf@x.values, length)))
    if (length(x.values) < 2 || length(unique(x.values)) < 2 || all(is.na(x.values)))
        return(NA)
    ## print(x.values)
    for (i in 1:length(perf@y.values)) {
        if ((length(perf.avg@y.values[[i]]) < 2 ||
             length(unique(perf.avg@y.values[[i]])) < 2 ||
             all(is.na(perf.avg@y.values[[i]])))){
            perf.avg@y.values[[i]] <- NA
            next
        }
        perf.avg@y.values[[i]] <-
            approxfun(perf@x.values[[i]], perf@y.values[[i]],
                      ties=mean, rule=2)(x.values)
    }
    perf.avg@y.values <- list(rowMeans( data.frame( perf.avg@y.values )))
    perf.avg@x.values <- list(x.values)
    perf.avg@alpha.values <- list()
    return(perf.avg)
}

perf_horiz_avg <- function(perf){
    perf <- prepare_perf(perf)
    perf.avg <- perf
    y.values <- seq(min(unlist(perf@y.values)), max(unlist(perf@y.values)),
                    length=max( sapply(perf@y.values, length)))
    if (length(y.values) < 2 || length(unique(y.values)) == 1 || all(is.na(y.values)))
        return(NA)
    for (i in 1:length(perf@x.values)) {
        perf.avg@x.values[[i]] <- approxfun(perf@y.values[[i]],
                                            perf@x.values[[i]],
                                            ties=mean, rule=2)(y.values)
    }
    perf.avg@x.values <- list(rowMeans( data.frame( perf.avg@x.values )))
    perf.avg@y.values <- list(y.values)
    perf.avg@alpha.values <- list()
    return(perf.avg)
}

perf_thresh_avg <- function(perf){
    perf <- prepare_perf(perf)
    perf.sampled <- perf
    alpha.values <- rev(seq(min(unlist(perf@alpha.values)),
                            max(unlist(perf@alpha.values)),
                            length=max( sapply(perf@alpha.values, length))))
    if (length(alpha.values) < 2 || length(unique(alpha.values)) == 1 || all(is.na(alpha.values)))
        return(NA)
    for (i in 1:length(perf.sampled@y.values)) {
        perf.sampled@x.values[[i]] <-
            approxfun(perf@alpha.values[[i]],perf@x.values[[i]],
                      rule=2, ties=mean)(alpha.values)
        perf.sampled@y.values[[i]] <-
            approxfun(perf@alpha.values[[i]], perf@y.values[[i]],
                      rule=2, ties=mean)(alpha.values)
    }

    ## compute average curve
    perf.avg <- perf.sampled
    perf.avg@x.values <- list( rowMeans( data.frame( perf.avg@x.values)))
    perf.avg@y.values <- list(rowMeans( data.frame( perf.avg@y.values)))
    perf.avg@alpha.values <- list( alpha.values )
    return(perf.avg)
}
