library(ggplot2)
source("src/rocr-helpers.r")
library(cowplot)

calc.vert.ci <- function(perf, ci.level){
    perf <- prepare_perf(perf)
    perf.avg <- perf
    x.values <- seq(min(unlist(perf@x.values)), max(unlist(perf@x.values)),
                    length=max( sapply(perf@x.values, length)))
    if (length(x.values) < 2 || length(unique(x.values)) < 2 || all(is.na(x.values)))
        return(NA)
    for (i in 1:length(perf@y.values)) {
        perf.avg@y.values[[i]] <-
            approxfun(perf@x.values[[i]], perf@y.values[[i]],
                      ties=mean, rule=2)(x.values)
    }
    y.ci.values <- list(apply(data.frame(perf.avg@y.values), 1,
                              function(row){
                                  row.sd <- sd(row, na.rm=TRUE)
                                  n.row <- length(row)
                                  err <- qt(0.025, df=(n.row-1))*row.sd/sqrt(n.row)
                                  return(err)
                              }))
    return(y.ci.values)
}

build.cvperf2d.df <- function(preds, x_measure, y_measure, feat_name, n, avg){
    perf <- performance(preds, measure=y_measure, x.measure=x_measure)
    if (avg){
        perf_avg <- perf_vert_avg(perf)
        y_ci_values <- calc.vert.ci(perf)
        if (is.na(perf_avg) || is.na(y_ci_values))
            return(NA)
        perf <- perf_avg
    }
    if (n > 1 && avg){
        perf_dat <- data.frame(x=perf@x.values[[1]],
                               y=perf@y.values[[1]],
                               y.ci=y_ci_values[[1]],
                               feature=rep(feat_name,
                                           length(perf@x.values[[1]])))
    }else if (n > 1 && !avg){
        perf_dat <- data.frame(x=unlist(perf@x.values),
                               y=unlist(perf@y.values),
                               feature=rep(feat_name,
                                           length(perf@x.values)))
    }else{
        perf_dat <- data.frame(x=perf@x.values[[1]],
                               y=perf@y.values[[1]],
                               feature=rep(feat_name,
                                           length(perf@x.values[[1]])))
    }
    return(perf_dat)
}

build.cvperf.sum.df <- function(preds, y_measure, feat_name, n, avg){
    perf <- performance(preds, measure=y_measure)
    if (avg){
        perf_avg <- perf_vert_avg(perf)
        y_ci_values <- calc.vert.ci(perf)
        if (is.na(perf_avg) || is.na(y_ci_values))
            return(NA)
        perf <- perf_avg
    }
    if (n > 1 && avg){
        perf_dat <- data.frame(y=perf@y.values[[1]],
                               y.ci=y_ci_values[[1]],
                               feature=rep(feat_name,
                                           length(perf@y.values[[1]])))
    }else if (n > 1 && !avg){
        perf_dat <- data.frame(y=unlist(perf@y.values),
                               feature=rep(feat_name,
                                           length(perf@y.values)))
    }else{
        perf_dat <- data.frame(y=perf@y.values[[1]],
                               feature=rep(feat_name,
                                           length(perf@y.values[[1]])))
    }
    return(perf_dat)
}

build.cvperf1d.df <- function(preds, y_measure, feat_name, n, avg){
    perf <- performance(preds, measure=y_measure)
    if (avg){
        perf_avg <- perf_vert_avg(perf)
        y_ci_values <- calc.vert.ci(perf)
        if (is.na(perf_avg) || is.na(y_ci_values))
            return(NA)
       perf <- perf_avg
    }
    if (n > 1 && avg){
        perf_dat <- data.frame(x=perf@x.values[[1]],
                               y=perf@y.values[[1]],
                               y.ci=y_ci_values[[1]],
                               feature=rep(feat_name,
                                           length(perf@y.values[[1]])))
    }else if (n > 1 && !avg){
        perf_dat <- data.frame(x=unlist(perf@x.values),
                               y=unlist(perf@y.values),
                               feature=rep(feat_name,
                                           length(perf@y.values)))
    }else{
        perf_dat <- data.frame(x=perf@x.values[[1]],
                               y=perf@y.values[[1]],
                               feature=rep(feat_name,
                                           length(perf@y.values[[1]])))
    }
    return(perf_dat)
}

build.perf2d.df <- function(preds, x_measure, y_measure, feat_names, n, avg){
    perf_dat <- NULL
    for (i in 1:length(preds)){
        if (is.na(preds[[i]]) || all(is.na(preds[[i]])))
            next
        perf <- performance(preds[[i]], measure=y_measure, x.measure=x_measure)
        if (is.na(perf))
            next
        if (avg){
            perf_avg <- perf_vert_avg(perf)
            y_ci_values <- calc.vert.ci(perf)
            if (is.na(perf_avg) || is.na(y_ci_values))
                next
            perf <- perf_avg
        }
        if (n > 1 && avg){
            new_dat <- data.frame(x=perf@x.values[[1]],
                                  y=perf@y.values[[1]],
                                  y.ci=y_ci_values[[1]],
                                  feature=rep(feat_names[i],
                                              length(perf@x.values[[1]])))
        }else if (n > 1 && !avg){
            new_dat <- data.frame(x=unlist(perf@x.values),
                                  y=unlist(perf@y.values),
                                  feature=rep(feat_names[i],
                                              length(perf@x.values)))
        }else{
            new_dat <- data.frame(x=perf@x.values[[1]],
                                  y=perf@y.values[[1]],
                                  feature=rep(feat_names[i],
                                              length(perf@x.values[[1]])))
        }
        if (is.null(perf_dat)){
            perf_dat <- new_dat
        }else{
            perf_dat <- rbind(perf_dat, new_dat)
        }
    }
    return(perf_dat)
}

build.perf1d.df <- function(preds, y_measure, feat_names, n, avg){
    perf_dat <- NULL
    for (i in 1:length(preds)){
        if (is.na(preds[[i]]))
            next
        perf <- performance(preds[[i]], measure=y_measure)
        if (avg){
            perf_avg <- perf_vert_avg(perf)
            y_ci_values <- calc.vert.ci(perf)
            if (is.na(perf_avg) || is.na(y_ci_values))
                next
            perf <- perf_avg
        }
        if (n > 1 && avg){
            new_dat <- data.frame(y=perf@y.values[[1]],
                                  y.ci=y_ci_values[[1]],
                                  feature=rep(feat_names[i],
                                              length(perf@y.values[[1]])))
        }else if (n > 1 && !avg){
            new_dat <- data.frame(y=unlist(perf@y.values),
                                  feature=rep(feat_names[i],
                                              length(perf@y.values)))
        }else{
            new_dat <- data.frame(y=perf@y.values[[1]],
                                  feature=rep(feat_names[i],
                                              length(perf@y.values[[1]])))
        }
        if (is.null(perf_dat)){
            perf_dat <- new_dat
        }else{
            perf_dat <- rbind(perf_dat, new_dat)
        }
    }
    return(perf_dat)
}

rocr2d.ggplot2 <- function(perf.data, x.name, y.name, n, abline, restrict_lims = TRUE){
    p <- ggplot(data=perf.data, aes(x=x, y=y, colour=feature)) +
        geom_line(size=1.0) +
        labs(x=x.name, y=y.name)
    if (restrict_lims){
        p <- p +
            xlim(0, 1) +
            ylim(0, 1)
    }
    if (n > 1){
        perf.data$y.min <- perf.data$y - perf.data$y.ci
        perf.data$y.max <- perf.data$y + perf.data$y.ci
        errorbar.rows <- c()
        for (f in levels(perf.data$feature)){
            perf.data.sub <- which(perf.data$feature==f)
            rows <- perf.data.sub[seq(1, length(perf.data.sub), length=6)[2:5]]
            errorbar.rows <- c(errorbar.rows, rows)
        }
        errorbar.df <- data.frame(perf.data[errorbar.rows,])
        p <- p + geom_linerange(mapping=aes(x=x, ymin=y.min, ymax=y.max),
                               data=errorbar.df,
                               inherit.aes=FALSE,
                               width=0.025)
    }
    if (length(unique(levels(perf.data$feature)))==1){
        p <- p + guides(colour=FALSE)##  +
            ## labs(title=unique(levels(perf.data$feature))[1])
    }
    if (!is.null(abline)){
        p <- p + geom_abline(intercept=abline[1], slope=abline[2], linetype=2,
                             colour="grey")
    }
    return(p)
}

rocr1d.ggplot2 <- function(perf.data, y.name, n, abline){
    p <- ggplot(data=perf.data, aes(x=alpha, y=y, colour=feature)) +
        geom_line(size=1.0) +
        xlim(0, 1) +
        ylim(0, 1) +
        labs(x="cutoff", y=y.name)
    if (n > 1){
        perf.data$y.min <- perf.data$y - perf.data$y.ci
        perf.data$y.max <- perf.data$y + perf.data$y.ci
        errorbar.rows <- c()
        for (f in levels(perf.data$feature)){
            perf.data.sub <- which(perf.data$feature==f)
            rows <- perf.data.sub[seq(1, length(perf.data.sub), length=6)[2:5]]
            errorbar.rows <- c(errorbar.rows, rows)
        }
        errorbar.df <- data.frame(perf.data[errorbar.rows,])
        p <- p + geom_linerange(mapping=aes(x=x, ymin=y.min, ymax=y.max),
                               data=errorbar.df,
                               inherit.aes=FALSE,
                               width=0.025)
    }
    if (length(unique(levels(perf.data$feature)))==1){
        p <- p + guides(colour=FALSE)##  +
            ## labs(title=unique(levels(perf.data$feature))[1])
    }
    if (!is.null(abline)){
        p <- p + geom_abline(intercept=abline[1], slope=abline[2], linetype=2,
                             colour="grey")
    }
    return(p)
}

rocr.sum.ggplot2 <- function(perf.data, y.name, n, abline){
    if (length(levels(perf.data$feature)) > 1){
        p <- ggplot(data=perf.data, aes(feature, y, colour=feature, fill=feature)) +
            ylim(0, 1) +
            ylab(y.name) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }else{
        p <- ggplot(data=perf.data, aes(y, colour=feature, fill=feature)) +
            xlim(0, 1) +
            xlab(y.name) +
            theme(legend.position="none")##  +
            ## labs(title=unique(levels(perf.data$feature))[1])
    }
    if (n > 1){
        if (length(levels(perf.data$feature)) > 1){
            p <- p + geom_violin(alpha=0.1, size=1.0)
            if (!is.null(abline)){
                p <- p + geom_abline(intercept=abline[1], slope=abline[2],
                                     linetype=2, colour="grey")
            }
        }else{
            p <- p + geom_density(alpha=0.1, size=1.0)
        }
    }else{
        p <- p + geom_point(size=3)
    }
    return(p)
}
