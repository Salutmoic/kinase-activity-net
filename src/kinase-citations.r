library(ggplot2)
library(classInt)
library(scales)

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    library(grid)

    ## Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)

    numPlots = length(plots)

    ## If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        ## Make the panel
        ## ncol: Number of columns of plots
        ## nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
    }

    if (numPlots==1) {
        print(plots[[1]])

    } else {
        ## Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

        ## Make each plot, in the correct location
        for (i in 1:numPlots) {
            ## Get the i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
        }
    }
}

kin.dat <- read.delim("out/kinase-citations.tsv", as.is=TRUE)
rel.dat <- read.delim("out/kinase-kinase-rels.tsv", as.is=TRUE)

num.cite.sub <- subset(kin.dat, count.type=="num.subs")$num.cites.filt
num.cite.ints <- classIntervals(num.cite.sub, 10, style="quantile")
num.cite.ints$brks <- round(num.cite.ints$brks)
num.cite.ints
kin.dat$cite.class <- cut(kin.dat$num.cites.filt, num.cite.ints$brks,
                          include.lowest=TRUE, dig.lab=4)

rel.dat$cite.class1 <- cut(rel.dat$num.cites.filt1, num.cite.ints$brks,
                           include.lowest=TRUE, dig.lab=4)
rel.dat$cite.class2 <- cut(rel.dat$num.cites.filt2, num.cite.ints$brks,
                           include.lowest=TRUE, dig.lab=4)

kin.dat$count.type <- factor(kin.dat$count.type,
                             levels=c("num.subs",
                                      "num.kin.subs",
                                      "num.phos.kins",
                                      "num.reg.sites"),
                             labels=c("substrates",
                                      "substrate kinases",
                                      "phosphorylating kinases",
                                      "regulatory sites"))

rel.dat$rel <- factor(rel.dat$rel,
                      levels=c("phos", "reg"),
                      labels=c("phosphorylation relationships", "regulatory relationships"))

pdf("img/kinase-citations-1.pdf", height=14, width=14)
ggplot(kin.dat, aes(x=cite.class, y=count, color=count.type)) +
    geom_count(shape=19, alpha=0.8) +
    scale_size(range=c(1, 12)) +
    facet_wrap(~count.type, scales="free") +
    theme_bw(base_size=20) +
    theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
    xlab("number of publications linked to protein kinase") +
    ylab("count") +
    scale_y_continuous(breaks=pretty_breaks()) +
    guides(colour=FALSE) +
    labs(size="number of\nkinases")
dev.off()

pdf("img/kinase-citations-2.pdf", height=7, width=14)
ggplot(subset(rel.dat, kinase1!=kinase2), aes(x=cite.class1, y=cite.class2)) +
    geom_bin2d() +
    facet_wrap(~rel, scales="free") +
    scale_fill_distiller(direction=1, palette="Blues", na.value="white",
                         guide_colourbar(title="number of\ninteracting\npairs")) +
    scale_x_discrete(drop=FALSE) +
    scale_y_discrete(drop=FALSE) +
    theme_bw(base_size=20) +
    theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    xlab("number of publications (kinase)") +
    ylab("number of publications (substrate kinase)")
dev.off()
