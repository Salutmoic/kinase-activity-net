library(reshape2)

argv <- commandArgs(TRUE)
if (length(argv) != 2){
    stop("USAGE: <script> DATA_FILE OUT_FILE")
}

kin.act.file <- argv[1]
out.file <- argv[2]

kin.act.tbl <- read.table(kin.act.file, as.is=TRUE)
names(kin.act.tbl) <- c("kinase", "condition", "activity")
kin.act.m <- as.data.frame(acast(kin.act.tbl, condition ~ kinase))

kin.cond.tbl <- read.delim("data/external/kinase_invivoconditions.csv",
                           as.is=TRUE, sep=",")

kin.act.cor <- cor(kin.act.m, method="spearman")
rownames(kin.act.cor) <- colnames(kin.act.m)
colnames(kin.act.cor) <- colnames(kin.act.m)

intxns <- NULL

for (kinase in colnames(kin.act.m)){
    message(paste(kinase, "-", paste0(match(kinase, colnames(kin.act.m)), "/", ncol(kin.act.m))))
    if (kinase %in% kin.cond.tbl$HGNC){
        kin.conds <- subset(kin.cond.tbl, HGNC==kinase)
        unperturbed <- setdiff(rownames(kin.act.m), kin.conds$Condition)
        kin.act.sub <- kin.act.m[unperturbed, order(kin.act.cor[kinase,], decreasing=TRUE)]
    }else{
        kin.act.sub <- kin.act.m[,order(kin.act.cor[kinase,], decreasing=TRUE)]
    }
    kin.model.formula <- as.formula(paste0(kinase, "~."))
    kin.model <- lm(kin.model.formula, data=kin.act.sub)
    #### kin.model.simple <- step(kin.model, k=log(nrow(kin.act.sub)), trace=0)
    ## kin.model.simple <- step(kin.model, trace=0)
    ## regulators <- attr(kin.model.simple$terms, "term.labels")
    ## kin.model.anova <- anova(kin.model.simple)
    regulators <- attr(kin.model$terms, "term.labels")
    kin.model.anova <- anova(kin.model)
    kinase1 <- colnames(kin.act.sub)
    kinase2 <- rep(kinase, ncol(kin.act.sub))
    f.values <- kin.model.anova$"F value"
    ## f.values <- f.values/max(f.values, na.rm=TRUE)
    f.value <- c()
    for (i in 1:length(kinase1)){
        if (kinase1[i] == kinase || !(kinase1[i] %in% regulators)){
            f.value <- c(f.value, 0.0)
        }else{
            kin.index <- match(kinase1[i], regulators)
            f.value <- c(f.value, f.values[kin.index])
        }
    }
    intxns.tmp <- cbind(kinase1, kinase2, f.value)
    if (is.null(intxns)){
        intxns <- intxns.tmp
    }else{
        intxns <- rbind(intxns, intxns.tmp)
    }
}

intxns <- as.data.frame(intxns)
names(intxns) <- c("node1", "node2", "f.value")
intxns <- intxns[order(intxns$node1, intxns$node2),]
write.table(intxns, out.file, row.names=FALSE, col.names=TRUE, sep="\t",
            quote=FALSE)
