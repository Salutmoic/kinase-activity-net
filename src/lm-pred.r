library(reshape2)

argv <- commandArgs(TRUE)
if (length(argv) != 3){
    stop("USAGE: <script> DATA_FILE PSSM_FILE OUT_FILE")
}

kin.act.file <- argv[1]
pssm.file <- argv[2]
out.file <- argv[3]

kin.act.tbl <- read.table(kin.act.file, as.is=TRUE)
names(kin.act.tbl) <- c("kinase", "condition", "activity")
kin.act.m <- acast(kin.act.tbl, condition ~ kinase)

kin.cond.tbl <- read.delim("data/external/kinase-condition-pairs.tsv",
                           as.is=TRUE, sep=",")

## kin.act.cor <- cor(kin.act.m, method="spearman")
## rownames(kin.act.cor) <- colnames(kin.act.m)
## colnames(kin.act.cor) <- colnames(kin.act.m)

pssm <- read.delim(pssm.file, as.is=TRUE)
pssm.m <- acast(pssm, node2 ~ node1)
pssm.m <- pssm.m[colnames(kin.act.m), colnames(kin.act.m)]

print(c(dim(kin.act.m), dim(pssm.m)))

kin1s <- c()
kin2s <- c()
coef.vals <- c()
p.vals <- c()
t.vals <- c()

for (kinase in colnames(kin.act.m)){
    kin.sub.pssm <- pssm.m[kinase,]
    if (all(is.na(pssm.m)) || all(pssm.m == 0))
        next
    kin.sub.pssm[kinase] <- 1.0
    kin.sub.pssm.m <- matrix(rep(kin.sub.pssm, nrow(kin.act.m)), nrow=nrow(kin.act.m),
                             byrow=TRUE)
    kin.act.mod <- kin.act.m * kin.sub.pssm.m
    good.cols <- which(apply(kin.act.mod, 2,
                             function(col){
                                 length(which(is.na(col)))/length(col) < 1/3
                             }))
    kin.sub.pssm <- kin.sub.pssm[good.cols]
    if (kinase %in% kin.cond.tbl$Kinase){
        kin.conds <- subset(kin.cond.tbl, Kinase==kinase)
        unperturbed <- setdiff(rownames(kin.act.mod), kin.conds$Condition)
        kin.act.sub <- kin.act.mod[unperturbed, good.cols]
        ## kin.act.sub <- kin.act.sub[,order(kin.sub.pssm, decreasing=TRUE)]
    }else{
        kin.act.sub <- kin.act.mod[,good.cols]
        ## kin.act.sub <- kin.act.sub[,order(kin.sub.pssm, decreasing=TRUE)]
    }
    if (!(kinase %in% colnames(kin.act.sub)) || ncol(kin.act.sub) < 3){
        next
    }
    kin.model.formula <- as.formula(paste0(kinase, "~ ."))
    kin.model <- lm(kin.model.formula, data=as.data.frame(kin.act.sub))
    ## BIC
    kin.model.simple <- step(kin.model, k=log(nrow(kin.act.sub)), trace=0)
    ## kin.model.simple <- step(kin.model, trace=0)
    regulators <- attr(kin.model.simple$terms, "term.labels")
    kin.model.anova <- anova(kin.model.simple)
    kin.model.sum <- summary(kin.model.simple)
    ## regulators <- attr(kin.model$terms, "term.labels")
    ## kin.model.anova <- anova(kin.model)
    kinase1 <- colnames(kin.act.sub)
    kinase2 <- rep(kinase, ncol(kin.act.sub))
    ## f.values <- kin.model.anova$"F value"
    coefs <- kin.model.simple$coefficients
    ## f.values <- f.values/max(f.values, na.rm=TRUE)
    ## f.value <- c()
    coef <- c()
    t.val <- c()
    p.val <- c()
    for (i in 1:length(kinase1)){
        if (kinase1[i] == kinase || !(kinase1[i] %in% regulators)){
            ## f.value <- c(f.value, 0.0)
            coef <- c(coef, NA)
            t.val <- c(t.val, NA)
            p.val <- c(p.val, NA)
        }else{
            kin.index <- match(kinase1[i], regulators)
            ## f.value <- c(f.value, f.values[kin.index])
            ## f.pval <- kin.model.anova$"Pr(>F)"[kin.index]
            t.val <- c(t.val, kin.model.sum$coefficients[kinase1[i], "t value"])
            p.val <- c(p.val, kin.model.sum$coefficients[kinase1[i], "Pr(>|t|)"])
            coef <- c(coef, coefs[kinase1[i]])
        }
    }
    kin1s <- c(kin1s, kinase1)
    kin2s <- c(kin2s, kinase2)
    coef.vals <- c(coef.vals, coef)
    t.vals <- c(t.vals, t.val)
    p.vals <- c(p.vals, p.val)
}
intxns <- data.frame(node1=kin1s, node2=kin2s, coefficient=coef.vals, t=t.vals,
                     p=p.vals, p.adj=p.adjust(p.vals, method="fdr"))
intxns <- intxns[complete.cases(intxns),]
intxns <- intxns[order(intxns$node1, intxns$node2),]
write.table(intxns, out.file, row.names=FALSE, col.names=TRUE, sep="\t",
            quote=FALSE)
