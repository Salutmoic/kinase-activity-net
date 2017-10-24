gen.possible.false.intxns <- function(all.kinases, prot.groups){
    possible.false.intxns <- c()
    for (k1 in all.kinases){
        for (k2 in all.kinases){
            if (k1 == k2)
                next
            if (!(k1 %in% prot.groups$protein) || !(k2 %in% prot.groups$protein))
                next
            k1.groups <- unique(subset(prot.groups, protein == k1)$group)
            k2.groups <- unique(subset(prot.groups, protein == k2)$group)
            if (length(k1.groups)>0 && length(k2.groups)>0 && length(intersect(k1.groups, k2.groups)) == 0)
                possible.false.intxns <- c(possible.false.intxns, paste(k1, k2, sep="-"))
        }
    }
    return (possible.false.intxns)
}

argv <- commandArgs(TRUE)
if (length(argv) != 4){
    stop("USAGE: <script> DATA_FILE TRUE_POS PROT_GROUPS OUT_FILE")
}

pred.score.file <- argv[1]
true.pos.file <- argv[2]
prot.groups.file <- argv[3]
out.file <- argv[4]

pred.score <- read.delim(pred.score.file, as.is=TRUE)
names(pred.score) <- c("prot1", "prot2", "pred.score")
pred.score <- subset(pred.score, prot1 != prot2)
rownames(pred.score) <- paste(pred.score$prot1, pred.score$prot2, sep="-")

true.intxns.tbl <- read.table(true.pos.file, as.is=TRUE)
true.intxns <- paste(true.intxns.tbl[,1], true.intxns.tbl[,2], sep="-")

prot.groups.full <- read.table(prot.groups.file, as.is=TRUE)
names(prot.groups.full) <- c("group", "protein")

all.kinases <- unique(c(pred.score$prot1, pred.score$prot2))
prot.groups <- subset(prot.groups.full, protein %in% all.kinases)
possible.false.intxns <- gen.possible.false.intxns(all.kinases, prot.groups)
## Remove true interactions from the set of possible false interactions
false.intxns <- setdiff(possible.false.intxns, true.intxns)

false.intxns.tbl <- as.data.frame(matrix(unlist(strsplit(false.intxns, split="-")),
                                         ncol=2, byrow=TRUE))

write.table(false.intxns.tbl, out.file, row.names=FALSE, col.names=FALSE,
            quote=FALSE, sep="\t")
