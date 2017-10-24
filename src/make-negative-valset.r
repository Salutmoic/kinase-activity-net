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
if (length(argv) != 3){
    stop("USAGE: <script> TRUE_POS PROT_GROUPS OUT_FILE")
}

true.pos.file <- argv[1]
prot.groups.file <- argv[2]
out.file <- argv[3]

true.intxns.tbl <- read.table(true.pos.file, as.is=TRUE)
true.intxns <- paste(true.intxns.tbl[,1], true.intxns.tbl[,2], sep="-")

prot.groups <- read.table(prot.groups.file, as.is=TRUE)
names(prot.groups) <- c("group", "protein")

all.kinases <- unique(c(true.intxns.tbl[,1], true.intxns.tbl[,2]))
possible.false.intxns <- gen.possible.false.intxns(all.kinases, prot.groups)
## Remove true interactions from the set of possible false interactions
false.intxns <- setdiff(possible.false.intxns, true.intxns)

false.intxns.tbl <- as.data.frame(matrix(unlist(strsplit(false.intxns, split="-")),
                                         ncol=2, byrow=TRUE))

write.table(false.intxns.tbl, out.file, row.names=FALSE, col.names=FALSE,
            quote=FALSE, sep="\t")
