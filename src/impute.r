library(reshape2)
suppressMessages(library(VIM))

impute.on.table <- function(tbl){
    tbl.imp <- data.frame(kNN(tbl, imp_var=FALSE))
    rownames(tbl.imp) <- rownames(tbl)
    colnames(tbl.imp) <- colnames(tbl)
    tbl.imp$kinase <- rownames(tbl.imp)
    return(tbl.imp)
}

argv <- commandArgs(TRUE)
if (length(argv) != 2){
    stop("USAGE: <script> DATA_FILE OUT_FILE")
}

data.tbl.file <- argv[1]
out.file <- argv[2]

data.tbl <- read.table(data.tbl.file, as.is=TRUE)
names(data.tbl) <- c("kinase", "condition", "activity")
data.tbl.m <- acast(data.tbl, kinase ~ condition)
data.tbl.imp <- impute.on.table(data.tbl.m)
data.tbl.imp.df <- melt(data.tbl.imp)

write.table(data.tbl.imp.df[order(data.tbl.imp.df[,1], data.tbl.imp.df[,2]),],
            out.file, col.names=FALSE, row.names=FALSE,
            quote=FALSE, sep="\t")
