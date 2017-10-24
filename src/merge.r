argv <- commandArgs(TRUE)
if (length(argv) < 3){
    stop("USAGE: <script> OUT_FILE PRED_FILES")
}

out.file <- argv[1]

started <- FALSE
for (f in argv[2:length(argv)]){
    if (!started){
        merged.dat <- read.delim(f, as.is=TRUE)
        started <- TRUE
        meth.match <- regexpr("([a-z0-9]*)\\.tsv", f, perl=TRUE)
        old.meth <- substring(f, attr(meth.match, "capture.start"),
                              attr(meth.match, "capture.start")
                              + attr(meth.match, "capture.length")-1)
        next
    }
    new.dat <- read.delim(f, as.is=TRUE)
    meth.match <- regexpr("([a-z0-9]*)\\.tsv", f, perl=TRUE)
    new.meth <- substring(f, attr(meth.match, "capture.start"),
                      attr(meth.match, "capture.start")
                      + attr(meth.match, "capture.length")-1)
    merged.dat <- merge(merged.dat, new.dat,
                        by.x=c("node1", "node2"),
                        by.y=c("node1", "node2"),
                        suffixes=c(paste0(".", old.meth), paste0(".", new.meth)),
                        all=TRUE)
    old.meth <- new.meth
}

write.table(merged.dat, out.file, row.names=FALSE, col.names=TRUE,
            sep="\t", quote=FALSE)
