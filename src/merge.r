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
        next
    }
    new.dat <- read.delim(f, as.is=TRUE)
    merged.dat <- merge(merged.dat, new.dat,
                        by.x=c("node1", "node2"),
                        by.y=c("node1", "node2"),
                        all=TRUE)
}

write.table(merged.dat, out.file, row.names=FALSE, col.names=TRUE,
            sep="\t", quote=FALSE)
