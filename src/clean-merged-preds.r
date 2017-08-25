argv <- commandArgs(TRUE)
if (length(argv) != 1){
    stop("USAGE: <script> MERGED_PREDS")
}

merged.preds.file <- argv[1]

merged.preds <- read.delim(merged.preds.file, as.is=TRUE)
num.feats <- ncol(merged.preds) - 2

main.feats <- setdiff(colnames(merged.preds[3:ncol(merged.preds)]),
                      c("node1", "node2", "prot1", "prot2", "kinase.type",
                        "sub.kinase.type", "go.bp.sem.sim", "go.cc.sem.sim",
                        "tissue.sel1", "tissue.sel2", "celline.sel1",
                        "celline.sel2"))
                        ## "tissue.coexpr", "tissue.sim", "prot.atlas.cell",
                        ## "prot.atlas.tissue"))

## Find rows where all of the phospho/pssm-based features are missing
bad.rows <- which(apply(merged.preds[,main.feats], 1,
                        function(row){
                            Reduce(`&&`, is.na(row))
                        }))

merged.preds[bad.rows, 3:ncol(merged.preds)] <- rep(NA, num.feats)

file.base <- strsplit(basename(merged.preds.file), split="\\.")[[1]][1]
out.file <- paste0("out/", file.base, "-clean.tsv")

write.table(merged.preds, out.file, quote=FALSE, sep="\t", row.names=FALSE,
            col.names=TRUE)
