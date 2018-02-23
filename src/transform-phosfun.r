argv <- commandArgs(TRUE)
if (length(argv) != 2){
    stop("USAGE: <script> PHOSFUN_FILE OUT_FILE")
}

phosfun_file <- argv[1]
out_file <- argv[2]

phosfun <- read.table(phosfun_file, sep="\t")
rownames(phosfun) <- paste(phosfun[,1], phosfun[,2], sep="-")
psites <- read.table("data/pride-phosphosites-kinome.tsv", as.is=TRUE, sep="\t")
rownames(psites) <- paste(psites[,1], psites[,2], sep="-")

phosfun <- phosfun[which(rownames(phosfun) %in% rownames(psites)),]

phosfun[,3] <- log(phosfun[,3])
phosfun[,3] <- (phosfun[,3]-min(phosfun[,3]))/(max(phosfun[,3])-min(phosfun[,3]))

write.table(phosfun, out_file, sep="\t", quote=FALSE, row.names=FALSE,
            col.names=FALSE)
