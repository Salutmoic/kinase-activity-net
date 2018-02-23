argv <- commandArgs(TRUE)
if (length(argv) != 2){
    stop("USAGE: <script> KIN_SUB_TBL OUT_FILE")
}

psite.plus.file <- argv[1]
out.file <- argv[2]

kinases <- read.table("data/human-kinome.txt", as.is=TRUE)[,1]

psites = read.table("data/pride-phosphosites.tsv", as.is=TRUE, sep="\t")
names(psites) <- c("protein", "position", "residue", "motif")
pride.sites <- paste(psites$protein, psites$position, sep="-")

psite.plus.tbl <- read.delim(psite.plus.file, as.is=TRUE, sep="\t")

psite.plus.filt <- subset(psite.plus.tbl, GENE %in% kinases)

if ("SUB_MOD_RSD" %in% colnames(psite.plus.filt)){
    psp.pos <- sub("[STY]", "", psite.plus.filt$SUB_MOD_RSD)
    psp.pos <- sub("-p", "", psp.pos)
    psp.sites <- paste(psite.plus.filt$SUB_GENE, psp.pos, sep="-")
    psite.plus.filt <- psite.plus.filt[which(psp.sites %in% pride.sites),]
}else{
    psp.pos <- sub("[STY]", "", psite.plus.filt$MOD_RSD)
    psp.pos <- sub("-p", "", psp.pos)
    psp.sites <- paste(psite.plus.filt$GENE, psp.pos, sep="-")
    psite.plus.filt <- psite.plus.filt[which(psp.sites %in% pride.sites),]
}

write.table(psite.plus.filt, out.file, sep="\t", quote=FALSE,
            row.names=FALSE, col.names=TRUE)
