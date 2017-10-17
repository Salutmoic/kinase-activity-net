load("data/log.wKSEA.kinase_condition.clean.Rdata")

argv <- commandArgs(TRUE)
if (length(argv) != 3){
    stop("USAGE: <script> KEGG_RELS ID_MAPPING OUT_FILE")
}

kegg.rels.file <- argv[1]
id.map.file <- argv[2]
out.file <- argv[3]

kegg.rels <- read.delim(kegg.rels.file, as.is=TRUE)
names(kegg.rels) <- c("uniprot1", "uniprot2")
id.map <- read.delim(id.map.file, as.is=TRUE)
names(id.map) <- c("uniprot", "gene.name")

merged.dat <- merge(kegg.rels, id.map,
                    by.x=c("uniprot1"),
                    by.y=c("uniprot"))
names(merged.dat) <- c("uniprot1", "uniprot2", "prot1")
merged.dat2 <- merge(merged.dat, id.map,
                    by.x=c("uniprot2"),
                    by.y=c("uniprot"))
names(merged.dat2) <- c("uniprot2", "uniprot1", "prot1", "prot2")

kin.sub.pairs <- unique(subset(merged.dat2, select=c("prot1", "prot2")))

kin.act <- log.wKSEA.kinase_condition.clean
kin.sub.filt <- subset(kin.sub.pairs, prot1 %in% rownames(kin.act) &
                                      prot2 %in% rownames(kin.act))

write.table(kin.sub.filt, out.file,
            quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
