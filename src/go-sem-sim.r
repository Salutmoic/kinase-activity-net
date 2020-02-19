library(GOSim)

argv <- commandArgs(TRUE)
if (length(argv) != 1){
    stop("USAGE: <script> N")
}
n <- as.integer(sub("^\\\\", "", argv[1]))

setEvidenceLevel(evidences = c("IMP", "IGI", "IPI", "ISS", "IDA", "IEP", "IEA"))

full.id.map <- read.table("data/external/HUMAN_9606_idmapping.dat", sep="\t",
                     as.is=TRUE)
colnames(full.id.map) <- c("uniprot.id", "id.type", "other.id")
uniprot.to.geneid <- subset(full.id.map,
                            id.type=="GeneID")[,c("uniprot.id", "other.id")]
colnames(uniprot.to.geneid) <- c("uniprot", "geneid")
uniprot.to.hgnc <- read.table("data/uniprot-id-map.tsv", as.is=TRUE)
colnames(uniprot.to.hgnc) <- c("uniprot", "hgnc")
kinases <- read.table("data/human-kinome.txt", as.is=TRUE)$V1

merged.map <- merge(uniprot.to.geneid, uniprot.to.hgnc)
hgnc.to.geneid <- unique(merged.map[c("geneid", "hgnc")])

kinases.geneid <- subset(hgnc.to.geneid, hgnc %in% kinases)$geneid
kinase.geneid <- kinases.geneid[n]

gene.sim <- getGeneSim(kinases.geneid[n],
                       kinases.geneid[n:length(kinases.geneid)],
                       similarity="OA",
                       similarityTerm="Lin", avg=TRUE)

prot1 <- rep(kinase, length(kinases.geneid)-n+1)
prot2 <- sapply(kinases.geneid[n:length(kinases.geneid)],
                function(geneid){
                    unique(subset(hgnc.to.geneid, geneid==geneid)$hgnc)
                })
gene.sim.hgnc.df <- data.frame(prot1=prot1, prot2=prot2,
                               go.sem.sim=t(gene.sim))

write.table(gene.sim.hgnc.df,
            paste0("out/go-sem-sim2/kinase-", n, "-sem-sim.tsv"),
            sep="\t", quote=FALSE, colnames=TRUE, rownames=FALSE)
