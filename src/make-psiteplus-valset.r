argv <- commandArgs(TRUE)
if (length(argv) != 2){
    stop("USAGE: <script> KINASE_SUBSTRATE REG_SITES")
}

kin.sub.file <- argv[1]
reg.sites.file <- argv[2]

kin.sub <- read.delim(kin.sub.file, as.is=TRUE)
reg.sites <- read.delim(reg.sites.file, as.is=TRUE)

merged.dat <- merge(kin.sub, reg.sites,
                    by.x=c("SUB_GENE", "SUB_MOD_RSD"),
                    by.y=c("GENE", "MOD_RSD"))

kin.sub.pairs <- unique(subset(merged.dat, select=c("GENE", "SUB_GENE")))

write.table(kin.sub.pairs, "data/validation-set-psiteplus.tsv",
            quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
