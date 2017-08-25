library(reshape2)

raw.dat <- read.delim("data/external/BaitPreyPairs_noFilters_BP2a.tsv", as.is=TRUE)

dat.filt <- raw.dat[which(!grepl("^tr", raw.dat$db_protein_id) &
                          !grepl("Fragment", raw.dat$prot_description) &
                          !grepl("Isoform", raw.dat$prot_description)),
                    c("bait_symbol", "symbol", "pInt")]
names(dat.filt) <- c("node1", "node2", "p.int")

kinome <- read.table("data/human-kinome.txt", as.is=TRUE)$V1

dat.kin <- subset(dat.filt, prot1 %in% kinome & prot2 %in% kinome &
                            prot1 != prot2)
rownames(dat.kin) <- paste(dat.kin$prot1, dat.kin$prot2, sep="-")

dat.kin.rev <- dat.kin
dat.kin.rev$prot1 <- dat.kin.rev$prot2
dat.kin.rev$prot2 <- dat.kin$prot1
rownames(dat.kin.rev) <- paste(dat.kin.rev$prot1, dat.kin.rev$prot2, sep="-")
dat.kin.rev <- dat.kin.rev[which(!(rownames(dat.kin.rev) %in% rownames(dat.kin))),]

dat.kin.supp <- rbind(dat.kin, dat.kin.rev)

write.table(dat.kin.supp, "out/bioplex.tsv", sep="\t", quote=FALSE,
            row.names=FALSE, col.names=TRUE)
