## phos stores kinase activiteis acorss conditions (Ochoa et al. 2016)
kinome = read.delim("data/human-kinome.txt", as.is=TRUE)

ksub = read.delim("data/psiteplus-kinase-substrates.tsv", as.is=TRUE)

kinases = kinome[,1]

overlap = c()

## see the overlap between each and every kinase pair
for(i in 1:(length(kinases)-1)){
    m =nrow(overlap)

    for(j in (i+1):length(kinases)){
        k1 = ksub[which(ksub$GENE == kinases[i]),]
        k2 = ksub[which(ksub$GENE == kinases[j]),]

        ## paste modified residue to substrate name
        subs1 = paste(k1$SUB_GENE, k1$SUB_MOD_RSD)
        subs2 = paste(k2$SUB_GENE, k2$SUB_MOD_RSD)

        subs1 = unique(subs1)
        subs2 = unique(subs2)

        inter = length(intersect(subs1, subs2))
        if (length(subs1) == 0 || length(subs2) == 0){
            perc.overlap <- 0.0
        }else{
            perc.overlap <- inter/min(length(subs1), length(subs2))
        }

        row = c(kinases[i], length(subs1), kinases[j], length(subs2), inter,
                perc.overlap)
        overlap = rbind(overlap, row)
    }

    if(is.null(m)){
        m = 0
    }
}

colnames(overlap) = c("kinase1", "no.substrates1", "kinase2", "no.substrates2",
                      "no.shared.substrates", "proportion.of.substrates.shared")

## We can do this filtering in gen-activity-table.r
## overlap.sub = overlap[which(as.numeric(overlap[,6]) !=0), ]
## overlap.sub = overlap.sub[which(as.numeric(overlap.sub[,2])  > 10 & as.numeric(overlap.sub[,4]) > 10), ]

## overlap.sub = overlap.sub[order(as.numeric(overlap.sub[,6]),decreasing = TRUE),]
write.table(overlap, "data/psiteplus-kinase-substrate-overlap.tsv", sep = "\t",
            quote=FALSE, row.names=FALSE, col.names=TRUE)
