## phos stores kinase activiteis acorss conditions (Ochoa et al. 2016)
kinome = read.delim("data/human-kinome.txt", as.is=TRUE)

ksub = read.delim("data/psiteplus-kinase-substrates.tsv", as.is=TRUE)

kinases = kinome[,1]

overlap = c()

## see the overlap between each and every kinase pair
for(i in 1:(length(kinases)-1)){
    m =nrow(overlap)

    for(j in (i+1):length(kinases)){
        k1 = ksub[which(toupper(ksub$GENE) == toupper(as.character(kinases[i]))),]
        k2 = ksub[which(toupper(ksub$GENE) == toupper(as.character(kinases[j]))),]

        ## paste modified residue to substrate name
        subs1 = paste(k1$SUBSTRATE, k1$SUB_MOD_RSD)
        subs2 = paste(k2$SUBSTRATE, k2$SUB_MOD_RSD)

        subs1 = unique(subs1)
        subs2 = unique(subs2)

        inter = length(intersect(toupper(subs1),toupper(subs2)))

        if(length(which(toupper(subs1) %in% toupper(subs2))) != inter){
            break
        }

        row = c(as.character(kinases[i]), as.character(length(subs1)),
                as.character(kinases[j]), as.character(length(subs2)),
                as.character(inter), as.character(inter/(min(length(subs2),length(subs1)))))
        overlap = rbind(overlap, row)

        ## print(c(i,j))
    }

    if(is.null(m)){
        m = 0
    }

    ## if(nrow(overlap) != m + (length(kinases)-i)){
    ##     break
    ## }

}

colnames(overlap) = c("kinase1", "no.substrates1", "kinase2", "no.substrates2",
                      "no.shared.substrates", "proportion.of.substrates.shared")

## We can do this filtering in gen-activity-table.r
## overlap.sub = overlap[which(as.numeric(overlap[,6]) !=0), ]
## overlap.sub = overlap.sub[which(as.numeric(overlap.sub[,2])  > 10 & as.numeric(overlap.sub[,4]) > 10), ]

## overlap.sub = overlap.sub[order(as.numeric(overlap.sub[,6]),decreasing = TRUE),]
write.table(overlap, "data/psiteplus-kinase-substrate-overlap.tsv", sep = "\t",
            quote=FALSE, row.names=FALSE, col.names=TRUE)
