## phos stores kinase activiteis acorss conditions (Ochoa et al. 2016)
load("data/external/log.wGSEA.kinase_condition.clean.Rdata")
kin.act <- log.wGSEA.kinase_condition.clean

ksub = read.delim("data/human_kinase_table.tsv", as.is=TRUE)

kinases = rownames(kin.act)

overlap = c()

v= 0

## see the overlap between each and every kinase pair
for(i in 1:length(kinases)){
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

        v = v+1
        ## print(c(i,j))
    }

    if(is.null(m)){
        m = 0
    }

    if(nrow(overlap) != m + (length(kinases)-i)){
        break
    }

}

colnames(overlap) = c("kinase1", "no.substrates1", "kinase2", "no.substrates2",
                      "no.shared.substrates", "proportion.of.substrates.shared")

## We can do this filtering in gen-activity-table.r
## overlap.sub = overlap[which(as.numeric(overlap[,6]) !=0), ]
## overlap.sub = overlap.sub[which(as.numeric(overlap.sub[,2])  > 10 & as.numeric(overlap.sub[,4]) > 10), ]

overlap.sub = overlap.sub[order(as.numeric(overlap.sub[,6]),decreasing = TRUE),]
write.table(overlap.sub, "data/kinase_substrate_overlap.tsv", sep = "\t",
            quote=FALSE, row.names=FALSE, col.names=TRUE)
