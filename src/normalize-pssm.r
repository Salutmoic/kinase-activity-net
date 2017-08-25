library(dplyr)

pssm_file <- "out/kinase-pssm.tsv"
pssm <- read.delim(pssm_file)

pssm_summ <- pssm %>% group_by(node1) %>%
    summarise(min.score=min(max.pssm.score, na.rm=TRUE),
              max.score=max(max.pssm.score, na.rm=TRUE))
rownames(pssm_summ) <- pssm_summ$node1

pssm$max.pssm.score <- sapply(1:nrow(pssm),
                              function(i){
                                  kin <- pssm$node1[i]
                                  min.score <- pssm_summ$min.score[kin]
                                  max.score <- pssm_summ$max.score[kin]
                                  score <- pssm$max.pssm.score[i]
                                  return((score-min.score)/(max.score-min.score))
                              })

write.table(pssm, pssm_file, sep="\t", quote=FALSE, row.names=FALSE,
            col.names=TRUE)
