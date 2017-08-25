library(reshape2)
library(dplyr)

kinbase.id.map.full <- read.table("data/kinbase-id-map.tsv", as.is=TRUE)
kinbase.id.map <- kinbase.id.map.full$V2
names(kinbase.id.map) <- kinbase.id.map.full$V1

map.id <- function(id, id.map){
    id.split <- strsplit(id, split="_")[[1]]
    id.split[1] <- id.map[id.split[1]]
    if (length(id.split)==1)
        return(id.split[1])
    return(paste(id.split[1], id.split[2], sep="_"))
}

sdr.data <- readRDS("data/external/SDR_score_mat.rds")
colnames(sdr.data) <- as.vector(sapply(colnames(sdr.data), map.id, kinbase.id.map))
rownames(sdr.data) <- colnames(sdr.data)

col.kinases <- as.vector(unlist(lapply(strsplit(colnames(sdr.data), split="_"),
                                       function(foo) foo[1])))
kinases <- unique(col.kinases)

kin.sub.data <- read.delim("data/psiteplus-kinase-substrates.tsv", as.is=TRUE)
kin.sub.data <- subset(kin.sub.data, GENE != SUB_GENE)

kins.sub.counts <- kin.sub.data %>%
    group_by(GENE) %>%
    summarise(n=n())
kins.sub.counts <- as.data.frame(kins.sub.counts)
rownames(kins.sub.counts) <- kins.sub.counts$GENE
kins.w.pssms <- kins.sub.counts %>%
    filter(n>=10) %>%
    `$`("GENE")

sdr.data.sub <- sdr.data[,which(col.kinases %in% kins.w.pssms)]
col.kinases.sub <- as.vector(unlist(lapply(strsplit(colnames(sdr.data.sub), split="_"),
                                           function(foo) foo[1])))

kins <- NULL
pssms <- NULL
scores <- NULL
for (kin in kinases){
    if (kin %in% kins.w.pssms){
        if (is.null(kins)){
            kins <- kin
            pssms <- kin
            scores <- 1.0
        }else{
            kins <- c(kins, kin)
            pssms <- c(pssms, kin)
            scores <- c(scores, 1.0)
        }
        next
    }
    kin.rows <- which(col.kinases==kin)
    if (length(kin.rows) == 1){
        row <- sdr.data.sub[kin.rows,]
        best <- which.max(row)
        best.val <- row[best]
        ## If more than one kinase is equally close, then choose the
        ## one with the most known substrates
        if (length(which(row==best.val)) > 1){
            matches <- which(row==best.val)
            kin.matches <- col.kinases.sub[matches]
            kin.m.sub.counts <- kins.sub.counts[kin.matches, "n"]
            best.sub.count <- which.max(kin.m.sub.counts)
            best <- matches[best.sub.count]
        }
        pssm <- col.kinases.sub[best]
        score <- sdr.data.sub[kin.rows, best]
    }else{
        rows <- sdr.data.sub[kin.rows,]
        best.cols <- apply(rows, 1, which.max)
        best.vals <- sapply(1:nrow(rows), function(n) rows[n, best.cols[n]])
        best.row <- which.max(best.vals)
        best.col <- best.cols[best.row]
        best.val <- max(best.vals)
        if (length(which(rows==best.val)) > 1){
            matches <- which(apply(rows, 2, function(col) any(col==best.val)))
            kin.matches <- col.kinases.sub[matches]
            kin.m.sub.counts <- kins.sub.counts[kin.matches, "n"]
            best.sub.count <- which.max(kin.m.sub.counts)
            best.col <- matches[best.sub.count]
            best.row <- which(rows[,best.col]==best.val)
        }
        pssm <- col.kinases.sub[best.col]
        score <- sdr.data.sub[kin.rows[best.row], best.col]
    }
    if (is.null(kins)){
        kins <- kin
        pssms <- pssm
        scores <- score
    }else{
        kins <- c(kins, kin)
        pssms <- c(pssms, pssm)
        scores <- c(scores, score)
    }
}

pssm.assigns <- data.frame(kinase=kins, pssm=pssms, score=scores)
pssm.assigns$kinase <- as.character(pssm.assigns$kinase)
pssm.assigns$pssm <- as.character(pssm.assigns$pssm)
pssm.assigns <- pssm.assigns[order(pssm.assigns$kinase),]

write.table(pssm.assigns, "data/sdr-pssm-assignments.tsv", quote=FALSE,
            col.names=TRUE, row.names=FALSE, sep="\t")
