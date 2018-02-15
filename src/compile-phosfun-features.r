get_pos_in_domain <- function(pos, dom_start, dom_end){
    if (pos < dom_start || pos > dom_end)
        return(NA)
    return((pos-dom_start)/(dom_end-dom_start))
}

phospho_table <- read.delim("data/phosphosites_reduced.tsv", as.is=TRUE)
phospho_table$residue <- sapply(phospho_table$MOD_RSD, substr, 1, 1)
phospho_table$MOD_RSD <- sub("^[STY]", "", phospho_table$MOD_RSD)
phospho_table$MOD_RSD <- sub("-p$", "", phospho_table$MOD_RSD)

sites <- phospho_table[,c("GENE", "ACC_ID", "MOD_RSD", "residue")]

diso <- read.delim("data/external/phosfun-features/feature_disopred",
                   as.is=TRUE, sep=",")
foldx <- read.delim("data/external/phosfun-features/feature_foldx",
                    as.is=TRUE, sep=",")
hotspots <- read.delim("data/external/phosfun-features/feature_hotspots",
                       as.is=TRUE, sep=",")
prot_length <- read.delim("data/external/phosfun-features/feature_proteinlength",
                          as.is=TRUE, sep=",")
pfam <- read.delim("data/human-pfam.tsv", as.is=TRUE, sep="\t")

diso <- diso[c("acc", "position", "disopred_score")]
foldx <- foldx[c("acc", "position", "exp3d_max_ddG_effect",
                 "exp3d_ala_ddG_effect")]
hotspots <- hotspots[c("acc", "position", "log10_hotspot_pval_min")]

features <- merge(sites, diso,
                  by.x=c("ACC_ID", "MOD_RSD"),
                  by.y=c("acc", "position"), all.x=TRUE, all.y=FALSE)
features <- merge(features, foldx,
                  by.x=c("ACC_ID", "MOD_RSD"),
                  by.y=c("acc", "position"), all.x=TRUE, all.y=FALSE)
features <- merge(features, hotspots,
                  by.x=c("ACC_ID", "MOD_RSD"),
                  by.y=c("acc", "position"), all.x=TRUE, all.y=FALSE)
features <- merge(features, prot_length,
                  by.x=c("ACC_ID"), by.y=c("acc"),
                  all.x=TRUE, all.y=FALSE)
features$pos_in_domain <- apply(features, 1,
                                  function(row){
                                      kinase <- row[["GENE"]]
                                      pos <- as.integer(row[["MOD_RSD"]])
                                      if (is.na(pos))
                                          print(row)
                                      if (!(kinase %in% pfam$protein)){
                                          return(NA)
                                      }
                                      domains <- subset(pfam,
                                                        protein==kinase &
                                                        domain %in% c("Pkinase",
                                                                      "Pkinase_Tyr"))
                                      if (nrow(domains) == 0)
                                          return(NA)
                                      pos_in_domain <-
                                          apply(domains, 1,
                                                function(domain){
                                                    dom_start <- as.integer(domain[["start"]])
                                                    dom_end <- as.integer(domain[["end"]])
                                                    get_pos_in_domain(pos, dom_start, dom_end)
                                                })
                                      in_dom <- which(!is.na(pos_in_domain) &
                                                      pos_in_domain >= 0 &
                                                      pos_in_domain <= 1)
                                      if (length(in_dom)>0)
                                          return(pos_in_domain[in_dom])
                                      dist_from_dom <- sapply(pos_in_domain,
                                                              function(x){
                                                                  if (is.na(x))
                                                                      return(NA)
                                                                  if (x < 0){
                                                                      return(abs(x))
                                                                  }else{
                                                                      return(x-1)
                                                                  }
                                                              })
                                      if (all(is.na(dist_from_dom)))
                                          return(NA)
                                      return(pos_in_domain[which.min(dist_from_dom)])
                                  })
features$pfam <- as.factor(apply(features, 1,
                                  function(row){
                                      kinase <- row[["GENE"]]
                                      pos <- as.integer(row[["MOD_RSD"]])
                                      if (!(kinase %in% pfam$protein)){
                                          return(NA)
                                      }
                                      domains <- subset(pfam,
                                                        protein==kinase &
                                                        domain %in% c("Pkinase",
                                                                      "Pkinase_Tyr"))
                                      if (nrow(domains) == 0)
                                          return(NA)
                                      pos_in_domain <-
                                          apply(domains, 1,
                                                function(domain){
                                                    dom_start <- as.integer(domain[["start"]])
                                                    dom_end <- as.integer(domain[["end"]])
                                                    return(pos >= dom_start && pos <= dom_end)
                                                })
                                      if (any(pos_in_domain)){
                                          return(domains[which(pos_in_domain)[1], "domain"])
                                      }
                                      return(NA)
                                  }))
features$is_tyr_kin <- apply(features, 1,
                               function(row){
                                   kinase <- row[["GENE"]]
                                   if (!(kinase %in% pfam$protein)){
                                       return(NA)
                                   }
                                   domains <- subset(pfam, protein==kinase)
                                   return ("Pkinase_Tyr" %in% domains$domain)
                               })
features$pos_perc <- sapply(1:nrow(features),
                              function(i){
                                  as.integer(features[i, "MOD_RSD"])/features[i, "prot_length"]
                              })

features <- features[, c(3, 2, 4:ncol(features))]

names(features)[1:2] <- c("kinase", "position")

good_cases <- which(apply(features[3:ncol(features)], 1,
                          function(row){
                              length(which(!is.na(row)))
                          }) >= 2)

final_features <- features[good_cases,]
write.table(final_features, "data/phosfun-features.tsv", quote=FALSE,
            sep="\t", row.names=FALSE, col.names=TRUE)
