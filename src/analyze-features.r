edge.feats = read.delim("kinase-merged-feats-clean-hinegs.tsv")
sign.feats = read.delim("sign-feats.tsv")
site.feats = read.delim("reg-site-sign-feats.tsv")

tissue.spec = abs(edge.feats$tissue.sel1 - edge.feats$tissue.sel2)
edge.feats = edge.feats[,-grep("tissue.sel1",colnames(edge.feats))]
edge.feats = cbind(edge.feats,tissue.spec)


feats = list(edge.feats,sign.feats,site.feats)
names(feats) = c("edges","signs","sites")

val.set = read.delim("validation-set-omnipath.tsv",header= F)
sign.set = read.delim("sign-validation-set-omnipath.tsv",header=F)
site.set = read.delim("site-sign-validation-set.tsv",header = F)

val.sets = list(val.set,sign.set,site.set)

for(i in 1:length(val.sets)){
	tbl = feats[[i]]
	if(i != 3){
		tbl = tbl[-which(tbl[,1] == tbl[,2]),]
	}
	val = val.sets[[i]]
	if(i != 1){
		val.act = val[which(val[,3] == "activates"),]
		val.inh = val[which(val[,3] == "inhibits"),]

		act.pairs = paste(val.act[,1],val.act[,2])
		inh.pairs = paste(val.inh[,1],val.inh[,2])
		
		tbl.pairs = paste(tbl[,1],tbl[,2])
		tbl.act = tbl[tbl.pairs %in% act.pairs,]
		tbl.inh = tbl[tbl.pairs %in% inh.pairs,]
		
		for(j in 3:ncol(tbl)){

			if(!(is.numeric(tbl[,j]))){
                                next
                        }

			if(median(tbl.act[,j],na.rm = T) < median(tbl.inh[,j],na.rm = T)){
				alt = "less"
			}else{
				alt = "greater"
			}	
			stats = wilcox.test(tbl.act[,j],tbl.inh[,j],alternative= alt)	
			row = data.frame(set= names(feats)[i],feat = colnames(tbl)[j],pval = stats[[3]],W = stats[[1]])
			write.table(row,"feature-stats.tsv",quote=F,sep= "\t",row.names = F,append = T)
		}
	}else{
		pos.pairs = paste(val[,1],val[,2])
		tbl.pairs = paste(tbl[,1],tbl[,2])
		neg.pairs = setdiff(tbl.pairs,pos.pairs)

		tbl.pos = tbl[tbl.pairs %in% pos.pairs,]
		tbl.neg = tbl[tbl.pairs %in% neg.pairs,]

		for(j in 3:ncol(tbl)){
			if(!(is.numeric(tbl[,j]))){
				next
			}

			 if(median(tbl.pos[,j],na.rm = T) < median(tbl.neg[,j],na.rm = T)){
                                alt = "less"
                        }else{
                                alt = "greater"
                        }

			stats = wilcox.test(tbl.pos[,j],tbl.neg[,j],alternative= alt)
			row = data.frame(set= names(feats)[i],feat = colnames(tbl)[j],pval = stats[[3]],W = stats[[1]])
			write.table(row,"feature-stats.tsv",quote=F,sep= "\t",row.names = F,append = T)	
		}
	}	
}	


