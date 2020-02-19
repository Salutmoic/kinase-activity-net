library(igraph)
library(ggplot2)
library(RColorBrewer)
library(easyGgplot2)
options(error=traceback)

# argv = [edge_cutoff, weighted = T OR F]
argv = commandArgs(TRUE)

cutoff = as.numeric(argv[1])
wtd = as.logical(argv[2])
#Data
#Kinase-kinase network, row: kin1	kin2	edge_probability
net = read.delim("../out/merged-predictor.tsv")
#Validation set kin1	kin2	
val.set = read.delim("../data/validation-set-omnipath.tsv",header = F)

#Results from ROTS analysis, Both MEK and PIK  

#PIK interactions added, as PIK is a lipid kinase
pairs = read.delim("data/PI3K-ints.tsv",header = F)
pairs = pairs[-which(pairs[,2] == "AKT3"),]

all.psite.mek = read.delim(paste("out/MEK-all-psites-8.tsv",sep=""))
all.psite.pik = read.delim("out/PIK-all-psites-8.tsv")

#Kinase substrate data from phosphosite plus

ksub = read.delim("data/kin-sub-psiteplus.tsv")
#ksub = read.delim("../data/psiteplus-kinase-substrates.tsv")
ksub = ksub[,c("GENE","SUB_GENE","position")]

#Edges from the haruna et al study/
vitro = read.delim("haruna-full-int-revised.tsv")
vitro = vitro[,c("kin_gene","gene","Position")]
vitro[,3] = gsub("S|Y|T","",vitro[,3])

#In vivo data
vivo = read.delim("cutillas-pairs.tsv")
vivo[,3] = gsub("Y|T|S","",vivo[,3])

#network processing
net = net[,1:3]
net = net[which(net[,3] >= cutoff),]

val.pairs = paste(val.set[,1],val.set[,2]) 

pairs = pairs[which(pairs[,2] %in% c(as.character(net[,1]),as.character(net[,2]))),]
pik.net = data.frame(prot1 = pairs[,1],prot2= pairs[,2],bart.pred.mean=1)

#Downregulated phosphosites added to list of length 3, 1 slot for each kinase
all.psite.mek[,2] = gsub("Y|T|S","",all.psite.mek[,2])
down.psite = all.psite.mek[which(all.psite.mek$adj.pval < 0.1 & all.psite.mek$logfc < -1),] 

print(dim(down.psite))

down.sites = list()
down.sites[["MAP2K1"]] = paste(down.psite[,1],down.psite[,2])
down.sites[["MAP2K2"]] = paste(down.psite[,1],down.psite[,2])
print(length(unique(down.sites[["MAP2K1"]])))

all.sites = paste(all.psite.mek[,1],all.psite.mek[,2])
 
all.psite.pik[,2] = gsub("Y|T|S","",all.psite.pik[,2])
down.psite = all.psite.pik[which(all.psite.pik$adj.pval < 0.1 & all.psite.pik$logfc < -1),] 
print(dim(down.psite))

down.sites[["PIK3"]] = paste(down.psite[,1],down.psite[,2])
print(length(unique(down.sites[["PIK3"]])))


#roots= nodes from which the distance to psites are calculated
roots=c("MAP2K1","MAP2K2","PIK3")

#######################################
#####Now add substrates PsitePlus######
#######################################

vivo.pairs = paste(vivo[,1],vivo[,2],vivo[,3])
vitro.pairs = paste(vitro[,1],vitro[,2],vitro[,3])

vivo.pairs.prot = paste(vivo[,1],vivo[,2])
vitro.pairs.prot = paste(vitro[,1],vitro[,2])

vivo = vivo[-which(duplicated(vivo.pairs)),]
vivo.pairs = vivo.pairs[-which(duplicated(vivo.pairs))]

set.name = "unknown"
if(length(argv) == 3){
	net.pairs = paste(net[,1],net[,2])
	net = net[which(net.pairs %in% union(union(val.pairs,vivo.pairs.prot),vitro.pairs.prot)),]		
	set.name = "known"
}

#take the intersaction of the vivo and vitro data set and add to the validation set. 
vivo= vivo[which(vivo.pairs %in% vitro.pairs),]
colnames(vivo) = colnames(ksub)

ksub = rbind(ksub,vivo)
ksub.sites=paste(ksub$SUB_GENE,ksub$position)
		
kinsite= paste(ksub$GENE,ksub.sites)
ksub.sites=ksub.sites[-which(duplicated(kinsite))]
ksub = ksub[-which(duplicated(kinsite)),]
ksub = ksub[which(ksub.sites %in% all.sites),]
ksub = ksub[which(ksub[,1] %in% c(as.character(net[,1]),as.character(net[,2]))),]

print(table(ksub[,1] %in% c(as.character(net[,1]),as.character(net[,2]))))
ksub.sites=paste(ksub$SUB_GENE,ksub$position)

print(table(ksub.sites %in% unlist(down.sites)))
print(ksub[which(ksub.sites %in% unlist(down.sites)),])
print(table(ksub.sites %in% all.sites))
print(length(unique(ksub.sites)))

#rows to be added to the network
add.rows = data.frame(prot1 = ksub$GENE,prot2= ksub.sites,bart.pred.mean =1)

#Calculate and compare distances between down regulated sites and the rest

if(length(which(ksub.sites %in% unlist(down.sites))) > 0){

down.dist= c()
all.dist= c()
j=0
#calculate down regulated dists

for(root in roots){
	print(j)
	#Here I gather phosphosites that ar down-regulated
	d.sites =  unique(ksub.sites[which(ksub.sites %in% down.sites[[root]])])
    print(length(d.sites))	
	for(site in d.sites){

		net.sub = net

		weights = net.sub[,3]
        weights = (weights-min(weights))/(max(weights)-min(weights))
		weights[which(weights ==1)] = 0.99
		net.sub= rbind(net.sub,add.rows)
		net.sub = rbind(net.sub,pik.net)
		weights = c(weights,add.rows[,3],pik.net[,3])
		weights = 1-weights
		


		if(root == "PIK3"){
			root.rows = grep(root,net.sub[,1])
		}else{
			root.genes = paste(roots[1],roots[2],sep ="|")
			root.rows = grep(root.genes,net.sub[,1])
		}

		if(site %in% net.sub[root.rows,2] ){
			next
		}

		graph = graph_from_edgelist(as.matrix(net.sub[,c(1,2)]), directed = TRUE)

		if(wtd == T){
			srtspths = all_shortest_paths(graph,root,site,mode = "out",weights = weights)$res
			 if(length(srtspths) > 1){
				print("Warning: paths may not the same length edge wise")
                print(srtspths)
            }
		}else{
			srtspths =all_shortest_paths(graph,root,site,mode = "out")$res
		}

		if(length(srtspths)==0){
            next
        }

		for(srp in srtspths){
			path = names(unlist((srp)))
			write.table(path,paste("out/",argv[1],"-",argv[2],"-",set.name,"-limma-reg-sites-PIK3-new.tsv",sep=""),quote=F,sep="\t",row.names = F,append = T)
		}

		if(wtd){
			print("here")
            d= 0
            for(no in 1:(length(path)-1)){
                d = d+weights[which(as.character(net.sub[,1]) == path[no] & as.character(net.sub[,2]) == path[no+1])]
            }
			print(d)
			print(path)
        }else{
            d = length(path) -1

			if(d==1){
                print(path)
                next
            }

            if(root == "PIK3" & d == 2){
				print(path)
                next
            }

        }
		if(d!=0){
			j=j+1
		}
		down.dist = c(down.dist,d)
		print(length(down.dist[which(down.dist>0)])) 
	
	}
}

for(root in roots){
	rest.sites = unique(ksub.sites[-which(ksub.sites %in% down.sites[[root]])])
	for(site in rest.sites){
			
		net.sub = net
			
		weights =net.sub[,3]
        weights = (weights-min(weights))/(max(weights)-min(weights))
		weights[which(weights ==1)] = 0.99
		net.sub = rbind(net.sub,add.rows)
		net.sub = rbind(net.sub,pik.net)
		weights = c(weights,add.rows[,3],pik.net[,3])

		weights = 1-weights

		if(root == "PIK3"){
            root.rows = grep(root,net.sub[,1])
        }else{
            root.genes = paste(roots[1],roots[2],sep ="|")
            root.rows = grep(root.genes,net.sub[,1])
        }

		if(site %in% net.sub[root.rows,2] ){
            next
        }

		graph = graph_from_edgelist(as.matrix(net.sub[,c(1,2)]), directed = TRUE)
	
		if(wtd == T){
			 srtspths = all_shortest_paths(graph,root,site,weights = weights,mode="out")$res
			 if(length(srtspths) > 1){
				 print("Warning: paths may not the same length edge wise")
				print(srtspths)
			}
        }else{
			srtspths =all_shortest_paths(graph,root,site,mode="out")$res
        }
		
		if(length(srtspths)==0){
            next
        }

		for(srp in srtspths){
            path = names(unlist((srp)))
		}
		if(wtd){
			d= 0
			for(no in 1:(length(path)-1)){
				d = d+weights[which(net.sub[,1] == path[no] & net.sub[,2] == path[no+1])]
			}
		}else{
			d = length(path) -1
			
			if(d==1){
				print(path)
				next
			}
			
			if(root == "PIK3" & d == 2){
				next
			}
		}

		if(d==0){
			print(path)
		}
		all.dist = c(all.dist,d)
	}
}

all.dist = all.dist[which(all.dist >0 )]
down.dist = down.dist[which(down.dist >0 )]
print(length(down.dist))
print(j)

print(wilcox.test(down.dist,all.dist,alternative="less"))

pdf(paste("img/",argv[1],"-",argv[2],"-",set.name,"-rots-down-reg-dist-4.pdf",sep=""))
boxplot(down.dist,all.dist,names = c("down-regulated","rest"))
dev.off()

#ggplot

ggplot.data = data.frame(distance= c(down.dist,all.dist),set = c(rep("down-regulated",length(down.dist)),rep("up/un-regulated",length(all.dist))))
print(head(ggplot.data))


pdf(paste("img/",argv[1],"-",argv[2],"-",set.name,"-ggplot-rots-down-reg-dist-4.pdf",sep=""))

p=ggplot(ggplot.data,aes(x = set,y = distance,fill = set))+xlab("set") + geom_boxplot() + theme_classic()+ scale_fill_brewer(palette="PRGn")
print(p)
dev.off()

pdf(paste("img/histo-",argv[1],"-",argv[2],"-",set.name,"-ggplot-rots-down-reg-dist-4.pdf",sep=""))

p=ggplot(ggplot.data,aes(x = distance,fill = set))+xlab("distance") + geom_histogram(aes(y=..density..),position="identity",alpha = 0.5,binwidth=1) + theme_classic()+ scale_fill_manual(values=c("#762A83","#1B7837"))

print(p)
dev.off()

pdf(paste("img/density-",argv[1],"-",argv[2],"-",set.name,"-ggplot-limma-down-reg-dist-4.pdf",sep=""))

p=ggplot(ggplot.data,aes(x = distance,fill = set))+xlab("distance") + geom_density(alpha = 0.5) + theme_classic()+ scale_fill_manual(values=c("#762A83","#1B7837"))

print(p)
dev.off()



}

