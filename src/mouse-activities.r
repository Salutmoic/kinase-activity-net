library("reshape2")
library("ksea")
require("gdata")
library(stringr)


find_psites<- function(cond,sites,reg.sites){

ksea = ksea(sites,cond,reg.sites,trial = 1000, significance = TRUE)

return(ksea)
}

map_to_protein <- function(peptides,proteins,modind,proteome){
    psites = c()

    for(i in 1:length(peptides)){
        seq = Proteome[[proteins[i]]][1]

        ind = unlist(gregexpr(peptides[i],seq))
        mods = unlist(modind[[i]])
        for(j in 1:length(mods)){
            aa = substr(seq,ind+mods[j]-(j+1),ind+mods[j]-(j+1))
            aa = paste(aa,ind+mods[j]-(j+1),sep = "")
            psites = c(psites,paste(proteins[i],aa,sep = "_"))

        }
    }


    return(psites)
} 


library(seqinr)
Proteome <- read.fasta("http://www.uniprot.org/uniprot/?query=taxonomy:10090&format=fasta", seqtype = "AA", as.string = T)
names(Proteome) <- gsub("sp\\||\\|.*$|tr\\||\\|.*$","",names(Proteome))

data = list("../control_FGF4_MEKinh_corrected.xlsx","../EScellsFGF_stimulation_corrected.xlsx","../FGF4_Srcinh_GSK3inh.xlsx")

ksub = read.delim("../data/mouse-kinome.tsv")
kinases = read.delim("../data/mouse-kinases.tsv")

reg.sites = paste(ksub$SUB_ACC_ID,ksub$SUB_MOD_RSD, sep = "_")

pepts = list()
conds = c()
vals= list()

for(tbls in data){

    tbl = read.xls (tbls, sheet = 1, header = TRUE)
    if(ncol(tbl) > 13){
        conds = c(conds,colnames(tbl)[13],colnames(tbl)[14])
        vals[[length(vals)+1]] = cbind(tbl[,13],tbl[,14])
    }else{
        conds = c(conds,colnames(tbl)[12],colnames(tbl)[13])    
        vals[[length(vals)+1]] = cbind(tbl[,12],tbl[,13])

    }
    proteins = tbl$Proteins
    proteins = str_split_fixed(proteins, "\\|", 3)[,2]

    peptides = tbl$Sequence
    modifications = tbl$Modified.Sequence
    modind = list()
    for(mods in modifications){
        modind[[length(modind) + 1]] = gregexpr("p",mods,perl = T)
    }
    pepts[[length(pepts)+1]] =  map_to_protein(peptides,proteins,modind)

}

#names(vals) = conds
kin.act = list()
for(kinase in kinases[,1]){
    regulons = reg.sites[which(ksub[,1] == kinase)]
    kin.act[[kinase]] = c()
    for(i in 1:length(pepts)){
        for(j in 1:2){
            kin.act[[kinase]] = c(kin.act[[kinase]],ksea(pepts[[i]],vals[[i]][,j],regulons))

        }

    }
print(kinase)
}

save(kin.act,file ="kinase-activities-mice.RData")

