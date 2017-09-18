checkifknown <- function(a,b,ksub,kinases){

t = 0

c = ksub[toupper(ksub[,1]) ==a,]

if(b %in% toupper(c[,2])){
t = 1
}

return(t)

}


find.indirect <-function(kinase,sub,ksub,d,kinases){
#kinase is a string that contains the name of the kinase while sub is a potential substrate (string)
#ksub is a data frame which contains a list of known kinase -> substrate interaction, d is the allowed distance
#between kinase and sub (if d = 1 then A-> C and C->B is would validate A-> B while A-> C, C->D and D-> B would not) 

ksub.sub = ksub[ksub[,1] == kinase,]


substrates = as.character(ksub.sub[,2])
substrates.next = c()

i = 0

while(i < d){

for(j in substrates){
#here I see if sub interacts with any of the substrates of kinase

ksub.sub = ksub[which(ksub[,1] == j),]

#if d > 1 I have to take a look and see if substrate's substrates interacts with sub
#therefore the substrates of the subsrates in substrates are collected.
substrates.next = c(substrates.next,as.character(ksub.sub[,2]))

if(sub %in%  ksub.sub[,2] & !(sub %in% kinases) ){

return(1)
}

}
print(substrates.next)
#substrates.next will be considered for the next foor loop
substrates = substrates.next

i = i+1
}

return(0)

}

ksub = rbind(c("A","B"),c("B","C"),c("B","D"),c("D","E"))
find.indirect("A","E",data.frame(ksub),3)

proport_known <- function(edges,ksub){
#takes identified edges and returns the proportion of empirically known edges
# that have been identified

known = paste(ksub[,1],ksub[,2])

inferred = paset(edges[,1],edges[,2])

n = itersect(inferred, known)

return(c(n,n/length(known)))



}



known = c()
all = c()
for(j in 1:nrow(scores)){

t = checkifknown(scores[j,1],scores[j,2],ksub)

if(t == 1){

known = c(known,as.numeric(scores[j,3]))
}

all = c(all,as.numeric(scores[j,3]))


}
