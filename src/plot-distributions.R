#known contains list of alledges found in the reactome dataset.
#all is a sample of length length(known) obtained from the list of all edges
# u is the clustering method used and file is the datasets' filename
# s is the name of the statistical method
# The function plots the score distribution of the known edges and the sample from all edges.


plot_dist<-function(known,all,s,u,file){

x = density(known[!is.na(known)])
y = density(sample(all[!is.na(all)],length(known)))

pdf(paste("/home/borgthor/Kinase_activities/Results/FunchiSq/",file,s, u,"matrix clustering .pdf",sep = " "))

hy =max(c(range(y[[2]])[2],range(x[[2]])[2]))
hx = max(c(range(y[[1]])[2],range(x[[1]])[2]))
lx = min(c(range(y[[1]])[1],range(x[[1]])[1]))
plot(x,main = "Comparison between known and unknown edges ",xlab = paste(s,"score"),xlim = c(lx,hx),ylim = c(0,hy + 0.2*hy),col = "red")

lines(y,col = "green")
legend(hx*0.6,hy*0.6,c("known edges","random edges"), text.col = c("red","green"), bty = "n")
dev.off()

}
