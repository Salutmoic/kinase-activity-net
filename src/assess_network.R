argv <- commandArgs(TRUE)
if (length(argv) != 3){
  stop("USAGE: <script> DATA_FILE OUT_FILE METHOD")
}


true.pos = argv[1]
true.neg = argv[2]
network = argv[3]

true.pos.pairs = paste(true.pos[,1],true.pos[,2])
true.false.pairs = paste(true.false[,1],true.false[,2])
network.pairs = paste(network[,1],network[,2])

net.true = network[which(network.pairs %in% true.pos.pairs),]
net.false = network[which(network.pairs %in% true.false.pairs),]

print("difference between true positives and true negatives")

print(wilcox.test(as.numeric(net.true[,3]),as.numeric(net.false[,3])))