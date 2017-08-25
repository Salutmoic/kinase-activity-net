library("glasso")

grlasso <-function(cor,p){

#input : cor is p * p weigthed, intergrated correlation matrix 
#outpur: wi is the inverse covariance matrix while w is the estimated covariance matrix
#  p is the regularization parameter, for instance 0.06 should remove ~ 60% of the edges

g =glasso(cor,p)



return( list(g$wi,g$w))
}



