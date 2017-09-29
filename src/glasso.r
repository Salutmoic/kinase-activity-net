library("glasso")

grlasso <-function(cor){

#input : cor is p * p weigthed, intergrated correlation matrix 
#outpur: wi is the inverse covariance matrix while w is the estimated covariance matrix
# I put the regularization parameter to 0.06 which should remove ~ 60% of the edges

g =glasso(cor,p = 0.06)



return( list(g$wi,g$w))
}



