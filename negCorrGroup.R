#Block diagonal model, blocks of highly correlated predictors
library(glmnet)
library(psych)
library(MASS)
library(gplots)

#some initialization
set.seed(09162016)
q = .01
n=1000
p=100
err.sigma = 3
s0 = seq(1:30) #true active set
true.beta <- c(rep(3,30),rep(0,70))
true.beta[seq(1:5)] = -3
true.beta[seq(11:15)] = -3
true.beta[seq(21:25)] = -3

varMat = matrix(1,10,10)

for(i in 1:10)
  for(j in 1:10)
  {
    if(i!=j)
      varMat[i,j] = .9
  }
#create block diagonal matrix
blockMat = superMatrix(rep(list(varMat), 10))

gen.data.BlockDiagNeg = function(p, n, true.beta, err.sigma)
{
  x.train = mvrnorm(n,rep(0,p),blockMat)
  
  for(i in 1:5)
  {
    x.train[,i] =  - x.train[,i]
  }
  
  for(i in 11:15)
  {
    x.train[,i] =  - x.train[,i]
  }
  for(i in 21:25)
  {
    x.train[,i] =  - x.train[,i]
  }
  
  y.train = drop (x.train %*% true.beta) + err.sigma*rnorm (n)
  
  x.valid = mvrnorm(n,rep(0,p), blockMat)
  for(i in 1:5)
  {
    x.valid[,i] =  - x.valid[,i]
  }
  
  for(i in 11:15)
  {
    x.valid[,i] =  - x.valid[,i]
  }
  for(i in 21:25)
  {
    x.valid[,i] =  - x.valid[,i]
  }
  y.valid = drop (x.valid %*% true.beta) + err.sigma*rnorm (n)
  
  result = list(x.train = x.train,  y.train = y.train, x.valid = x.valid, y.valid=y.valid)
  return(result)
}

gen.heatmap = function(d)
{
  d = gen.data.BlockDiagNeg(p, n, true.beta, err.sigma)
  x = d$x.train
  sig = t(x)%*%x/n
  #heatmap(sig,   Rowv=NA, keep.dendro = FALSE,Colv = "Rowv") 
  heatmap.2(sig, Rowv=NA, Colv=NA, dendrogram="none", col=heat.colors(64),
            trace="none", labRow=NA, labCol=NA, margin=c(6, 6) )
  
}


