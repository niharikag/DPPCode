#Block diagonal model, blocks of highly correlated predictors
library(glmnet)
library(psych)
library(MASS)
library(gplots)

#some initialization
set.seed(09162016)
q = .01
n=100
p=50
err.sigma = 1.5
S0 = seq(1:15) #true active set
true.beta <- c(rep(2,15),rep(0,35))


gen.data.linDep = function(p, n, true.beta, err.sigma)
{
  x.train=matrix(0,n,p)
  x.valid =matrix(0,n,p)
  
  for(i in 1:50)
  {
    x.train[,i] = rnorm(n)
    x.valid[,i] = rnorm(n)
  }
  
  
  x.train[,1] = x.train[,2]+x.train[,3]+x.train[,4]+x.train[,5] + rnorm(n,0,q)
  x.train[,6] = x.train[,7]+x.train[,8]+x.train[,9]+x.train[,10] + rnorm(n,0,q)
  x.train[,11] = x.train[,12]+x.train[,13]+x.train[,14]+x.train[,15] + rnorm(n,0,q)
  
  x.valid[,1] = x.valid[,2]+x.valid[,3]+x.valid[,4]+x.valid[,5] + rnorm(n,0,q)
  x.valid[,6] = x.valid[,7]+x.valid[,8]+x.valid[,9]+x.valid[,10] + rnorm(n,0,q)
  x.valid[,11] = x.valid[,12]+x.valid[,13]+x.valid[,14]+x.valid[,15] + rnorm(n,0,q)
  
  y.train = drop (x.train %*% true.beta) + err.sigma*rnorm (n)
  
  y.valid = drop (x.valid %*% true.beta) + err.sigma*rnorm (n)
  
  result = list(x.train = x.train,  y.train = y.train, x.valid = x.valid, y.valid=y.valid)
  return(result)
}

iterate.linDep = function(iter)
{
  predErr.Lasso = rep(0,iter)
  predErr.LDSRE = rep(0,iter)
  predErr.enet = rep(0,iter)
  
  TPR.Lasso = rep(0,iter)
  TPR.LDSRE = rep(0,iter)
  TPR.enet = rep(0,iter)
  
  for(index in 1:iter)
  {
    d = gen.data.linDep(p, n, true.beta, err.sigma)
    
    result = fit.Lasso(d)
    predErr.Lasso[index] = result$predErr
    S = which(result$beta != 0)
    TPR.Lasso[index] = length(intersect(S, S0))/30
    
    result = fit.LDSRE(d)
    predErr.LDSRE[index] = result$predErr
    S = which(result$beta != 0)
    TPR.LDSRE[index] = length(intersect(S, S0))/30
    
    result = fit.enet(d)
    predErr.enet[index] = result$predErr
    S = which(result$beta != 0)
    TPR.enet[index] = length(intersect(S, S0))/30
    
  }
  
  print("Lasso")
  print(min(predErr.Lasso))
  print(sd(predErr.Lasso))
  ind = which.min(predErr.Lasso)
  print(TPR.Lasso[ind])
  
  print("LDSRE")
  print(min(predErr.LDSRE))
  print(sd(predErr.LDSRE))
  ind = which.min(predErr.LDSRE)
  print(TPR.LDSRE[ind])
  
  print("ENET")
  print(min(predErr.enet))
  print(sd(predErr.enet))
  ind = which.min(predErr.enet)
  print(TPR.enet[ind])
  
  
}