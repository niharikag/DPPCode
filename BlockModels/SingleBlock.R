library(glmnet)
library(MASS)
#some initialization
set.seed(09092016)
q = .01
p=100
n=100
err.sigma = 1.5
s = 30 #number of active features
S0 = seq(1:30) #true active set
true.beta = c(rep(2,s),rep(0,p-s))

########generate data
gen.data.singleBlock = function(p, n, true.beta, err.sigma)
{
  x.train=matrix(0,n,p)
  x.valid =matrix(0,n,p)
  
  Z1 = rnorm(n)
  
  #first 10 features are correlated
  for(i in 1:10)
  {
    x.train[,i] = Z1 + rnorm(n,0,q)
    x.valid[,i]= Z1 + rnorm(n,0,q)
  }
  
  for(i in 11:100)
  {
    x.train[,i] = rnorm(n)
    x.valid[,i] = rnorm(n)
  }

  y.train = drop (x.train %*% true.beta) + err.sigma*rnorm (n)
  y.valid = drop (x.valid %*% true.beta) + err.sigma*rnorm (n)
  
  result = list(x.train = x.train,  y.train = y.train, x.valid = x.valid, y.valid=y.valid)
  return(result)
}


iterate.singleBlock = function(iter)
{
  predErr.Lasso = rep(0,iter)
  predErr.LDSRE = rep(0,iter)
  predErr.enet = rep(0,iter)
  
  TPR.Lasso = rep(0,iter)
  TPR.LDSRE = rep(0,iter)
  TPR.enet = rep(0,iter)
  
  for(index in 1:iter)
  {
    d = gen.data.singleBlock(p, n, true.beta, err.sigma)
    
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

