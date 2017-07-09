#to check GIC for correlated variables
#adding variables correlated to active variables  
library(glmnet)
library(MASS)
library(lars)

set.seed(10102016)
q = .001
n=400
p=500
err.sigma = 1
s = 100 #number of active features
true.beta = c(rep(1,s),rep(0,p-s))
n.valid = 500
S0 = c(1:s)
########generate orthonormal data
generate.data.corr = function(p, n, true.beta, err.sigma)
{
  x.train=matrix(0,n,p)
  x.valid =matrix(0,n.valid,p)
  
  #independent features
  i=1
  while(i <= p)
  {
    x.train[,i] = rnorm(n)
    x.valid[,i] = rnorm(n.valid)
    if(i<100)
    {
      for(j in 1:4)
      {
        x.train[,i+j] = x.train[,i] + rnorm(n,0,q)
        x.valid[,i+j] = x.valid[,i] +  rnorm(n,0,q)
      }
      i = i+5
    }
    else{i=i+1}
  }
  
  y.train = drop (x.train %*% true.beta) + err.sigma*rnorm (n)
  y.valid = drop (x.valid %*% true.beta) + err.sigma*rnorm (n.valid)
  
  result = list(x.train = x.train,  y.train = y.train, x.valid = x.valid, y.valid=y.valid)
  return(result)
}

iterate.corr = function(iter)
{
  predErr.Lasso = rep(0,iter)
  predErr.LDSRE = rep(0,iter)
  predErr.enet = rep(0,iter)
  
  TPR.Lasso = rep(0,iter)
  TPR.LDSRE = rep(0,iter)
  TPR.enet = rep(0,iter)
  
  for(index in 1:iter)
  {
    d = generate.data.corr(p, n, true.beta, err.sigma)
    
    result = fit.Lasso(d)
    predErr.Lasso[index] = result$predErr
    S = which(result$beta != 0)
    TPR.Lasso[index] = length(intersect(S, S0))/100
    
    result = fit.LDSRE(d)
    predErr.LDSRE[index] = result$predErr
    S = which(result$beta != 0)
    TPR.LDSRE[index] = length(intersect(S, S0))/100
    
    result = fit.enet(d)
    predErr.enet[index] = result$predErr
    S = which(result$beta != 0)
    TPR.enet[index] = length(intersect(S, S0))/100
    
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

