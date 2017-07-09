#orthogonal design
library(glmnet)
library(psych)
library(MASS)
library(gplots)

#some initialization
set.seed(09162016)
q = .01
n=100
p=100
err.sigma = 1.5
S0 = seq(1:30) #true active set
true.beta <- c(rep(2,30),rep(0,70))


#create Identity matrix
blockMat = diag(1,100)

gen.data.Ortho = function(p, n, true.beta, err.sigma)
{
  x.train = mvrnorm(n,rep(0,p),blockMat)
  y.train = drop (x.train %*% true.beta) + err.sigma*rnorm (n)
  
  x.valid = mvrnorm(n,rep(0,p), blockMat)
  y.valid = drop (x.valid %*% true.beta) + err.sigma*rnorm (n)
  
  result = list(x.train = x.train,  y.train = y.train, x.valid = x.valid, y.valid=y.valid)
  return(result)
}

iterate.Ortho = function(iter)
{
    predErr.Lasso = rep(0,iter)
    predErr.LDSRE = rep(0,iter)
    predErr.enet = rep(0,iter)
    
    TPR.Lasso = rep(0,iter)
    TPR.LDSRE = rep(0,iter)
    TPR.enet = rep(0,iter)
    
    for(index in 1:iter)
    {
      d = gen.data.Ortho(p, n, true.beta, err.sigma)
      
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

gen.heatmap.Ortho = function(d)
{
  #d = gen.data.Ortho(p, n, true.beta, err.sigma)
  x = d$x.train
  sig = t(x)%*%x/n
  #heatmap(sig,   Rowv=NA, keep.dendro = FALSE,Colv = "Rowv") 
  heatmap.2(sig, Rowv=NA, Colv=NA, dendrogram="none", col=heat.colors(64),
            trace="none", labRow=NA, labCol=NA, margin=c(6, 6) )
  
}


