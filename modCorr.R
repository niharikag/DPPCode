#moderately correlated groups,not highly correlated
library(glmnet)
library(psych)
library(MASS)
#some initialization
set.seed(09112016)
q = .01
set.seed (333)
n=100
p=100
err.sigma = 3
s0 = seq(1:20) #true active set
seq1 = seq(.1,2,.1)

true.beta <- c(rep(1,20),rep(0,p-20))
true.beta[1:20] = sample(seq(.1,2,.1),20, replace=FALSE)    
true.beta[1:20] = rep(1,20)
varMat = matrix(1,10,10)

for(i in 1:10)
  for(j in 1:10)
  {
    if(i!=j)
      varMat[i,j] = .6
  }
#create block diagonal matrix
blockMat1 = superMatrix(rep(list(varMat),1))
blockMat1 = varMat

for(i in 1:10)
  for(j in 1:10)
  {
    if(i!=j)
      varMat[i,j] = .9
  }
#create block diagonal matrix
blockMat2 = superMatrix(rep(list(varMat),9))

blockMat = superMatrix(blockMat1,blockMat2)
generate.data = function(p, n, true.beta, err.sigma)
{
  x.train = mvrnorm(n,rep(0,p),blockMat)
  y.train = drop (x.train %*% true.beta) + err.sigma*rnorm (n)
  
  x.valid = mvrnorm(n,rep(0,p), blockMat)
  y.valid = drop (x.valid %*% true.beta) + err.sigma*rnorm (n)
  
  result = list(x.train = x.train,  y.train = y.train, x.valid = x.valid, y.valid=y.valid)
  return(result)
}

compareMethods = function(observations = 100)
{
  n  = observations
  #true.beta[1:s] = sample(c(1,2,3), s, replace = TRUE)
  d = generate.data(p, n = n, true.beta, err.sigma)
  x = d$x.train
  y = d$y.train
  
  #sig = t(x)%*%x
  ##heatmap(sig,Rowv=NA)
  #heatmap(sig,   Rowv=NA, keep.dendro = FALSE,Colv = "Rowv")   
  
  #x = (scale(d$x.train)) #standardize the predictors)
  #y = scale(d$y.train,scale=FALSE) #centralize the response
  
  #Find lamda_max
  temp = rep(0,p)
  for(i in 1:p){
    temp[i] = abs(t(x[,i])%*%y)
  }
  lambda_max  = max(temp)/n
  print(lambda_max)
  
  fit.LDSRE(d, n, lambda_max)
  fit.enet(d, n, lambda_max)
  
}

fit.LDSRE = function(d, n, lambda_max)
{
  #d = generate.data(p, n, true.beta, err.sigma)
  x = d$x.train
  y = d$y.train
  x.valid = d$x.valid
  y.valid = d$y.valid
  
  res.lasso = glmnet(x, y, nlambda = 100, alpha = 1, intercept = FALSE)
  y.estimate = predict.glmnet(res.lasso, x.valid)
  gridLen = dim(y.estimate)[2]
  
  #compute residuals
  y.residuals = y.valid - y.estimate
  y.residuals = (y.residuals)^2
  y.residuals = colMeans(y.residuals)
  
  #minimum Prediction Error
  indx = which.min(y.residuals)
  predErr =min(y.residuals)
  lambda.opt = res.lasso$lambda[indx]
  
  
  res.lasso = glmnet(x,y,lambda = lambda.opt, alpha = 1, intercept = FALSE)
  s.lasso = (which(res.lasso$beta!=0))
  
  #compute the lasso in-sample prediction error
  y.hat = predict.glmnet(res.lasso, x, lambda = lambda.opt, alpha = 1, intercept = FALSE)
  
  print("Lasso Error")
  print(predErr)
  print("Lasso selection")
  print(s.lasso)
  
  #Step2(a): compute theta, for each lamda
  #theta = (y - (x%*%res.lasso$beta))
  theta = (y - y.hat)
  
  dual.feasible = rep(0, p)
  
  for(j in 1:p){
    dual.feasible[j] = (t(x[,j])%*%theta)
  }
  #take the absolute of dual feasible
  dual.feasible = round(abs(dual.feasible)/n,1)
  
  #find dual active variables
  s.dual = which(dual.feasible >= lambda.opt)
  
  #print the L1 norm
  print(sum(res.lasso$beta))
  
    res.ridge = glmnet(x[,s.dual], y, lambda = seq(0, lambda_max, .1), alpha = 0, intercept = FALSE)
    y.estimate = predict.glmnet(res.ridge, x.valid[,s.dual])
    #compute residuals
    y.residuals = y.valid - y.estimate
    y.residuals = (y.residuals)^2
    y.residuals = colMeans(y.residuals)
    
    #minimum Prediction Error
    indx = which.min(y.residuals)
    predErr =min(y.residuals)
    
    print(sum(res.ridge$beta[,indx]))
    
    print("Ridge Error")
    print(predErr)
    print("Ridge selection")
    print(s.dual)
}


fit.enet = function(d, n, lambda_max)
{
  #d = generate.data(p, n, true.beta, err.sigma)
  x = d$x.train
  y = d$y.train
  x.valid = d$x.valid
  y.valid = d$y.valid
  
  #res.enet = glmnet(x, y, lambda = lam_range, alpha = 0.5, intercept = FALSE)
  res.enet = glmnet(x, y, lambda = seq(0, lambda_max, .1), alpha = 0.5, intercept = FALSE)
  y.estimate = predict.glmnet(res.enet, x.valid)
  
  
  y.residuals = y.valid - y.estimate
  y.residuals = (y.residuals)^2
  y.residuals = colMeans(y.residuals)
  
  #minimum Prediction Error
  indx = which.min(y.residuals)
  predErr = min(y.residuals)
  s.enet = which(res.enet$beta[,indx]!=0)
  print("Enet Error")
  print(predErr)
  print("Enet selection")
  print(s.enet)
}

