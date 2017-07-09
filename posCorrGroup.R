library(glmnet)
library(MASS)
#some initialization
set.seed(09092016)
q = .01
p=100
err.sigma = 1
s = 30 #number of active features
true.beta = c(rep(3,s),rep(0,p-s))
true.beta[1:5] = rep(-1,5)

########generate data
gen.data = function(p, n, true.beta, err.sigma)
{
  x.train=matrix(0,n,p)
  x.valid =matrix(0,n,p)
  
  Z1 = rnorm(n)
  Z2 = rnorm(n)
  Z3 = rnorm(n)
  
  #first 10 features are correlated
  for(i in 1:10)
  {
    x.train[,i] = Z1 + rnorm(n,0,q)
    x.valid[,i]= Z1 + rnorm(n,0,q)
  }
  
  
  #make the next five correlated within group and between next group as well
  for(i in 11:20)
  {
    x.train[,i] = Z2 + rnorm(n,0,q)
    x.valid[,i]= Z2 + rnorm(n,0,q)
  }
  
  #make the next five correlated within group and between next group as well
  for(i in 21:30)
  {
    x.train[,i] = Z3 + rnorm(n,0,q)
    x.valid[,i]= Z3 + rnorm(n,0,q)
  }
  
  for(i in 31:100)
  {
    x.train[,i] = rnorm(n)
    x.valid[,i] = rnorm(n)
  }
  
  
  y.train = drop (x.train %*% true.beta) + err.sigma*rnorm (n)
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
  
  #value of beta for range (0, lamda_max)
  from = lambda_max/20
  to = lambda_max - .1
  
  lam_range = seq(from, to, by = ((to - from)/20))
  
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
  
  #for(i in 1:length(lam_range))
  
  #Step1: fit lasso for the range of lamda
  #res.lasso = glmnet(x, y, lambda = lam_range, alpha = 1, intercept = FALSE)
  #y.estimate = predict.glmnet(res.lasso, x.valid, lambda = lam_range, alpha = 0.5, intercept = FALSE)
  #gridLen = length(lam_range)
  
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
  
  #Step1: fit lasso for the range of lamda
  res.lasso = glmnet(x,y,lambda = lambda.opt, alpha = 1, intercept = FALSE)
  
  s.lasso = (which(res.lasso$beta!=0))
  
  #compute the lasso prediction error
  y.hat = predict.glmnet(res.lasso, x, lambda = lambda.opt, alpha = 1, intercept = FALSE)
  #pred.err.lasso = sum((y-y.hat)^2)/n
  
  print("Lasso Error")
  print(predErr)
  print("Lasso selection")
  print(s.lasso)
  
  #Step2(a): compute theta, for each lamda
  #theta = (y - (x%*%res.lasso$beta))
  theta = (y - y.hat)
  
  dual.feasible = rep(0, p)
  
  for(j in 1:p){
    #dual.feasible[j] = ((abs(t(x[,j])%*%theta)/n)[1][1])/lam_range[i]
    dual.feasible[j] = (t(x[,j])%*%theta)/lambda.opt
  }
  #take the absolute of dual feasible
  dual.feasible = round(abs(dual.feasible)/n,1)
  
  #find the minimum slope
  min.theta = min(dual.feasible[s.lasso])
  
  #find dual active variables
  s.dual = which(dual.feasible >= min.theta)
  
  
  #just to check the intercept
  #y.hat - (x%*%res.lasso$beta)
  
  
  #for computing the sum of absolute
  #beta.sum = (t(t(x)%*%theta)) %*% res.lasso$beta/n
  #print(beta.sum)
  #print(sum(res.lasso$beta))
  
  
  if(length(s.dual) >0)
  {
    #res.ridge = glmnet(x[,s.dual], y, nlambda = 100, alpha = 0, intercept = FALSE)
    res.ridge = glmnet(x[,s.dual], y, lambda = seq(0, lambda_max, .1), alpha = 0, intercept = FALSE)
    y.estimate = predict.glmnet(res.ridge, x.valid[,s.dual])
    #compute residuals
    y.residuals = y.valid - y.estimate
    y.residuals = (y.residuals)^2
    y.residuals = colMeans(y.residuals)
    
    #minimum Prediction Error
    indx = which.min(y.residuals)
    predErr =min(y.residuals)
    
    print("Ridge Error")
    print(predErr)
    print("Ridge selection")
    print(s.dual)
  }  
  
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

