#one small experiment
q = .01
p=3
err.sigma = 1
true.beta = rep(3,3)
n=100
true.beta[3] = -1

genData = function()
{
  x.train=matrix(0,n,p)
  x.valid =matrix(0,n,p)
  
  Z1 = rnorm(n)
  
  for(i in 1:3)
  {
    x.train[,i] = Z1 + rnorm(n,0,q)
    x.valid[,i]= Z1 + rnorm(n,0,q)
  }
  
  #x.train[,3]= -x.train[,3]
  #x.valid[,3]= -x.valid[,3]
  
  x.train[,2]= -x.train[,1]
  x.valid[,2]= -x.valid[,1]
  
  y.train = drop (x.train %*% true.beta) + err.sigma*rnorm (n)
  y.valid = drop (x.valid %*% true.beta) + err.sigma*rnorm (n)
  
  result = list(x.train = x.train,  y.train = y.train, x.valid = x.valid, y.valid=y.valid)
  return(result)
}


LassoMethod = function()
{
  d = genData()
  x = d$x.train
  y = d$y.train
  

  res.lasso = cv.glmnet(x, y, nlambda = 20, alpha = 1, intercept = FALSE)
  min.lambda = res.lasso$lambda.min
  indx = which(res.lasso$lambda == min.lambda)
  res.lasso$glmnet.fit$beta[,min.lambda]
  
  #ridge
  res.ridge = glmnet(x, y, nlambda = 20, alpha = 0, intercept = FALSE)
  res.ridge$beta
  
  res.enet = cv.glmnet(x, y, nlambda = 20, alpha = 0.5, intercept = FALSE)
  res.enet$glmnet.fit$beta
  
}