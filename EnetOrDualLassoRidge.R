#Comparision of Enet and Dual Lasso Ridge, for fixed lambda_1 and lambda_2
library(glmnet)

set.seed(5102016)
q = .01
n=100
p=10
err.sigma = 1
s = 8 #number of active features
true.beta = c(rep(1,s),rep(0,p-8))
true.beta = rep(1,10)
#Correlated design
########generate data
generate.data = function(p, n, true.beta, err.sigma)
{
  x.train=matrix(0,n,p)
  x.valid =matrix(0,n,p)
  
  Z1 = rnorm(n)
  
  #correlated features
  for(i in 1:10)
  {
    x.train[,i] = Z1 + rnorm(n,0,q)
    x.valid[,i] = Z1 + rnorm(n,0,q)
    
  }
  
  #independent features
  #for(i in 9:p)
  #{
   # x.train[,i] = rnorm(n)
    #x.valid[,i]=rnorm(n)
  #}
  
  y.train = drop (x.train %*% true.beta) + err.sigma*rnorm (n)
  y.valid = drop (x.valid %*% true.beta) + err.sigma*rnorm (n)
  
  result = list(x.train = x.train,  y.train = y.train, x.valid = x.valid, y.valid=y.valid)
  return(result)
}


#Correlated design
########generate data
generate.data2 = function(p, n, true.beta, err.sigma)
{
  x.train=matrix(0,n,p)
  x.valid =matrix(0,n,p)
  
  #independent features
  for(i in 1:p)
  {
    x.train[,i] = rnorm(n)
    x.valid[,i]=rnorm(n)
  }
  
  y.train = drop (x.train %*% true.beta) + err.sigma*rnorm (n)
  y.valid = drop (x.valid %*% true.beta) + err.sigma*rnorm (n)
  
  result = list(x.train = x.train,  y.train = y.train, x.valid = x.valid, y.valid=y.valid)
  return(result)
}

d = generate.data(p, n, true.beta, err.sigma)

x = d$x.train
y = d$y.train

#one set of lambda, alpha
lambda_1 = 2
lambda_2 = .01

alpha = lambda_1/(lambda_1+lambda_2)
lambda = lambda_1+lambda_2

#elastic net, using direct glmnet
res.enet = glmnet(x, y, alpha = alpha, intercept = FALSE)
enet.beta = coef(res.enet, lambda/n)

#Apply ridge directly on active set.
x.red = x #x[,-c(9,10)]
#elastic net, using direct glmnet
res.rr = glmnet(x.red, y, alpha = 0, intercept = FALSE)
rr.beta = coef(res.rr, lambda_2/n)

enet.beta
rr.beta

################################################################
#for orthonormal design
d = generate.data2(p, n, true.beta, err.sigma)

x = d$x.train
y = d$y.train

#one set of lambda, alpha
lambda_1 = 0.5
lambda_2 = 2

alpha = lambda_1/(lambda_1+lambda_2)
lambda = lambda_1+lambda_2

#elastic net, using direct glmnet
res.enet = glmnet(x, y, alpha = alpha, intercept = FALSE)
enet.beta = coef(res.enet, lambda/n)

#Apply ridge directly on active set.
x.red = x #x[,-c(9,10)]
#elastic net, using direct glmnet
res.rr = glmnet(x.red, y, alpha = 0, intercept = FALSE)
rr.beta = coef(res.rr, lambda_2/n)

enet.beta
rr.beta