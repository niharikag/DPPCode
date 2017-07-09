#get Enet estimate using different approaches
library(glmnet)
library(MASS)
library(lars)

set.seed(27092016)
q = .001
n=100
p=10
err.sigma = 1
s = 5 #number of active features
true.beta = c(rep(1,s),rep(0,p-s))
n.valid = 1000

########generate data
generate.data = function(p, n, true.beta, err.sigma)
{
  x.train=matrix(0,n,p)
  x.valid =matrix(0,n.valid,p)
  
  Z1 = rnorm(n)
 
  #correlated features
  for(i in 1:5)
  {
    x.train[,i] = Z1 + rnorm(n,0,q)
    x.valid[,i] = Z1 + rnorm(n.valid,0,q)
    
  }
  
  #independent features
  for(i in 6:p)
  {
    x.train[,i] = rnorm(n)
    x.valid[,i]=rnorm(n.valid)
  }
  
  y.train = drop (x.train %*% true.beta) + err.sigma*rnorm (n)
  y.valid = drop (x.valid %*% true.beta) + err.sigma*rnorm (n.valid)
  
  result = list(x.train = x.train,  y.train = y.train, x.valid = x.valid, y.valid=y.valid)
  return(result)
}

d = generate.data(p, n, true.beta, err.sigma)

x = d$x.train
y = d$y.train

x = scale(x)
y = scale(y)

#one set of lambda, alpha
lambda_1 = 1
lambda_2 = 3

alpha = lambda_1/(lambda_1+lambda_2)
lambda = lambda_1+lambda_2

#elastic net, using direct glmnet
#res.enet = glmnet(x, y, alpha = alpha, lambda = lambda/n)
res.enet = glmnet(x, y, alpha = alpha, intercept = FALSE, standardize = FALSE)
enet.beta = coef(res.enet, lambda/n)

#compute enet applying Lasso on the modified design matrix
dg = diag(sqrt(lambda_2),p) 
x1 = rbind(x, dg) #Augment the design matrix with d
x1 = x1/sqrt(1+lambda_2)
y1 =  c(y, rep(0,p)) #Augment the response vector with 0
lmd = lambda_1/sqrt(1+lambda_2)
res.lasso =  lars(x1, y1, type = "lasso", intercept = FALSE, normalize = FALSE)
lasso.beta = coef(res.lasso, s=lmd, mode="lambda")

#Compare both the estimates
enet.beta
#lasso.beta*sqrt((1+lambda_2))
lasso.beta/sqrt((1+lambda_2))#Naive elastic net
lasso.beta*sqrt((1+lambda_2))#Enet solution

#to test both are giving the same predicted fir/resposne, so LARS can be used
#with scaled X and Y , since LARS response is more interpretable
y.hat = predict(res.enet, x, s= lambda/n, type="response")
y1.hat =predict.lars(res.lasso, x1, s= lmd, type="fit", mode="lambda")

#predict out of sample
d.test = generate.data(p, n, true.beta, err.sigma)
x.test = d.test$x.train
y.test = d.test$y.train
y.hat = predict(res.enet, x.test, s= lambda/n, type="response")
y1.hat = predict.lars(res.lasso, x.test/sqrt(1+lambda_2), s= lmd, type="fit", mode="lambda")
y1.hat = y1.hat$fit*sqrt(1+lambda_2)
mean((y.test-y1.hat)^2)
mean((y.test-y.hat)^2)

RRandOLS = function(){
  
  #try ridge, and augmented lasso
  d = generate.data(p, n, true.beta, err.sigma)
  
  x = d$x.train
  y = d$y.train
  
  x = scale(x)
  y = scale(y)
  
  #one set of lambda, alpha
  lambda_1 = 1
  lambda_2 = 4
  
  #elastic ridge, using direct glmnet
  #res.enet = glmnet(x, y, alpha = 0, lambda = lambda_2/n, intercept = FALSE, standardize=FALSE)
  res.enet = glmnet(x, y, alpha = 0, intercept = FALSE, standardize=FALSE)
  enet.beta = coef(res.enet, s = lambda_2/n)
  #compute ridge applying Lasso on the modified design matrix OLS
  dg = diag(sqrt(lambda_2),p) 
  x1 = rbind(x, dg) #Augment the design matrix with d
  y1 =  c(y, rep(0,p)) #Augment the response vector with 0
  #res.rr =  glmnet(x1, y1, alpha = 1, lambda = 0, intercept = FALSE, standardize=FALSE)
  res.rr =  glmnet(x1, y1, alpha = 1, intercept = FALSE, standardize=FALSE)
  rr.beta = coef(res.rr, 0)
  
  
  ridgec <- lm.ridge (y ~ x, lambda = lambda_2)
  
  # LM ridge and modified OLS works fine, they are the same
  rr.beta
  ridgec$coef
  #but ENet evaluated through glmnet is not fine
  enet.beta
  
  #using LARS algorithm
  res.lars =  lars(x1, y1, type = "lasso", intercept = FALSE, normalize = FALSE)
  lars.beta = coef(res.lars, s=0, mode="lambda")
}

