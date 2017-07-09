#try with independent variables, (NOT correlated), we can show that
# the theta is always less than 1 and primal and dual optimal coincide
library(glmnet)

set.seed(01092016)
q = .01
n=100
p=10
err.sigma = 1
s = 8 #number of active features
true.beta = c(rep(1,s),rep(0,p-8))

########generate data
generate.data.indep = function(p, n, true.beta, err.sigma)
{
  x.train=matrix(0,n,p)
  x.valid =matrix(0,n,p)
  
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


d = generate.data.indep(p, n, true.beta, err.sigma)

x = d$x.train
y = d$y.train

#Find lamda_max
temp = rep(0,p)
for(i in 1:p){
  temp[i] = abs(t(x[,i])%*%y)
}
lambda_max  = max(temp)/n
lambda_max

#value of beta for range (0, lamda_max)
lam_range = seq(0.05, lambda_max-.05, .05)

dual.feasible = rep(0, p)

for(i in 1:length(lam_range))
{
  #Step1: fit lasso for the range of lamda
  res.lasso = glmnet(x,y,lambda = lam_range[i])
  
  res.lasso$beta
  
  #compute the lasso prediction error
  y.hat = predict.glmnet(res.lasso, x, lambda = lam_range[i])
  pred.err.lasso = sum((y-y.hat)^2)/n
  pred.err.lasso
  
  #Step2(a): compute theta, for each lamda
  #theta = (y - (x%*%res.lasso$beta))/lam_range[i]
  theta = (y - y.hat)/lam_range[i]
  
  for(j in 1:p){
    dual.feasible[j] = (abs(t(x[,j])%*%theta)/n)[1][1]
  }
  
  round(dual.feasible,1)
  
  
  #Step2(b): find the feature index satisfying the feasibility at the dual optimal
  s.dual = which(dual.feasible >=0.9)
  print(s.dual)
  res.lasso$beta
  print((t(t(x)%*%theta)) %*% res.lasso$beta/n)
  print(sum(res.lasso$beta))
  if(length(s.dual) >0)
  {
    #Step 3: fit ridge regression for the set of features selected by the Lasso dual
    #Step1: fit lasso for the range of lamda
    res.ridge = glmnet(x[,s.dual],y, alpha = lam_range[i], lambda = 0)
    res.ridge$beta
    
    #compute the prediction error
    y.hat = predict.glmnet(res.ridge, x[,s.dual])
    pred.err.ridge = sum((y-y.hat)^2)/n
    pred.err.ridge
  }
}

