#get Enet estimate using different approaches
library(glmnet)
library(MASS)
library(gplots)

set.seed(27092016)
q = .01
n=100
p=3
err.sigma = 1
s = 8 #number of active features
true.beta = c(3,3,0)

#generate heat map
gen.heatmap = function(x)
{
  sig = t(x)%*%x/n
  sig = cor(x)
  my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 1000)
  heatmap.2(sig, Rowv=NA, Colv=NA, dendrogram="none", col=my_palette,
            trace="none", labRow=NA, labCol=NA, margin=c(6, 6) )
  
}

########generate data
generate.data = function(p, n, true.beta, err.sigma)
{
  x.train=matrix(0,n,p)
  x.valid =matrix(0,n,p)
  
  Z1 = rnorm(n)
  Z2 = rnorm(n)
  Z3 = rnorm(n)
  
  #correlated features
  for(i in 1:2)
  {
    x.train[,i] = Z1 + rnorm(n,0,q)
    x.valid[,i] = Z1 + rnorm(n,0,q)
    
  }
  
  for(i in 2:2)
  {
    x.train[,i] = Z2 + rnorm(n,0,q)
    x.valid[,i] = Z2 + rnorm(n,0,q)
    
  }
  
  #independent features
  for(i in 3:3)
  {
    x.train[,i] = Z1+Z2 + rnorm(n,0,q)
    x.valid[,i] = Z1+Z2 + rnorm(n,0,q)
  }
  
  y.train = drop (x.train %*% true.beta) + err.sigma*rnorm (n)
  y.valid = drop (x.valid %*% true.beta) + err.sigma*rnorm (n)
  
  result = list(x.train = x.train,  y.train = y.train, x.valid = x.valid, y.valid=y.valid)
  return(result)
}

d = generate.data(p, n, true.beta, err.sigma)

x = d$x.train
y = d$y.train

res.enet = glmnet(x, y, alpha = 1, intercept = FALSE)
enet.beta = coef(res.enet, .1)

gen.data.mat = function(p, n, true.beta, err.sigma)
{
  x.train = mvrnorm(n,rep(0,p),blockMat)
  y.train = drop (x.train %*% true.beta) + err.sigma*rnorm (n)
  
  x.valid = mvrnorm(n,rep(0,p), blockMat)
  y.valid = drop (x.valid %*% true.beta) + err.sigma*rnorm (n)
  
  result = list(x.train = x.train,  y.train = y.train, x.valid = x.valid, y.valid=y.valid)
  return(result)
}

rho = -.4
Mtr = matrix(c(1,rho,rho,1),2,2)
Mcol = c(.5,.5)
blockMat = cbind(Mtr,Mcol)
Mrow = c(.5,.5,1)
blockMat = rbind(blockMat,Mrow)

d = gen.data.mat(p, n, true.beta, err.sigma)

x = d$x.train
y = d$y.train

res.enet = glmnet(x, y, alpha = 1, intercept = FALSE)
enet.beta = coef(res.enet, .1)
enet.beta

Mcol%*%solve(Mtr) #to check if IC holds are not.
Mtr

