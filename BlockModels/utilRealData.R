#Real Data Riboflavin
library(hdi)
library(glmnet)

rib = data(riboflavin)
dim(riboflavin)
x.ribo = riboflavin$x #4088 gene predictor
y.ribo = riboflavin$y #response variable (riboflavin production)
dim(x.ribo)
length(y.ribo)
#select first 1000 variables having highest variances

p= dim(x.ribo)[2]
n= dim(x.ribo)[1]

x = x.ribo
y = y.ribo

performLassoRibo = function(x, y)
{
  
  res = cv.glmnet(x, y)
  y.estimate = predict.cv.glmnet(res, x)
  
  y.residuals = y - y.estimate
  y.residuals = y.residuals^2
  
  predErr =  colMeans(y.residuals)
  
  lassoResult = coef(res,s="lambda.min")
  cfs = lassoResult[-1]
  lasso.index = which(cfs!=0)
}

performEnetRibo = function(x, y)
{
  res.enet = cv.glmnet(x, y, alpha = 0.5)
  y.estimate = predict.cv.glmnet(res.enet, x)
  
  y.residuals = y - y.estimate
  y.residuals = y.residuals^2
  
  predErr =  colMeans(y.residuals)
  indx = which.min(predErr)
  enetResult = coef(res.enet,s="lambda.min")
  cfs = enetResult[-1]
  enet.index = which(cfs!=0)
}


performLDSRERibo = function(x, y)
{
  res.lasso = cv.glmnet(x, y, standardize = FALSE)
  lambda.opt = res.lasso$lambda.min
  
  res.lasso = glmnet(x, y, lambda = lambda.opt, standardize = FALSE)
  y.estimate = predict.glmnet(res.lasso, x)
  
  #compute residuals
  residuals = y - y.estimate
  y.residuals = (residuals)^2
  predErr = mean(y.residuals)
  
  #minimum Prediction Error
  s.lasso = which(res.lasso$beta != 0)
  
  y.hat = res.lasso$a0+x[,s.lasso]%*% res.lasso$beta[s.lasso,]
  theta = y - y.hat
  
  dual.feasible = rep(0, p)
  
  for(j in 1:p){
    dual.feasible[j] = (t(x[,j])%*%theta)
  }
  #take the absolute of dual feasible
  dual.feasible = abs(dual.feasible)/n
  
  #truncate to second decimal.
  dual.feasible = trunc(dual.feasible*100,2)/100
  new.data = dual.feasible/(lambda.opt)
  feas.min = min(dual.feasible[s.lasso])
  s.dual = which(dual.feasible >= feas.min)
  
  
  #res.ridge = glmnet(x[,s.dual], y, lambda = seq(0, lambda_max, .1), alpha = 0, intercept = FALSE)
  res.ridge = glmnet(x[,s.dual], y, alpha = 0)
  y.estimate = predict.glmnet(res.ridge, x[,s.dual])
  #compute residuals
  y.residuals = y - y.estimate
  y.residuals = (y.residuals)^2
  y.residuals = colMeans(y.residuals)
  
  #minimum Prediction Error
  indx = which.min(y.residuals)
  predErr =min(y.residuals)
  beta =  res.ridge$beta[,indx]
  
  #print(sum(res.ridge$beta[,indx]))
  
  result = list(predErr = predErr,  beta = beta)
  return(result)
}


#generate heat map
gen.heatmap = function(x)
{
  sig = t(x)%*%x/n
  my_palette <- colorRampPalette(c("green", "grey", "red"))(n = 1000)
  heatmap.2(sig, Rowv=NA, Colv=NA, dendrogram="none", col=my_palette,
            trace="none", labRow=NA, labCol=NA, margin=c(6, 6) )
  
}