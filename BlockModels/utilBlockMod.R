#Utility file
library(gplots)
#generate heat map
gen.heatmap = function(d)
{
  #d = gen.data.noiseOverlap(p, n, true.beta, err.sigma)
  x = scale(d$x.train)
  sig = t(x)%*%x/n
  my_palette <- colorRampPalette(c("green", "grey", "red"))(n = 1000)
  heatmap.2(sig, Rowv=NA, Colv=NA, dendrogram="none", col=my_palette,
            trace="none", labRow=NA, labCol=NA, margin=c(6, 6) )
  
}

#####fit Lasso using LARS algorithms
fit.Lasso = function(d)
{
    x = scale(d$x.train)
    y = scale(d$y.train)
    x.valid = d$x.valid
    y.valid = d$y.valid
    
    #compute Lambda range 
    n = dim(x)[1]
    p = dim(x)[2]
    lambda.min = sqrt(log(p)/n)
    temp = rep(0,p)
    for(i in 1:p){
      temp[i] = abs(t(x[,i])%*%y)
    }
    lambda_max  = max(temp)/n
    lambda.seq = seq(lambda_max, lambda.min, length.out= 100)
    
    ##Using LARS algorithm
    res.lasso =  lars(x, y, type = "lasso", intercept = FALSE, normalize = FALSE)
    fits.lasso = predict.lars(res.lasso, x.valid, s = lambda.seq, mode="lambda", type="fit")
    y.estimate = fits.lasso$fit
    
    #compute residuals
    y.residuals = y.valid - y.estimate
    y.residuals = (y.residuals)^2
    y.residuals = colMeans(y.residuals)
    
    #minimum Prediction Error
    indx = which.min(y.residuals)
    predErr =min(y.residuals)
    beta =  coef(res.lasso, s=lambda.seq[indx], mode="lambda")
    #which(beta!=0)
    #print(res.lasso$lambda[indx])
    result = list(predErr = predErr,  beta = beta)
    return(result)
    
    
    ###to verfiy : glmnet is also giving almost the SAME values
    res.lasso = glmnet(x, y, alpha = 1, intercept = FALSE, standardize = FALSE)
    y.estimate = predict(res.lasso, x.valid, s = lambda.seq/n, mode="lambda", type="response")
    
    #compute residuals
    y.residuals = y.valid - y.estimate
    y.residuals = (y.residuals)^2
    y.residuals = colMeans(y.residuals)
    
    #minimum Prediction Error
    indx = which.min(y.residuals)
    predErr =min(y.residuals)
    beta =  coef(res.lasso, s=lambda.seq[indx]/n, mode="lambda")
}

fit.LDSRE = function(d, n, lambda_max)
{
  x = scale(d$x.train)
  y = d$y.train
  x.valid = d$x.valid
  y.valid = d$y.valid
  
  #compute Lambda range 
  n = dim(x)[1]
  p = dim(x)[2]
  lambda.min = sqrt(log(p)/n)
  temp = rep(0,p)
  for(i in 1:p){
    temp[i] = abs(t(x[,i])%*%y)
  }
  lambda_max  = max(temp)/n
  lambda.seq = seq(lambda_max, lambda.min, length.out= 100)
  
  #using LARS
  res.lasso =  lars(x, y, type = "lasso", intercept = FALSE, normalize = FALSE)
  fits.lasso = predict.lars(res.lasso, x.valid, s = lambda.seq, mode="lambda", type="fit")
  y.estimate = fits.lasso$fit
  
  #compute residuals
  y.residuals = y.valid - y.estimate
  y.residuals = (y.residuals)^2
  y.residuals = colMeans(y.residuals)
  
  #minimum Prediction Error
  indx = which.min(y.residuals)
  predErr =min(y.residuals)
  lasso.beta =  coef(res.lasso, s=lambda.seq[indx], mode="lambda")
  
  s.lasso = which(lasso.beta!=0)
  lasso.pred = predict.lars(res.lasso, x, s= lambda.seq[indx], type="fit", mode="lambda")
  theta = y - lasso.pred$fit
    
  dual.feasible = rep(0, p)
  for(j in 1:p){
    dual.feasible[j] = (t(x[,j])%*%theta)
  }
  #take the absolute of dual feasible
  dual.feasible = abs(dual.feasible)
  dual.feasible = dual.feasible/lambda.seq[indx]
  dual.feasible = round(dual.feasible,1)
  s.dual = which(dual.feasible == 1)
  
  
  lambda2.seq =  seq(0,100,0.01)
  ridgec <- lm.ridge (y ~ x[,s.dual], lambda = 1/n)
  #ridge.coef = coef(ridgec, lambda = 1/n)
  #ridge.coef = coef(ridgec, lambda = lambda2.seq)
  ridge.est = x.valid[,s.dual]%*% ridgec$coef
  
  #compute residuals
  y.residuals = y.valid - ridge.est
  y.residuals = (y.residuals)^2
  #y.residuals = colMeans(y.residuals)
  predErr = mean(y.residuals)
  
  #minimum Prediction Error
  beta = rep(0,p)
  beta[s.dual] =  ridgec$coef
  
  result = list(predErr = predErr,  beta = beta)
  return(result)
  
  
  ####################ignore for the time being
  res.ridge = glmnet(x[,s.dual], y, alpha = 0, intercept = FALSE, standardize = FALSE)
  y.estimate = predict.glmnet(res.ridge, x.valid[,s.dual], s = lambda.seq/n, mode="lambda", type="response")
  #compute residuals
  y.residuals = y.valid - y.estimate
  y.residuals = (y.residuals)^2
  y.residuals = colMeans(y.residuals)
  
  #minimum Prediction Error
  indx = which.min(y.residuals)
  predErr =min(y.residuals)
  lmd = res.ridge$lambda[indx]
  beta = coef(res.ridge, s=lmd, mode="lambda")
  
  result = list(predErr = predErr,  beta = beta)
  return(result)
  
  #using glmnet
  res.lasso = glmnet(x, y, intercept = FALSE, standardize = FALSE)
  y.estimate = predict(res.lasso, x.valid, s = lambda.seq/n, mode="lambda", type="response")
  
}

fit.enet = function(d, n, lambda_max)
{
  x = scale(d$x.train)
  y = d$y.train
  x.valid = d$x.valid
  y.valid = d$y.valid
  
  #compute Lambda range 
  n = dim(x)[1]
  p = dim(x)[2]
  lambda.min = sqrt(log(p)/n)
  temp = rep(0,p)
  for(i in 1:p){
    temp[i] = abs(t(x[,i])%*%y)
  }
  lambda_max  = max(temp)/n
  #print(lambda.min)
  #print(lambda_max)
  lambda.seq = seq(lambda_max-1, lambda.min/2, length.out= 100)
  
  ##using LARS algorithm
  #lambda_2 = lambda.min
  lambda_2 = 1
  #compute enet applying Lasso on the modified design matrix
  dg = diag(sqrt(lambda_2),p) 
  x1 = rbind(x, dg) #Augment the design matrix with d
  x1 = x1/sqrt(1+lambda_2)
  y1 =  c(y, rep(0,p)) #Augment the response vector with 0
  lmd = lambda.seq/sqrt(1+lambda_2)
  res.lasso =  lars(x1, y1, type = "lasso", intercept = FALSE, normalize = FALSE)
  fits.lasso = predict.lars(res.lasso, x.valid, s= lmd, type="fit", mode="lambda")
  y.estimate = fits.lasso$fit
  
  #compute residuals
  y.residuals = y.valid - y.estimate
  y.residuals = (y.residuals)^2
  y.residuals = colMeans(y.residuals)
  
  #minimum Prediction Error
  indx = which.min(y.residuals)
  predErr =min(y.residuals)
  beta =  coef(res.lasso, s=lmd[indx], mode="lambda")
  result = list(predErr = predErr,  beta = beta)
  
  return(result)
  
  
  #############################
  #d = generate.data(p, n, true.beta, err.sigma)
  
  
  
  alpha = 1/n
  
  ####using glmnet, to check if we are getting the same value
  res.enet = glmnet(x, y,  alpha = alpha, intercept = FALSE, standardize = FALSE)
  y.estimate = predict(res.enet, x.valid, s = lambda.seq/n, mode="lambda", type="response")
  
  
  y.residuals = y.valid - y.estimate
  y.residuals = (y.residuals)^2
  y.residuals = colMeans(y.residuals)
  
  #minimum Prediction Error
  indx = which.min(y.residuals)
  predErr = min(y.residuals)
  s.enet = which(res.enet$beta[,indx]!=0)
  beta =  coef(res.enet, s=lambda.seq[indx]/n, mode="lambda")
  
  result = list(predErr = predErr,  beta = beta)
}


TBD = function()
{
  ###########to be deleted
  
  res.ridge = glmnet(x[,s.dual], y, alpha = 0, intercept = FALSE, standardize = FALSE)
  y.estimate = predict.glmnet(res.ridge, x.valid[,s.dual], 
                              s=lambda.seq, mode="lambda", type="response")
  #compute residuals
  y.residuals = y.valid - y.estimate
  y.residuals = (y.residuals)^2
  y.residuals = colMeans(y.residuals)
  
  #minimum Prediction Error
  indx = which.min(y.residuals)
  predErr =min(y.residuals)
  beta = coef(res.ridge, s=lambda.seq[indx]/n, mode="lambda")
}