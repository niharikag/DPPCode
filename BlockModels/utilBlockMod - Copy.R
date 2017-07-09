#Utility file

#generate heat map
gen.heatmap = function(d)
{
  #d = gen.data.noiseOverlap(p, n, true.beta, err.sigma)
  x = d$x.train
  sig = t(x)%*%x/n
  my_palette <- colorRampPalette(c("green", "grey", "red"))(n = 1000)
  heatmap.2(sig, Rowv=NA, Colv=NA, dendrogram="none", col=my_palette,
            trace="none", labRow=NA, labCol=NA, margin=c(6, 6) )
  
}


fit.Lasso = function(d)
{
    x = d$x.train
    y = d$y.train
    x.valid = d$x.valid
    y.valid = d$y.valid
    n = dim(x)[1]
    p = dim(x)[2]
    #compute Lambda range 
    lambda.min = 1.5*sqrt(log(p)/n)
    temp = rep(0,p)
    for(i in 1:p){
      temp[i] = abs(t(x[,i])%*%y)
    }
    lambda_max  = max(temp)/n
    #print(lambda.min)
    #print(lambda_max)
    lambda.seq = seq(lambda.min, lambda_max-1, length.out= 100)
    #lambda.seq = seq(0.01, lambda.min, length.out= 30)
    #res.lasso = glmnet(x, y, lambda = lambda.seq, alpha = 1, intercept = FALSE)
    res.lasso = glmnet(x, y, intercept = FALSE)
    y.estimate = predict.glmnet(res.lasso, x.valid)
    gridLen = dim(y.estimate)[2]
    
    #compute residuals
    y.residuals = y.valid - y.estimate
    y.residuals = (y.residuals)^2
    y.residuals = colMeans(y.residuals)
    
    #minimum Prediction Error
    indx = which.min(y.residuals)
    predErr =min(y.residuals)
    beta =  res.lasso$beta[,indx]
    #which(beta!=0)
    #print(res.lasso$lambda[indx])
    result = list(predErr = predErr,  beta = beta)
    
    return(result)
}

fit.LDSRE = function(d, n, lambda_max)
{
  x = d$x.train
  y = d$y.train
  x.valid = d$x.valid
  y.valid = d$y.valid
  n = dim(x)[1]
  p = dim(x)[2]
  #compute Lambda range 
  lambda.min = 1.5*sqrt(log(p)/n)
  temp = rep(0,p)
  for(i in 1:p){
    temp[i] = abs(t(x[,i])%*%y)
  }
  lambda_max  = max(temp)/n
  #print(lambda.min)
  #print(lambda_max)
  lambda.seq = seq(lambda.min, lambda_max-1, length.out= 100)
  #lambda.seq = seq(lambda.min/20, lambda.min, length.out= 50)
  
  res.lasso = glmnet(x, y, lambda = lambda.seq, alpha = 1, intercept = FALSE)
  #res.lasso = glmnet(x, y, intercept = FALSE)
  y.estimate = predict.glmnet(res.lasso, x)
  gridLen = dim(y.estimate)[2]
  
  #compute residuals
  residuals = y - y.estimate
  y.residuals = (residuals)^2
  y.residuals = colMeans(y.residuals)
  
  #minimum Prediction Error
  indx = which.min(y.residuals)
  lambda.opt = res.lasso$lambda[indx]
  s.lasso = which(res.lasso$beta[,indx]!=0)
  theta = residuals[,indx]
  
    
  dual.feasible = rep(0, p)
  
  for(j in 1:p){
    dual.feasible[j] = (t(x[,j])%*%theta)
  }
  #take the absolute of dual feasible
  dual.feasible = abs(dual.feasible)/n
  
  #truncate to second decimal.
  dual.feasible = trunc(dual.feasible*1000,3)/1000
  
  #find dual active variables
  feas.min = min(dual.feasible[s.lasso])
  s.dual = which(dual.feasible >= feas.min)
  
  x.valid = d$x.valid
  y.valid = d$y.valid
  
  #res.ridge = glmnet(x[,s.dual], y, lambda = seq(0, lambda_max, .1), alpha = 0, intercept = FALSE)
  res.ridge = glmnet(x[,s.dual], y, alpha = 0, intercept = FALSE)
  y.estimate = predict.glmnet(res.ridge, x.valid[,s.dual])
  #compute residuals
  y.residuals = y.valid - y.estimate
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


fit.enet = function(d, n, lambda_max)
{
  #d = generate.data(p, n, true.beta, err.sigma)
  x = d$x.train
  y = d$y.train
  x.valid = d$x.valid
  y.valid = d$y.valid
  
  #res.enet = glmnet(x, y, lambda = seq(0, lambda_max, .1), alpha = 0.5, intercept = FALSE)
  res.enet = glmnet(x, y,  alpha = 0.6, intercept = FALSE)
  y.estimate = predict.glmnet(res.enet, x.valid)
  
  y.residuals = y.valid - y.estimate
  y.residuals = (y.residuals)^2
  y.residuals = colMeans(y.residuals)
  
  #minimum Prediction Error
  indx = which.min(y.residuals)
  predErr = min(y.residuals)
  s.enet = which(res.enet$beta[,indx]!=0)
  beta =  res.enet$beta[,indx]
  
  result = list(predErr = predErr,  beta = beta)
  return(result)
}

someExperiment= function()
{
  t.beta = sum(abs())
  t.beta = sum(abs(res.lasso$beta[,indx]))
  
  beta.ridge = coef(res.ridge,s=(t.beta)^2/length(s.dual), exact = TRUE)
  beta.ridge = coef(res.ridge,s=.01, exact = TRUE)
  which(beta.ridge!=0)
  s.lasso
  sum(abs(beta.ridge))
  
  length(beta.ridge)
  y.hat = beta.ridge[1]+ x[,s.dual] %*%beta.ridge[2:31]
  mean((y - y.hat)^2)
  
  
  s1 = colSums(abs(res.ridge$beta))
  which(s1 >= t.beta)
  dim(res.ridge$beta)
  