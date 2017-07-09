#####
##Backup of the compare.R file
############
library(glmnet)

#some initialization
set.seed(07092016)
q = .01
n=50
p=100
err.sigma = 1
s = 30 #number of active features
true.beta = c(rep(3,s),rep(0,p-s))
true.beta[1:s] = sample(c(1,2,3), s, replace = TRUE)

########generate data
generate.data = function(p, n, true.beta, err.sigma)
{
  x.train=matrix(0,n,p)
  x.valid =matrix(0,n,p)
  
  Z1 = rnorm(n)
  Z2 = rnorm(n)
  Z3 = rnorm(n)
  Z4 = rnorm(n)
  Z5 = rnorm(n)
  
  #first 10 features are correlated
  for(i in 1:10)
  {
    x.train[,i] = Z1 + rnorm(n,0,q)
    x.valid[,i]= Z1 + rnorm(n,0,q)
  }
  
  #make the last two  negatively correlated
  x.train[,5] = -x.train[,5]
  x.valid[,5] = -x.valid[,5]
  x.train[,9] = -x.train[,9]
  x.valid[,9] = -x.valid[,9]
  x.train[,10] = -x.train[,10]
  x.valid[,10] = -x.valid[,10]
  
  #make the next five correlated within group and between next group as well
  for(i in 11:15)
  {
    x.train[,i] = Z2 + rnorm(n,0,q)
    x.valid[,i]= Z2 + rnorm(n,0,q)
  }
  
  #make the next five correlated
  for(i in 16:20)
  {
    x.train[,i] = Z3 + Z4+rnorm(n,0,q)
    x.valid[,i]= Z3 + Z4+rnorm(n,0,q)
  }
  
  #make the next six un-correlated
  for(i in 21:26)
  {
    x.train[,i] = rnorm(n)
    x.valid[,i]= rnorm(n)
  }
  
  #make the next four linear dependent
  for(i in 27:30)
  {
    x.train[,i] = rnorm(n)
    x.valid[,i]= rnorm(n)
  }
  
  x.train[,27] = x.train[,28]+x.train[,29]+x.train[,30]+rnorm(n,0,q)
  x.valid[,27] = x.valid[,28]+x.valid[,29]+x.valid[,30]+rnorm(n,0,q)
  #
  for(i in 31:100)
  {
    x.train[,i] = rnorm(n)
    x.valid[,i]=rnorm(n)
  }
  
  
  y.train = drop (x.train %*% true.beta) + err.sigma*rnorm (n)
  y.valid = drop (x.valid %*% true.beta) + err.sigma*rnorm (n)
  
  result = list(x.train = x.train,  y.train = y.train, x.valid = x.valid, y.valid=y.valid)
  return(result)
}

compareMethods = function()
{
  d = generate.data(p, n, true.beta, err.sigma)
  x = d$x.train
  y = d$y.train
  
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
  from = .1
  to = lambda_max/10
  
  lam_range = seq(from, to, by = ((to - from)/20))
  
  for(i in 1:length(lam_range))
  {
    #Step1: fit lasso for the range of lamda
    res.lasso = glmnet(x,y,lambda = lam_range[i], alpha = 1, intercept = FALSE)
    
    s.lasso = (which(res.lasso$beta!=0))
    
    #compute the lasso prediction error
    y.hat = predict.glmnet(res.lasso, x, lambda = lam_range[i], intercept = FALSE)
    pred.err.lasso = sum((y-y.hat)^2)/n
    
    print("Lasso Error")
    print(pred.err.lasso)
    print("Lasso selection")
    print(s.lasso)
    
    #plot lasso
    
    
    #Step2(a): compute theta, for each lamda
    #theta = (y - (x%*%res.lasso$beta))
    theta = (y - y.hat)
    
    dual.feasible = rep(0, p)
    
    for(j in 1:p){
      #dual.feasible[j] = ((abs(t(x[,j])%*%theta)/n)[1][1])/lam_range[i]
      dual.feasible[j] = (t(x[,j])%*%theta)/lam_range[i]
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
      #Step 3: fit ridge regression for the set of features selected by the Lasso dual
      #Step1: fit lasso for the range of lamda
      res.ridge = glmnet(x[,s.dual],y, alpha = 0, lambda = lam_range[i])
      res.ridge$beta
      
      #compute the prediction error
      y.hat = predict.glmnet(res.ridge, x[,s.dual])
      pred.err.ridge = sum((y-y.hat)^2)/n
      print("Ridge Error")
      print(pred.err.ridge)
      print("Ridge selection")
      print(s.dual)
    }
    
    
    #fit Elastic net for comparison
    res.enet = glmnet(x, y, lambda = lam_range[i], alpha = 0.5, intercept = FALSE)
    s.enet = (which(res.enet$beta!=0))
    #compute the enet prediction error
    y.hat = predict.glmnet(res.enet, x, lambda = lam_range[i],  alpha = 0.5, intercept = FALSE)
    pred.err.enet = sum((y-y.hat)^2)/n
    print("Enet Error")
    print(pred.err.enet)
    print("Enet selection")
    print(s.enet)
    
    # plot(res.lasso, "lambda", label=FALSE)
    #  text (6, 0.4, "A", cex=1.8, font=1)
    
  }
}

run.lars = function ()
{
  #LARS fit, difficult to understand
  res.lars <- lars(x,y,type="lar")
  plot(res.lars)
}

###bakcup
fit.enet.backup = function(d, lam_range)
{
  d = generate.data(p, n, true.beta, err.sigma)
  x = d$x.train
  y = d$y.train
  x.valid = d$x.valid
  y.valid = d$y.valid
  #fit Elastic net for comparison
  res.enet = glmnet(x, y, lambda = lam_range, alpha = 0.5, intercept = FALSE)
  
  y.estimate = predict.glmnet(res.enet, x.valid, lambda = lam_range, alpha = 0.5, intercept = FALSE)
  
  gridLen = length(lam_range)
  
  #compute residuals
  y.residuals = y.valid - y.estimate
  MSE = rep(0,gridLen)
  for(i in 1:gridLen){
    MSE[i] = mean(y.residuals[,i]^2)
  }
  
  #mean Prediction Error
  indx = which.min(MSE)
  predErr = MSE[indx]
  fit = coef(res.enet , res.enet$lambda[indx]) 
  fit = fit[-1]
  # result = list(predErr=predErr, fit=fit, lamda_fit = fit1$lambda[indx])
  #return(result)
  
  s.enet = (which(fit!=0))
  #compute the enet prediction error
  
  print("Enet Error")
  print(predErr)
  print("Enet selection")
  print(s.enet)
  
  # plot(res.lasso, "lambda", label=FALSE)
  #  text (6, 0.4, "A", cex=1.8, font=1)
  
}

heatMapRelated = function()
{
  if (!require("gplots")) {
    install.packages("gplots", dependencies = TRUE)
    library(gplots)
  }
  if (!require("RColorBrewer")) {
    install.packages("RColorBrewer", dependencies = TRUE)
    library(RColorBrewer)
  }
  # creates a own color palette from red to green
  my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)
  
  heatmap(sig,   keep.dendro = FALSE,Colv = "Rowv")   
  
  heatmap.2(sig,
            #cellnote = mat_data,  # same data set for cell labels
            main = "Correlation", # heat map title
            notecol="black",      # change font color of cell labels to black
            density.info="none",  # turns off density plot inside color legend
            trace="none",         # turns off trace lines inside the heat map
            margins =c(9,9),     # widens margins around plot
            #col=my_palette,       # use on color palette defined earlier
            #breaks=col_breaks,    # enable color transition at specified limits
            dendrogram="row",     # only draw a row dendrogram
            Colv="NA")            # turn off column clustering
}

Backup.ridge = function()
{
  #Step 3: fit ridge regression for the set of features selected by the Lasso dual
  #Step1: fit lasso for the range of lamda
  #res.ridge = glmnet(x[,s.dual],y, alpha = 0, lambda = lambda.opt)
  res.ridge = glmnet(x[,s.dual],y, alpha = 0)
  
  #compute the prediction error
  y.estimate = predict.glmnet(res.ridge, x.valid[,s.dual], alpha = 0)
  
  gridLen = dim(y.estimate)[2] #length(lam_range)
  
  #compute residuals
  y.residuals = y.valid - y.estimate
  MSE = rep(0,gridLen)
  for(i in 1:gridLen){
    MSE[i] = mean(y.residuals[,i]^2)
  }
  
  #minimum Prediction Error
  indx = which.min(MSE)
  predErr = MSE[indx]
  
  
  x.red = x[,s.dual]
  res.ridge = lm.ridge(y~x.red,lambda = seq(0,0.1,0.001))  
  select(res.ridge)
  x.valid.red = x.valid[,s.dual]
  pred.err.ridge = (y.valid - mean(y.valid) - x.valid.red%*%res.ridge$coef)^2
  pred.err.ridge = mean(pred.err.ridge)
  
  x.red = x[,s.dual]
  res.ridge = lm.ridge(y~x.red, lambda = .1)#lambda = lambda.opt)  
  x.valid.red = x.valid[,s.dual]
  pred.err.ridge = (y.valid - mean(y) - x.valid.red%*%res.ridge$coef)^2
  pred.err.ridge = mean(pred.err.ridge)
  
  
  
  x.red = x[,s.dual]
  res.ridge = lm.ridge(y~x.red,lambda = seq(0,2,0.1))  
  x.valid.red = x.valid[,s.dual]
  y.valid = y.valid - mean(y.valid)
  pred.err.ridge = ( y.valid - x.valid.red%*%res.ridge$coef)
  pred.err.ridge = (pred.err.ridge)^2
  pred.err.ridge = colMeans(pred.err.ridge)
  pred.err.ridge = min(pred.err.ridge)
  
}

for(i in 29:30)
{
  x.train[,i] = -x.train[,i]
  x.valid[,i]= -x.valid[,i]
}

for(i in 19:20)
{
  x.train[,i] = -x.train[,i]
  x.valid[,i]= -x.valid[,i]
}