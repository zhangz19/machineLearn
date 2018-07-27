
rm(list=ls())
options(scipen=500)

# library(lme4)
library(R.matlab)
library(gplots) 
require(RColorBrewer)
require(randomForest)
require(xgboost)
require(rms)
# require(asbio)   #for partial.R2 function
# source('util.R')
evalPred <- function(y, yhat){
  # rmse <- sqrt(mean((y-yhat)^2)) # rmad=mean(abs(yhat/y-1)), 
  r2 <- summary(lm(y~yhat))$r.squared
  out <- r2  #c(rmse, r2)
  return(round(out,3))
}

dat <- read.csv('VS.csv')    #59 obs. of  20 variables:

table(dat$Year, dat$VS)

# exclude year 1999
inds <- which(dat$Year==1999)
if(length(inds)) dat <- droplevels(dat[-inds,])

dat$Year <- factor(dat$Year, levels=sort(unique(dat$Year)))
dat$VS <- factor(dat$VS, levels=sort(unique(dat$VS)))
dat$LakeInVS <- factor(dat$LakeInVS, levels=sort(unique(dat$LakeInVS)))
str(dat)
table(dat$Year, dat$VS)


# model as of 2017-8-31: drop D_DSDam_m + D_USDam_m +
X0 <- model.matrix( ~ Year + MRW_m + Depth_avg + Ln_VS_area + N_DS_barr +  D_LkOnt_m + Ag_pct + For_pct + Urb_pct + Wetl_pct + Depth_sd + Dep_div + Vel_div + Sub_div + WatTemp + Hab_div, data=dat)[,-1] #exclude the intercept

# standardize it
X0 <- scale(X0)

y0 <- dat$SR_jack1  #hist(y) shows the normality is good
mat <- as.data.frame(X0); mat$y <- y0

# check collinearity
tmp <- mat[,ncol(mat):1]; names(tmp)[1] <- 'SR_jack1'
corMat <- cor(tmp)  # use the correlation matrix
# png(paste0('./heat.png'), units='in', height=5, width=6, pointsize=11, res=600)
par(mfrow=c(1,1), mar=c(0,0,0,0)+.4, cex.axis=.9)
wrapper <- function(x) rev(brewer.pal(x,'RdBu'))
heatmap.2(corMat, cexRow=.8, col=wrapper, breaks=10, trace='none', key.title=NA, revC=TRUE, key.xlab='Correlation')

n <- length(y0)
nmod <- length(mods <- c('Multiple regression', 'Random forests', 
                         'Extreme gradient boosting'))
yhat <- matrix(NA, n, nmod)
nfold <- 10
binSize <- ceiling(n/nfold)

set.seed(8)
ind0 <- sample(n, n) #random permutation
for(k in 1:nfold){
  indv <- ind0[(k-1)*binSize + c(1:binSize)]
  vd <- mat[indv,] #validate
  td <- mat[-indv,] #train
  
  m0 <- lm(y~., data=td) #linear model
  yhat[indv,1] <- predict(m0, newdata=vd)
  
  set.seed(15)
  m1 <- randomForest(y~., data=td, ntree=500, importance=TRUE)
  yhat[indv,2] <- predict(m1, newdata=vd)
  
  myparam <- list(objective="reg:linear", metrics="rmse", max_depth=6, 
                  eta=.004, gamma=0.144, subsample=0.78, min_child_weight=0.18, 
                  colsample_bytree=0.54, nthread=1)
  X <- data.matrix(td[,-which(names(td)=='y')])
  train <- xgb.DMatrix(data=X, label=td$y)
  # ptm <- proc.time()[3]
  set.seed(15)
  m2 <- xgb.train(data=train, params=myparam, nrounds=500)
  # cputime <- as.numeric(proc.time()[3]-ptm)
  # cat('XGBoost took', round(cputime/60,3), 'minutes to complete.','\n')
  yhat[indv,3] <- predict(m2, data.matrix(vd[, -which(names(td)=='y')]))
}

# cbind(y0, yhat)
for(i in 1:nmod){
  if(i==1) cat('R-squared based on the',nfold,'fold cross validation:\n')
  cat(mods[i],': ', evalPred(y0, yhat[,i]),'\n', sep='')
}; cat('\n')

# ranking by refitting the models to the entire data
m0 <- ols(y~., data=mat)
cat(mods[1], 'R-squared overall fit: ', evalPred(mat$y, predict(m0)),'\n')
v0 <- rev(plot(anova(m0), what='partial R2', pl=FALSE))

set.seed(15)
m1 <- randomForest(y~., data=mat, ntree=500, importance=TRUE, keep.forest=TRUE)
cat(mods[2],': ', evalPred(mat$y, predict(m1, newdata=mat)),'\n', sep='') 
v1 <- importance(m1)[,'%IncMSE'] #for regression tree
v1 <- sort(v1)

myparam <- list(objective="reg:linear", metrics="rmse", max_depth=6, eta=.004, gamma=0.144, subsample=0.78, min_child_weight=0.18, colsample_bytree=0.54, nthread=1)
X <- data.matrix(mat[,-which(names(mat)=='y')])
train <- xgb.DMatrix(data=X, label=mat$y)
# ptm <- proc.time()[3]
set.seed(15)
m2 <- xgb.train(data=train, params=myparam, nrounds=500)
cat(mods[3],': ', evalPred(mat$y, predict(m2, data.matrix(mat[, -which(names(td)=='y')]))), sep='')
tmp <- xgb.importance(model=m2)[,c('Feature','Gain')]
v2 <- tmp$Gain; names(v2) <- tmp$Feature
v2 <- sort(v2)

par(mfrow=c(1,3), mar=c(2,.5,1.5,1.2), mgp=c(.9,.1,0), cex.axis=.8, tck=-0.01, cex.main=.86, cex=.9)
dotchart(v0, pch=16, ylab='', xlab="Partial R-squared", main='Multiple regression')
dotchart(v1, pch=16, ylab='', xlab="Prediction accuracy", main='Random Forests')
dotchart(v2, pch=16, ylab='', xlab="Information gain", main='eXtreme Gradient Boosting')


# partial influence for RF
# par(mfrow=c(4,4), mar=c(1,2,1.2,.5), mgp=c(1,.2,0), cex.axis=.8, tck=-0.02, cex.main=.86, cex=.9)
lmat <- matrix(1:16, nrow=4, ncol=4, byrow=TRUE)
layout(mat=lmat, heights=rep(.25,4))
par(oma=c(0,2,2,0), font.lab=1, mar=c(2,.1,.1,.7)+.4,mgp=c(1.2,.2,0), tck=-0.02, cex.axis=1, cex.lab=1.1, cex.main=1.2 )
xnam <- rev(names(v1))
for(i in 1:length(xnam)){
  partialPlot(x=m1, pred.data=mat, x.var=xnam[i], main='', ylab='', xlab=xnam[i], col='blue', lwd=2) 
  grid(col='gray70')
  partialPlot(x=m1, pred.data=mat, x.var=xnam[i], col='blue', lwd=2, add=TRUE)
}
title('Random forests', outer=TRUE, font.main=1, cex.main=1.5)
mtext("Partial dependence", side=2, line=.6, outer=TRUE, cex=1, font=1)




