
rm(list=ls())
dat <- read.csv('HGMI-1.csv', check=F)


saveFig <- FALSE
saveDir <- './'

# exclude missing values
mat <- na.omit(dat[,-1])  #remove watershed
row.names(mat) <- 1:nrow(mat)


# now extract the matrix with categorical variables such as year
myFormula <- as.formula(paste(' ~', paste(names(mat), collapse='+'),'-1' ))
mat <- as.data.frame( model.matrix(myFormula , data=mat))


####### util functions

.panelhist <- function(x, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE, breaks=8)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col="cyan4", ...)
    if(length(unique(x)) == 2){
     text(breaks[c(1,nB-1)]+.05, y[c(1,nB-1)]+.2, c("F","M"))
    }
}


.panelcor <- function(x, y, digits=2, prefix="", cex.cor)
{
    if(length(unique(y)) == 2) return(NULL)
    else{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- (cor(x, y, use="complete.obs", method='spearman'))
    #tmp <- cor.test(x,y, exact=F, method="spearman")
    #r <- tmp$estimate; pvalue <- tmp$p.value
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep=""); tcol <- 'black'
    # if(tmp$p.value < 0.05) tcol <- 'red'
    #else if(r < .15) tcol <- 'blue2'
    text(0.5, 0.5, txt, cex = .5/strwidth(txt), col = tcol)
    }
}

library(car)
.panellm <- function(x,y)
{ 
  #points(x,y,pch=1,col=cols0); cfs <- coef(lm(y~x)); abline(a=cfs[1],b=cfs[2],col='cyan4', lwd=1) 
  # dataEllipse(x,y,levels=c(0.975),add=T, plot.points=FALSE, center.cex=0, col='blue', lwd=1)
  dataEllipse(x,y,levels=c(0.95),add=T, plot.points=FALSE, center.cex=0, col='blue', lwd=1)
  points(x,y,pch=1,col='black'); cfs <- coef(lm(y~x)); abline(a=cfs[1],b=cfs[2],col='red', lwd=1) 
}          #lmcols0


# red = gain, blue = loss


if(saveFig) jpeg(paste0(saveDir,'pairs.jpg'), quality=100, height=1800,width=3000, pointsize=12, res=300) 
par(mar=c(0,0,0,0)+.0,mgp=c(.6,.3,0), tck=-0.04, cex.axis=1.2, cex.lab=.8, cex.main=1)
cols0 <- c('blue','red')[as.integer(dat$X.Change.in.Groundwater.Recharge > 0)+1]
pairs(HGMI  ~ ., data=mat,
  main=paste('Scatterplot matrix',sep=''), 
        diag.panel=.panelhist, upper.panel=.panellm,
        lower.panel=.panelcor, #labels=nams,
        cex.labels=.85, gap=0.0, font.labels=1)
if(saveFig) dev.off()


# scatterplot ends above



# calculate the prediction measure
### calculating measures for prediction
evalPred <- function(y,x){ #y=observed, x=fitted
  ym <- mean(y)
  xm <- mean(x)
  # R2 <- (sum((y-ym)*(x-xm)))^2 / sum((y-ym)^2*(x-xm)^2)
  RMSE <- sqrt(mean((y-x)^2))
  NSE <- 1- sum((y-x)^2)/sum((y-ym)^2)
  PBIAS <- 100 * sum(y-x)/sum(y)
  z <- c(RMSE,NSE,PBIAS); names(z) <- c('RMSE','NSE','PBIAS')
  return(z)
}

# for random forests
require(randomForest)

set.seed(18) # for reproducibility
# m1 <- randomForest(HGMI~., data=mat, ntree=100, mtry=2, importance=TRUE, keep.forest=TRUE) #without Watershed
m1 <- randomForest(HGMI~., data=mat, ntree=1000,importance=TRUE, keep.forest=TRUE)

print(m1)
# plot(m1)

# variable importance
if(saveFig) jpeg(paste0(saveDir,'RelativeInfluence.jpg'), quality=100, height=1800,width=3000, pointsize=14, res=300) 
varImpPlot(m1, sort = T, main="Variable Importance", pch=16)
dev.off()

imp <- importance(m1)
impvar <- rownames(imp)[order(imp[, 1], decreasing=TRUE)]
impvarLab <- paste0(impvar, ' (',sort(round(imp[, 1],2), decreasing=TRUE), ')')

# if(saveFig) jpeg(paste0(saveDir,'PartialDependence.jpg'), quality=100, height=1800,width=3000, pointsize=14, res=300)
# par(mfrow=c(2, 4), mar=c(1.5,1.5,1.2,.1)+.4,mgp=c(1.4,.3,0), tck=-0.01, cex.axis=1.1, cex.lab=1.3, cex.main=1.2)
# for (i in seq_along(impvar)){
#   partialPlot(m1, mat, impvar[i], xlab='', col='blue', lwd=2, main=impvar[i], ylim=c(36, 52))
#   grid(col='gray70'); 
#   partialPlot(m1, mat, impvar[i], xlab=impvar[i], col='blue', lwd=2, add=T)
# }
# if(saveFig) dev.off()

if(saveFig) jpeg(paste0(saveDir,'PartialDependence.jpg'), quality=100, height=1800,width=3000, pointsize=14, res=300)

lmat <- matrix(1:8,nrow = 2,ncol = 4,byrow = TRUE)
layout(mat = lmat, heights = rep(.5,2))
par(oma=c(0,2,2,0), font.lab=1, mar=c(2,.1,.1,.7)+.4,mgp=c(1.2,.2,0), tck=-0.02, cex.axis=1, cex.lab=1.1, cex.main=1.2 )
for (i in seq_along(impvar)){
  partialPlot(m1, mat, impvar[i], xlab=impvarLab[i], col='blue', lwd=2, main='', ylab='', ylim=c(36, 52))
  grid(col='gray70');
  partialPlot(m1, mat, impvar[i], xlab=impvarLab[i], col='blue', lwd=2, add=T)
}

par(mar=c(0,0,0,0)); plot(1, type = "n", axes=FALSE, xlab="", ylab="")  #fill the 8th panel since we only have 7 vars
title('Random forests', outer=T, font.main=2, cex.main=1.5)
mtext("Partial dependence", side=2, line=.6, outer=TRUE, cex=1, font=2)

if(saveFig) dev.off()


# get predictive measures
# plot(m1$predicted ~ mat$HGMI, pch=16); abline(0,1)
evalPred(mat$HGMI, m1$predicted)


# for cross validation
require(ipred) #wrapper of random forest for cross validation
# k=5 fold cross validation
set.seed(20) # for reproducibility
errR <- errorest(HGMI~., data=mat, model=randomForest, estimator='cv', est.para=control.errorest(k=5, predictions=TRUE, getmodels=TRUE), ntree=500, mtry=2, importance=T) 

print(errR)

# get predictive measures for 5-fold cross validation
n <- nrow(mat) # sample size
vec <- numeric(); K <- 5
y0 <- errR$predictions
for(k in 1:K){
  inds <- 1:n
  model <- errR$models[[k]]
  inds <- inds[- as.integer(row.names(as.data.frame(model$predicted)))]
  vec <- rbind(vec, evalPred(mat$HGMI[inds], y0[inds]))
}
(eval1CV <- round(colMeans(vec), 5))  # average prediction based on the 5-fold cross validation




#*********************************************** compare with boosted regression tree

# GBT: Boosted regression tree
require(dismo)

set.seed(10)
m2 <- gbm.step(data=mat, gbm.x=2:ncol(mat), gbm.y=1, family='gaussian', n.folds = 10, learning.rate = 0.01, tree.complexity = 1,  bag.fraction = 0.75, n.trees = 50, step.size = 50, max.trees = 10000) 

m2$n.trees

# # setting of Yang et al. 2016
# set.seed(10)
# m2 <- gbm.step(data=mat, gbm.x=2:ncol(mat), gbm.y=1, family='gaussian', n.folds = 10, learning.rate = 0.0025, tree.complexity = 9,  bag.fraction = 0.75, n.trees = 50, step.size = 50, max.trees = 10000) 

if(saveFig) jpeg(paste0(saveDir,'VI_BRT.jpg'), quality=100, height=1800,width=3000, pointsize=14, res=300)
par(mfrow=c(1,1), mar=c(2,12,1.5,0)+.4,mgp=c(1.3,.3,0), tck=-0.01, cex.axis=1, cex.lab=1, cex.main=1, las=2) 
par(xaxt='n')
a <- summary(m2, plotit=T)
par(las=1,xaxt='s'); axis(1,at=seq(0,30,by=5))
dev.off()

row.names(a) <- NULL
print(a)  # this is variable importance

(eval2 <- evalPred(mat$HGMI, fitted(m2)))


# get predictive measures for 5-fold cross validation
n <- nrow(mat) # sample size
vec <- numeric(); K <- 5
set.seed(20)
for(k in 1:K){
  print(k)
  inds <- 1:n
  modelRF <- errR$models[[k]]
  indTrain <- as.integer(row.names(as.data.frame(modelRF$predicted)))
  a <- droplevels(mat[indTrain,])
  
  model <- gbm.step(data=a, gbm.x=2:ncol(mat), gbm.y=1, family='gaussian', n.folds = 10, learning.rate = 0.01, tree.complexity = 1,  bag.fraction = 0.75)
  
  # # setting of Yang et al. 2016
  # model <- gbm.step(data=droplevels(mat[indTrain,]), gbm.x=2:ncol(mat), gbm.y=1, family='gaussian', n.folds = 10, learning.rate = 0.01, tree.complexity = 1,  bag.fraction = 0.75)
  
  if(is.null(model)){
    model <- gbm.step(data=a, gbm.x=2:ncol(mat), gbm.y=1, family='gaussian', n.folds = 10, learning.rate = 0.01, tree.complexity = 1,  bag.fraction = 0.75)
    
    # # setting of Yang et al. 2016
    # model <- gbm.step(data=droplevels(mat[indTrain,]), gbm.x=2:ncol(mat), gbm.y=1, family='gaussian', n.folds = 10, learning.rate = 0.01, tree.complexity = 1,  bag.fraction = 0.75)
  }
  
  indTest <- inds[- indTrain]
  b <- droplevels(mat[indTest,])
  y2 <- predict(model, b, n.trees=model$n.trees)
  vec <- rbind(vec, evalPred(mat$HGMI[indTest], y2))
}
(eval2CV <- round(colMeans(vec), 5) ) # average prediction based on the 5-fold cross validation



if(saveFig) jpeg(paste0(saveDir,'PP_BRT.jpg'), quality=100, height=1800,width=3000, pointsize=14, res=300)

lmat <- matrix(1:8,nrow = 2,ncol = 4,byrow = TRUE)
layout(mat = lmat, heights = rep(.5,2))
par(oma=c(0,2,2,0), font.lab=1, mar=c(2,.1,.1,.7)+.4,mgp=c(1.2,.2,0), tck=-0.02, cex.axis=1, cex.lab=1.1, cex.main=1.2 )
source('gbm.plot2.R')
Xnam <- names(mat)[-1]
gbm.plot2(m2, plot.layout=c(2,4), variable.no=0, write.title=FALSE, x.label=Xnam, y.label="")

par(mar=c(0,0,0,0)); plot(1, type = "n", axes=FALSE, xlab="", ylab="")  #fill the 8th panel since we only have 7 vars
title('Boosted regression trees', outer=T, font.main=2, cex.main=1.5)
mtext("Fitted function", side=2, line=.6, outer=TRUE, cex=1, font=2)

if(saveFig) dev.off()




#*********************************************** XGBoost

require('xgboost')
nrounds <- 1000; #1e3
a <- data.matrix(mat[,-1])
b <- as.numeric(mat$NJIS)

set.seed(8)
m3 <- xgboost(data = a, label = b, nrounds = nrounds, objective = "reg:linear")
y3 <- predict(m3, a)
(eval3 <- evalPred(mat$NJIS, y3) )


# get predictive measures for 5-fold cross validation
n <- nrow(mat) # sample size
vec <- numeric(); K <- 5
for(k in 1:K){
  inds <- 1:n
  modelRF <- errR$models[[k]]
  indTrain <- as.integer(row.names(as.data.frame(modelRF$predicted)))
  tmp <- droplevels(mat[indTrain,])
  d0 <- data.matrix(tmp[,-1])
  e0 <- data.matrix(tmp[,1])
  model <- xgboost(data = d0, label = e0, nrounds = nrounds, objective = "reg:linear")
  indTest <- inds[- indTrain]
  tmp <- droplevels(mat[indTest,])
  d1 <- data.matrix(tmp[,-1])
  e1 <- tmp[,1]
  y3tmp <- predict(model, d1)
  vec <- rbind(vec, evalPred(e1, y3tmp))
}
(eval3CV <- round(colMeans(vec), 5) )  # average prediction based on the 5-fold cross validation






#++++++++++++++++++++ combine GBM and RF to plot the variable importance
foo <- m2$contributions; row.names(foo) <- NULL
foo1 <- data.frame(var=row.names(imp), IncMSE=imp[,1,drop=F]); row.names(foo1) <- NULL
foo <- merge(foo, foo1)

xlabs <- c('Relative influence', 'Mean decrease in accuracy')
mains <- c('Boosted regression trees', 'Random forests')

if(saveFig) jpeg(paste0(saveDir,'importanceBRTnRF.png'), quality=100, height=1800,width=3000, pointsize=14, res=300)

par(mfrow=c(1, 2), mar=c(2,1.3,.8,.1)+.4,mgp=c(1.2,.3,0), tck=-0.03, cex.axis=1, cex.lab=1, cex.main=1)
for(i in 1:2){
  inds <- order(foo[,i+1])
  dotchart(foo[inds,i+1],labels=foo[inds,1], pch=16, pt.cex=1.5, xlab=xlabs[i], ylab = "", main=mains[i])
}

if(saveFig) dev.off()




# check GDD_diff

pdf(file=paste('./results/raw/simu/GDD',case,'.pdf',sep=''),pointsize=13,width=6.7,height=3)

par(mfrow=c(1, 3), mar=c(2,2,.8,.1)+.4,mgp=c(1,.1,0), tck=-0.01, cex.axis=1, cex.lab=1.1, cex.main=1.15)
source('gbm.plot2.R')
gbm.plot2(m2, plot.layout=c(1,2), variable.no=33, write.title=FALSE, y.label="Fitted function", x.label='GDD_Diff',main=mains[1], show.contrib=F)
partialPlot(x=m1, pred.data=mat, x.var='V33', col='blue', lwd=2, main=mains[2], xlab='GDD_Diff', ylab='Partial dependence') 
grid(col='gray40')
plot(y~V33, data=a, cex=1.8, xlab='GDD_Diff', ylab='Net distance (log)', main='Data')
# abline(lm(y~V33, data=a), col='blue', lwd=2 )
grid(col='gray40')

dev.off()


# pdf(file=paste('./tex/gdd1.pdf',sep=''),pointsize=13,width=6,height=3)
# a <- mat
# par(mfrow=c(1,1), mar=c(1.8,1.8, .2,.1)+.4,mgp=c(1.1,.3,0), tck=-0.03, cex.axis=.8, cex.lab=1, cex.main=1)
# plot(y~V33, data=a, ylab='log-Distance', xlab='GDD_Diff')
# dev.off()




#******************************************* get table
evalFit <- rbind(eval1, eval2, eval3)
row.names(evalFit) <- c('Random Forests', 'Boosted Regression Tree', 'Extreme Gradient Boosting')
print(evalFit)

evalCV <- rbind(eval1CV, eval2CV, eval3CV)
row.names(evalCV) <- c('Random Forests', 'Boosted Regression Tree', 'Extreme Gradient Boosting')
print(evalCV)


save(file=paste0('out',datId), evalFit, evalCV)

# not run














