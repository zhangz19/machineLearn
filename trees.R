
rm(list=ls())
options(scipen = 500)
library(R.matlab)

# calculate the prediction measure
### calculating measures for prediction
evalPred <- function(y,x){ #y=observed, x=fitted
  ym <- mean(y)
  xm <- mean(x)
  RMSE <- sqrt(mean((y-x)^2))
  NSE <- 1- sum((y-x)^2)/sum((y-ym)^2)
  PBIAS <- 100 * sum(y-x)/sum(y)
  z <- c(RMSE,NSE,PBIAS); names(z) <- c('RMSE','NSE','PBIAS')
  return(z)
}

# import your data object called mat

# RF: random forests approach
require('randomForest')

set.seed(12)
m1 <- randomForest(y~., data=mat, importance=T)

evalPred(m1$y, predict(m1))

foo1 <- importance(m1, scale=TRUE)  #mean decrease in accuracy,  mean decrease in node impurity

# # random forests with variable selection: optional, can be slow
# require(caret)
# con <- rfeControl(functions=rfFuncs, method='CV', number=10)
# set.seed(12)
# m11 <- rfe(y~., data=mat, rfeControl=con)
# myForm <- as.formula(paste('y~',paste(m11$optVariables, collapse='+')))
# set.seed(12)
# m12 <- randomForest(myForm, data=mat, importance=T)
# evalPred(m1$y, predict(m12))


# # save importance
pdf(file=paste('./RFimp.pdf',sep=''), pointsize=9,width=7,height=5,paper='special')

source('myvarImpPlot.R')
myvarImpPlot(m12, sort = T, main="Variable Importance", pch=16, n.var=p)

dev.off()



# save partial dependence
imp <- importance(m1)
impvar <- rownames(imp)[Is <- order(imp[, 'IncNodePurity'], decreasing=TRUE)]
labs <- Xnam[Is]

pdf(file=paste('./RFpd.pdf',sep=''),pointsize=13,width=7,height=8,paper='special')

par(mfrow=c(5, 4), mar=c(1.3,1.3,1.2,.1)+.4,mgp=c(1.4,.3,0), tck=-0.03, cex.axis=1, cex.lab=1.3, cex.main=1.3)
for (i in 1:20){ # top 20
  print(i)
  partialPlot(x=m1, pred.data=mat, x.var=impvar[i], xlab='', col='blue', lwd=2, main=labs[i])
  grid(col='gray70')
}

dev.off()



# GBT: Boosted regression tree
require(dismo)

set.seed(10)
m2 <- gbm.step(data=mat, gbm.x=1:p, gbm.y=p+1, family='gaussian')  # this is too slow!

a <- summary(m2, plotit=T)

row.names(a) <- NULL
print(a)  # this is variable importance

evalPred(m1$y, fitted(m2))


pdf(file=paste('./BRT.pdf',sep=''),pointsize=13,width=7.4,height=9,paper='special')

source('gbm.plot2.R')
gbm.plot2(m2, plot.layout=c(9,4), variable.no=0, write.title=FALSE, x.label=Xnam)

dev.off()




