
rm(list=ls())
options(scipen = 500)
library(R.matlab)

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



o <- readMat('yan.mat')
load('Xnam')

case <- 2

mat <- as.data.frame(o$X); names(mat) <- Xnam
if(case==1) mat <- mat[, 1:33] else mat <- mat[, -33] 
mat$y <- o$Y

Xnam[Xnam=='tag_length'] <- 'length'
Xnam[Xnam=='years_diff'] <- 'years_lag'
# Xnam[Xnam=='d5_tag'] <- 'GDD'
Xnam[Xnam=='Dip'] <- 'Diporeia'
Xnam[Xnam=='sexFemale'] <- 'sex: Female'
Xnam[Xnam=='GDD1'] <- 'GDD_Diff'
Xnam[Xnam=='GDD2_M9'] <- 'Sep*GDDlake'
Xnam[Xnam=='GDD2_M10'] <- 'Oct*GDDlake'
Xnam[Xnam=='GDD2_M11'] <- 'Nov*GDDlake'
Xnam[Xnam=='GDD2_M12'] <- 'Dec*GDDlake'
Xnam <- gsub('rec_Y','rec_Y: ',Xnam)
Xnam <- gsub('tag_Y','tag_Y: ',Xnam)
Xnam <- gsub('rec_M','rec_M: ',Xnam)
Xnam <- gsub('tag_site','tag_site: ',Xnam)
Xnam <- paste(Xnam, " ") # add some gap from axis

if(case==1) Xnam <- Xnam[1:33] else Xnam <- Xnam[-33] 

p <- length(Xnam)



names(mat) <- c( paste0('V',1:p),  'y')


# RF
require('randomForest')

set.seed(12)
m1 <- randomForest(y~., data=mat, importance=T)

evalPred(m1$y, predict(m1))

foo1 <- importance(m1, scale=TRUE)  #mean decrease in accuracy,  mean decrease in node impurity

# # random forests with variable selection
# require(caret)
# con <- rfeControl(functions=rfFuncs, method='CV', number=10)
# set.seed(12)
# m11 <- rfe(y~., data=mat, rfeControl=con)
# 
# myForm <- as.formula(paste('y~',paste(m11$optVariables, collapse='+')))
# set.seed(12)
# m12 <- randomForest(myForm, data=mat, importance=T)
# evalPred(m1$y, predict(m12))




# aa <- readMat(paste0('./results/raw/paras',case,'.mat') )
# evalPred(m1$y, aa$yhat0)
# evalPred(m1$y, aa$yhat1)



# load BVS results for comparison
load('BVSout')  #vs
a <- vs[[case]]
all(a[,1] == Xnam)  #must be true
XnamP <- gsub(" ","",Xnam)
XnamP <- paste0(Xnam, '(',a$Selectivity,')')


# # save importance
# pdf(file=paste('./tex/RFimp',case,'.pdf',sep=''),pointsize=9,width=7,height=5,paper='special')
# 
# source('myvarImpPlot.R')
# myvarImpPlot(m12, sort = T, main="Variable Importance", pch=16, n.var=p)
# 
# dev.off()



# # save partial dependence
# imp <- importance(m1)
# # impvar <- rownames(imp)[Is <- order(imp[, '%IncMSE'], decreasing=TRUE)]
# impvar <- rownames(imp)[Is <- order(imp[, 'IncNodePurity'], decreasing=TRUE)]
# labs <- Xnam[Is]
# 
# pdf(file=paste('./tex/RFpd',case,'.pdf',sep=''),pointsize=13,width=7,height=8,paper='special')
# 
# # for (i in seq_along(impvar)){
# par(mfrow=c(5, 4), mar=c(1.3,1.3,1.2,.1)+.4,mgp=c(1.4,.3,0), tck=-0.03, cex.axis=1, cex.lab=1.3, cex.main=1.3)
# for (i in 1:20){ # top 
#   print(i)
#   partialPlot(x=m1, pred.data=mat, x.var=impvar[i], xlab='', col='blue', lwd=2, main=labs[i]) #, ylim=c(36, 52)
#   grid(col='gray70'); 
#   # partialPlot(m1, mat, impvar[i], col='blue', lwd=2, add=T)
# }
# 
# # partialPlot(m1, mat, impvar[1],col='blue', lwd=2)
# 
# dev.off()



# GBT
require(dismo)

set.seed(10)
m2 <- gbm.step(data=mat, gbm.x=1:p, gbm.y=p+1, family='gaussian')  # this is too slow!

a <- summary(m2, plotit=T)

# # dev.off()
# row.names(a) <- NULL
# print(a)  # this is variable importance
# 
# evalPred(m1$y, fitted(m2))


# pdf(file=paste('./tex/BRT',case,'.pdf',sep=''),pointsize=13,width=7.4,height=9,paper='special')
# 
# source('gbm.plot2.R')
# gbm.plot2(m2, plot.layout=c(9,4), variable.no=0, write.title=FALSE, x.label=Xnam)
# 
# dev.off()


# combine GBM and RF
foo <- m2$contributions
foo$var <- as.numeric(substr(foo$var, 2, 100L))
foo <- foo[order(foo$var),]
foo$var <- XnamP
# row.names(foo) <- XnamP
foo$impRF <- foo1[,1]  #%IncMSE
# foo <- foo[,-1]


mains <- c('Boosted regression trees', 'Random forests')


pdf(file=paste('./results/raw/simu/tree',case,'.pdf',sep=''),pointsize=9,width=7,height=4.8)

par(mfrow=c(1, 2), mar=c(2,1.3,.8,.1)+.4,mgp=c(1.2,.3,0), tck=-0.03, cex.axis=1, cex.lab=1, cex.main=1)
xlabs <- c('Relative influence', 'Prediction accuracy')
for(i in 1:2){
  inds <- order(foo[,i+1])
  dotchart(foo[inds,i+1],labels=foo[inds,1], pch=16, pt.cex=1.5, xlab=xlabs[i], ylab = "", main=mains[i])
}

dev.off()




# check GDD_diff
if(case==1){
  
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
  
}

save(file='tree2', m1, m2)



# pdf("out_gbm.pdf",  width = 7, height = 9, pointsize=11)
# 
# source('gbm.plot2.R')
# gbm.plot2(m2, plot.layout=c(7,6), variable.no=0, write.title=FALSE, x.label=Xnam)
# 
# save(file='tree2', m1, m2)
# 
# dev.off()
# 
# plot(y~GDD1, data=mat)
# plot(y~GDD1, data=mat, xlim=c(-.2,2))
# 
# mat[mat$GDD1>9, ]
# 
# # set.seed(20)
# # m2 <- gbm.step(data=mat,  gbm.x=1:p, gbm.y=p+1, family='gaussian', tree.complexity = 5, learning.rate = 0.01, bag.fraction = 0.5, n.folds=5, keep.fold.fit=TRUE, keep.fold.vector=TRUE)



