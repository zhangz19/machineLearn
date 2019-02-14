

require(unbalanced)
require(randomForest)
require(xgboost)

#### load data:  tdat, tlab, tnam, indDW, indMOE, propYes
load('demo.rda')

names(tdat) <- paste0('x', 1:ncol(tdat))  #some names have sep which may not work well for RF in R

# partition jobs that explore all optimal combinations of tuning parameters
ns <- length(svec <- seq(0,500,by=100))
p0 <- 100*mean(tlab=='Yes') #6.5 #%  mean(mat$y)
nr <- length(rvec <- seq(p0, 20,length=9))
K <- 10  #K-fold cross validation
nb <- length(bvec <- 1:K)  #each bin of the K-fold cross validation
nd <- length(dvec <- 1:3)  #feature sources: 3 data sets

Svec <- rep( rep(1:ns, nr), nb*nd)
Rvec <- rep( rep(1:nr, each=ns), nb*nd)
Bvec <- rep( rep(bvec, each=nr*ns), nd)
Dvec <- rep( rep(dvec, each=nr*ns), each=nb)
# can be more decent if using merge function

# repl is the job ID which presumably runs over all possible combinations of tuning parameters for the pre-balancing + model.
repl <- as.numeric(commandArgs(TRUE))

sid <- Svec[repl]
ps <- svec[sid]  #percentage for SMOTE
rid <- Rvec[repl]
pr <- rvec[rid]  #percentage for random subsampling
binID <- Bvec[repl]  #bin ID for K-CV
datID <- Dvec[repl]  #data ID

cat('datID=',datID,', binID=',binID,', sid=',sid,', rid=',rid,':\n',sep='')

# function to evaluate the performance
stab <- function(y, yhat, returnObj=FALSE){
  # yhat <- predict(m)
  tab <- table(y, yhat); acc <- sum(diag(tab))/sum(tab)
  n <- sum(tab); Pd <- diag(tab)/n; Pc <- colSums(tab)/n; Pr <- rowSums(tab)/n
  P0 <- sum(Pd); Pe <- sum(Pc * Pr)
  kappa <- (P0 - Pe) / (1 - Pe) ; #cat("Cohen's kappa =",round(kappa,3),"\n")
  tab <- cbind(tab, rowSums(tab)); tab <- rbind(tab, colSums(tab))
  TP <- tab[2,2]; FP <- tab[1,2]; TN <- tab[1,1]; FN <- tab[2,1]
  phi <- exp(log(TP*TN - FP*FN) - 0.5*(log(TP+FN)+log(TN+FP)+log(TP+FP)+log(TN+FN))) # this is Matthews correlation coefficient
  f1 <- 2*TP/(2*TP + FP + FN)
  # pred <- prediction(m$votes[,2], y) #require(ROCR)
  # perf_AUC <- performance(pred,"auc") #Calculate the AUC value
  # auc <- perf_AUC@y.values[[1]]
  # perf_ROC <- performance(pred,"tpr","fpr")
  # plot(perf_ROC, main="ROC plot")
  # text(0.5,0.5,paste("AUC = ",format(AUC, digits=5, scientific=FALSE)))
  # print(tab)
  # out <- data.frame(acc=acc, kappa=kappa, f1=f1, phi=phi) #, auc=auc
  out <- c(acc, kappa, f1, phi)
  if(returnObj) return(out) else print(round(out,3))
}


#++++++++++++++ for RF classifier, if > nRF: for xgboost
nt <- length(ntrees <- c(50,100,200, 500, 1000, 1, 2, 3) )  
nRF <- 5  #match with nt


# 2. standardization?
standardization <- TRUE
if(standardization){
  # if(any(ls()=='vdat')) tmp <- scale(rbind(tdat, vdat)) else tmp <- scale(tdat)
  #Scale the training data even there is validation data. Make the training process not involving any validation data. 
  tmp <- scale(tdat)
  
  inds <- 1:nrow(tdat)
  tdat <- tmp[inds,]   #if(any(ls()=='vdat')) vdat <- tmp[-inds, ]
}

mat <- tdat 
p <- ncol(mat)
mat <- as.data.frame(mat); mat$y <- tlab


# stratified validation to guarantee we have enough samples for both classes
# same random seeds for separating the data
set.seed(20)

n0 <- length(ind0 <- which(mat$y=="No"))
n1 <- length(ind1 <- which(mat$y=="Yes"))

# stratified cross validation
m0 <- ceiling(n0/K)
inds0 <- rep(NA, m0*K); inds0[1:n0] <- sample(ind0) #random permutation
inds0All <- t(matrix(inds0, K, m0))

m1 <- ceiling(n1/K)
inds1 <- rep(NA, m1*K); inds1[1:n1] <- sample(ind1) #random permutation
inds1All <- t(matrix(inds1, K, m1))


# random seeds now vary
set.seed(myseed <-  repl*10)

ptm <- proc.time()[3]
print('CV begins...')

metricsV <- metricsT <- array(NA, dim=c(K,4,nd,nt), 
                              dimnames=list(paste0('CV',1:K), c('acc', 'kappa', 'f1', 'phi'), 
                                            c('DW','MOE','Combined'), ntrees))

for(i in datID:datID){  #1:nd: now it fits the data with descriptor set indicated by datID
  xinds <- c(indDW, indMOE)  #1:ncol(mat)  #for i==3: combined features
  if(i==1) xinds <- indDW  #use Data Warrior descriptors only
  if(i==2) xinds <- indMOE  #use MOE descriptors only
  xinds <- c(xinds, which(names(mat)=='y'))
  
  for(k in binID:binID ){ #1:K: now it fits the data in bin id indicated by binID
    # print(k)
    
    vinds <- na.omit(c(inds0All[,k], inds1All[,k]))
    vd <- mat[vinds, xinds]; td <- mat[-vinds, xinds]
    X <- td[,-which(names(td)=='y')]
    tmp <- list(X=as.data.frame(X), Y=as.factor(c(0,1)[as.integer(td$y)]) )
    
    # U-step
    if(rid!=1) tmp <- ubUnder(tmp$X, tmp$Y, perc=pr, method="percPos")
    # S-step
    if(sid!=1){
      tmp1 <- ubSMOTE(tmp$X, tmp$Y, k=5, perc.over=ps, perc.under=0)
      tmp$X <- rbind(tmp$X, tmp1$X)
      tmp$Y <- factor(levels(tmp$Y)[c(tmp$Y, tmp1$Y)], levels=levels(tmp$Y))
    }
    # E-step
    tmp <- ubENN(tmp$X, tmp$Y, k=3)  #removed from majority class 
    # tmp <- ubTomek(tmp$X, tmp$Y)  #from both classes
    
    mat2 <- as.data.frame(tmp$X)
    mat2$y <- as.factor(tmp$Y); levels(mat2$y) <- c('No','Yes')
    
    # require(e1071)
    # m.nb <- naiveBayes(y ~ ., data=mat2)
    # table(vd$y, predict(m.nb, vd))
    
    for(j in 1:nt){  # for all models (e.g. different numbers of trees for RF)
      cat(j,"")
      
      if(j<=nRF){  #random forests
        m4 <- randomForest(y~., data=mat2, ntree=ntrees[j], importance=T)
        y1 <- predict(m4, newdata=vd )
        # table(vd$y, y1)
        metricsT[k,,i,j] <- stab(mat2$y, predict(m4), returnObj=TRUE) # this is OOB
        metricsV[k,,i,j] <- stab(vd$y, y1, returnObj=TRUE)
      }else{  #j>nRF: for xgboost
        if(ntrees[j]==1){  #default setup
          # I think nthread does not matter as we run in single thread without invoking openMP
          myparam <-list(objective = "binary:logistic", metrics="error", 
                         max_depth=6, eta=.3, gamma=0.144, subsample=1, 
                         min_child_weight=1, colsample_bytree=1, nthread=2)
        }
        if(ntrees[j]==2){  #setup for NTP data challenge
          myparam <-list(objective = "binary:logistic", metrics="error", 
                         max_depth=10, eta=.02, gamma=0.144, subsample=0.78, 
                         min_child_weight=0.18, colsample_bytree=0.54, nthread=4)
        }
        if(ntrees[j]==3){  #setup inbetween
          myparam <-list(objective = "binary:logistic", metrics="error", 
                         max_depth=8, eta=.16, gamma=0.2, subsample=0.9, 
                         min_child_weight=0.5, colsample_bytree=.8, nthread=2)
        }
        train_x <- data.matrix(mat2[, -which(names(mat2)=='y')])
        test_x <- data.matrix(vd[, -which(names(vd)=='y')])
        train_y <- as.numeric(mat2$y) - 1  #[1,2] -> [0,1]
        test_y <- as.numeric(vd$y) - 1 
        dtrain <- xgb.DMatrix(data=train_x, label=train_y)
        mtmp <- xgb.train(data=dtrain, params=myparam, nrounds=1000)
        metricsT[k,,i,j] <- stab(train_y, predict(mtmp, train_x)>0.5, returnObj=TRUE)
        metricsV[k,,i,j] <- stab(test_y, predict(mtmp, test_x)>0.5, returnObj=TRUE)
      } 
      
    }; cat('\n')
    
  }
}

cputime <- as.numeric(proc.time()[3]-ptm)
cputime <- cputime/60
cat('\nCPUtime', round(cputime,3), 'minutes: completed!','\n')


save(file=paste('out',repl,sep=''), metricsT, metricsV, ntrees, myseed, cputime, datID, binID, rid, sid)



combinIt <- FALSE
if(combinIt){
  
  load('demo.rda')
  ns <- length(svec <- seq(0,500,by=100))
  p0 <- 100*mean(tlab=='Yes')
  nr <- length(rvec <- seq(p0, 20,length=9))
  K <- 10  #K-fold cross validation
  nb <- length(bvec <- 1:K)  #each bin
  nd <- length(dvec <- 1:3)  #feature sources: 3 data sets
  
  Svec <- rep( rep(1:ns, nr), nb*nd)
  Rvec <- rep( rep(1:nr, each=ns), nb*nd)
  Bvec <- rep( rep(bvec, each=nr*ns), nd)
  Dvec <- rep( rep(dvec, each=nr*ns), each=nb)
  
  nt <- length(ntrees <- c(50, 100, 200, 500, 1000, 1, 2, 3) ) 
  
  metricsV_all <- metricsT_all <- array(NA, dim=c(nr, ns, K, 4, nd, nt))
  nsim <- nr*ns*nb*nd
  for(repl in 1:nsim){
    cat(repl,' '); if(!repl%%20) cat('\n')
    if(file.exists(fnam <- paste('out',repl,sep=''))){
      load(fnam)
      metricsV_all[rid, sid, binID,,datID,] <- metricsV[binID,,datID,]  #3=f, 4=phi(mcc)
      metricsT_all[rid, sid, binID,,datID,] <- metricsT[binID,,datID,] 
    }
  }
  
  dimnames(metricsV_all) <- list(rvec, svec, 1:K, c('acc', 'kappa', 'f1', 'phi'), 
                                 c('DW','MOE','Combined'), ntrees)
  
  metricsV_all_F1 <- metricsV_all[,,,3,,]
  metricsV_all_F1 <- apply(metricsV_all_F1, c(1,2,4,5), mean, na.rm=T) #average over CV bins
  
  nRF <- 5
  # tmp <- metricsV_all_F1  #best overall 
  tmp <- metricsV_all_F1[,,,1:nRF]  #best for RF
  # tmp <- metricsV_all_F1[,,,-c(1:nRF)] #best for xgboost
  f1max <- tmp[which.max(tmp)] 
  oind <- arrayInd(which.max(tmp), dim(tmp))
  optcase <- mapply(`[[`, dimnames(tmp), oind) 
  names(optcase) <- c('RU','SMOTE','Data','ntree')
  
  print(c(optcase, f1max=f1max))
  
  save(file='sout', metricsV_all_F1, optcase, f1max)
  
}


# not run

