
# codes for NTP data challege on early-phase toxicity prediction
rm(list=ls())
require(Metrics)
library(moments)
library(MASS)
require(xgboost)
require(randomForest)
windowsFonts(Times=windowsFont("Times New Roman")) #some journals ask for the font in plot


#########################################
# step1. preprocessing data
# extract featureID, delete observations with y=na


# mat <- read.csv('./trainingset_171127_DW_MOE_More_Complete.csv')  #8994  264
# 
# X <- read.csv('./trainingset_171127_DW_moe.csv')  #8994 2597, new features start from col 258
# 
# mat$Molecule.Name <- as.character(mat$Molecule.Name)
# X$Molecule.Name <- as.character(X$Molecule.Name)
# length(unique(mat$Molecule.Name)) == nrow(mat)  #TRUE
# length(unique(X$Molecule.Name)) == nrow(X)  #TRUE
# all(mat$zagreb==X$zagreb)
# all(mat$Molecule.Name %in% X$Molecule.Name)
# all(X$Molecule.Name %in% mat$Molecule.Name)
# mat$Molecule.Name[!mat$Molecule.Name%in%X$Molecule.Name]
# # X[!X$Molecule.Name%in%mat$Molecule.Name, c('Molecule.Name','Molweight')]
# # mat[!mat$Molecule.Name%in%X$Molecule.Name, c('Molecule.Name','Molweight')]
# 
# inds <- which(!mat$Molecule.Name%in%X$Molecule.Name)
# 
# recode <- function(x){  #2/1/1918 -> 1918-02-1
#   ind1 <- unlist(gregexpr('/',x))
#   y <- substr(x,1,ind1[1]-1); if(nchar(y)==1) y <- paste0('0',y)
#   y <- paste0(substr(x,ind1[2]+1,100L),'-',y,'-',substr(x,ind1[1]+1,ind1[2]-1))
#   return(y)
# }
# 
# tmp <-  sapply(mat$Molecule.Name[inds], recode)
# all(tmp %in% X$Molecule.Name[!X$Molecule.Name%in%mat$Molecule.Name])
# mat$Molecule.Name[inds] <- tmp
# all(mat$Molecule.Name %in% X$Molecule.Name)  #TRUE
# 
# mat <- merge(mat, X[,c(2, 258:ncol(X))], by='Molecule.Name')  # 8994 2604
# 
# head(cbind(mat$LD50_mgkg.x, mat$LD50_mgkg.y))
# 
# save(file='mat', mat)
load('mat')

# exclude some redundant variables due to merge
mat <- mat[,-which(names(mat) %in% c('GHS_category.y','ADMET_BBB','EPA_category.y',
                                     'Nasty.Functions.y','LD50_mgkg.y') )]  #now 8994 2599
names(mat)[which(names(mat)=='LD50_mgkg.x')] <- 'LD50_mgkg'


mat[1:5,c('TOPKAT_Rat_Oral_LD50','LD50_mgkg')]
plot(mat[,c('TOPKAT_Rat_Oral_LD50','LD50_mgkg')])  #TOPKAT_Rat_Inhalational_LC50


# After comparing the two dataset, there are two diffs:
# 1) row number decreased
# 2) a new column named "Stereo.Configuration"
dim(mat)
mat <- mat[, c(which(names(mat) %in% c('CASRN','InChI.Key_QSARr','LD50_mgkg')), 17:ncol(mat))]
mat <- mat[, -c(which(names(mat) %in% c('LE.from.LD50_mgkg','LELP.from.LD50_mgkg','LLE.from.LD50_mgkg')))]
LD50_mgkg.missing.inds <- which(is.na(mat$LD50_mgkg)); 

# delete features with missing values (no such feature)
if(length(LD50_mgkg.missing.inds)) mat <- mat[-LD50_mgkg.missing.inds,]
LD50_mgkg.missing.data = mat[LD50_mgkg.missing.inds,]
(p <- length(featureInd <- 4:ncol(mat)))
(tab <- table(nmiss <- apply(mat[,featureInd], 2, function(x) sum(is.na(x)))))
names(which(nmiss==as.numeric(names(tab[2]))))  
#only Relative.PSA has two missing values
mat[is.na(mat$Relative.PSA), 1:10]

# Need to remove some suspicious data points: 
#          CASRN             InChI.Key_QSARr LD50_mgkg Molweight Monoisotopic.Mass cLogP cLogS H.Acceptors
# 7421   74-82-8 VNWKTOKETHGBQD-UHFFFAOYSA-N      2000   16.0428           16.0313     0 -0.53           0
# 7460 7440-44-0 VNWKTOKETHGBQD-UHFFFAOYSA-N      2000   16.0428           16.0313     0 -0.53           0
mat <- droplevels(mat[-which(is.na(mat$Relative.PSA)), ])


#remove redundant features which have only one unique value
noVariation <- apply(mat[,featureInd], 2, function(x) length(unique(x))==1)
# mat <- mat[,-featureInd[which( nmiss>0 |  noVariation )]]
mat <- mat[,-featureInd[which( noVariation )]]


table(mat$ADMET_BBB_Level)


#deal with categorical variables
(p <- length(featureInd <- 4:ncol(mat)))
table(types <- sapply(mat[,featureInd], class))
isFactor <- names((which(types=='factor')))
mat = mat[,-which(names(mat) %in% isFactor[5:6])]
mat[,isFactor[1]] = as.numeric(mat[,isFactor[1]] )
mat[,isFactor[2]] = as.numeric(mat[,isFactor[2]] )
mat[,isFactor[3]] = as.numeric(mat[,isFactor[3]] )
mat[,isFactor[4]] = as.numeric(mat[,isFactor[4]] )

(p <- length(featureInd <- 4:ncol(mat))) # update number of features and index
tab <- table(mat$InChI.Key_QSARr)
(repeatedInChI <- names(tab)[tab>1]) #check those with repeated InChIKey
# mat <- droplevels(mat[-which(mat$InChI.Key_QSARr %in% repeatedInChI), ])

#+++ since now we don't use InChI for prediction, no need to drop

# (p <- length(featureInd <- 3:ncol(mat))) # update number of features and index
# table(types <- sapply(mat[,featureInd], class))
# isNumeric <- names((which(types=='numeric')))

n <- nrow(mat)
row.names(mat) <- 1:n



##################################################
myparam <- list(objective = "reg:linear",metrics =  "rmse",max_depth = 10,
                eta = .02, gamma = 0.144, subsample = 0.78,
                min_child_weight = 0.18, colsample_bytree = 0.54,nthread=4)
myseed <- 692
set.seed(123)
inds <- sample(1:10,size=n,prob=rep(0.1,each=10),replace=T)
LD50 <- mat$LD50_mgkg

nams <- names(mat)
# write.csv(file='tmp.csv', mat[1:5,c('LD50_mgkg','LD50_mgkg.y')])
# mat[1:5,c('LD50_mgkg','LD50_mgkg.y')]

best_fold_loss <- Inf
best_id <- 0
lossvec <- rep(NA, 10)
for(i in 1:10){
  
  test_id <- which(inds==i)
  id <- c(1:n)
  train_id <- id[-test_id]
  train.x <- mat[train_id,featureInd]
  train.x <- train.x[!names(train.x) %in% c('LD50_mgkg')]
  train.x <- data.matrix(train.x)
  train.y <- log(1+LD50[train_id])
  
  test.x <- mat[test_id,featureInd]
  test.x <- test.x[!names(test.x) %in% c('LD50_mgkg')]
  test.x <- data.matrix(test.x)
  test.y <- log(1+LD50[test_id])
  
  dtrain <- xgb.DMatrix(data = train.x, label = train.y)
  set.seed(myseed)
  model <- xgb.train(data=dtrain, params = myparam, nrounds = 1000)
  test.pred <- predict(model,test.x)
  loss <- rmse(test.pred,test.y)
  
  lossvec[i] <- loss
  
  if(loss<best_fold_loss){
    best_fold_loss <- loss
    best_id <- i
    best_model <- model
  }
  
  cputime <- as.numeric(proc.time()[3]-ptm)/60 
  cat('\nCPUtime', round(cputime,3), 'minutes: completed!','\n')
  
  # RF CPUtime 118.531 minutes: completed!  loss=1.336219
  ptm <- proc.time()[3] 
  set.seed(20)
  tmp <- mat[train_id,featureInd]; tmp$y <- train.y
  m <- randomForest(y~., data=tmp)
  cputime <- as.numeric(proc.time()[3]-ptm)/60 
  cat('\nCPUtime', round(cputime,3), 'minutes: completed!','\n')
  
  yhat <- predict(m, mat[test_id,featureInd])
  lossRF <- rmse(yhat, test.y)
  
  head(cbind(yhat, test.y, test.pred))
  summary(lm(yhat~test.y))$r.squared
  summary(lm(test.pred~test.y))$r.squared
  
}

model <- best_model

save(file='model', model,loss)

save(file='modelRF', m,loss)

# not run

