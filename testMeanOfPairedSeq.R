

#+++++ codes for paper: 
# "Evaluation of the Impacts of Hydrologic Model Calibration Methods on Predictability of Ecologically-Relevant Hydrologic Indices"
# published in Journal of Hydrology, 2018: https://www.sciencedirect.com/science/article/pii/S0022169418305705

# functions to test equivalence of paired sequences (can be irregularly spaced time series)
# The data was extracted from the time series but specific dates were selected that are not necessarily continuous.

rm(list=ls())
# require(heavy)  # for heavyLM
require(MASS)  #for fitdist
# require(actuar)  #for Pareto when data>0
# require(fitdistrplus)  #fitdist
require(mgcv)
require(nlme)

runIt <- TRUE
runItLocally <- TRUE
runParallel <- FALSE
jobID <- 1; if(runParallel) jobID <- as.numeric(commandArgs(TRUE))

evalPred <- function(y,x){ #y=observed, x=fitted
  ym <- mean(y)
  xm <- mean(x)
  RMSE <- sqrt(mean((y-x)^2))
  NSE <- 1- sum((y-x)^2)/sum((y-ym)^2)
  PBIAS <- 100 * sum(y-x)/sum(y)
  z <- c(RMSE,NSE,PBIAS); names(z) <- c('RMSE','NSE','PBIAS')
  return(z)
}

# fit distributions via assuming the samples are independent
fitD <- function(y, den, prob=c(0.25, 0.75), plotIt=FALSE){
  paras <- NA
  if(den=='normal'){
    m <- fitdistr(y, densfun='normal')
    if(plotIt){
      qqplot(qnorm(ppoints(length(y)), mean=coef(m)[['mean']], sd=coef(m)[['sd']]), y)
      qqline(y, distribution = function(p) qnorm(p, mean=coef(m)[['mean']], sd=coef(m)[['sd']]), prob=prob, col = 2)
    }
  }
  if(den=='t'){
    m <- try(fitdistr(y, densfun='t'))
    if(class(m)=='try-error'){
      # m <- fitdist(y, "t", method="mle",start=list(df=10))
      m <- fitdistr(y, densfun='t',method='SANN')
    }
    if(plotIt){
      qqplot(qt(ppoints(length(y)), df=coef(m)[['df']]), y)
      qqline(y, distribution = function(p) qt(p, df=coef(m)[['df']]), prob = prob, col = 2)
    }
    paras <- coef(m)[['df']]
  }
  if(den=='cauchy'){
    m <- fitdistr(y, densfun='cauchy')
    qqplot(qcauchy(ppoints(length(y)), location=coef(m)[['location']], scale=coef(m)[['scale']]), y)
    qqline(y, distribution = function(p) qcauchy(p, location=coef(m)[['location']], scale=coef(m)[['scale']]), prob = prob, col = 2)
  }
  if(den=='pareto'){  # or generalized pareto in Package 'gPdtest', but note, data should > 0
    #empirical raw moment
    memp <- function(x, order) mean(x^order)
    m <-  fitdist(y, "pareto", method="mme", order=c(1, 2), memp="memp", start=list(shape=10, scale=10), lower=1, upper=Inf)
    qqplot(qpareto(ppoints(length(y)), shape=coef(m)[['shape']], scale=coef(m)[['scale']]), y)
    qqline(y, distribution = function(p) qpareto(p, shape=coef(m)[['shape']], scale=coef(m)[['scale']]), prob = prob, col = 2)
  }
  
  noDifference <- ifelse(prod(confint(m)[1,])>0, NA, 'TRUE')
  return(data.frame( df=paras, mean=coef(m)[1], lower=confint(m)[1,1], upper=confint(m)[1,2], AIC=AIC(m), noDifference=noDifference ))  #assume the first row is mean or location
}


loop <- 0 
for(nopt in 1:2){
  dir0 <- paste0('./Data/Option ',nopt,'/')
  dirs <- dir(path=dir0, pattern='SIM*')
  cate <- substr(dirs, 17, 17)
  (nf <- length(dirs))
  
  out <- numeric()
  for(i in 1:nf){
    dat <- t(read.csv(fnam <- paste0(dir0,dirs[i]), h=F))
    
    dates <- dates0 <- read.csv(fnam <- paste0(dir0, 'Time_output_Opt',nopt,'_',cate[i],'.dat'), h=F)
    dates[,1] <- substr(dates[,1], 3,4)  # to fit the formula for as.Date
    dates <- apply(dates[, c(2,3,1)], 1, paste, collapse='/')
    length(unique(dates))
    dates <- sapply(dates, as.Date, format="%m/%d/%y")
    length(unique(dates))
    # scale to 0-1? Not necessary
    # dates <- (dates-min(dates))/diff(range(dates))
    dates <- dates - min(dates) + 1
    
    # print(dim(dat))
    # matplot(dat[,1:2], type='b', col='black', pch=c(16,1))
    cat(i,'')
    y <- dat[,1]
    
    for(k in 2:ncol(dat)){
      loop <- loop + 1
      
      if( runIt && (!runParallel || loop==jobID) ){
        if(runParallel || runItLocally){
          x <- dat[,k]
          # t.test(x,y, paired=T)$p.value
          mat <- data.frame( diffs=as.numeric(y-x), time=dates, month=as.factor(dates0$V2), year=as.factor(dates0$V1))
          
          # plot(mat$diffs, pch=16)
          # hist(mat$diffs)
          
          m0 <- gls(diffs ~ 1, data=mat, correlation = corCAR1(form = ~ time), method='ML')
          tmp <- as.numeric(summary(m0)$tTable[,c('Value','p-value')])
          foo <- data.frame(mean=tmp[1], lower=confint(m0)[1], upper=confint(m0)[2], AIC=AIC(m0), noDifference=ifelse(prod(confint(m0))>0, NA, 'TRUE'), Pvalue=tmp[2])
          
          # m01 <- lme(diffs ~ 1, random=~1|month, data=mat, correlation = corCAR1(form = ~ time), method='ML')
          # m02 <- lme(diffs ~ 1, random=~1|year, data=mat, correlation = corCAR1(form = ~ time), method='ML')#consider year random effect and CAR(1) within each year
          # anova(m01, m02, test=F)
          
          (Pvalue <- summary(m1 <- lm(diffs~1, data=mat))$coefficients[,"Pr(>|t|)"])  #same as above
          # plot(m1, which = 2): this is equivalent to the following
          # z <- rstandard(m1)  #rstudent  residuals
          # qqnorm(z);  qqline(z, col = 2)
          # evalPred(mat$diffs, fitted(m1))[1]
          
          # m2 <- heavyLm(diffs~1, data=mat, family = Student(df=3.579880 ))
          # # m2 <- heavyLm(diffs~1, data=mat, family = Cauchy())
          # a <- summary(m2)
          # a$coefficients[,"p-value"]
          # evalPred(mat$diffs, fitted(m2))[1]  # does not see improve over lm results 
          # z <- residuals(m2)  #rstudent  residuals
          
          filename <- ifelse(k==2, fnam, NA)
          difference <- paste0('row 1 - ',k)
          tmp <- cbind(filename, difference, fitD(mat$diffs, 'normal')[-1], 
                       Pvalue, 
                       fitD(mat$diffs, 't'), 
                       foo)  #note fitD with normal should produce same results as m1 the lm fit, check AIC(m1)
        }
        
        if(runParallel) save(file=paste0('out',jobID), tmp)  else  load(paste0('out',loop))
        out <- rbind(out, tmp)
      }
    }
    
  }
  cat('\n')
  
  if(!runParallel){
    inds <- which(names(out)=='AIC')
    out$BetterT <- ifelse(out[,inds[2]] < out[,inds[1]], 'TRUE', NA)
    out$BestGLS <- ifelse(out[,inds[3]]<apply(out[,inds[1:2]],1,min), "Yes", NA)
    write.csv(file=paste0('output_Option',nopt,'.csv'), out, row.names=F, na='')
  }
}



# # check dates
# a <- t(read.csv('./Data/Option 1/SIM_output_Opt1_H_L.dat', h=F))
# d1 <- read.csv('./Data/Option 1/Time_output_Opt1_H.dat', h=F)
# https://stackoverflow.com/questions/12623027/how-to-analyse-irregular-time-series-in-r
# mod <- gamm(response ~ s(dayOfYear, bs = "cc") + s(timeOfSampling), data = foo,
#             correlation = corCAR1(form = ~ timeOfSampling))

# not run


