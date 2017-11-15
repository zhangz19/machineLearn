

# functions to test equivalence of paired sequences (can be irregularly spaced time series)
# The data was extracted from the time series but specific dates were selected that are not necessarily continuous.

rm(list=ls())

# require(heavy)  # for heavyLM
require(MASS)  #for fitdist
# require(actuar)  #for Pareto when data>0
# require(fitdistrplus)  #fitdist

evalPred <- function(y,x){ #y=observed, x=fitted
  ym <- mean(y)
  xm <- mean(x)
  RMSE <- sqrt(mean((y-x)^2))
  NSE <- 1- sum((y-x)^2)/sum((y-ym)^2)
  PBIAS <- 100 * sum(y-x)/sum(y)
  z <- c(RMSE,NSE,PBIAS); names(z) <- c('RMSE','NSE','PBIAS')
  return(z)
}

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


for(nopt in 1:2){
  dir0 <- paste0('./Data/Option ',nopt,'/')
  dirs <- dir(path=dir0, pattern='*.dat')
  (nf <- length(dirs))
  
  out <- numeric()
  for(i in 1:nf){
    dat <- t(read.csv(fnam <- paste0(dir0,dirs[i]), h=F))
    # print(dim(dat))
    # matplot(dat[,1:2], type='b', col='black', pch=c(16,1))
    cat(i,'')
    y <- dat[,1]
    
    for(k in 2:ncol(dat)){
      x <- dat[,k]
      # t.test(x,y, paired=T)$p.value
      mat <- data.frame( diffs=as.numeric(y-x) )
      
      # hist(mat$diffs)
      
      (pval <- summary(m1 <- lm(diffs~1, data=mat))$coefficients[,"Pr(>|t|)"])  #same as above
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
      tmp <- cbind(filename, difference, fitD(mat$diffs, 'normal')[-1], pval, fitD(mat$diffs, 't'))  #note fitD with normal should produce same results as m1 the lm fit, check AIC(m1)
      
      out <- rbind(out, tmp)
    }
    
  }
  cat('\n')
  
  inds <- which(names(out)=='AIC')
  out$BetterT <- (out[,inds[2]] < out[,inds[1]])
  write.csv(file=paste0('output_Option',nopt,'.csv'), out, row.names=F, na='')
}

# not run


