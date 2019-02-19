library(glmpath)
library(mice)

########################################################################
## lasso.lm    2017/11/30
##
## LASSO for linear model, minimim AIC is used for model selection
## LASSO will be used to variable selection purpose and the selection
## Variables will be fitted using ordinary least square method.
##
## datain is the imputed data set
## yname is the name of the outcome
## xname is the name of the predictors
########################################################################

lasso.lm <- function(datain, yname, xname)
{ 
flm1 <- as.formula(paste(yname, "~", paste(xname, collapse= "+")))
A <- model.matrix(flm1, data=datain)
data1 <- cbind(y=datain[, yname], A)    # this data does not carry ptid

fit1 <- glmpath(as.matrix(data1[, -c(1,2)]), data1[, yname], family='gaussian')

b.pred.aic <- cbind(fit1$b.predictor, aic = fit1$aic)     # estimated coefficients in each step
stepo <- which(b.pred.aic[, ncol(b.pred.aic)] == min(fit1$aic))
t <-  length(xname)+1
b.pred.aicm <- b.pred.aic[stepo,][-c(1, length(xname)+2)]
x.index <- names(b.pred.aicm[b.pred.aicm !=0])
flm2 <- as.formula(paste(yname, "~", paste(x.index, collapse= "+")))
b.pred.lm <- summary(lm(flm2, data=datain))$coefficients[,1:2]
b.pred.lm <- data.frame(varselect=rownames(b.pred.lm), b.pred.lm)
b.pred.lm <- merge(datanull, b.pred.lm, by='varselect', all.x=T)

b.pred.lm[is.na(b.pred.lm)] <- 0
lb <- b.pred.lm[, 2]-1.96*b.pred.lm[, 3]
ub <- b.pred.lm[, 2]+1.96*b.pred.lm[, 3]
ciw <- ub-lb
covp <- rep(0, length(beta))

b.ints.lb <- NULL
b.ints.ub <- NULL
for (index in 1:length(beta))
{ if (beta[index]<=ub[index]&beta[index]>=lb[index]) covp[index] <- 1
  b.ints.lb[index] <- ifelse(beta[index]<lb[index], lb[index]-beta[index], 0)
  b.ints.ub[index] <- ifelse(beta[index]>ub[index], beta[index]-ub[index], 0)
}

b.ints <- ciw + 2/0.05*(b.ints.lb+b.ints.ub)

b.pred.lm <- cbind(b.pred.lm, lb, ub, covp, ciw, b.ints)
return(b.pred.lm)
}

########################################################################
## lasso.lm2    2018/10/30
##
## Linear model without model selection
##
## datain is the imputed data set
## yname is the name of the outcome
## xname is the name of the predictors
########################################################################

lasso.lm2 <- function(datain, yname, xname)
{ 
  flm2 <- as.formula(paste(yname, "~", paste(xname, collapse= "+")))
  b.pred.lm <- summary(lm(flm2, data=datain))$coefficients[,1:2]
  b.pred.lm <- data.frame(varselect=rownames(b.pred.lm), b.pred.lm)
  b.pred.lm <- merge(datanull, b.pred.lm, by='varselect', all.x=T)
  
  b.pred.lm[is.na(b.pred.lm)] <- 0
  lb <- b.pred.lm[, 2]-1.96*b.pred.lm[, 3]
  ub <- b.pred.lm[, 2]+1.96*b.pred.lm[, 3]
  ciw <- ub-lb
  covp <- rep(0, length(beta))
  
  b.ints.lb <- NULL
  b.ints.ub <- NULL
  for (index in 1:length(beta))
  { if (beta[index]<=ub[index]&beta[index]>=lb[index]) covp[index] <- 1
  b.ints.lb[index] <- ifelse(beta[index]<lb[index], lb[index]-beta[index], 0)
  b.ints.ub[index] <- ifelse(beta[index]>ub[index], beta[index]-ub[index], 0)
  }
  
  b.ints <- ciw + 2/0.05*(b.ints.lb+b.ints.ub)
   
  b.pred.lm <- cbind(b.pred.lm, lb, ub, covp, ciw, b.ints)
  
  return(b.pred.lm)
}



#############################################################################
##  MIS.lm():    2017/11/30
##
##  Multiple Imputations + model Selection
##  Impute data 5 times, for each imputed data set, perform model selection
##  Called Rubin's Rules approach in paper
############################################################################

MIS.lm <- function(dataimp)
{
beta.1all <- datanull
sd.1all <- datanull

for (i in 1:M)
{
	impi <- complete(dataimp, i)
	outtmp <- lasso.lm(impi, yname, xname)  
    beta.1all <- cbind(beta.1all, outtmp[,2])
    sd.1all <- cbind(sd.1all, outtmp[,3])
}

beta.1all <- beta.1all[, -1]
beta.1 <- apply(beta.1all, 1, mean)   # point estimate from MI

# Overall variance is calculated using Rubin's Rules
sd.1all <- sd.1all[, -1]
meanvar <- apply(sd.1all^2, 1, mean)   # mean of variances for each imputation
varbeta <- apply(beta.1all, 1, sd)^2   # variance of estimate

sd.1 <- sqrt(meanvar+ (1+1/M)*varbeta)

lb.1 <- beta.1-1.96*sd.1
ub.1 <- beta.1+1.96*sd.1 
ciw.1 <- ub.1-lb.1
covp.1 <- rep(0, length(beta))

beta.ints.lb1 <- NULL
beta.ints.ub1 <- NULL
for (index in 1:length(beta))
{if (beta[index]<=ub.1[index]&beta[index]>=lb.1[index]) covp.1[index] <- 1
  beta.ints.lb1[index]<-ifelse(beta[index]<lb.1[index], lb.1[index]-beta[index], 0)
  beta.ints.ub1[index]<-ifelse(beta[index]>ub.1[index], beta[index]-ub.1[index], 0)
}
beta.ints.1 <- ciw.1 + 2/0.05*(beta.ints.lb1+beta.ints.ub1)
                       
return(cbind(beta.1, sd.1, lb.1, ub.1, covp.1, ciw.1, beta.ints.1))
}



######################################################################
##  MI.lm():    2018/10/30
##
##  Multiple Imputations without model Selection
##  Impute data 5 times
######################################################################

MI.lm <- function(dataimp)
{ 
  beta.1all <- datanull
  sd.1all <- datanull
  
  for (i in 1:M)
  {
    impi <- complete(dataimp, i)
    outtmp <- lasso.lm2(impi, yname, xname)  
    beta.1all <- cbind(beta.1all, outtmp[,2])
    sd.1all <- cbind(sd.1all, outtmp[,3])
  }
  
  beta.1all <- beta.1all[, -1]
  beta.1 <- apply(beta.1all, 1, mean)   # point estimate from MI
   
  # Overall variance is calculated using Rubin's Rules
  sd.1all <- sd.1all[, -1]
  meanvar <- apply(sd.1all^2, 1, mean)   # mean of variances for each imputation
  varbeta <- apply(beta.1all, 1, sd)^2   # variance of estimate
  
  sd.1 <- sqrt(meanvar+ (1+1/M)*varbeta)
  
  lb.1 <- beta.1-1.96*sd.1
  ub.1 <- beta.1+1.96*sd.1 
  ciw.1 <- ub.1-lb.1
  covp.1 <- rep(0, length(beta))
  
  beta.ints.lb1 <- NULL
  beta.ints.ub1 <- NULL
  for (index in 1:length(beta))
  {if (beta[index]<=ub.1[index]&beta[index]>=lb.1[index]) covp.1[index] <- 1
  beta.ints.lb1[index]<-ifelse(beta[index]<lb.1[index], lb.1[index]-beta[index], 0)
  beta.ints.ub1[index]<-ifelse(beta[index]>ub.1[index], beta[index]-ub.1[index], 0)
  }   
  beta.ints.1 <- ciw.1 + 2/0.05*(beta.ints.lb1+beta.ints.ub1)
   
  return(cbind(beta.1, sd.1, lb.1, ub.1, covp.1, ciw.1, beta.ints.1))
}



#######################################################################################################
##  SBISv2.lm():  2017/11/30
##  
##  (Smoothed) Bootstrap + single Imputation (SI) + model Selection (MS) for linear model
## 
##  This funciton will calculate the estimate and variance using both 
##  bootstrap SI MS (BIS, called Impute-Select approach in paper) 
##  and smoothed bootstrap SI MS (SBIS, called Efron's Rules approach in paper)
##  
##  input:  
##  dataimp1 -  singly imputed data
##  datasimu -  simulated data set with missing data
##  
##  Default paramters from main program: B, yname, xname, n, datanull
##
##  output:
##  beta.3    -  estimated coefficient from a single imputation (BIS - Impute-Select approach)
##  sd.3      -  estimated variance for beta.3
##  beta.4    -  smoothed estimated coefficient 
##               (average of estimates from B boostrap samples, SBIS - Efron's Rules apprach)
##  sd.4      -  smoothed variance estimate for beta.4
##  ystar.all -  n rows by B+1 columns, includes all boostrap samples
##               1st column: ptid, 
##               2nd to last columns: the count that each record is selected in each bootstrap sample
##  beta.3all -  all estimated coefficient using boot SI MS 
##               (B rows and q columns, q is the number of predictors) 
########################################################################################################

SBISv2.lm <- function(dataimp1, datasimu)
{ 
ql <- function(x){quantile(x, 0.025)}
qu <- function(x){quantile(x, 0.975)}
	
# Estimate from original data (point estimate of method BIS), NA is also included
beta.3 <- lasso.lm(dataimp1, yname, xname)	

idall <- datasimu$ptid
# To save estimated coefficient for all predictors (will be 0 if not selected)
beta.3all <- NULL 
ystar.all <- data.frame(ptid=idall)

for (j in 1:B)
{
	id <- sample(1:length(idall), n, replace=TRUE, prob=NULL)
	datasimua <- datasimu[id,]   
    imp <- mice(datasimua[, c(xname, yname)], m=1, method=mein)
    imp <- complete(imp, 1)
    outtmp<- lasso.lm(imp, yname, xname)
    beta.3all <- cbind(beta.3all, outtmp[,2])  #save estimated coefficient from each bootstrap sample
    
    ystar <- data.frame(table(datasimua$ptid))
    colnames(ystar) <- c('ptid', paste('count', j, sep='.'))
    ystar.all <- merge(ystar.all, ystar, by='ptid', all.x=T)    
}

# BIS variance estimate
sd.3 <- apply(beta.3all, 1, sd)     

# SBIS variance estimate
beta.4 <- apply(beta.3all, 1, mean)  # t*, bootstrap estimates (an average of bootstrap estimates, smoothed point estimate)

ystar.all[is.na(ystar.all)] <- 0
ystar.all2 <- ystar.all[,-1]         # Y*_{ij} in Efron
zstar.all2 <- ystar.all2-1           # n by B matrix

ym <- apply(ystar.all2, 1, mean)     # Y*_{.j} in Efron, 1 by n vector
ystar.allm <- ystar.all2-ym          # n by B dimension, Y*_{ij} - Y*_{.j}

beta.3allm <- t(beta.3all-beta.4)    # t*_i - t*, B by p matrix, p=the number of predictors+1 (intercept)

cov <- as.matrix(ystar.allm)%*%beta.3allm/B  # product of two matrices, get n by p matrix
cov2 <- cov^2

# Variance bias correction
crct.all <- NULL
for (p in 1:length(beta.4))
{
		zstar <- as.matrix(t(zstar.all2))*beta.3allm[,p]  # B by n
		# not a product of matrices, ith row of the first matrix will be multiplied by the ith element in the second vector. 
         crct <- t(t(zstar) - cov[,p])^2   # B by n
crct.all <- c(crct.all, sum(crct)/B^2) # a vection with length of p
}

sd.4 <- sqrt(apply(cov2, 2, sum)-crct.all)

lb.3 <- apply(beta.3all, 1, ql)
ub.3 <- apply(beta.3all, 1, qu)
ciw.3 <- ub.3-lb.3

covp.3 <- rep(0, length(beta))
beta.ints.lb3 <- NULL
beta.ints.ub3 <- NULL

for (index in 1:length(beta))
{ 
if (beta[index]<=ub.3[index]&beta[index]>=lb.3[index]) covp.3[index] <- 1
beta.ints.lb3[index]<-ifelse(beta[index]<lb.3[index], lb.3[index]-beta[index],0)
beta.ints.ub3[index]<-ifelse(beta[index]>ub.3[index], beta[index]-ub.3[index], 0)
}

beta.ints.3 <- ciw.3 + 2/0.05*(beta.ints.lb3+beta.ints.ub3)

lb.4 <- beta.4-1.96*sd.4
ub.4 <- beta.4+1.96*sd.4
ciw.4 <- ub.4-lb.4

covp.4 <- rep(0, length(beta))

beta.ints.lb4 <- NULL
beta.ints.ub4 <- NULL
for (index in 1:length(beta))
{
if (beta[index]<=ub.4[index]&beta[index]>=lb.4[index]) covp.4[index] <- 1
beta.ints.lb4[index]<-ifelse(beta[index]<lb.4[index], lb.4[index]-beta[index],0)
beta.ints.ub4[index]<-ifelse(beta[index]>ub.4[index], beta[index]-ub.4[index], 0)
}
beta.ints.4 <- ciw.4 + 2/0.05*(beta.ints.lb4+beta.ints.ub4)

return(list(estimate=cbind(beta.3[,2], sd.3, beta.4, sd.4, beta.ints.3, beta.ints.4), ystar.all=ystar.all, beta.3all=beta.3all, ci=cbind(lb.3, ub.3, lb.4, ub.4), covp=cbind(covp.3, covp.4), ciw=cbind(ciw.3, ciw.4)))
}


##################################################
##  simucopula.R - 2016/08/01
##  Data simulation
##
##  Use data in Schomaker and Heumann (2014), 
##  Use copula to add dependence of predictors
##################################################

simucopula <- function(seedin)
{
library(copula)
set.seed(seedin)	
mycop.clayton <- archmCopula(family='clayton', dim=6, param=1)
mymvd <- mvdc(mycop.clayton, margins=c('norm', 'lnorm', 'weibull', 'exp', 'gamma', 'norm'), paramMargins=list(list(mean=0.5, sd=1), list(meanlog=0.5, sdlog=0.5), list(shape=1.75, scale=1.9), 
                  list(rate=1), list(shape=0.25, scale=2), list(mean=0.25, sd=1)))
x <- rMvdc(n, mymvd)

x <- data.frame(1, x)
colnames(x) <- c('1', 'x1', 'x2', 'x3', 'x4', 'x5', 'x6')

x1 <- x[, 'x1']
x2 <- x[, 'x2']
x3 <- x[, 'x3']
x4 <- x[, 'x4']
x5 <- x[, 'x5']
x6 <- x[, 'x6']

beta <- c(2.5, -3, -0.25, 0, -1.5, 0, 0.35)
mu <- as.vector(beta%*%t(x))
y<- NULL
for (i in 1:length(mu)) y[i] <- rnorm(1, mu[i], exp(1.25))  #exp(1.25)=3.49, check this step!!
ptid <- 1:n
dt <- data.frame(ptid, mu, y, x[, 2:7])

# Missing probability
px1 <- 1-((0.15*y)^2+1)^-1
px4 <- 1-(1+0.02*x2^3)^-1
px5 <- 1-(1+exp(1-2*x3))^-1

# Generate missing index
a1 <- runif(n, 0, 1)
a4 <- runif(n, 0, 1)
a5 <- runif(n, 0, 1)

datamis <- data.frame(ptid, mu, y, x[, 2:7], a1, a4, a5)
datamis <- transform(datamis, x1.mis=ifelse(a1<px1, NA, x1),
						      x4.mis=ifelse(a4<px4, NA, x4),
						      x5.mis=ifelse(a5<px5, NA, x5))

# Original complete data
datacomp <- datamis[, c('ptid', 'y', 'x1', 'x2', 'x3', 'x4', 'x5', 'x6')]
# simulated missing data			        
datasimu <- datamis[, c('ptid', 'y', 'x1.mis', 'x2', 'x3', 'x4.mis', 'x5.mis', 'x6')]					  
colnames(datasimu) <- c('ptid', 'y', 'x1', 'x2', 'x3', 'x4', 'x5', 'x6')	

x1.misrate <- prop.table(table(is.na(datasimu$x1)))[2]
x4.misrate <- prop.table(table(is.na(datasimu$x4)))[2]
x5.misrate <- prop.table(table(is.na(datasimu$x5)))[2]
datacompcor <- cor(datacomp[, -c(1,2)])

return(list(datacomp=datacomp, datasimu=datasimu, misrate=c(x1.misrate, x4.misrate, x5.misrate), datacompcor))
}

