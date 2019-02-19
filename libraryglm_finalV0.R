library(glmpath)
library(mice)

########################################################################
## lasso.glm    2017/11/29
##
## LASSO for generlized liner model, 
##       minimim AIC is used for model selection
## LASSO will be used to variable selection purpose and the selection
## variables will be fitted using ordinary glm method.
##
## datain is the imputed data set
## yname is the name of the outcome
## xname is the name of the predictors
########################################################################

lasso.glm <- function(datain, yname, xname)
{ 
library(glmpath)

# x1=age, x2=lwt, x3=black, x4=other
datain <- transform(datain, x3=as.numeric(x3)-1,   #x3 and x4 are factors in original data
                            x4=as.numeric(x4)-1)
  
flm1 <- as.formula(paste(yname, "~", paste(xname, collapse= "+")))
A <- model.matrix(flm1, data=datain)
data1 <- cbind(y=datain[, yname], A)                     # this data does not carry ptid

fit1 <- glmpath(as.matrix(data1[, -c(1,2)]), data1[, yname], family='binomial')

b.pred.aic <- cbind(fit1$b.predictor, aic = fit1$aic)      # estimated coefficients in each step
stepo <- which(b.pred.aic[, ncol(b.pred.aic)] == min(fit1$aic))  

t <-  length(xname)+1
b.pred.aicm <- b.pred.aic[stepo,][-(length(xname)+2)]
x.index <- names(b.pred.aicm[b.pred.aicm !=0])[-1]
if(length(x.index)!=0) flm2 <- as.formula(paste(yname, "~", paste(x.index, collapse= "+")))
if(length(x.index)==0) flm2 <- as.formula(paste(yname, "~", paste(1, collapse= "+")))  # fit an intercept model

coef <- summary(fit2 <- glm(flm2, data=datain, family='binomial'))$coefficients
varselect <- rownames(coef)
b.pred.glm <- data.frame(varselect, Estimate=coef[, 1], Std.error=coef[,2])
b.pred.glm <- merge(datanull, b.pred.glm, by='varselect', all.x=T)

b.pred.glm[is.na(b.pred.glm)] <- 0
lb <- b.pred.glm[, 2]-1.96*b.pred.glm[, 3]
ub <- b.pred.glm[, 2]+1.96*b.pred.glm[, 3]
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

b.pred.glm <- cbind(b.pred.glm, lb, ub, covp, ciw, b.ints)

datain <- transform(datain, ypred=predict(fit2))
datain <- transform(datain, p=exp(ypred)/(1+exp(ypred)))

# The following results match Hjort and Claeskens (2003) Table 1 (mode 345), SE used delta method
pwhitefit <- predict(fit2, newdata=data.frame(x1=agemw, x2=lwtmw, x3=0, x4=0), type='response', se.fit=T)
pblackfit <- predict(fit2, newdata=data.frame(x1=agemb, x2=lwtmb, x3=1, x4=0), type='response', se.fit=T)
pwhite <- pwhitefit[[1]]
pblack <- pblackfit[[1]]
sewhite <- pwhitefit[[2]]
seblack <- pblackfit[[2]]

pwb.glm <- cbind(c(pwhite, pblack), c(sewhite, seblack))

lbpwb <- pwb.glm[, 1]-1.96*pwb.glm[, 2]
ubpwb <- pwb.glm[, 1]+1.96*pwb.glm[, 2]
ciw.pwb <- ubpwb-lbpwb

covp.pwb <- rep(0, length(pwb))

pwb.ints.lb <- NULL
pwb.ints.ub <- NULL
for (index in 1:length(pwb))
{  if (pwb[index]<=ubpwb[index]&pwb[index]>=lbpwb[index]) covp.pwb[index] <- 1
   pwb.ints.lb[index] <- ifelse (pwb[index]<lbpwb[index], lbpwb[index]-pwb[index], 0)
   pwb.ints.ub[index] <- ifelse (pwb[index]>ubpwb[index], pwb[index]-ubpwb[index],0)
}

pwb.ints <- ciw.pwb + 2/0.05*(pwb.ints.lb+pwb.ints.ub)

pwb.glm <- cbind(pwb.glm, lbpwb, ubpwb, covp.pwb, ciw.pwb, pwb.ints)
b.pred.glm <- b.pred.glm[order(b.pred.glm[,1]),]

return(list(b.pred.glm, pwb.glm))
}



############################################################################
##  MIS.glm():    2017/11/29
##
##  Multiple Imputations + model Selection
##  Impute data 5 times, for each imputed data set, perform model selection
##  Both beta coefficient and predicted probabity will be estimated
##  Called Rubin's Rules approach in paper
############################################################################

MIS.glm <- function(dataimp)
{
	
beta.1all <- datanull
sd.1all <- datanull
pwb.1all <- NULL
sdpwb.1all <- NULL

for (i in 1:M)
{
	impi <- complete(dataimp, i)
	outtmp <- lasso.glm(impi, yname, xname)  
	beta.1all <- merge(beta.1all, outtmp[[1]][,1:2], by='varselect') # this is revised from earlier version 
	sd.1all <- merge(sd.1all, outtmp[[1]][,c(1,3)],by='varselect')
    pwb.1all <- cbind(pwb.1all, outtmp[[2]][,1])
    sdpwb.1all <- cbind(sdpwb.1all, outtmp[[2]][,2])
}

varselect <- beta.1all[, 1]
beta.1all <- beta.1all[, -1]          
beta.1 <- apply(beta.1all, 1, mean)   # point estimate from MI
pwb.1 <- apply(pwb.1all, 1, mean)

# Overall variance is calculated using Rubin's Rule
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
     
# Estimate predicted probability             
meanvarpwb <- apply(sdpwb.1all^2, 1, mean)   # mean of variances for each imputation
varpwb <- apply(pwb.1all, 1, sd)^2   # variance of estimate

sdpwb.1 <- sqrt(meanvarpwb+ (1+1/M)*varpwb)

lbpwb.1 <- pwb.1-1.96*sdpwb.1
ubpwb.1 <- pwb.1+1.96*sdpwb.1 
ciwpwb.1 <- ubpwb.1-lbpwb.1

covppwb.1 <- rep(0, length(pwb))

pwb.ints.lb1 <- NULL
pwb.ints.ub1 <- NULL
for (index in 1:length(pwb))
{  
if (pwb[index]<=ubpwb.1[index]&pwb[index]>=lbpwb.1[index]) covppwb.1[index] <- 1
pwb.ints.lb1[index]<-ifelse(pwb[index]<lbpwb.1[index],lbpwb.1[index]-pwb[index], 0)
pwb.ints.ub1[index]<-ifelse(pwb[index]>ubpwb.1[index], pwb[index]-ubpwb.1[index],0)
}

pwb.ints.1 <- ciwpwb.1 + 2/0.05*(pwb.ints.lb1+pwb.ints.ub1)

outmis <- data.frame(varselect, beta.1, sd.1, lb.1, ub.1, covp.1, ciw.1, beta.ints.1)
outmis <- outmis[order(outmis[,1]), ]
return(list(outmis, outp=cbind(pwb.1, sdpwb.1, lbpwb.1, ubpwb.1, covppwb.1, ciwpwb.1, pwb.ints.1)))
}



###########################################################################################################
##  SBISv2.glm():  2017/11/29
##  
##  (Smoothed) Bootstrap + single Imputation + model Selection for linear model
## 
##  This funciton will calculate the estimate and variance using both 
##  bootstrap SI MS (BIS, called Impute-Select approach in paper)
##  and smoothed bootstrap SI MS (SBIS, called Efron's Rules approach in paper)
##  
##  input:  
##  dataimp1 -  singly imputed data
##  datasimu -  simulated data set with missing data
##  
##  Default parameters (defined in main program): B, yname, xname, n, datanull
##
##  output:
##  beta.3    -  estimated coefficient from a single imputation (BIS)
##  sd.3      -  estimated variance for beta.3
##  beta.4    -  smoothed estimated coefficient (average of estimates from B bootstrap samples, SBIS)
##  sd.4      -  smoothed variance estimate for beta.4
##  ystar.all -  n rows by B+1 columns, includes all bootstrap samples
##               1st column: ptid, 
##               2nd to last columns: the count that each record is selected in each bootstrap sample
##  beta.3all -  all estimated coefficient using bootstrap SI MS (B rows and q columns, 
##               q is the number of predictors) 
###########################################################################################################

SBISv2.glm <- function(dataimp1, datasimu)
{ 
	
ql <- function(x){quantile(x, 0.025)}
qu <- function(x){quantile(x, 0.975)}
	
# Estimate from original data (point estimate of method BIS), NA is also included
templasso <- lasso.glm(dataimp1, yname, xname)
varselect.3 <- 	templasso[[1]][,1]	
beta.3 <- templasso[[1]][,2]
pwb.3 <- templasso[[2]][,1]

idall <- datasimu$ptid
# To save estimated coefficient for all predictors (will be 0 if not selected)
beta.3all <- datanull
ystar.all <- data.frame(ptid=idall)
pwb.3all  <- NULL

for (j in 1:B)
{
#	set.seed(20160731+j)
	id <- sample(1:length(idall), n, replace=TRUE, prob=NULL)
	datasimua <- datasimu[id,]   
    imp <- mice(datasimua[, c(xname, yname)], m=1, method=mein)  
    imp <- complete(imp, 1)
    outtmp<- lasso.glm(imp, yname, xname)
    beta.3all <- merge(beta.3all, outtmp[[1]][,1:2], by='varselect')  # save estimated coefficient from each bootstrap samples
    
    ystar <- data.frame(table(datasimua$ptid))
    colnames(ystar) <- c('ptid', paste('count', j, sep='.'))
    ystar.all <- merge(ystar.all, ystar, by='ptid', all.x=T)
    
    pwb.3all <- cbind(pwb.3all, outtmp[[2]][,1])
    
}

# BIS variance estimate
sd.3 <- apply(beta.3all[,-1], 1, sd)     
sdpwb.3 <- apply(pwb.3all, 1, sd) 

# SBIS variance estimate
varselect.4 <- beta.3all[,1]
beta.4 <- apply(beta.3all[,-1], 1, mean)  # t*, bootstrap estimates (an average of bootstrapped estimates, smoothed point estimate)
pwb.4 <- apply(pwb.3all, 1, mean)

ystar.all[is.na(ystar.all)] <- 0
ystar.all2 <- ystar.all[,-1]         # Y*_{ij} in Efron
zstar.all2 <- ystar.all2-1           # n by B matrix

ym <- apply(ystar.all2, 1, mean)     # Y*_{.j} in Efron, 1 by n vector
ystar.allm <- ystar.all2-ym          # n by B dimension, Y*_{ij} - Y*_{.j}

beta.3allm <- t(beta.3all[,-1]-beta.4)    # t*_i - t*, B by p matrix, p=the number of predictors+1 (intercept)
pwb.3allm <- t(pwb.3all-pwb.4)

cov <- as.matrix(ystar.allm)%*%beta.3allm/B  # product of two matrices, get n by p matrix
covpwb <- as.matrix(ystar.allm)%*%pwb.3allm/B

cov2 <- cov^2
covpwb2 <- covpwb^2

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

crctpwb.all <- NULL
for (p in 1:length(pwb.4))
{
  zstarpwb <- as.matrix(t(zstar.all2))*pwb.3allm[,p]  # B by n
  # not a product of matrices, ith row of the first matrix will be multiplied by the ith element in the second vector. 
  crctpwb <- t(t(zstarpwb) - covpwb[,p])^2   # B by n
  crctpwb.all <- c(crctpwb.all, sum(crctpwb)/B^2) # a vection with length of p
}
sdpwb.4 <- sqrt(apply(covpwb2, 2, sum)-crctpwb.all)

lb.3 <- apply(beta.3all[,-1], 1, ql)
ub.3 <- apply(beta.3all[,-1], 1, qu)
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

lbpwb.3 <- apply(pwb.3all, 1, ql)
ubpwb.3 <- apply(pwb.3all, 1, qu)
ciwpwb.3 <- ubpwb.3-lbpwb.3

covppwb.3 <- rep(0, length(pwb))

pwb.ints.lb3 <-NULL
pwb.ints.ub3 <-NULL
for (index in 1:length(pwb))
{  
if (pwb[index]<=ubpwb.3[index]&pwb[index]>=lbpwb.3[index]) covppwb.3[index] <- 1
pwb.ints.lb3[index]<-ifelse(pwb[index]<lbpwb.3[index],lbpwb.3[index]-pwb[index], 0)
pwb.ints.ub3[index] <-ifelse(pwb[index]>ubpwb.3[index],pwb[index]-ubpwb.3[index],0)
}

pwb.ints.3 <- ciwpwb.3 + 2/0.05*(pwb.ints.lb3+pwb.ints.ub3)

lbpwb.4 <- pwb.4-1.96*sdpwb.4
ubpwb.4 <- pwb.4+1.96*sdpwb.4
ciwpwb.4 <- ubpwb.4-lbpwb.4

covppwb.4 <- rep(0, length(pwb))

pwb.ints.lb4 <-NULL
pwb.ints.ub4 <-NULL
for (index in 1:length(pwb))
{  
if (pwb[index]<=ubpwb.4[index]&pwb[index]>=lbpwb.4[index]) covppwb.4[index] <- 1
pwb.ints.lb4[index]<-ifelse(pwb[index]<lbpwb.4[index],lbpwb.4[index]-pwb[index], 0)
pwb.ints.ub4[index] <-ifelse(pwb[index]>ubpwb.4[index],pwb[index]-ubpwb.4[index],0)
}

pwb.ints.4 <- ciwpwb.4 + 2/0.05*(pwb.ints.lb4+pwb.ints.ub4)

estimate.3 <- data.frame(varselect=varselect.3, beta.3, sd.3, lb.3, ub.3, covp.3, ciw.3, beta.ints.3)
estimate.4 <- data.frame(varselect=varselect.4, beta.4, sd.4, lb.4, ub.4, covp.4, ciw.4, beta.ints.4)
outp<- cbind(pwb.3, sdpwb.3, pwb.4, sdpwb.4)
cip <- cbind(lbpwb.3, ubpwb.3, lbpwb.4, ubpwb.4, ciwpwb.3, ciwpwb.4, pwb.ints.3, pwb.ints.4)

estimate.3 <- estimate.3[order(estimate.3[,1]), ]
estimate.4 <- estimate.4[order(estimate.4[,1]), ]

return(list(estimate.3, estimate.4,  
            outp=cbind(pwb.3, sdpwb.3, pwb.4, sdpwb.4), 
            covppwb=cbind(covppwb.3, covppwb.4), 
            ciwpwb=cbind(ciwpwb.3, ciwpwb.4),  
            cip, 
            pwb.ints=cbind(pwb.ints.3, pwb.ints.4)))
}



####################################################
##  simu.glm - 2016/08/08
##  Data simulation
##
####################################################

simu.glm <- function(datain, seedin)
{
set.seed(seedin)
datain <- transform(datain, ptid=1:n)
idall <- datain$ptid
id <- sample(1:length(idall), n, replace=TRUE, prob=NULL)
dataina <- datain[id, c('age', 'lwt', 'race')]
yp <- predict(fitwt, dataina)
dataina <- data.frame(dataina, yp)
dataina <- transform(dataina, x1=age,
                              x2=lwt,
                              x3=ifelse(race==2, 1, 0),   #black
                              x4=ifelse(race==3, 1, 0),   #other
                              p=exp(yp)/(1+exp(yp)))
datainb <- dataina[, c('yp', 'x1', 'x2', 'x3', 'x4', 'p')]
y <- rbinom(length(datainb$p), 1, datainb$p)

datainb <- data.frame(datainb, y=y) 	
datainc <- datainb[, c('y', 'x1', 'x2', 'x3', 'x4')]

x1 <- datainc$x1
x2 <- datainc$x2
x3 <- datainc$x3
x4 <- datainc$x4
                      
# Missing probability
px2 <- 1/(y+0.007*x1^2)
px3 <- 1-1/(1+0.008*x1)
px4 <- 1-1/(1+0.005*x1)

# Generate missing index
a2 <- runif(n, 0, 1)
a3 <- runif(n, 0, 1)
a4 <- runif(n, 0, 1)

datamis <- data.frame(ptid=1:n, y, x1, x2, x3, x4, a2, a3, a4)
datamis <- transform(datamis, x2.mis=ifelse(a2<px2, NA, x2),
						                  x3.mis=as.factor(ifelse(a3<px3, NA, x3)),
						                  x4.mis=as.factor(ifelse(a4<px4, NA, x4)))
datacomp <- datamis[, c('ptid', 'y', 'x1', 'x2', 'x3', 'x4')]	
datacomp <- transform(datacomp, 
                      x3=as.factor(x3),
                      x4=as.factor(x4))
datasimu <- datamis[, c('ptid','y', 'x1', 'x2.mis', 'x3.mis', 'x4.mis')]						  
colnames(datasimu) <- c( 'ptid','y', 'x1', 'x2', 'x3', 'x4')	

x2.misrate <- prop.table(table(is.na(datasimu$x2)))[2]
x3.misrate <- prop.table(table(is.na(datasimu$x3)))[2]
x4.misrate <- prop.table(table(is.na(datasimu$x4)))[2]

return(list(datacomp=datacomp, datasimu=datasimu, misrate=c(x2.misrate, x3.misrate, x4.misrate)))
}
