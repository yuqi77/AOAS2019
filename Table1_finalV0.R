#####################################################################################
#                                                                                   #
#                              AOAS Table 1 R code                                  #
#  Comparison of post imputation-selection estimators in a linear regression model  #
#                                                                                   #
#                                                                                   #
#####################################################################################

source(paste(path, 'librarylm_finalV0.R', sep=''))

M <- 5        # multiple imputation size
R <- 500      # number of replications
B <- 200      # bootstrap size

# Sample size in each replication, 
# reduced to 250 from original 500 in Schomaker and Heumann (2014)
n <- 250 

beta <- c(2.5, -3, -0.25, 0, -1.5, 0, 0.35)
yname <- 'y'
xname <- c('x1', 'x2', 'x3', 'x4', 'x5', 'x6')

# All missing predictors are continuous, using 
# predictive mean matching for imputation 
mein <- c('pmm', '', '', 'pmm', 'pmm', '','')    
datanull <- data.frame(varselect=c('(Intercept)', xname))

# Initialize the output variables

# Estimated beta
beta.mso.all <- NULL
beta.msc.all <- NULL
beta.bis.all <- NULL
beta.mis.all <- NULL
beta.sbis.all <- NULL
beta.mi.all <- NULL
beta.wid.all <- NULL

# Estimated SD for beta
sd.mso.all <- NULL
sd.msc.all <- NULL
sd.bis.all <- NULL
sd.mis.all <- NULL
sd.sbis.all <- NULL
sd.mi.all <- NULL
sd.wid.all <- NULL

# False positive (FP) rate
fp.mso.all <- NULL
fp.msc.all <- NULL
fp.bis.all <- NULL
fp.mis.all <- NULL
fp.sbis.all <- NULL
fp.mi.all <- NULL
fp.wid.all <- NULL

# False negative (FN) rate
fn.mso.all <- NULL
fn.msc.all <- NULL
fn.bis.all <- NULL
fn.mis.all <- NULL
fn.sbis.all <- NULL
fn.mi.all <- NULL
fn.wid.all <- NULL

# Coverage probability
covp.mso.all <- NULL
covp.msc.all <- NULL
covp.bis.all <- NULL
covp.mis.all <- NULL
covp.sbis.all <- NULL
covp.mi.all <- NULL
covp.wid.all <- NULL

# Width of CI
ciw.mso.all <- NULL
ciw.msc.all <- NULL
ciw.bis.all <- NULL
ciw.mis.all <- NULL
ciw.sbis.all <- NULL
ciw.mi.all <- NULL
ciw.wid.all <- NULL

# Interval score
ints.mso.all <- NULL
ints.msc.all <- NULL
ints.bis.all <- NULL
ints.mis.all <- NULL
ints.sbis.all <- NULL
ints.mi.all <- NULL
ints.wid.all <- NULL

# CPU time
time.mso.all <- NULL
time.msc.all <- NULL
time.bis.all <- NULL
time.mis.all <- NULL
time.sbis.all <- NULL
time.mi.all <- NULL
time.wid.all <- NULL

# Missing rate
misrate.all <- NULL

# Correlation coefficient
corlist <- 0

for(r in 1:R)
{
  cat('Replication =', r, '\n')
  
  dataall <- simucopula(20170728+r)
  datacomp <- dataall[[1]]
  datasimu <- dataall[[2]]
  misrate.all <- rbind(misrate.all, dataall[[3]])
  
  corlist <- corlist+dataall[[4]]
  # Multiple imputation
  dataimp <- mice(datasimu[, c(xname, yname)], m=M, method = mein)  
  # Single imputation
  dataimp1 <- mice(datasimu[, c(xname, yname)], m=1, method = mein) 
  dataimp1 <- complete(dataimp1, 1)
  
  # MS (model selection) with complete data
  time.mso.all <- c(time.mso.all, system.time(outtmp.mso <- lasso.lm(datacomp, yname, xname))[1])
  beta.mso.all <- cbind(beta.mso.all, outtmp.mso[,2])
  sd.mso.all <- cbind(sd.mso.all, outtmp.mso[,3])
  # count the number of non-zeros in b3 and b5
  fp.mso.all <- c(fp.mso.all, sum(beta.mso.all[c(4,6),r]!=0))
  # count the number of non-zeros in b1, b2, b4 and b6   
  fn.mso.all <- c(fn.mso.all, sum(beta.mso.all[c(2, 3, 5, 7),r]==0))
  covp.mso.all <- cbind(covp.mso.all, outtmp.mso[,6])
  ciw.mso.all <- cbind(ciw.mso.all, outtmp.mso[,7])
  ints.mso.all <- cbind(ints.mso.all, outtmp.mso[,8])
  
  # MS with complete cases
  idc <- complete.cases(datasimu)
  datasimuc <- datasimu[idc,]
  time.msc.all <- c(time.msc.all, system.time(outtmp.msc <- lasso.lm(datasimuc, yname, xname))[1])
  beta.msc.all <- cbind(beta.msc.all, outtmp.msc[,2])
  sd.msc.all <- cbind(sd.msc.all, outtmp.msc[,3])
  fp.msc.all <- c(fp.msc.all, sum(beta.msc.all[c(4,6),r]!=0))   
  fn.msc.all <- c(fn.msc.all, sum(beta.msc.all[c(2, 3, 5, 7),r]==0))
  covp.msc.all <- cbind(covp.msc.all, outtmp.msc[,6])
  ciw.msc.all <- cbind(ciw.msc.all, outtmp.msc[,7])
  ints.msc.all <- cbind(ints.msc.all, outtmp.msc[,8])
  
  # Wide model (no model selection) with single imputation
  time.wid.all <- c(time.wid.all, system.time(outtmp.wid <- lasso.lm2(dataimp1, yname, xname))[1])
  beta.wid.all <- cbind(beta.wid.all, outtmp.wid[,2])
  sd.wid.all <- cbind(sd.wid.all, outtmp.wid[,3])
  # count the number of non-zeros in b3 and b5
  fp.wid.all <- c(fp.wid.all, sum(beta.wid.all[c(4,6),r]!=0))
  # count the number of non-zeros in b1, b2, b4 and b6   
  fn.wid.all <- c(fn.wid.all, sum(beta.wid.all[c(2, 3, 5, 7),r]==0))
  covp.wid.all <- cbind(covp.wid.all, outtmp.wid[,6])
  ciw.wid.all <- cbind(ciw.wid.all, outtmp.wid[,7])
  ints.wid.all <- cbind(ints.wid.all, outtmp.wid[,8])
  
  # Wide model with multiple imputation
  time.mi.all <- c(time.mi.all, system.time(outtmp.mi <- MI.lm(dataimp))[1])
  beta.mi.all <- cbind(beta.mi.all, outtmp.mi[,1])
  sd.mi.all <- cbind(sd.mi.all, outtmp.mi[,2])
  fp.mi.all <- c(fp.mi.all, sum(outtmp.mi[c(4,6),1]!=0))   
  fn.mi.all <- c(fn.mi.all, sum(outtmp.mi[c(2, 3, 5, 7),1]==0))
  covp.mi.all <- cbind(covp.mi.all, outtmp.mi[,5])
  ciw.mi.all <- cbind(ciw.mi.all, outtmp.mi[, 6])
  ints.mi.all <- cbind(ints.mi.all, outtmp.mi[,7])
  
  # Multiple Imputation and model selection (MIS) - Rubin's Rules approach
  time.mis.all <- c(time.mis.all, system.time(outtmp.mis <- MIS.lm(dataimp))[1])
  beta.mis.all <- cbind(beta.mis.all, outtmp.mis[,1])
  sd.mis.all <- cbind(sd.mis.all, outtmp.mis[,2])
  fp.mis.all <- c(fp.mis.all, sum(outtmp.mis[c(4,6),1]!=0))   
  fn.mis.all <- c(fn.mis.all, sum(outtmp.mis[c(2, 3, 5, 7),1]==0))
  covp.mis.all <- cbind(covp.mis.all, outtmp.mis[,5])
  ciw.mis.all <- cbind(ciw.mis.all, outtmp.mis[, 6])
  ints.mis.all <- cbind(ints.mis.all, outtmp.mis[,7])
  
  # Bootstrap, single imputation then model selection (BIS) - Impute-Select approach
  time.bis.all <- c(time.bis.all, system.time(outtmp <- SBISv2.lm(dataimp1, datasimu))[1])
  outtmp.bis <- outtmp[[1]][,1:2]
  beta.bis.all <- cbind(beta.bis.all, outtmp.bis[,1])
  sd.bis.all <- cbind(sd.bis.all, outtmp.bis[,2])
  fp.bis.all <- c(fp.bis.all, sum(outtmp.bis[c(4,6),1]!=0))   
  fn.bis.all <- c(fn.bis.all, sum(outtmp.bis[c(2, 3, 5, 7),1]==0))
  covp.bis.all <- cbind(covp.bis.all, outtmp[[5]][,1])
  ciw.bis.all <- cbind(ciw.bis.all, outtmp[[6]][,1])
  ints.bis.all <- cbind(ints.bis.all, outtmp[[1]][,5])
  
  # Smoothed Bootstrap, single imputation then model selection (SBIS) - Efron's Rules approach
  outtmp.sbis <- outtmp[[1]][,3:4]
  beta.sbis.all <- cbind(beta.sbis.all, outtmp.sbis[,1])
  sd.sbis.all <- cbind(sd.sbis.all, outtmp.sbis[,2])
  fp.sbis.all <- c(fp.sbis.all, sum(outtmp.sbis[c(4,6),1]!=0))   
  fn.sbis.all <- c(fn.sbis.all, sum(outtmp.sbis[c(2, 3, 5, 7),1]==0)) 
  covp.sbis.all <- cbind(covp.sbis.all, outtmp[[5]][,2])
  ciw.sbis.all <- cbind(ciw.sbis.all, outtmp[[6]][,2])
  ints.sbis.all <- cbind(ints.sbis.all, outtmp[[1]][,6])
}

corlist <- corlist/R
mincor <- min(corlist)
maxcor <- max(corlist[corlist!=1])
corrange <- c(mincor, maxcor)

time.sbis.all <- time.bis.all
colnames(misrate.all) <- c('Missing.rate(b1)', 'Missing.rate(b4)', 'Missing.rate(b5)')

sd.mso.emp<- apply(beta.mso.all, 1, sd)
sd.mso1 <- apply(sd.mso.all, 1, mean)  

sd.msc.emp<- apply(beta.msc.all, 1, sd)
sd.msc1 <- apply(sd.msc.all, 1, mean)

sd.wid.emp<- apply(beta.wid.all, 1, sd)
sd.wid1 <- apply(sd.wid.all, 1, mean)

sd.bis.emp<- apply(beta.bis.all, 1, sd)
sd.bis1 <- apply(sd.bis.all, 1, mean)

sd.mis.emp<- apply(beta.mis.all, 1, sd)
sd.mis1 <- apply(sd.mis.all, 1, mean)

sd.sbis.emp<- apply(beta.sbis.all, 1, sd)
sd.sbis1 <- apply(sd.sbis.all, 1, mean)

sd.mi.emp<- apply(beta.mi.all, 1, sd)
sd.mi1 <- apply(sd.mi.all, 1, mean)

# Missing rate
misrate.avg <- apply(misrate.all, 2, mean)

# Bias (averaged squared loss)
bias.mso <- c(apply((beta.mso.all-beta)^2, 1, mean),  sum((beta.mso.all-beta)^2)/R)
bias.msc <- c(apply((beta.msc.all-beta)^2, 1, mean),  sum((beta.msc.all-beta)^2)/R)
bias.wid <- c(apply((beta.wid.all-beta)^2, 1, mean),  sum((beta.wid.all-beta)^2)/R)
bias.bis <- c(apply((beta.bis.all-beta)^2, 1, mean),  sum((beta.bis.all-beta)^2)/R)
bias.mis <- c(apply((beta.mis.all-beta)^2, 1, mean),  sum((beta.mis.all-beta)^2)/R)
bias.sbis<- c(apply((beta.sbis.all-beta)^2, 1, mean), sum((beta.sbis.all-beta)^2)/R)
bias.mi <- c(apply((beta.mi.all-beta)^2, 1, mean),  sum((beta.mi.all-beta)^2)/R)

# CPU time
timeavg.mso <- c(mean(time.mso.all), sd(time.mso.all))
timeavg.msc <- c(mean(time.msc.all), sd(time.msc.all))
timeavg.wid <- c(mean(time.wid.all), sd(time.wid.all))
timeavg.bis <- c(mean(time.bis.all), sd(time.bis.all))
timeavg.mis <- c(mean(time.mis.all), sd(time.mis.all))
timeavg.sbis <- timeavg.bis
timeavg.mi <- c(mean(time.mi.all), sd(time.mi.all))

# FP and FN, does not count beta0
fpavg.mso <-c(mean(fp.mso.all), sd(fp.mso.all))
fpavg.msc <-c(mean(fp.msc.all), sd(fp.msc.all))
fpavg.wid <-c(mean(fp.wid.all), sd(fp.wid.all))
fpavg.bis <-c(mean(fp.bis.all), sd(fp.bis.all))
fpavg.mis <-c(mean(fp.mis.all), sd(fp.mis.all))
fpavg.sbis <-c(mean(fp.sbis.all), sd(fp.sbis.all))
fpavg.mi <-c(mean(fp.mi.all), sd(fp.mi.all))

fnavg.mso <-c(mean(fn.mso.all), sd(fn.mso.all))
fnavg.msc <-c(mean(fn.msc.all), sd(fn.msc.all))
fnavg.wid <-c(mean(fn.wid.all), sd(fn.wid.all))
fnavg.bis <-c(mean(fn.bis.all), sd(fn.bis.all))
fnavg.mis <-c(mean(fn.mis.all), sd(fn.mis.all))
fnavg.sbis <-c(mean(fn.sbis.all), sd(fn.sbis.all))
fnavg.mi <-c(mean(fn.mi.all), sd(fn.mi.all))

# Coverage probability
covp.mso <- apply(covp.mso.all, 1, mean)
covp.msc <- apply(covp.msc.all, 1, mean)
covp.wid <- apply(covp.wid.all, 1, mean)
covp.bis <- apply(covp.bis.all, 1, mean)
covp.mis <- apply(covp.mis.all, 1, mean)
covp.sbis <- apply(covp.sbis.all, 1, mean)
covp.mi <- apply(covp.mi.all, 1, mean)

# Median width of 95% CI
ciw.mso <- apply(ciw.mso.all, 1, median)
ciw.msc <- apply(ciw.msc.all, 1, median)
ciw.wid <- apply(ciw.wid.all, 1, median)
ciw.bis <- apply(ciw.bis.all, 1, median)
ciw.mis <- apply(ciw.mis.all, 1, median)
ciw.sbis <- apply(ciw.sbis.all, 1, median)
ciw.mi <- apply(ciw.mi.all, 1, median)

# Mean interval score
ints.mso <- apply(ints.mso.all, 1, mean)
ints.msc <- apply(ints.msc.all, 1, mean)
ints.wid <- apply(ints.wid.all, 1, mean)
ints.bis <- apply(ints.bis.all, 1, mean)
ints.mis <- apply(ints.mis.all, 1, mean)
ints.sbis <- apply(ints.sbis.all, 1, mean)
ints.mi <- apply(ints.mi.all, 1, mean)

# Median interval score
ints.mso.med <- apply(ints.mso.all, 1, median)
ints.msc.med <- apply(ints.msc.all, 1, median)
ints.wid.med <- apply(ints.wid.all, 1, median)
ints.bis.med <- apply(ints.bis.all, 1, median)
ints.mis.med <- apply(ints.mis.all, 1, median)
ints.mi.med <- apply(ints.mi.all, 1, median)
ints.sbis.med <- apply(ints.sbis.all, 1, median)

# Output labels
# CD - Complete Data
# CC - Complete Cases
# IS - Impute-Select
# RR - Rubin's Rules
# ER - Efron's Rules
# SI - wide model with Single Imputation
# MI - wide model with Multiple Imputation

tabints <- rbind(ints.mso, ints.msc, ints.bis, ints.mis, ints.sbis, ints.wid, ints.mi)
tabints.med <- rbind(ints.mso.med, ints.msc.med, ints.bis.med, ints.mis.med, ints.sbis.med, ints.wid.med, ints.mi.med)
row.names(tabints) <- c('CD', 'CC', 'IS', 'RR', 'ER', 'SI', 'MI')
colnames(tabints) <- c('b0', 'b1', 'b2', 'b3', 'b4', 'b5', 'b6')
row.names(tabints.med) <- c('CD', 'CC', 'IS', 'RR', 'ER', 'SI', 'MI')
colnames(tabints.med) <- c('b0', 'b1', 'b2', 'b3', 'b4', 'b5', 'b6')

tabbeta <- rbind(bias.mso, bias.msc, bias.bis, bias.mis, bias.sbis, bias.wid, bias.mi)
row.names(tabbeta) <- c('CD', 'CC', 'IS', 'RR', 'ER', 'SI', 'MI')
colnames(tabbeta) <- c('b0', 'b1', 'b2', 'b3', 'b4', 'b5', 'b6', 'Overall')

tabtime <- rbind(timeavg.mso, timeavg.msc, timeavg.bis, timeavg.mis, timeavg.sbis, timeavg.wid, timeavg.mi)
row.names(tabtime) <- c('CD', 'CC', 'IS', 'RR', 'ER', 'SI', 'MI')
colnames(tabtime) <- c("Mean CPU time", 'SD')

tabfp <- rbind(fpavg.mso, fpavg.msc, fpavg.bis, fpavg.mis, fpavg.sbis, fpavg.wid, fpavg.mi)
row.names(tabfp) <- c('CD', 'CC', 'IS', 'RR', 'ER', 'SI', 'MI')
colnames(tabfp) <- c("Mean FP", 'SD')

tabfn <- rbind(fnavg.mso, fnavg.msc, fnavg.bis, fnavg.mis, fnavg.sbis, fnavg.wid, fnavg.mi)
row.names(tabfn) <- c('CD', 'CC', 'IS', 'RR', 'ER', 'SI', 'MI')
colnames(tabfn) <- c("Mean FN", 'SD')

tabcovp <- rbind(covp.mso, covp.msc, covp.bis, covp.mis, covp.sbis, covp.wid, covp.mi)
row.names(tabcovp) <- c('CD', 'CC', 'IS', 'RR', 'ER', 'SI', 'MI')
colnames(tabcovp) <- c('b0', 'b1', 'b2', 'b3', 'b4', 'b5', 'b6')

tabciw <- rbind(ciw.mso, ciw.msc, ciw.bis, ciw.mis, ciw.sbis, ciw.wid, ciw.mi)
row.names(tabciw) <- c('CD', 'CC', 'IS', 'RR', 'ER', 'SI', 'MI')
colnames(tabciw) <- c('b0', 'b1', 'b2', 'b3', 'b4', 'b5', 'b6')

tabse <- rbind(c(sd.mso1[1], sd.mso.emp[1], sd.mso1[2], sd.mso.emp[2], sd.mso1[3], sd.mso.emp[3], sd.mso1[4], sd.mso.emp[4],  sd.mso1[5], sd.mso.emp[5], sd.mso1[6], sd.mso.emp[6], sd.mso1[7], sd.mso.emp[7]),
               
               c(sd.msc1[1], sd.msc.emp[1], sd.msc1[2], sd.msc.emp[2], sd.msc1[3], sd.msc.emp[3], sd.msc1[4], sd.msc.emp[4], sd.msc1[5], sd.msc.emp[5], sd.msc1[6], sd.msc.emp[6], sd.msc1[7], sd.msc.emp[7]),
               
               c(sd.bis1[1], sd.bis.emp[1], sd.bis1[2], sd.bis.emp[2], sd.bis1[3], sd.bis.emp[3], sd.bis1[4], sd.bis.emp[4], sd.bis1[5], sd.bis.emp[5], sd.bis1[6], sd.bis.emp[6], sd.bis1[7], sd.bis.emp[7]),
               
               c(sd.mis1[1], sd.mis.emp[1], sd.mis1[2], sd.mis.emp[2], sd.mis1[3], sd.mis.emp[3], sd.mis1[4], sd.mis.emp[4], sd.mis1[5], sd.mis.emp[5], sd.mis1[6], sd.mis.emp[6], sd.mis1[7], sd.mis.emp[7]),
               
               c(sd.sbis1[1], sd.sbis.emp[1], sd.sbis1[2], sd.sbis.emp[2], sd.sbis1[3], sd.sbis.emp[3], sd.sbis1[4], sd.sbis.emp[4], sd.sbis1[5], sd.sbis.emp[5], sd.sbis1[6], sd.sbis.emp[6], sd.sbis1[7], sd.sbis.emp[7]),
               
               c(sd.wid1[1], sd.wid.emp[1], sd.wid1[2], sd.wid.emp[2], sd.wid1[3], sd.wid.emp[3], sd.wid1[4], sd.wid.emp[4], sd.wid1[5], sd.wid.emp[5], sd.wid1[6], sd.wid.emp[6], sd.wid1[7], sd.wid.emp[7]),
               
               c(sd.mi1[1], sd.mi.emp[1], sd.mi1[2], sd.mi.emp[2], sd.mi1[3], sd.mi.emp[3], sd.mi1[4], sd.mi.emp[4], sd.mi1[5], sd.mi.emp[5], sd.mi1[6], sd.mi.emp[6], sd.mi1[7], sd.mi.emp[7]))

row.names(tabse) <- c('CD', 'CC', 'IS', 'RR', 'ER', 'SI', 'MI') 
colnames(tabse) <- c('b0-est', 'b0-emp', 'b1-est', 'b1-emp','b2-est', 'b2-emp','b3-est', 'b3-emp','b4-est', 'b4-emp','b5-est', 'b5-emp','b6-est', 'b6-emp')

# Percent bias in SD
tabsep <- matrix(0, nrow=7, ncol=7)
for (i in 0:(length(beta)-1))
  tabsep[,i+1] <- tabse[,2*i+1]/tabse[,2*i+2]
tabsep<- tabsep-1  # >0 if overestimating, <0 if underestimating
row.names(tabsep) <- c('CD', 'CC', 'IS', 'RR', 'ER', 'SI', 'MI')
colnames(tabsep) <- c('b0', 'b1', 'b2', 'b3', 'b4', 'b5', 'b6')

# Relative MSE of beta
tabbetap <- matrix(0, nrow=7, ncol=8)
for (i in 2:7)
  for (j in 1:ncol(tabbeta))
    tabbetap[i,j] <- tabbeta[i,j]/tabbeta[1,j]

for (j in 1:ncol(tabbeta))
  tabbetap[1,j] <- tabbeta[1,j]
row.names(tabbetap) <- c('CD', 'CC', 'IS', 'RR', 'ER', 'SI', 'MI')
colnames(tabbetap) <- c('b0', 'b1', 'b2', 'b3', 'b4', 'b5', 'b6', 'Overall')



##
## Table 1 results 
cat('Sample size n=', n, '\n')
cat('Replication size R=', R, '\n')
cat('Multiple imputation size M=', M, '\n')
cat('Bootstrap size B=', B, '\n')
cat('Corrleation coefficient range', corrange, '\n')
cat('Missing rate \n')
round(misrate.avg, 4)
cat('CPU time \n')
round(tabtime, 4)

cat('Relative MSE to complete data \n')
round(tabbetap, 4)
cat('Percentage bias of SE \n')
round(tabsep, 4)
cat('Coverage probability \n')
round(tabcovp, 4)
cat('Median width of 95% CI \n')
round(tabciw, 4)
cat('Mean interval score of 95% CI for coefficient \n')
round(tabints, 4)


