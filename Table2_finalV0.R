######################################################################################
#                                                                                    #
#                              AOAS Table 2 R code                                   #
#  Comparison of post-imputation-selection estimators in a generalized linear model  #
#                                                                                    #
#                                                                                    #
######################################################################################

source(paste(path, 'libraryglm_finalV0.R', sep=''))

library(mice)

M <- 5        # multiple imputation size
R <- 500      # number of replications
B <- 200      # boostrap size
n <- 189      # sample size in each replication

yname <- 'y'
xname <- c('x1', 'x2', 'x3', 'x4')

mein <- c('', 'pmm', 'logreg', 'logreg', '')    
datanull <- data.frame(varselect=c('(Intercept)', xname))

o <- datanull$varselect

# Estimated beta
beta.mso.all <- NULL
beta.msc.all <- NULL
beta.bis.all <- NULL
beta.mis.all <- NULL
beta.sbis.all <- NULL

# Estimated SD for beta
sd.mso.all <- NULL
sd.msc.all <- NULL
sd.bis.all <- NULL
sd.mis.all <- NULL
sd.sbis.all <- NULL

# False negative rate 
# no false positive since all variables are in the true model
fn.mso.all <- NULL
fn.msc.all <- NULL
fn.bis.all <- NULL
fn.mis.all <- NULL
fn.sbis.all <- NULL

# Coverage probability
covp.mso.all <- NULL
covp.msc.all <- NULL
covp.bis.all <- NULL
covp.mis.all <- NULL
covp.sbis.all <- NULL

# Width of CI
ciw.mso.all <- NULL
ciw.msc.all <- NULL
ciw.bis.all <- NULL
ciw.mis.all <- NULL
ciw.sbis.all <- NULL

# Interval score
ints.mso.all <- NULL
ints.msc.all <- NULL
ints.bis.all <- NULL
ints.mis.all <- NULL
ints.sbis.all <- NULL

#probability of outcomes p(black) and p(white)
pwb.mso.all <- NULL
pwb.msc.all <- NULL
pwb.bis.all <- NULL
pwb.mis.all <- NULL
pwb.sbis.all <- NULL

# SD for the outcome probability
sdpwb.mso.all <- NULL
sdpwb.msc.all <- NULL
sdpwb.bis.all <- NULL
sdpwb.mis.all <- NULL
sdpwb.sbis.all <- NULL

# Coverage probability for outcome
covppwb.mso.all <- NULL
covppwb.msc.all <- NULL
covppwb.bis.all <- NULL
covppwb.mis.all <- NULL
covppwb.sbis.all <- NULL

# Width of CI for outcome
ciwpwb.mso.all <- NULL
ciwpwb.msc.all <- NULL
ciwpwb.bis.all <- NULL
ciwpwb.mis.all <- NULL
ciwpwb.sbis.all <- NULL

# Interval score for outcome
intspwb.mso.all <- NULL
intspwb.msc.all <- NULL
intspwb.bis.all <- NULL
intspwb.mis.all <- NULL
intspwb.sbis.all <- NULL


# CPU time
time.mso.all <- NULL
time.msc.all <- NULL
time.bis.all <- NULL
time.mis.all <- NULL
time.sbis.all <- NULL

misrate.all <- NULL

library(MASS)
fitwt<- glm(low~age+lwt+as.factor(race), data=birthwt, family='binomial')

# True beta
beta <- as.vector(fitwt$coefficient)
birthwt2 <- transform(birthwt, ypred=predict(fitwt))
birthwt2 <- transform(birthwt2, p=exp(ypred)/(1+exp(ypred)))

# True prob(white) is the prob for average whites (average age and weight)
agemw <- mean(birthwt$age[birthwt$race==1])
agemb <- mean(birthwt$age[birthwt$race==2])
lwtmw <- mean(birthwt$lwt[birthwt$race==1])
lwtmb <- mean(birthwt$lwt[birthwt$race==2])


# The following results match Hjort and Claeskens (2003) Table 1 (model 345)
pwhitefit <- predict(fitwt, newdata=data.frame(age=agemw, lwt=lwtmw, race=1), type='response', se.fit=T)
pblackfit <- predict(fitwt, newdata=data.frame(age=agemb, lwt=lwtmb, race=2), type='response', se.fit=T)
pwhite <- pwhitefit[[1]]
pblack <- pblackfit[[1]]
sewhite <- pwhitefit[[2]]
seblack <- pblackfit[[2]]

pwb <- c(pwhite, pblack)  #true p(white) and p(black)

for(r in 1:R)
{
  cat('Replication =', r, '\n')
  
  dataall <- simu.glm(birthwt, 20160808+r)
  datacomp <- dataall[[1]]
  datasimu <- dataall[[2]]
  misrate.all <- rbind(misrate.all, dataall[[3]])
  
  # Multiple imputation
  dataimp <- mice(datasimu[, c(xname, yname)], m=M, method = mein)
  # Single imputation
  dataimp1 <- mice(datasimu[, c(xname, yname)], m=1, method = mein)
  dataimp1 <- complete(dataimp1, 1)  
  
  # MS (model selection) original (complate data)
  time.mso.all <- c(time.mso.all, system.time(outtmp.mso <- lasso.glm(datacomp, yname, xname))[1])
  beta.mso.all <- cbind(beta.mso.all, outtmp.mso[[1]][,2])
  sd.mso.all <- cbind(sd.mso.all, outtmp.mso[[1]][,3])
  #count the number of non-zeros, not all true coefficients are nonzeroes.  
  fn.mso.all <- c(fn.mso.all, sum(beta.mso.all[-1,r]==0))
  covp.mso.all <- cbind(covp.mso.all, outtmp.mso[[1]][,6])
  ciw.mso.all <- cbind(ciw.mso.all, outtmp.mso[[1]][,7])
  ints.mso.all <- cbind(ints.mso.all, outtmp.mso[[1]][,8])
  
  pwb.mso.all <- cbind(pwb.mso.all, outtmp.mso[[2]][,1])
  sdpwb.mso.all <- cbind(sdpwb.mso.all, outtmp.mso[[2]][,2])
  covppwb.mso.all <- cbind(covppwb.mso.all, outtmp.mso[[2]][,5])
  ciwpwb.mso.all <- cbind(ciwpwb.mso.all, outtmp.mso[[2]][,6])
  intspwb.mso.all <- cbind(intspwb.mso.all, outtmp.mso[[2]][,7])
  
  # MS complete case
  idc <- complete.cases(datasimu)
  datasimuc <- datasimu[idc,]
  time.msc.all <- c(time.msc.all, system.time(outtmp.msc <- lasso.glm(datasimuc, yname, xname))[1])
  
  beta.msc.all <- cbind(beta.msc.all, outtmp.msc[[1]][,2])
  sd.msc.all <- cbind(sd.msc.all, outtmp.msc[[1]][,3])  
  #exclude b0 (actually this does not affect results)
  fn.msc.all <- c(fn.msc.all, sum(beta.msc.all[-1,r]==0))  
  covp.msc.all <- cbind(covp.msc.all, outtmp.msc[[1]][,6])
  ciw.msc.all <- cbind(ciw.msc.all, outtmp.msc[[1]][,7])
  ints.msc.all <- cbind(ints.msc.all, outtmp.msc[[1]][,8])  
  pwb.msc.all <- cbind(pwb.msc.all, outtmp.msc[[2]][,1])
  sdpwb.msc.all <- cbind(sdpwb.msc.all, outtmp.msc[[2]][,2])
  covppwb.msc.all <- cbind(covppwb.msc.all, outtmp.msc[[2]][,5])
  ciwpwb.msc.all <- cbind(ciwpwb.msc.all, outtmp.msc[[2]][,6])
  intspwb.msc.all <- cbind(intspwb.msc.all, outtmp.msc[[2]][,7])
  
  
  # Multiple imputation and model selection (MIS) - Rubin's Rules approach
  time.mis.all <- c(time.mis.all, system.time(outtmp.mis <- MIS.glm(dataimp))[1])
  beta.mis.all <- cbind(beta.mis.all, outtmp.mis[[1]][,2])
  sd.mis.all <- cbind(sd.mis.all, outtmp.mis[[1]][,3])  
  fn.mis.all <- c(fn.mis.all, sum(beta.mis.all[-1,r]==0))
  covp.mis.all <- cbind(covp.mis.all, outtmp.mis[[1]][,6])
  ciw.mis.all <- cbind(ciw.mis.all, outtmp.mis[[1]][, 7])
  ints.mis.all <- cbind(ints.mis.all, outtmp.mis[[1]][,8])
  pwb.mis.all <- cbind(pwb.mis.all, outtmp.mis[[2]][,1])
  sdpwb.mis.all <- cbind(sdpwb.mis.all, outtmp.mis[[2]][,2])
  covppwb.mis.all <- cbind(covppwb.mis.all, outtmp.mis[[2]][,5])
  ciwpwb.mis.all <- cbind(ciwpwb.mis.all, outtmp.mis[[2]][,6])
  intspwb.mis.all <- cbind(intspwb.mis.all, outtmp.mis[[2]][,7])
  

  # Bootstrap single imputation then model selection (BIS) - Impute-Select approach
  time.bis.all <- c(time.bis.all, system.time(outtmp <- SBISv2.glm(dataimp1, datasimu))[1])
  outtmp.bis <- outtmp[[1]]
  beta.bis.all <- cbind(beta.bis.all, outtmp.bis[,2])
  sd.bis.all <- cbind(sd.bis.all, outtmp.bis[,3])
  fn.bis.all <- c(fn.bis.all, sum(beta.bis.all[-1,r]==0))
  covp.bis.all <- cbind(covp.bis.all, outtmp.bis[,6])
  ciw.bis.all <- cbind(ciw.bis.all, outtmp.bis[,7])
  ints.bis.all <- cbind(ints.bis.all, outtmp.bis[,8])  
  pwb.bis.all <- cbind(pwb.bis.all, outtmp$outp[,1])
  sdpwb.bis.all <- cbind(sdpwb.bis.all, outtmp$outp[,2])
  covppwb.bis.all <- cbind(covppwb.bis.all, outtmp$covppwb[,1])
  ciwpwb.bis.all <- cbind(ciwpwb.bis.all, outtmp$ciwpwb[,1])
  intspwb.bis.all <- cbind(intspwb.bis.all, outtmp$pwb.ints[,1])
  
  
  ## Smooth bootstrap single imputation then model selection (SBIS) - Efron's Rules approach
  outtmp.sbis <- outtmp[[2]]
  beta.sbis.all <- cbind(beta.sbis.all, outtmp.sbis[,2])
  sd.sbis.all <- cbind(sd.sbis.all, outtmp.sbis[,3])   
  fn.sbis.all <- c(fn.sbis.all, sum(beta.sbis.all[-1,r]==0)) 
  covp.sbis.all <- cbind(covp.sbis.all, outtmp.sbis[,6])
  ciw.sbis.all <- cbind(ciw.sbis.all, outtmp.sbis[,7])
  ints.sbis.all <- cbind(ints.sbis.all, outtmp.sbis[,8])
  pwb.sbis.all <- cbind(pwb.sbis.all, outtmp$outp[,3])
  sdpwb.sbis.all <- cbind(sdpwb.sbis.all, outtmp$outp[,4])
  covppwb.sbis.all <- cbind(covppwb.sbis.all, outtmp$covppwb[,2])
  ciwpwb.sbis.all <- cbind(ciwpwb.sbis.all, outtmp$ciwpwb[,2])
  intspwb.sbis.all <- cbind(intspwb.sbis.all, outtmp$pwb.ints[,2])
}


time.sbis.all <- time.bis.all
colnames(misrate.all) <- c('Missing.rate(LWT)', 'Missing.rate(Race.black)', 'Missing.rate(Race.other)')

sd.mso.emp<- apply(beta.mso.all, 1, sd)
sd.mso1 <- apply(sd.mso.all, 1, mean)  

sd.msc.emp<- apply(beta.msc.all, 1, sd)
sd.msc1 <- apply(sd.msc.all, 1, mean)

sd.bis.emp<- apply(beta.bis.all, 1, sd)
sd.bis1 <- apply(sd.bis.all, 1, mean)

sd.mis.emp<- apply(beta.mis.all, 1, sd)
sd.mis1 <- apply(sd.mis.all, 1, mean)

sd.sbis.emp<- apply(beta.sbis.all, 1, sd)
sd.sbis1 <- apply(sd.sbis.all, 1, mean)

sdpwb.mso.emp<- apply(pwb.mso.all, 1, sd)
sdpwb.mso1 <- apply(sdpwb.mso.all, 1, mean)  

sdpwb.msc.emp<- apply(pwb.msc.all, 1, sd)
sdpwb.msc1 <- apply(sdpwb.msc.all, 1, mean)

sdpwb.bis.emp<- apply(pwb.bis.all, 1, sd)
sdpwb.bis1 <- apply(sdpwb.bis.all, 1, mean)

sdpwb.mis.emp<- apply(pwb.mis.all, 1, sd)
sdpwb.mis1 <- apply(sdpwb.mis.all, 1, mean)

sdpwb.sbis.emp<- apply(pwb.sbis.all, 1, sd)
sdpwb.sbis1 <- apply(sdpwb.sbis.all, 1, mean)


# Missing rate
misrate.avg <- apply(misrate.all, 2, mean)

# Bias (averaged squared loss)
bias.mso <- c(apply((beta.mso.all-beta)^2, 1, mean),  sum((beta.mso.all-beta)^2)/R)
bias.msc <- c(apply((beta.msc.all-beta)^2, 1, mean),  sum((beta.msc.all-beta)^2)/R)
bias.bis <- c(apply((beta.bis.all-beta)^2, 1, mean),  sum((beta.bis.all-beta)^2)/R)
bias.mis <- c(apply((beta.mis.all-beta)^2, 1, mean),  sum((beta.mis.all-beta)^2)/R)
bias.sbis<- c(apply((beta.sbis.all-beta)^2, 1, mean), sum((beta.sbis.all-beta)^2)/R)

biaspwb.mso <- c(apply((pwb.mso.all-pwb)^2, 1, mean))
biaspwb.msc <- c(apply((pwb.msc.all-pwb)^2, 1, mean))
biaspwb.bis <- c(apply((pwb.bis.all-pwb)^2, 1, mean))
biaspwb.mis <- c(apply((pwb.mis.all-pwb)^2, 1, mean))
biaspwb.sbis<- c(apply((pwb.sbis.all-pwb)^2, 1, mean))

# CPU time
timeavg.mso <- c(mean(time.mso.all), sd(time.mso.all))
timeavg.msc <- c(mean(time.msc.all), sd(time.msc.all))
timeavg.bis <- c(mean(time.bis.all), sd(time.bis.all))
timeavg.mis <- c(mean(time.mis.all), sd(time.mis.all))
timeavg.sbis <- timeavg.bis

# FN, does not count beta0
fnavg.mso <-c(mean(fn.mso.all), sd(fn.mso.all))
fnavg.msc <-c(mean(fn.msc.all), sd(fn.msc.all))
fnavg.bis <-c(mean(fn.bis.all), sd(fn.bis.all))
fnavg.mis <-c(mean(fn.mis.all), sd(fn.mis.all))
fnavg.sbis <-c(mean(fn.sbis.all), sd(fn.sbis.all))

# Coverage probability
covp.mso <- apply(covp.mso.all, 1, mean)
covp.msc <- apply(covp.msc.all, 1, mean)
covp.bis <- apply(covp.bis.all, 1, mean)
covp.mis <- apply(covp.mis.all, 1, mean)
covp.sbis <- apply(covp.sbis.all, 1, mean)

# Median width of 95% CI
ciw.mso <- apply(ciw.mso.all, 1, median)
ciw.msc <- apply(ciw.msc.all, 1, median)
ciw.bis <- apply(ciw.bis.all, 1, median)
ciw.mis <- apply(ciw.mis.all, 1, median)
ciw.sbis <- apply(ciw.sbis.all, 1, median)

# Mean interval score
ints.mso <- apply(ints.mso.all, 1, mean)
ints.msc <- apply(ints.msc.all, 1, mean)
ints.bis <- apply(ints.bis.all, 1, mean)
ints.mis <- apply(ints.mis.all, 1, mean)
ints.sbis <- apply(ints.sbis.all, 1, mean)

# Median interval score
ints.mso.med <- apply(ints.mso.all, 1, median)
ints.msc.med <- apply(ints.msc.all, 1, median)
ints.bis.med <- apply(ints.bis.all, 1, median)
ints.mis.med <- apply(ints.mis.all, 1, median)
ints.sbis.med <- apply(ints.sbis.all, 1, median)

# Coverage probability
covppwb.mso <- apply(covppwb.mso.all, 1, mean)
covppwb.msc <- apply(covppwb.msc.all, 1, mean)
covppwb.bis <- apply(covppwb.bis.all, 1, mean)
covppwb.mis <- apply(covppwb.mis.all, 1, mean)
covppwb.sbis <- apply(covppwb.sbis.all, 1, mean)

# Median width of 95% CI
ciwpwb.mso <- apply(ciwpwb.mso.all, 1, median)
ciwpwb.msc <- apply(ciwpwb.msc.all, 1, median)
ciwpwb.bis <- apply(ciwpwb.bis.all, 1, median)
ciwpwb.mis <- apply(ciwpwb.mis.all, 1, median)
ciwpwb.sbis <- apply(ciwpwb.sbis.all, 1, median)

# Mean interval score of 95% CI
intspwb.mso <- apply(intspwb.mso.all, 1, mean)
intspwb.msc <- apply(intspwb.msc.all, 1, mean)
intspwb.bis <- apply(intspwb.bis.all, 1, mean)
intspwb.mis <- apply(intspwb.mis.all, 1, mean)
intspwb.sbis <- apply(intspwb.sbis.all, 1, mean)

# Mean interval score of 95% CI
intspwb.mso.med <- apply(intspwb.mso.all, 1, median)
intspwb.msc.med <- apply(intspwb.msc.all, 1, median)
intspwb.bis.med <- apply(intspwb.bis.all, 1, median)
intspwb.mis.med <- apply(intspwb.mis.all, 1, median)
intspwb.sbis.med <- apply(intspwb.sbis.all, 1, median)

tabints <- rbind(ints.mso, ints.msc, ints.bis, ints.mis, ints.sbis)
tabintspwb <- rbind(intspwb.mso, intspwb.msc, intspwb.bis, intspwb.mis, intspwb.sbis)

# Output labels
# CD - Complete Data
# CC - Complete Cases
# IS - Impute-Select
# RR - Rubin's Rules
# ER - Efron's Rules

row.names(tabints) <- c('CD', 'CC', 'IS', 'RR', 'ER')
colnames(tabints) <- c('b0', 'Age', 'LWT', 'Race.black', 'Race.other')
row.names(tabintspwb) <- c('CD', 'CC', 'IS', 'RR', 'ER')
colnames(tabintspwb) <- c('White', 'Black')

tabints.med <- rbind(ints.mso.med, ints.msc.med, ints.bis.med, ints.mis.med, ints.sbis.med)
tabintspwb.med <- rbind(intspwb.mso.med, intspwb.msc.med, intspwb.bis.med, intspwb.mis.med, intspwb.sbis.med)
row.names(tabints.med) <- c('CD', 'CC', 'IS', 'RR', 'ER')
colnames(tabints.med) <- c('b0', 'Age', 'LWT', 'Race.black', 'Race.other')
row.names(tabintspwb.med) <- c('CD', 'CC', 'IS', 'RR', 'ER')
colnames(tabintspwb.med) <- c('White', 'Black')


tabbeta <- rbind(bias.mso, bias.msc, bias.bis, bias.mis, bias.sbis)
row.names(tabbeta) <- c('CD', 'CC', 'IS', 'RR', 'ER')
colnames(tabbeta) <- c('b0', 'Age', 'LWT', 'Race.black', 'Race.other', 'Overall')


tabtime <- rbind(timeavg.mso, timeavg.msc, timeavg.bis, timeavg.mis, timeavg.sbis)
row.names(tabtime) <- c('CD', 'CC', 'IS', 'RR', 'ER')
colnames(tabtime) <- c("Mean CPU time", 'SD')

tabfn <- rbind(fnavg.mso, fnavg.msc, fnavg.bis, fnavg.mis, fnavg.sbis)
row.names(tabfn) <- c('CD', 'CC', 'IS', 'RR', 'ER')
colnames(tabfn) <- c("Mean FN", 'SD')


tabcovp <- rbind(covp.mso, covp.msc, covp.bis, covp.mis, covp.sbis)
row.names(tabcovp) <- c('CD', 'CC', 'IS', 'RR', 'ER')
colnames(tabcovp) <- c('b0', 'Age', 'LWT', 'Race.black', 'Race.other')


tabciw <- rbind(ciw.mso, ciw.msc, ciw.bis, ciw.mis, ciw.sbis)
row.names(tabciw) <- c('CD', 'CC', 'IS', 'RR', 'ER')
colnames(tabciw) <- c('b0', 'Age', 'LWT', 'Race.black', 'Race.other')


tabse <- rbind(c(sd.mso1[1], sd.mso.emp[1], sd.mso1[2], sd.mso.emp[2], sd.mso1[3], sd.mso.emp[3], sd.mso1[4], sd.mso.emp[4],  sd.mso1[5], sd.mso.emp[5]),
               
               c(sd.msc1[1], sd.msc.emp[1], sd.msc1[2], sd.msc.emp[2], sd.msc1[3], sd.msc.emp[3], sd.msc1[4], sd.msc.emp[4], sd.msc1[5], sd.msc.emp[5]),
               
               c(sd.bis1[1], sd.bis.emp[1], sd.bis1[2], sd.bis.emp[2], sd.bis1[3], sd.bis.emp[3], sd.bis1[4], sd.bis.emp[4], sd.bis1[5], sd.bis.emp[5]),
               
               c(sd.mis1[1], sd.mis.emp[1], sd.mis1[2], sd.mis.emp[2], sd.mis1[3], sd.mis.emp[3], sd.mis1[4], sd.mis.emp[4], sd.mis1[5], sd.mis.emp[5]),
               
               c(sd.sbis1[1], sd.sbis.emp[1], sd.sbis1[2], sd.sbis.emp[2], sd.sbis1[3], sd.sbis.emp[3], sd.sbis1[4], sd.sbis.emp[4], sd.sbis1[5], sd.sbis.emp[5]))

row.names(tabse) <- c('CD', 'CC', 'IS', 'RR', 'ER')
colnames(tabse) <- c('b0-est', 'b0-emp', 'Age-est', 'Age-emp','LWT-est', 'LWT-emp','Race.black-est', 'Race.black-emp','Race.other-est', 'Race.other-emp')

tabsep <- matrix(0, nrow=5, ncol=5)
for (i in 0:(length(beta)-1))
  tabsep[,i+1] <- tabse[,2*i+1]/tabse[,2*i+2]
tabsep<- tabsep-1  
row.names(tabsep) <- c('CD', 'CC', 'IS', 'RR', 'ER')
colnames(tabsep) <- c('b0', 'Age', 'LWT', 'Race.black', 'Race.other')

tabbetap <- matrix(0, nrow=5, ncol=6)
for (i in 2:5)
  for (j in 1:ncol(tabbeta))
    tabbetap[i,j] <- tabbeta[i,j]/tabbeta[1,j]

for (j in 1:ncol(tabbeta))
  tabbetap[1,j] <- tabbeta[1,j]
row.names(tabbetap) <- c('CD', 'CC', 'IS', 'RR', 'ER')
colnames(tabbetap) <- c('b0', 'Age', 'LWT', 'Race.black', 'Race.other', 'Overall')



tabpwb <- rbind(biaspwb.mso, biaspwb.msc, biaspwb.bis, biaspwb.mis, biaspwb.sbis)
row.names(tabpwb) <- c('CD', 'CC', 'IS', 'RR', 'ER')
colnames(tabpwb) <- c('White', 'Black')


tabcovppwb <- rbind(covppwb.mso, covppwb.msc, covppwb.bis, covppwb.mis, covppwb.sbis)
row.names(tabcovppwb) <- c('CD', 'CC', 'IS', 'RR', 'ER')
colnames(tabcovppwb) <- c('White', 'Black')


tabciwpwb <- rbind(ciwpwb.mso, ciwpwb.msc, ciwpwb.bis, ciwpwb.mis, ciwpwb.sbis)
row.names(tabciwpwb) <- c('CD', 'CC', 'IS', 'RR', 'ER')
colnames(tabciwpwb) <- c('White', 'Black')


tabsepwb <- rbind(c(sdpwb.mso1[1], sdpwb.mso.emp[1], sdpwb.mso1[2], sdpwb.mso.emp[2]),
                  
                  c(sd.msc1[1], sd.msc.emp[1], sd.msc1[2], sd.msc.emp[2]),
                  
                  c(sd.bis1[1], sd.bis.emp[1], sd.bis1[2], sd.bis.emp[2]),
                  
                  c(sd.mis1[1], sd.mis.emp[1], sd.mis1[2], sd.mis.emp[2]),
                  
                  c(sd.sbis1[1], sd.sbis.emp[1], sd.sbis1[2], sd.sbis.emp[2]))

row.names(tabsepwb) <- c('CD', 'CC', 'IS', 'RR', 'ER')
colnames(tabsepwb) <- c('White-est', 'White-emp', 'Black-est', 'Black-emp')

tabsepwbp <- matrix(0, nrow=5, ncol=2)
for (i in 0:(length(pwb)-1))
  tabsepwbp[,i+1] <- tabsepwb[,2*i+1]/tabsepwb[,2*i+2]
tabsepwbp<- tabsepwbp-1
row.names(tabsepwbp) <- c('FU', 'CC', 'IS ', 'MI', 'BI')
colnames(tabsepwbp) <- c('White', 'Black')

tabpwbp <- matrix(0, nrow=5, ncol=2)
for (i in 2:5)
  for (j in 1:ncol(tabpwb))
    tabpwbp[i,j] <- tabpwb[i,j]/tabpwb[1,j]

for (j in 1:ncol(tabpwb))
  tabpwbp[1,j] <- tabpwb[1,j]
row.names(tabpwbp) <- c('FU', 'CC', 'IS ', 'MI', 'BI')
colnames(tabpwbp) <- c('White', 'Black')


# Table 2 results
cat('Sample size n=', n, '\n')
cat('Replication size R=', R, '\n')
cat('Multiple imputation size M=', M, '\n')
cat('Bootstrap size B (smoothed) =', B, '\n')
cat('Missing rate \n')
round(misrate.avg, 4)
cat('CPU time \n')
round(tabtime, 4)

cat('Relative MSE to complete data for beta \n')
round(tabbetap, 4)
cat('Percentage bias of SE for beta \n')
round(tabsep, 4)
cat('Mean interval score of 95% CI for coefficient \n')
round(tabints, 4)
cat('Coverage probability for beta \n')
round(tabcovp, 4)
cat('Median width of 95% CI for beta \n')
round(tabciw, 4)


cat('Relative MSE to complete data for predicted probability \n')
round(tabpwbp, 4)
cat('Percentage bias of SE for predicted probability \n')
round(tabsepwbp, 4)
cat('Mean interval score of 95% CI for predicted probability \n')
round(tabintspwb, 4)
cat('Coverage probability  for predicted probability \n')
round(tabcovppwb, 4)
cat('Median width of 95% CI for predicted probability \n')
round(tabciwpwb, 4)