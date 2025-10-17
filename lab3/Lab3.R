######################################################################################
# BIOL 570 Lab 3: Generalized linear mixed effects models for salmon populations
# Written by Mark Lewis February 2021, based on code from Stephanie Peacock.
######################################################################################
install.packages("lme4")
library(lme4)

install.packages("gplots")
library(gplots)

install.packages('Rcpp')
library(Rcpp)

rm(list=ls())
data<-read.csv("lab3/Pink_S-R_data.csv")
head(data)
data$Survival<-log(data$Recruits/data$Spawners)
data$Population<-paste(data$River,data$Odd_Even)

# Adding louse exposure category:
data$exposure<-factor("unexposed", levels=c("exposed", "unexposed"))
data$exposure[which(data$Area==12)]<-factor("exposed")

# Infestations occured in spring 2000-2005; better treatment plans in place by spring 2006
# Fallow management intervention in 2003 (affected fish returning in 2004)
data$time_period<-factor("before", levels=c("before", "during", "after", "fallow"))
data$time_period[which(data$Return_Yr>=2002&data$Return_Yr <=2006)]<-factor("during")
data$time_period[which(data$Return_Yr>=2007)]<-factor("after")
data$time_period[which(data$Return_Yr ==2004)]<-factor("fallow")

#Exposure category
data$exposure_category<-as.factor(paste(data$exposure, data$time_period, sep=" "))

#Fit model 1
fit1a<-lm(Survival~Population:Spawners,data=data)
fit1b<-lm(Survival~0+exposure_category+Population:Spawners, data=data)

#Bar chart for model 1
bp1<-barplot2(summary(fit1b)$coefficients[1:8,1], plot.ci=TRUE, ci.l=summary(fit1b)$coefficients[1:8,1]-1.96*summary(fit1b)$coefficients[1:8,2], ci.u=summary(fit1b)$coefficients[1:8,1]+1.96*summary(fit1b)$coefficients[1:8,2], col=grey(c(0.4,0.4,0.4,0.4,0.8,0.8,0.8,0.8)), names.arg=FALSE,ylab="Population growth rate")
factorlabels<-c("after","before","during","fallow","after","before","during","fallow")
for(i in 1:8) text(bp1[i], -0.5, factorlabels[i], xpd=TRUE)
legend("topright", fill=grey(c(0.4, 0.8)), c("Exposed", "Unexposed"), bty="n")


# Log likelihood LRT and AIC for model 1
ll1a<-logLik(fit1a)
ll1b<-logLik(fit1b)
teststat1 <- -2 * (as.numeric(ll1a)-as.numeric(ll1b))
p.val1 <- pchisq(teststat1, df = 7, lower.tail = FALSE)
AIC1a <- AIC(fit1a)
AIC1b <- AIC(fit1b)

#Fit model 3 REML=TRUE
fit3a<-lmer(Survival~Population:Spawners+(1|Return_Yr), REML=TRUE, data=data)
fit3b<-lmer(Survival~0+exposure_category+Population:Spawners+(1|Return_Yr), REML=TRUE, data=data)

#AIC for model 3 REML=TRUE
ll3a<-logLik(fit3a)
ll3b<-logLik(fit3b)
AIC3a <- AIC(fit3a)
AIC3b <- AIC(fit3b)

#Fit model 4 REML=TRUE
fit4a<-lmer(Survival~Population:Spawners+(1|Area), REML=TRUE, data=data)
fit4b<-lmer(Survival~0+exposure_category+Population:Spawners+(1|Area), REML=TRUE, data=data)

#AIC for model 4 REML=TRUE
ll4a<-logLik(fit4a)
ll4b<-logLik(fit4b)
AIC4a <- AIC(fit4a)
AIC4b <- AIC(fit4b)

#Fit model 7 REML=TRUE
fit7a<-lmer(Survival~Population:Spawners+(1|Return_Yr/Area), REML=TRUE, data=data)
fit7b<-lmer(Survival~0+exposure_category+Population:Spawners+(1|Return_Yr/Area), REML=TRUE, data=data)

#AIC for model 7 REML=TRUE
ll7a<-logLik(fit7a)
ll7b<-logLik(fit7b)
AIC7a <- AIC(fit7a)
AIC7b <- AIC(fit7b)

#Fit model 3 REML=FALSE
fit3a_fixed<-lmer(Survival~Population:Spawners+(1|Return_Yr), REML=FALSE, data=data)
fit3b_fixed<-lmer(Survival~0+exposure_category+Population:Spawners+(1|Return_Yr), REML=FALSE, data=data)

# Log likelihood, AIC and LRT for model 3 REML=FALSE
ll3_fixeda<-logLik(fit3a_fixed)
ll3b_fixed<-logLik(fit3b_fixed)
AIC3a_fixed <- AIC(fit3a_fixed)
AIC3b_fixed <- AIC(fit3b_fixed)
teststat3 <- -2 * (as.numeric(ll3a_fixed)-as.numeric(ll3b_fixed))
p.val3 <- pchisq(teststat3, df = 7, lower.tail = FALSE)

#Fit model 4 REML=FALSE
fit4a_fixed<-lmer(Survival~Population:Spawners+(1|Area), REML=FALSE, data=data)
fit4b_fixed<-lmer(Survival~0+exposure_category+Population:Spawners+(1|Area), REML=FALSE, data=data)

# Log likelihood, AIC and LRT for model 4 REML=FALSE
ll4a_fixed<-logLik(fit4a_fixed)
ll4b_fixed<-logLik(fit4b_fixed)
AIC4a_fixed <- AIC(fit4a_fixed)
AIC4b_fixed <- AIC(fit4b_fixed)
teststat4 <- -2 * (as.numeric(ll4a_fixed)-as.numeric(ll4b_fixed))
p.val4 <- pchisq(teststat7_fixed, df = 7, lower.tail = FALSE)

#Fit model 7 REML=FALSE
fit7a_fixed<-lmer(Survival~Population:Spawners+(1|Return_Yr/Area), REML=FALSE, data=data)
fit7b_fixed<-lmer(Survival~0+exposure_category+Population:Spawners+(1|Return_Yr/Area), REML=FALSE, data=data)

# Log likelihood, AIC and LRT for model 7 REML=FALSE
ll7a_fixed<-logLik(fit7a_fixed)
ll7b_fixed<-logLik(fit7b_fixed)
AIC7a_fixed <- AIC(fit7a_fixed)
AIC7b_fixed <- AIC(fit7b_fixed)
teststat7_fixed <- -2 * (as.numeric(ll7a_fixed)-as.numeric(ll7b_fixed))
p.val7 <- pchisq(teststat7_fixed, df = 7, lower.tail = FALSE)

#Bar chart for model 7 REML=FALSE
fixef.7b<-summary(fit7b)$coefficients[1:8,1:2]
bp7<-barplot2(summary(fit7b)$coefficients[1:8,1], plot.ci=TRUE, ci.l=summary(fit7b)$coefficients[1:8,1]-1.96*summary(fit7b)$coefficients[1:8,2],ci.u=summary(fit7b)$coefficients[1:8,1]+1.96*summary(fit7b)$coefficients[1:8,2],col=grey(c(0.4,0.4,0.4,0.4,0.8,0.8,0.8,0.8)),names.arg=FALSE,ylab="Population growth rate")
for(i in 1:8) text(bp1[i], -0.5, factorlabels[i], xpd=TRUE)
legend("topright", fill=grey(c(0.4, 0.8)), c("Exposed", "Unexposed"), bty="n")

NPopulation<-length(unique(data$Population))
NReturn_Yr<-length(unique(data$Return_Yr))
NArea<-length(unique(data$Area))

              