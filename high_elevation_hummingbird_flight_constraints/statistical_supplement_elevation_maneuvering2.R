
# Statistical Supplement
##########################################
# In support of the following manuscript submitted to Current Biology:
# Title: Mechanical constraints on flight at high elevation decrease maneuvering performance
# Manuscript Authors: Paolo S. Segre, Roslyn Dakin, Tyson J. G. Read, Andrew D. Straw, and Douglas L. Altshuler

# Analysis and Script by: Roslyn Dakin
# Last Modified: August 25, 2016
##########################################

# the following packages are used to perform the analyses:
# nlme
# MuMIn
# multcomp
# qvalue
# visreg

# note that qvalue can be installed by running the following two lines: 
# source("https://bioconductor.org/biocLite.R")
# biocLite("qvalue")

library(nlme)
library(MuMIn)
library(multcomp)
library(qvalue)
library(visreg)

gas <- read.csv('gas_substitution_experiment.csv')
trans <- read.csv('translocation_experiment.csv')
tableS3 <- read.csv('bird_traits.csv')

head(trans) # data from elevation translocation experiment
head(gas) # data from gas substitution experiment
head(tableS3) # individual traits



# Table S3
##########################################

tableS3 # individual morphology & load-lifting parameters for 16 Anna's hummingbirds

# comparing birds studied at two sites
ttests <- data.frame(trait=colnames(tableS3)[3:11], tstat=NA, pvalue=NA)
for(i in 1:9){
	ttests$tstat[i] <- t.test(as.formula(paste(ttests$trait[i],'~site',sep='')), data=tableS3)$statistic
	ttests$pvalue[i] <- round(t.test(as.formula(paste(ttests$trait[i],'~site',sep='')), data=tableS3)$p.value,5)
}
ttests # birds at the two sites (California & British Columbia) differ in several morphological and wingbeat parameters



# Elevation translocation experiment
##########################################

metrics <- c('vel', 'hor_accel', 'hor_decel', 'pitch_up', 'pitch_down', 'yaw', 'PRT_deg_turned', 'PRT_time', 'arc_rad', 'arc_vel', 'arc_cent_accel', 'PRT_percent') 
trans$treatment <- factor(trans$treatment, levels=c('sea level', 'elev')) # sea level is the control

# take the absolute value of the negative metrics
# such that higher values = higher performance, for all performance metrics 
trans[trans$performance.metric=='hor_decel', 'avg.performance'] <- abs(trans[trans$performance.metric=='hor_decel', 'avg.performance'])
trans[trans$performance.metric=='pitch_down', 'avg.performance'] <- abs(trans[trans$performance.metric=='pitch_down', 'avg.performance'])

# for PRT%, the no. maneuvers is the total no. complex turns
trans[trans$performance.metric=='PRT_percent', 'no.maneuvers'] <- trans[trans$performance.metric=='PRT_time', 'no.maneuvers'] + trans[trans$performance.metric=='arc_vel', 'no.maneuvers']

# Table S1: elevation translocation effect
# sample sizes, model coefficients, and signficance tests for elevation

tableS1 <- data.frame(metric=metrics, low.elev.traj=NA, hi.elev.traj=NA, treatment.coef=NA, coef.lower=NA, coef.upper=NA, test.stat=NA, p.value=NA)

# effect of no. maneuvers
no.man1 <- data.frame(metric=metrics, deltaAICc=NA, coef=NA, coef.lower=NA, coef.upper=NA)

# model residuals
mod.resid <- list()

# final models
mod.finalTRANS <- list()

# fit all 12 models and save the results:
for(i in 1:12){
	tableS1[i, 2:3] <- tapply(subset(trans, performance.metric==metrics[i])$no.maneuvers, subset(trans, performance.metric==metrics[i])$treatment, 'sum')
	
	# each model includes fixed effects for competitor presence, body mass, load lifting, capture group, and elevation
	# we use AICc to evaluate whether including no.maneuvers as a covariate improves model fit enough to justify the additional parameter
	modA <- lme(avg.performance ~ solo.competition + body.mass + mass.lifted + capture.group + treatment, random=~1|bird, data=subset(trans, performance.metric==metrics[i]), method='ML')
	modB <- lme(avg.performance ~ solo.competition + body.mass + mass.lifted + capture.group + no.maneuvers + treatment, random=~1|bird, data=subset(trans, performance.metric==metrics[i]), method='ML')
	
	no.man1$deltaAICc[i] <- AICc(modB) - AICc(modA) 
			
	if( (AICc(modB) - AICc(modA)) <= -2 ){
		mod <- update(modB, method='REML')
		no.man1$coef[i] <- summary(mod)$tTable['no.maneuvers', 1]
		no.man1$coef.lower[i] <- summary(mod)$tTable['no.maneuvers', 1] - 1.96*summary(mod)$tTable['no.maneuvers', 2]
		no.man1$coef.upper[i] <- summary(mod)$tTable['no.maneuvers', 1] + 1.96*summary(mod)$tTable['no.maneuvers', 2]
			
	} else {
		mod <- update(modA, method='REML')
	}
		
	mod.resid[[i]] <- residuals(mod)
	mod.finalTRANS[[i]] <- mod
	
	tableS1$treatment.coef[i] <- summary(mod)$tTable['treatmentelev', 1]
	tableS1$coef.lower[i] <- summary(mod)$tTable['treatmentelev', 1] - 1.96*summary(mod)$tTable['treatmentelev', 2]
	tableS1$coef.upper[i] <- summary(mod)$tTable['treatmentelev', 1] + 1.96*summary(mod)$tTable['treatmentelev', 2]
	tableS1$test.stat[i] <- summary(mod)$tTable['treatmentelev', 4]
	tableS1$p.value[i] <- summary(mod)$tTable['treatmentelev', 5]
}
	
par(mfrow=c(4,4)) # check normality of model residuals
for(i in 1:12) {hist(mod.resid[[i]], main=metrics[i])}
for(i in 1:12) {print(shapiro.test(mod.resid[[i]]))}
 
tableS1 # compare high and low elevation

no.man1 # no. of maneuvers was retained in the model for PRT time only
summary(mod.finalTRANS[[8]])
visreg(mod.finalTRANS[[8]], 'no.maneuvers', main='residual PRT time (s)') # birds that performed more PRTs executed them in less time



# Gas substitution experiment
##########################################

gas$treatment <- factor(gas$treatment, levels=c('air', 'nit', 'hel')) # air is the control for both comparisons

# take the absolute value of the negative metrics
gas[gas$performance.metric=='hor_decel', 'avg.performance'] <- abs(gas[gas$performance.metric=='hor_decel', 'avg.performance'])
gas[gas$performance.metric=='pitch_down', 'avg.performance'] <- abs(gas[gas$performance.metric=='pitch_down', 'avg.performance'])

# for PRT%, the no. maneuvers is the total no. complex turns
gas[gas$performance.metric=='PRT_percent', 'no.maneuvers'] <- gas[gas$performance.metric=='PRT_time', 'no.maneuvers'] + gas[gas$performance.metric=='arc_vel', 'no.maneuvers']

# Table S2a: air density
# sample sizes, model coefficients, and signficance tests for air density

tableS2a <- data.frame(metric=metrics, air.traj=NA, hypodense.traj=NA, treatment.coef=NA, coef.lower=NA, coef.upper=NA, test.stat=NA, p.value=NA)

# Table S2b: oxygen level
# sample sizes, model coefficients, and signficance tests for oxygen

tableS2b <- data.frame(metric=metrics, air.traj=NA, hypoxia.traj=NA, treatment.coef=NA, coef.lower=NA, coef.upper=NA, test.stat=NA, p.value=NA)

# Table S2c: time in captivity
# sample sizes, model coefficients, and significance tests for days in captivity

tableS2c <- data.frame(metric=metrics, traj=NA, days.captive.coef=NA, coef.lower=NA, coef.upper=NA, test.stat=NA, p.value=NA)

# effect of no. maneuvers
no.man2 <- data.frame(metric=metrics, deltaAICc=NA, coef=NA, coef.lower=NA, coef.upper=NA)

# model residuals
mod.resid <- list()

# final models
mod.finalGAS <- list()

# fit 12 models and enter the results:

for(i in 1:12){
	tableS2a[i, 2:3] <- tapply(subset(gas, performance.metric==metrics[i])$no.maneuvers, subset(gas, performance.metric==metrics[i])$treatment, 'sum')[c(1,3)]
	tableS2b[i, 2:3] <- tapply(subset(gas, performance.metric==metrics[i])$no.maneuvers, subset(gas, performance.metric==metrics[i])$treatment, 'sum')[c(1,2)]
	tableS2c[i, 2] <- sum(tapply(subset(gas, performance.metric==metrics[i])$no.maneuvers, subset(gas, performance.metric==metrics[i])$treatment, 'sum'))
	
	# each model includes fixed effects for competitor presence, body mass, load lifting, days in captivity, and treatment
	# we use AICc to evaluate whether including no.maneuvers as a covariate improves model fit enough to justify the additional parameter
	modA <- lme(avg.performance ~ solo.competition + body.mass + mass.lifted + days.after + treatment, random=~1|bird, data=subset(gas, performance.metric==metrics[i]), method='ML')
	modB <- lme(avg.performance ~ solo.competition + body.mass + mass.lifted + days.after + no.maneuvers + treatment, random=~1|bird, data=subset(gas, performance.metric==metrics[i]), method='ML')
					
	no.man2$deltaAICc[i] <- AICc(modB) - AICc(modA)
	
	if( (AICc(modB) - AICc(modA)) <= -2 ){
		mod <- update(modB, method='REML')	
		no.man2$coef[i] <- summary(mod)$tTable['no.maneuvers', 1]
		no.man2$coef.lower[i] <- summary(mod)$tTable['no.maneuvers', 1] - 1.96*summary(mod)$tTable['no.maneuvers', 2]
		no.man2$coef.upper[i] <- summary(mod)$tTable['no.maneuvers', 1] + 1.96*summary(mod)$tTable['no.maneuvers', 2]
	} else {
		mod <- update(modA, method='REML')
	}

	mod.resid[[i]] <- residuals(mod)
	mod.finalGAS[[i]] <- mod
	
	tableS2b$treatment.coef[i] <- summary(mod)$tTable['treatmentnit', 1] 
	tableS2b$coef.lower[i] <- summary(mod)$tTable['treatmentnit', 1] - 1.96*summary(mod)$tTable['treatmentnit', 2]
	tableS2b$coef.upper[i] <- summary(mod)$tTable['treatmentnit', 1] + 1.96*summary(mod)$tTable['treatmentnit', 2]
	
	tableS2a$treatment.coef[i] <- summary(mod)$tTable['treatmenthel', 1] 
	tableS2a$coef.lower[i] <- summary(mod)$tTable['treatmenthel', 1] - 1.96*summary(mod)$tTable['treatmentnit', 2]
	tableS2a$coef.upper[i] <- summary(mod)$tTable['treatmenthel', 1] + 1.96*summary(mod)$tTable['treatmentnit', 2]
	
	tableS2c$days.captive.coef[i] <- summary(mod)$tTable['days.after', 1]
	tableS2c$coef.lower[i] <- summary(mod)$tTable['days.after', 1] - 1.96*summary(mod)$tTable['days.after', 2]
	tableS2c$coef.upper[i] <- summary(mod)$tTable['days.after', 1] + 1.96*summary(mod)$tTable['days.after', 2]
	
	tableS2b$test.stat[i] <- summary(glht(mod, linfct=mcp(treatment='Dunnet')))$test$tstat[1]
	tableS2b$p.value[i] <- summary(glht(mod, linfct=mcp(treatment='Dunnet')))$test$pvalues[1]	
	
	tableS2a$test.stat[i] <- summary(glht(mod, linfct=mcp(treatment='Dunnet')))$test$tstat[2]
	tableS2a$p.value[i] <- summary(glht(mod, linfct=mcp(treatment='Dunnet')))$test$pvalues[2]
	
	tableS2c$test.stat[i] <- summary(mod)$tTable['days.after', 4]
	tableS2c$p.value[i] <- summary(mod)$tTable['days.after', 5]
	
}

par(mfrow=c(4,4)) # check normality of model residuals
for(i in 1:12) {hist(mod.resid[[i]], main=metrics[i])}
for(i in 1:12) {print(shapiro.test(mod.resid[[i]]))}

tableS2a # effect of air density
tableS2b # effect of oxygen availability
tableS2c # effect of time in captivity

no.man2 # no. of maneuvers retained for arc cent accel only
summary(mod.finalGAS[[11]])
visreg(mod.finalGAS[[11]], 'no.maneuvers', main='residual Arc cent, max (m/s2)') # birds that performed more arcing turns had greater centripetal acceleration during these turns



# Check FDR correction
##########################################

alltests <- qvalue(c(tableS1$p.value, tableS2a$p.value, tableS2b$p.value, tableS2c$p.value), fdr.level=0.05, pi0.method="bootstrap")$significant

tableS1$test <- ifelse(alltests[1:12]==T, '*', 'ns')
tableS2a$test <- ifelse(alltests[13:24]==T, '*', 'ns')
tableS2b$test <- ifelse(alltests[25:36]==T, '*', 'ns')
tableS2c$test <- ifelse(alltests[37:48]==T, '*', 'ns')

tableS1
tableS2a
tableS2b
tableS2c

##########################################

# residual plots for time in captivity (Figure S2):
visreg(mod.finalGAS[[1]], 'days.after', main='residual Vel max (m/s)')
visreg(mod.finalGAS[[7]], 'days.after', main='residual PRT deg (º)')
visreg(mod.finalGAS[[8]], 'days.after', main='residual PRT time (s)')




