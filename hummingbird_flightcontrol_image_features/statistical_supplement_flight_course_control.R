
# Statistical Supplement
################################################################################
# Analysis and Script by: Roslyn Dakin
# Last Modified: May 14, 2016
# In support of the following study:
# Title: Visual guidance of forward flight in hummingbirds reveals control based on image features instead of pattern velocity
# Authors: Roslyn Dakin, Tyee K. Fellows, and Douglas L. Altshuler


################################################################################

library(dplyr)
library(nlme)
library(multcomp)

exp1 <- read.csv('experiment1_flightmeans.csv')
exp2 <- read.csv('experiment2_flightmeans.csv')
exp3 <- read.csv('experiment3_flightmeans.csv')
exp4 <- read.csv('experiment4_flightmeans.csv') # delete extra columns??
exp5a <- read.csv('experiment5_differentorient_flightmeans.csv') # ditto
exp5b <- read.csv('experiment5_sameorient_flightmeans.csv') # ditto
samewalls <- read.csv('experiment4-5_identicalgratings_flightmeans.csv')

# Experiment 1: Nasal-to-temporal motion of vertical gratings (Fig 2A, Table S2)
################################################################################

dim(exp1) # 297 flights

# model (Table S2)
exp1$to.from.feeder <- factor(exp1$to.from.feeder, levels=c('to', 'from'))
exp1$avg.x.velocitySTD <- scale(exp1$avg.x.velocity*100, center=T, scale=F) # abs. value of forward speed in cm/s, mean-centered
# convert to bird's frame of reference:
exp1$left.wall.speedBIRD <- ifelse(exp1$to.from.feeder=='from', -exp1$left.wall.speed, exp1$left.wall.speed) 
exp1$right.wall.speedBIRD <- ifelse(exp1$to.from.feeder=='from', -exp1$right.wall.speed, exp1$right.wall.speed)
exp1$avgyBIRD <- ifelse(exp1$to.from.feeder=='from', -exp1$avg.y.position, exp1$avg.y.position)
exp1$left.wall.speedBIRD <- factor(exp1$left.wall.speedBIRD, levels=c('0', '-0.34', '0.34'))
exp1$right.wall.speedBIRD <- factor(exp1$right.wall.speedBIRD, levels=c('0', '-0.34', '0.34'))
exp1mod <- lme(avgyBIRD*100 ~ left.wall.speedBIRD * avg.x.velocitySTD + right.wall.speedBIRD * avg.x.velocitySTD + replicate + to.from.feeder, random=~1|birdID/bird.treatment.block, data=exp1)
plot(exp1mod)
hist(residuals(exp1mod))
summary(exp1mod)$tTable # estimates in cm
head(model.matrix(exp1mod, data=exp1))
tests <- rbind('L wall back'=c(0,1,0,0,0,0,0,0,0,0,0,0), 'L wall fwd'=c(0,0,1,0,0,0,0,0,0,0,0,0), 'R wall back'=c(0,0,0,0,1,0,0,0,0,0,0,0), 'R wall fwd'=c(0,0,0,0,0,1,0,0,0,0,0,0))  
exp1tests <- glht(exp1mod, tests)
summary(exp1tests) # post-hoc contrasts, estimates in cm

# effect sizes comparing each motion treatment to the stationary control (Fig 2A)
exp1summary <- group_by(exp1, birdID, left.wall.speedBIRD, right.wall.speedBIRD)
exp1summary <- summarize(exp1summary, avgy=mean(avgyBIRD))
exp1summary <- group_by(exp1summary, birdID)
(exp1summary <- summarize(exp1summary, ST=(avgy[left.wall.speedBIRD=='0'&right.wall.speedBIRD=='0'])*100, Lpos=(avgy[left.wall.speedBIRD=='0.34'&right.wall.speedBIRD=='0'])*100, Lneg=(avgy[left.wall.speedBIRD=='-0.34'&right.wall.speedBIRD=='0'])*100, Rpos=(avgy[left.wall.speedBIRD=='0'&right.wall.speedBIRD=='0.34'])*100, Rneg=(avgy[left.wall.speedBIRD=='0'&right.wall.speedBIRD=='-0.34'])*100, effect1=(Lpos-ST), effect2=(Lneg-ST), effect3=(Rpos-ST), effect4=(Rneg-ST))) # effect sizes for individual birds, in cm.
# L fwd = 1; L back = 2; R fwd = 3; R back = 4
(exp1summary.all <- summarize(exp1summary, mean1=mean(effect1), mean2=mean(effect2), mean3=mean(effect3), mean4=mean(effect4), sd1=sd(effect1), sd2=sd(effect2), sd3=sd(effect3), sd4=sd(effect4), sample=6, lower1=mean1-1.96*sd1/sqrt(6), upper1=mean1+1.96*sd1/sqrt(6), lower2=mean2-1.96*sd2/sqrt(6), upper2=mean2+1.96*sd2/sqrt(6), lower3=mean3-1.96*sd3/sqrt(6), upper3=mean3+1.96*sd3/sqrt(6), lower4=mean4-1.96*sd4/sqrt(6), upper4=mean4+1.96*sd4/sqrt(6))) # overall effect sizes in cm and 95% confidence intervals

# average velocities in experiment 1, with stationary vertical gratings on both walls
exp1birdvel <- aggregate(subset(exp1, wall.speeds=='L 0 R 0')$avg.x.velocity, by=list(subset(exp1, wall.speeds=='L 0 R 0')$birdID), 'mean')
mean(exp1birdvel$x) # grand mean flight speed 2.0 m/s
mean(exp1birdvel$x) - 1.96*(sd(exp1birdvel$x)/sqrt(6)) 
mean(exp1birdvel$x) + 1.96*(sd(exp1birdvel$x)/sqrt(6)) # 95% CI [1.8, 2.2]


# Experiment 2: Nasal-to-temporal motion of dotfields (Fig 2B, Table S3)
################################################################################

dim(exp2) # 860 flights

# model (Table S3)
exp2$to.from.feeder <- factor(exp2$to.from.feeder, levels=c('to', 'from'))
exp2$dot.density <- factor(ifelse(exp2$dot.density==3200, 'low', ifelse(exp2$dot.density==6400, 'med', 'hi')), levels=c('low', 'med', 'hi'))
exp2$avg.x.velocitySTD <- scale(exp2$avg.x.velocity*100, center=T, scale=F) # abs. value of forward speed in cm/s, mean-centered
# convert to bird's frame of reference:
exp2$wall.diff.L.RBIRD <- ifelse(exp2$to.from.feeder=='from', -exp2$wall.diff.L.R, exp2$wall.diff.L.R)
exp2$wall.motion <- ifelse(exp2$wall.diff.L.RBIRD==-0.5144, 'Lneg', ifelse(exp2$wall.diff.L.RBIRD==0.5144, 'Lpos', 'none'))
exp2$wall.motion <- factor(exp2$wall.motion, levels=c('none', 'Lneg', 'Lpos'))
exp2$avgyBIRD <- ifelse(exp2$to.from.feeder=='from', -exp2$avg.y.position, exp2$avg.y.position)
exp2mod <- lme(avgyBIRD*100 ~ wall.motion * dot.density + replicate + wall.motion * avg.x.velocitySTD + to.from.feeder, random=~1|birdID/bird.treatment.block, data=exp2)
plot(exp2mod) # check residuals
hist(residuals(exp2mod))
summary(exp2mod)$tTable # estimates in cm
head(model.matrix(exp2mod, data=exp2))
tests <- rbind('low.den L back'=c(0,1,rep(0,12)), 'low.den L fwd'=c(0,0,1,rep(0,11)), 'med.den L back'=c(0,1,rep(0,6),1,rep(0,5)), 'med.den L fwd'=c(0,0,1,rep(0,6),1,rep(0,4)), 'high.den L back'=c(0,1,rep(0,8),1,rep(0,3)), 'high.den L fwd'=c(0,0,1,rep(0,8),1,rep(0,2)) )
exp2tests <- glht(exp2mod, tests)
summary(exp2tests) # post-hoc contrasts, estimates in cm

# effect sizes comparing each motion treatment to the stationary control (Fig 2B)
exp2summary <- group_by(exp2, birdID, dot.density, wall.motion)
exp2summary <- summarize(exp2summary, avgy=mean(avgyBIRD))
(exp2summary <- summarize(exp2summary, ST=(avgy[wall.motion=='none'])*100, Lpos=(avgy[wall.motion=='Lpos'])*100, Lneg=(avgy[wall.motion=='Lneg'])*100, effect1=(Lpos-ST), effect2=(Lneg-ST))) # effect sizes for individual birds, in cm
# L fwd, R back = 1; L back, R fwd = 2
exp2summary <- group_by(exp2summary, dot.density)
(exp2summary.all <- summarize(exp2summary, mean1=mean(effect1), mean2=mean(effect2), sd1=sd(effect1), sd2=sd(effect2), sample=10, lower1=mean1-1.96*sd1/sqrt(10), upper1=mean1+1.96*sd1/sqrt(10), lower2=mean2-1.96*sd2/sqrt(10), upper2=mean2+1.96*sd2/sqrt(10))) # overall effect sizes in cm and 95% confidence intervals

# average velocities in experiment 2, with stationary dotfields on both walls
exp2birdvel <- aggregate(subset(exp2, wall.motion=='none')$avg.x.velocity, by=list(subset(exp2, wall.motion=='none')$birdID), 'mean')
mean(exp2birdvel$x) # grand mean flight speed 2.2 m/s
mean(exp2birdvel$x) - 1.96*(sd(exp2birdvel$x)/sqrt(10)) 
mean(exp2birdvel$x) + 1.96*(sd(exp2birdvel$x)/sqrt(10)) # 95% CI [1.9, 2.5]


# Experiment 3: Up/downward motion of horizontal gratings (Fig 2C, Table S4)
################################################################################

dim(exp3) # 148 flights

# model (Table S4)
exp3$to.from.feeder <- factor(exp3$to.from.feeder, levels=c('to', 'from'))
exp3$grating.up.down.speed <- factor(exp3$grating.up.down.speed, levels=c('0', '-0.34', '0.34'))
exp3$avg.x.velocitySTD <- scale(abs(exp3$avg.x.velocity)*100, center=T, scale=F) # abs. value of forward speed in cm/s, mean-centered
exp3mod <- lme(avg.z.position*100 ~ to.from.feeder + grating.up.down.speed * avg.x.velocitySTD + replicate, random=~1|birdID/bird.treatment.block, data=exp3)
plot(exp3mod) # treatments have unequal variances, therefore adjust model variance structure to acccount f or this
exp3mod <- lme(avg.z.position*100 ~ to.from.feeder + grating.up.down.speed * avg.x.velocitySTD + replicate, random=~1|birdID/bird.treatment.block, weights=varIdent(form=~1|grating.up.down.speed), data=exp3)
plot(exp3mod) # better
summary(exp3mod)$tTable # estimates in cm
head(model.matrix(exp3mod, data=exp3))
tests <- rbind('walls down'=c(0,0,1,0,0,0,0,0), 'walls up'=c(0,0,0,1,0,0,0,0))  
exp3tests <- glht(exp3mod, tests)
summary(exp3tests) # post-hoc contrasts, estimates in cm

# effect sizes comparing each motion treatment to the stationary control (Fig 2C)
exp3summary <- group_by(exp3, birdID, bird.treatment.block, grating.up.down.speed)
exp3summary <- summarize(exp3summary, avgz=mean(avg.z.position))
exp3summary <- group_by(exp3summary, birdID)
(exp3summary <- summarize(exp3summary, ST=(avgz[grating.up.down.speed=='0'])*100, UP=(avgz[grating.up.down.speed=='0.34'])*100, DN=(avgz[grating.up.down.speed=='-0.34'])*100, effect1=(UP-ST), effect2=(DN-ST))) # effect sizes for individual birds, in cm
# walls up = 1; walls down = 2
(exp3summary.all <- summarize(exp3summary, mean1=mean(effect1), mean2=mean(effect2), sd1=sd(effect1), sd2=sd(effect2), sample=5, lower1=mean1-1.96*sd1/sqrt(5), upper1=mean1+1.96*sd1/sqrt(5), lower2=mean2-1.96*sd2/sqrt(5), upper2=mean2+1.96*sd2/sqrt(5))) # overall effect sizes in cm and 95% confidence intervals


# Experiment 4: Effect of grating orientation (Fig 3, Table S5)
################################################################################

dim(exp4) # 1180 flights with different-orientation gratings on each wall

# check whether color of large top stripe matters
checktop <- lme(avg.y.position*100 ~ left.grating.orientation * top.black.red + replicate + to.from.feeder, random=~1|birdID/bird.treatment.block, data=subset(exp4, grating.period==18.4))
test <- rbind(interaction=c(rep(0,5),1))
checktoptest <- glht(checktop, test)
summary(checktoptest) # interaction is not significant, indicating steering was not sig. affected by top stripe color

# model (Table S5)
exp4$grating.period <- factor(exp4$grating.period)
exp4$to.from.feeder <- factor(exp4$to.from.feeder, levels=c('to', 'from'))
exp4mod <- lme(avg.y.position*100 ~ left.grating.orientation * grating.period + replicate + to.from.feeder, data=exp4, random=~1|birdID/bird.treatment.block)
plot(exp4mod)
hist(residuals(exp4mod)) # check residuals
summary(exp4mod)$tTable # estimates in cm
head(model.matrix(exp4mod, data=exp4))
tests <- rbind(orient0.58=c(0,1,rep(0,12)), orient1.15=c(0,1,rep(0,7),1,rep(0,4)), orient2.3=c(0,1,rep(0,8),1,rep(0,3)), orient4.6=c(0,1,rep(0,9),1,rep(0,2)), orient9.2=c(0,1,rep(0,10),1,rep(0,1)), orient18.4=c(0,1,rep(0,11),1))
exp4tests <- glht(exp4mod, tests)
summary(exp4tests) # post-hoc contrasts, estimates in cm

# effect sizes comparing mirrored treatments (Fig 3)
exp4summary <- group_by(exp4, birdID, bird.treatment.block, left.grating.orientation, grating.period)
exp4summary <- summarize(exp4summary, fusion=mean(percent.vstripes.fused), avgy=mean(avg.y.position))
exp4summary <- group_by(exp4summary, birdID, grating.period)
(exp4summary <- summarize(exp4summary, LH=(avgy[left.grating.orientation=='hor'])*100, LV=(avgy[left.grating.orientation=='ver'])*100, fusion=mean(fusion), effect.size=(LH-LV))) # effect sizes for individual birds, in cm
exp4summary.all <- group_by(exp4summary, grating.period)
summarize(exp4summary.all, effect=mean(effect.size), effectSD=sd(effect.size), sample.size=10, lower=effect-1.96*effectSD/sqrt(10), upper=effect+1.96*effectSD/sqrt(10)) # overall effect sizes in cm and 95% confidence intervals

# percent of vertical bars fused (Fig 3)
summarize(exp4summary.all, grand.mean.fusion=round(100*mean(fusion),1)) # grand mean % of vertical bars fused


# Experiment 5: Effect of horizontal feature size (Fig 4, Tables S6, S7, S10)
################################################################################

# treatments with different grating orientations (Table S6)
dim(exp5a) # 474 flights with different-orientation gratings on each wall

# check whether color of large top stripe matters
checktop <- lme(avg.y.position*100 ~ top.black.red * left.grating.orientation + V.grating.period + replicate + to.from.feeder, random=~1|birdID/bird.treatment.block, data=subset(exp5a, H.grating.period==18.4))
test <- rbind('top by size'=c(rep(0,6),1))
checktoptest <- glht(checktop, test)
summary(checktoptest) # interaction is not significant, indicating steering was not sig. affected by top stripe color

# model (Table S6)
exp5a$H.grating.period <- factor(exp5a$H.grating.period)
exp5a$V.grating.period <- factor(exp5a$V.grating.period)
exp5a$to.from.feeder <- factor(exp5a$to.from.feeder, levels=c('to', 'from'))
exp5amod1 <- lme(avg.y.position*100 ~ left.grating.orientation * walls.stripe.size + replicate + to.from.feeder, random=~1|birdID/bird.treatment.block, data=exp5a)
plot(exp5amod1)
hist(residuals(exp5amod1)) # check residuals
summary(exp5amod1)$tTable # estimates in cm
head(model.matrix(exp5amod1, data=exp5a))
tests <- rbind('Hsm Vlg'=c(0,1,rep(0,3),0,0,1,0,0), 'Hsm Vsm'=c(0,1,rep(0,3),0,0,0,0,0), 'Hlg Vlg'=c(0,1,rep(0,3),0,0,0,0,1), 'Hlg Vsm'=c(0,1,rep(0,3),0,0,0,1,0))  
exp5atests <- glht(exp5amod1, tests)
summary(exp5atests) # post-hoc contrasts, in cm

# additional contrasts to evaluate the effect of horizontal grating period sizes, and vertical grating period size, respectively
exp5amod2 <- lme(avg.y.position*100 ~ left.grating.orientation*H.grating.period + left.grating.orientation*V.grating.period + replicate + to.from.feeder, random=~1|birdID/bird.treatment.block, data=exp5a)
summary(exp5amod2)$tTable
plot(exp5amod2)
hist(residuals(exp5amod2)) # check residuals
head(model.matrix(exp5amod2, data=exp5a))
tests <- rbind('orient:hsize'=c(0,0,0,0,0,0,1,0), 'orient:vsize'=c(0,0,0,0,0,0,0,1) )  
exp5atests2 <- glht(exp5amod2, tests)
summary(exp5atests2) # post-hoc contrasts, in cm

# treatments with the same grating orientation (Table S7)
dim(exp5b) # 236 flights with same-orientation gratings on both walls

# model (Table S7)
exp5b$to.from.feeder <- factor(exp5b$to.from.feeder, levels=c('to', 'from'))
exp5bmod <- lme(avg.y.position*100 ~ side.small.period * left.grating.orientation + replicate + to.from.feeder, random=~1|birdID/bird.treatment.block, data=exp5b)
plot(exp5bmod)
hist(residuals(exp5bmod)) # check residuals
summary(exp5bmod)$tTable # estimates in cm
head(model.matrix(exp5bmod, data=exp5b))
tests <- rbind('small wall IF walls H'=c(0,1,0,0,0,0), 'small wall IF walls V'=c(0,1,0,0,0,1), 'small wall H vs small wall V'=c(rep(0,5),1))  
exp5btests <- glht(exp5bmod, tests)
summary(exp5btests) # post-hoc contrasts, in cm

# effect sizes comparing mirrored treatments (Fig 4)
exp5asummary <- group_by(exp5a, birdID, bird.treatment.block, left.grating.orientation, walls.stripe.size)
exp5asummary <- summarize(exp5asummary, avgy=mean(avg.y.position))
exp5asummary <- group_by(exp5asummary, birdID, walls.stripe.size)
(exp5asummary <- summarize(exp5asummary, LH=(avgy[left.grating.orientation=='hor'])*100, LV=(avgy[left.grating.orientation=='ver'])*100, effect.size=(LH-LV))) # effect sizes for individual birds, in cm
exp5asummary.all <- group_by(exp5asummary, walls.stripe.size)
summarize(exp5asummary.all, effect=mean(effect.size), effectSD=sd(effect.size), sample.size=6, lower=effect-1.96*effectSD/sqrt(6), upper=effect+1.96*effectSD/sqrt(6)) # overall effect sizes in cm and 95% confidence intervals

exp5bsummary <- group_by(exp5b, birdID, bird.treatment.block, left.grating.orientation, side.small.period)
exp5bsummary <- summarize(exp5bsummary, avgy=mean(avg.y.position))
exp5bsummary <- group_by(exp5bsummary, birdID, left.grating.orientation)
(exp5bsummary <- summarize(exp5bsummary, Lsmall=(avgy[side.small.period=='L'])*100, Rsmall=(avgy[side.small.period=='R'])*100, effect.size=(Lsmall-Rsmall))) # effect sizes for individual birds, in cm
exp5bsummary.all <- group_by(exp5bsummary, left.grating.orientation)
summarize(exp5bsummary.all, effect=mean(effect.size), effectSD=sd(effect.size), sample.size=6, lower=effect-1.96*effectSD/sqrt(6), upper=effect+1.96*effectSD/sqrt(6)) # overall effect sizes in cm and 95% confidence intervals


# Additional treatments with identical gratings on the left and right walls (Tables S8, S9)
################################################################################

dim(samewalls) # 628 flights with the same grating on the left and right
samewalls$to.from.feeder <- factor(samewalls$to.from.feeder, levels=c('to', 'from'))

# descriptive statistics (Table S8)
samewalls.summary <- group_by(samewalls, birdID, treatment.name)
samewalls.summary <- summarize(samewalls.summary, avgy=mean(avg.y.position), avgabsy=mean(avg.abs.y.position)) # bird means
(samewalls.summary.all <- group_by(samewalls.summary, treatment.name))
(samewalls.summary.all <- summarize(samewalls.summary.all, mean.y=mean(avgy)*100, min.y=min(avgy)*100, max.y=max(avgy)*100, mean.abs.y=mean(avgabsy)*100, min.abs.y=min(avgabsy)*100, max.abs.y=max(avgabsy)*100)) # descriptive statistics, in cm
table(samewalls$treatment.name) # sample sizes

# analyses (Table S9)
samewallsmod <- lme(avg.abs.y.position*100 ~ treatment.name + replicate + to.from.feeder + avg.x.velocity, random=~1|birdID/bird.treatment.block, data=samewalls)
plot(samewallsmod)
hist(residuals(samewallsmod))
summary(samewallsmod)
samewallstest <- glht(samewallsmod, mcp(treatment.name='Tukey'))
summary(samewallstest)


# Pattern velocity with stationary vertical gratings (Table S1)
################################################################################

samewallsV <- subset(samewalls, grating.orientation=='ver')
table(samewallsV$grating.period)
samewallsV <- group_by(samewallsV, birdID, grating.period)
samewalls.summary <- summarize(samewallsV, L=mean(left.eye.pattern.velocity), R=mean(right.eye.pattern.velocity))
samewalls.summary <- group_by(samewalls.summary, grating.period) # individual bird means
summarize(samewalls.summary, L.mean=mean(L), L.lower=L.mean-1.96*(sd(L)/sqrt(16)), L.upper=L.mean+1.96*(sd(L)/sqrt(16)), R.mean=mean(R), R.lower=R.mean-1.96*(sd(R)/sqrt(16)), R.upper=R.mean+1.96*(sd(R)/sqrt(16))) # grand means and 95% CI








