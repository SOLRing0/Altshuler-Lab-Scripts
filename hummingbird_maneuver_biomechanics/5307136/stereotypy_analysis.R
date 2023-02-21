

# Stereotypic behaviour (in the sense of stereotypy resulting from captivity) is defined as a repetitive, invariant pattern of behavior.
# Examples include: pacing in cats, feather pecking in chickens, grooming in rodents, swimming in circles in fish
# Are the hummingbird behaviors repetitive and invariant?

library(stringr)
library(dplyr)
library(tidyverse)
library(nlme)

mdata4 <- read.csv('mdata4.csv')[,-1]
head(mdata4)

mdata4$id <- paste(mdata4$species, ifelse(nchar(mdata4$id)==1, paste('00', mdata4$id,sep=''), paste('0', mdata4$id,sep='')), sep='')
mdata4$mid.frame <- rowMeans(mdata4[,c('framee','frames')])
mdata4$id <- factor(mdata4$id)
mdata4$video <- factor(mdata4$video)
mdata4$cat <- factor(mdata4$cat)

vid.table <- table(factor(mdata4$id), mdata4$video)!=0
rownames(vid.table)[rowSums(vid.table)>1] # these 7 individuals have multi-part videos

mdata4 <- mdata4[with(mdata4, order(id, video, frames)),]
head(mdata4)
mdataX <- mdata4
mdataX <- subset(mdataX, cat %in% c('total_vel','hor_acel','hor_decel','pitch_up','pitch_down','yaw','arc_force','PRT_time'))
dim(mdataX) # 333,117 limits to one metric per maneuver, for transition probability analyses
mdataX$cat <- factor(mdataX$cat, levels=c("total_vel","hor_acel","hor_decel","pitch_up","pitch_down","yaw","arc_force","PRT_time"))  
mdata4$cat <- factor(mdata4$cat, levels=c("total_vel","hor_acel","hor_decel","pitch_up","pitch_down","yaw","arc_radius","arc_avg_vel","arc_force","PRT_degrees","PRT_time"))

# 1. How variable is the performance of maneuvers within a trial?
# If the behaviors are fixed (invariant), we expect low variation within trials.
# Quantify variation within trials, for each metric, as: (i) variance component %, attributed to variation within individual trials, and (ii) coefficient of variation (scale-free)

(varresults <- data.frame(metric=levels(mdata4$cat), percentvaramong=NA, percentvarwithin=NA, CV=NA))

coefv <- function(x){ # coefficient of variation
	sd(x,na.rm=T)/mean(abs(x),na.rm=T)
}

for(i in 1:11){
	mod <- lme(val ~ 1, random=~1|id, subset(mdata4, cat==levels(mdata4$cat)[i]))
	varresults$percentvarwithin[i] <- as.numeric(VarCorr(mod)[,1])[2]/sum(as.numeric(VarCorr(mod)[,1]))*100
	varresults$percentvaramong[i] <- as.numeric(VarCorr(mod)[,1])[1]/sum(as.numeric(VarCorr(mod)[,1]))*100
	cvs <- rep(NA, 207)
	for(j in 1:207){
		temp <- subset(mdata4, cat==levels(mdata4$cat)[i] & id==levels(mdata4$id)[j])
		cvs[j] <- coefv(temp$val)
	}
	varresults$CV[i] <- mean(cvs, na.rm=T)
	print(i)
}
varresults
summary(varresults) # CVs (a scale-free measure of within-trial dispersion) for each maneuver range from 0.25-0.68 on average, and 70-89% of the variation in total (raw) dataset can be attributed to within-trial variation

# 2. Transition probabilities: the frequency with which a given maneuver is followed by another specific maneuver (defined as being within a bout of flight; 1 second interval)
# If the behavior patterns are fixed (invariant), we expect very high transition probabilities (i.e., maneuver X is always followed by maneuver Y for individual A, etc.)
# Quantify transition probabilities for each bird-trial-metric and examine the max and the # of option for each metric with probability >10% of occurring after that maneuver type

dim(mdataX)
mdataX$frames_next <- lead(mdataX$frames)
mdataX$cat_next <- lead(mdataX$cat)
mdataX$tdiff <- (mdataX$frames_next - mdataX$frames)/200 # time difference in seconds
mdataX$cat_next[(mdataX$id!=lead(mdataX$id)) | (mdataX$video!=lead(mdataX$video)) | (mdataX$tdiff > 1) ] <-NA
mdataX$tdiff[(mdataX$id!=lead(mdataX$id)) | (mdataX$video!=lead(mdataX$video)) | (mdataX$tdiff > 1) ] <-NA
head(mdataX)
hist(mdataX$tdiff, breaks=1000)

# maybe max transition probabilities?
tprob <- data.frame(matrix(NA, nrow=207, ncol=8))
t10 <- data.frame(matrix(NA, nrow=207, ncol=8))
length(levels(factor(mdataX$id)))
for(i in 1:207){
	mytab <- table(subset(mdataX, id==levels(factor(mdataX$id))[i])[,c('cat_next','cat')])
	mytab <- scale(mytab, center = FALSE, scale = colSums(mytab))
	tprob[i,] <- apply(mytab, 2, 'max')
	t10[i,] <- colSums(mytab>0.1)
	print(i)
}
names(tprob) <- names(apply(mytab, 2, 'max'))
names(t10) <- names(tprob)
head(tprob)
head(t10)
range(colMeans(t10, na.rm=T), na.rm=T) # on average, each maneuver has 2.3-3.9 other maneuvers with transition probabilities >10%
range(colMeans(tprob, na.rm=T), na.rm=T) # the maximum transition probability for a maneuver ranges between 35-58%

varresults$nchoices <- NA
varresults$maxTP <- NA

varresults$nchoices[c(1:6,7,10)] <- as.numeric(colMeans(t10, na.rm=T))
varresults$maxTP[c(1:6,7,10)] <- as.numeric(colMeans(tprob, na.rm=T))














