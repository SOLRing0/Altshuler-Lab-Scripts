
# Evolution reveals the biomechanical organization of maneuvering flight in hummingbirds
# Dakin, Segre, Straw and Altshuler
# Figure 4, Figure S1

library(tidyverse)
library(nlme)
library(ape)

load("comparative_maneuverability_models.RData") # load the models RData file

means <- aggregate(maneuv[,c(16:27,36,39,40,50)], by=list(species=maneuv$species), 'mean', na.rm=T)
sds <- aggregate(maneuv[,c(16:27)], by=list(species=maneuv$species), 'sd', na.rm=T)
names(sds)[2:13] <- paste('sd', names(sds)[2:13], sep='.')
ns <- summarize(group_by(maneuv, species), ss=n())

# ensure that the working directory has "ll_morphology_20170517.csv"

# set up the data:
LLsample <- read.csv('ll_morphology_20170517.csv')
head(LLsample); summary(LLsample); dim(LLsample)
mass.sds <- aggregate(LLsample$mass, by=list(LLsample$species), 'sd', na.rm=T)
names(mass.sds)[2] <- 'sd.mass'
mass.ss <- summarize(group_by(LLsample, species), ss.mass=n())

forscaleplot <- cbind(means, sds[,2:13], ns[,2], mass.sds$sd.mass, mass.ss$ss.mass)
names(forscaleplot)[30:31] <- c('sd.mass', 'ss.mass')

forscaleplot$hor_decel <- -forscaleplot$hor_decel
forscaleplot$pitch_down <- -forscaleplot$pitch_down
forscaleplot$PRT_time <- -forscaleplot$PRT_time

LLsample$wing.load <- LLsample$mass/LLsample$wing.area
allom <- summarize(group_by(LLsample, species), mass=mean(mass, na.rm=T), WL=mean(wing.load, na.rm=T), area=mean(wing.area, na.rm=T))

# re-fit models to raw (unstandardized) data for plots:
smod.listraw <- list()
predictlist <- list()
smod.listraw[[1]] <- lme(total_vel ~ hi.elev + sp.mass + rel.mass + ss.V, random=~1|species, data=maneuv, na.action=na.omit, control=lmeControl(opt='optim'))
predictlist[[1]] <- predict(smod.listraw[[1]], newdata=data.frame(hi.elev=0.5, rel.mass=0, ss.V=mean(maneuv$ss.V), sp.mass=seq(2, 9.5, by=0.01)), level=0)
smod.listraw[[2]] <- lme(hor_acel ~ hi.elev + sp.mass + rel.mass + ss.HA, random=~1|species, data=maneuv, na.action=na.omit, control=lmeControl(opt='optim'))
predictlist[[2]] <- predict(smod.listraw[[2]], newdata=data.frame(hi.elev=0.5, rel.mass=0, ss.HA=mean(maneuv$ss.HA), sp.mass=seq(2, 9.5, by=0.01)), level=0)
smod.listraw[[3]] <- lme(-hor_decel ~ hi.elev + sp.mass + rel.mass + ss.HD, random=~1|species, data=maneuv, na.action=na.omit, control=lmeControl(opt='optim'))
predictlist[[3]] <- predict(smod.listraw[[3]], newdata=data.frame(hi.elev=0.5, rel.mass=0, ss.HD=mean(maneuv$ss.HD), sp.mass=seq(2, 9.5, by=0.01)), level=0)
smod.listraw[[4]] <- lme(pitch_up ~ hi.elev + sp.mass + rel.mass + ss.PU, random=~1|species, data=maneuv, na.action=na.omit, control=lmeControl(opt='optim'))
predictlist[[4]] <- predict(smod.listraw[[4]], newdata=data.frame(hi.elev=0.5, rel.mass=0, ss.PU=mean(maneuv$ss.PU), sp.mass=seq(2, 9.5, by=0.01)), level=0)
smod.listraw[[5]] <- lme(-pitch_down ~ hi.elev + sp.mass + rel.mass + ss.PD, random=~1|species, data=maneuv, na.action=na.omit, control=lmeControl(opt='optim'))
predictlist[[5]] <- predict(smod.listraw[[5]], newdata=data.frame(hi.elev=0.5, rel.mass=0, ss.PD=mean(maneuv$ss.PD), sp.mass=seq(2, 9.5, by=0.01)), level=0)
smod.listraw[[6]] <- lme(yaw ~ hi.elev + sp.mass + rel.mass + ss.Y, random=~1|species, data=maneuv, na.action=na.omit, control=lmeControl(opt='optim'))
predictlist[[6]] <- predict(smod.listraw[[6]], newdata=data.frame(hi.elev=0.5, rel.mass=0, ss.Y=mean(maneuv$ss.Y), sp.mass=seq(2, 9.5, by=0.01)), level=0)
smod.listraw[[7]] <- lme(arc_radius ~ hi.elev + sp.mass + rel.mass + ss.ARC, random=~1|species, data=maneuv, na.action=na.omit, control=lmeControl(opt='optim'))
predictlist[[7]] <- predict(smod.listraw[[7]], newdata=data.frame(hi.elev=0.5, rel.mass=0, ss.ARC=mean(maneuv$ss.ARC), sp.mass=seq(2, 9.5, by=0.01)), level=0)
smod.listraw[[8]] <- lme(arc_avg_vel ~ hi.elev + sp.mass + rel.mass + ss.ARC, random=~1|species, data=maneuv, na.action=na.omit, control=lmeControl(opt='optim'))
predictlist[[8]] <- predict(smod.listraw[[8]], newdata=data.frame(hi.elev=0.5, rel.mass=0, ss.ARC=mean(maneuv$ss.ARC), sp.mass=seq(2, 9.5, by=0.01)), level=0)
smod.listraw[[9]] <- lme(arc_force ~ hi.elev + sp.mass + rel.mass + ss.ARC, random=~1|species, data=maneuv, na.action=na.omit, control=lmeControl(opt='optim'))
predictlist[[9]] <- predict(smod.listraw[[9]], newdata=data.frame(hi.elev=0.5, rel.mass=0, ss.ARC=mean(maneuv$ss.ARC), sp.mass=seq(2, 9.5, by=0.01)), level=0)
smod.listraw[[10]] <- lme(PRT_degrees ~ hi.elev + sp.mass + rel.mass + ss.PRT, random=~1|species, data=maneuv, na.action=na.omit, control=lmeControl(opt='optim'))
predictlist[[10]] <- predict(smod.listraw[[10]], newdata=data.frame(hi.elev=0.5, rel.mass=0, ss.PRT=mean(maneuv$ss.PRT), sp.mass=seq(2, 9.5, by=0.01)), level=0)
smod.listraw[[11]] <- lme(-PRT_time ~ hi.elev + sp.mass + rel.mass + ss.PRT, random=~1|species, data=maneuv, na.action=na.omit, control=lmeControl(opt='optim'))
predictlist[[11]] <- predict(smod.listraw[[11]], newdata=data.frame(hi.elev=0.5, rel.mass=0, ss.PRT=mean(maneuv$ss.PRT), sp.mass=seq(2, 9.5, by=0.01)), level=0)
smod.listraw[[12]] <- lme(pPRT ~ hi.elev + sp.mass + rel.mass, random=~1|species, data=maneuv, na.action=na.omit, control=lmeControl(opt='optim'))
predictlist[[12]] <- predict(smod.listraw[[12]], newdata=data.frame(hi.elev=0.5, rel.mass=0, sp.mass=seq(2, 9.5, by=0.01)), level=0)

ymins <- c(1.3, 3, 3, 0.8, 0.7, 1.5, 0.3, 0.8, 3, 65, -0.54, 0.32)
ymax <- c(2.5, 8, 8, 1.6, 1.5, 2.3, 0.8, 2, 8, 130, -0.385, 1)

# Figure 4 among-species scaling:

dev.new(width=0.7, height=6); par(mfrow=c(8,1), mar=c(2,2,2,0.25), bty='l', las=1, mgp=c(1.25,0.5,0), cex=0.5)
for(i in c(10,2:4,5,7,6,13)){ #13
	plot(forscaleplot[,i] ~ forscaleplot$sp.mass, pch=plottable$pch, col=as.character(plottable$color), main=names(forscaleplot)[i], xlab='', ylab='', xlim=c(2.5,9.5), ylim=c(ymins[i-1], ymax[i-1]), xaxt='n', yaxt='n')
	axis(1, at=c(3,9))
	axis(2, at=c(ymins[i-1], ymax[i-1]))
	points(predictlist[[i-1]] ~ seq(2, 9.5, by=0.01), type='l', lty=3, lwd=0.94)
}

mynorm <- function(x){
	if(any(x<0)){
		return(  (x + abs(min(x)))/max((x + abs(min(x))))  )
	} else {
		return(  (x - abs(min(x)))/max((x - abs(min(x))))  )
	}
}

dev.new(width=4, height=5.5)
par(mfrow=c(3,2), bty='l', las=1, mar=c(1.5,1.5,2.5,0.25), mgp=c(1.25,0.45,0), tck=-0.035, cex.lab=0.75, cex.axis=0.5)
i = 3
plot(forscaleplot[,i] ~ forscaleplot$sp.mass, pch=16, col=rgb(0,0,0,mynorm(forscaleplot$sp.muscle)), main=names(forscaleplot)[i], xlab='', ylab='', xlim=c(2.5,9.5), ylim=c(ymins[i-1], ymax[i-1]))
text(x=forscaleplot$sp.mass, y=forscaleplot$hor_acel, labels= round(forscaleplot$sp.muscle, 1), cex=0.5, adj=c(1,2))
points(predictlist[[i-1]] ~ seq(2, 9.5, by=0.01), type='l', lty=3, lwd=0.94)
i = 4
plot(forscaleplot[,i] ~ forscaleplot$sp.mass, pch=16, col=rgb(0,0,0,mynorm(forscaleplot$sp.muscle)), main=names(forscaleplot)[i], xlab='', ylab='', xlim=c(2.5,9.5), ylim=c(ymins[i-1], ymax[i-1]))
text(x=forscaleplot$sp.mass, y=forscaleplot$hor_decel, labels= round(forscaleplot$sp.muscle, 1), cex=0.5, adj=c(1,2))
points(predictlist[[i-1]] ~ seq(2, 9.5, by=0.01), type='l', lty=3, lwd=0.94)
i = 5
plot(forscaleplot[,i] ~ forscaleplot$sp.mass, pch=16, col=rgb(0,0,0,mynorm(forscaleplot$sp.muscle)), main=names(forscaleplot)[i], xlab='', ylab='', xlim=c(2.5,9.5), ylim=c(ymins[i-1], ymax[i-1]))
text(x=forscaleplot$sp.mass, y=forscaleplot$pitch_up, labels= round(forscaleplot$sp.muscle, 1), cex=0.5, adj=c(1,2))
points(predictlist[[i-1]] ~ seq(2, 9.5, by=0.01), type='l', lty=3, lwd=0.94)
i = 10
plot(forscaleplot[,i] ~ forscaleplot$sp.mass, pch=16, col=rgb(0,0,0,mynorm(forscaleplot$sp.muscle)), main=names(forscaleplot)[i], xlab='', ylab='', xlim=c(2.5,9.5), ylim=c(ymins[i-1], ymax[i-1]))
text(x=forscaleplot$sp.mass, y=forscaleplot$arc_force, labels= round(forscaleplot$sp.muscle, 1), cex=0.5, adj=c(1,2))
points(predictlist[[i-1]] ~ seq(2, 9.5, by=0.01), type='l', lty=3, lwd=0.94)
i = 6
plot(forscaleplot[,i] ~ forscaleplot$sp.mass, pch=16, col=rgb(0,0,0,mynorm(allom$WL)), main=names(forscaleplot)[i], xlab='', ylab='', xlim=c(2.5,9.5), ylim=c(ymins[i-1], ymax[i-1]))
text(x=forscaleplot$sp.mass, y=forscaleplot$pitch_down, labels= round(allom$WL, 2), cex=0.5, adj=c(1,2))
points(predictlist[[i-1]] ~ seq(2, 9.5, by=0.01), type='l', lty=3, lwd=0.94)
i = 7
plot(forscaleplot[,i] ~ forscaleplot$sp.mass, pch=16, col=rgb(0,0,0,mynorm(allom$WL)), main=names(forscaleplot)[i], xlab='', ylab='', xlim=c(2.5,9.5), ylim=c(ymins[i-1], ymax[i-1]))
text(x=forscaleplot$sp.mass, y=forscaleplot$yaw, labels= round(allom$WL, 2), cex=0.5, adj=c(1,2))
points(predictlist[[i-1]] ~ seq(2, 9.5, by=0.01), type='l', lty=3, lwd=0.94)


# Figure S1 wing area and wing loading:

dev.new(width=4, height=3.75)
par(mfrow=c(2,2), bty='l', las=1, mar=c(2.5,2.5,0.25,0.25), mgp=c(1.25,0.45,0), tck=-0.025, cex.lab=0.75, cex.axis=0.5)
plot(allom$area ~ allom$mass, xlim=c(2,10), ylim=c(8,40), xlab='sp body mass (g)', ylab='sp wing area (cm2)', col=as.character(plottable$color), pch=plottable$pch)
abline(lm(allom$area ~ allom$mass), lty=1, lwd=1.3)
plot(log(allom$area, 10) ~ log(allom$mass, 10), xlim=log(c(2,10), 10), ylim=log(c(8,40), 10), xaxt='n', yaxt='n', xlab='log body mass (g)', ylab='log wing area (cm2)', col=as.character(plottable$color), pch=plottable$pch)
axis(1, at=log(c(2,4,6,8,10), 10), labels=c(2,4,6,8,10))
axis(2, at=log(c(10,20,30,40), 10), labels=c(10,20,30,40))
abline(lm(log(allom$area, 10) ~ log(allom$mass, 10)), lty=1, lwd=1.3)
plot(allom$WL ~ allom$mass, xlim=c(2,10), ylim=c(0.15,0.35), xlab='sp body mass (g)', ylab='sp wing load (g/cm2)', col=as.character(plottable$color), pch=plottable$pch)
abline(lm(allom$WL ~ allom$mass), lty=3)
plot(log(allom$WL, 10) ~ log(allom$mass, 10), xaxt='n', yaxt='n', xlim=log(c(2,10), 10), ylim=log(c(0.15,0.35), 10), xlab='log body mass (g)', ylab='log wing load (g/cm2)', col=as.character(plottable$color), pch=plottable$pch)
axis(1, at=log(c(2,4,6,8,10), 10), labels=c(2,4,6,8,10))
axis(2, at=log(c(0.15,0.2,0.25,0.3,0.35), 10), labels=c(0.15,0.20,0.25,0.30,0.35))
abline(lm(log(allom$WL, 10) ~ log(allom$mass, 10)), lty=3)
summary(lm(log(allom$WL, 10) ~ log(allom$mass, 10)))










