
# Evolution reveals the biomechanical organization of maneuvering flight in hummingbirds
# Dakin, Segre, Straw and Altshuler
# data, Figure 1, Figure S1, and Figure S2

library(nlme)
library(visreg)
library(dplyr)
library(ape)
library(phytools)

# get load-lifting data and derive sp. muscle capacity:
LLsample <- read.csv('ll_morphology_20170517.csv')
head(LLsample); summary(LLsample); dim(LLsample) # 263 birds in 25 species
# the 18 individuals with maneuv.id have separate measurements from the free flight assay of maneuverability

# wing loading in g/cm2
LLsample$wing.load <- LLsample$mass/LLsample$wing.area; hist(LLsample$wing.load)

# make 'Male' the default level for sex
LLsample$Sex <- factor(LLsample$Sex, levels=c('Male','Unknown','Female'))

# derive sp. muscle capacity using load-lifting, adjusted for body mass, elevation of the test site, and sex
mod1 <- lme(lifted.beads ~ elev + mass + Sex, random=~Sex|species, data=LLsample, na.action=na.omit, control=lmeControl(opt='optim'))
summary(mod1)
plot(mod1)
par(mfrow=c(3,3)); visreg(mod1); hist(residuals(mod1))
ranef(mod1) # use the random intercept (for males) as the adjusted sp. value of muscle capacity 
mean(ranef(mod1)[,1]) # it is centered on 0

# species average traits, raw and standardized:
LLsampleg <- group_by(LLsample, species)
sptraits <- cbind(data.frame(summarize(LLsampleg, mass=mean(mass, na.rm=T), wing.area=mean(wing.area, na.rm=T), AR=mean(AR, na.rm=T), wing.load=mean(wing.load, na.rm=T))), ranef(mod1)[,1])
names(sptraits)[6] <- 'muscle'
names(sptraits)[2:6] <- paste('sp', names(sptraits)[2:6], sep='.')
for(i in 2:6){
	sptraits[,paste(names(sptraits)[i], 'STD', sep='.')] <- as.numeric(scale(sptraits[,i]))
}
summary(sptraits) # STD values = centered on 0 and scaled to have SD of 1

# get maneuverability assay of individual birds:
maneuv <- read.csv('SA_maneuver_data_170724.csv')
head(maneuv); dim(maneuv) # 213 individuals
table(maneuv$species, maneuv$sex) # only 1 sp. have variable sex, S. flammula only one where females used
maneuv <- subset(maneuv, !(sex=='F'))
dim(maneuv) # 207 individuals remaining

# add species traits & calculate individual values relative to the raw sp. average:
maneuv <- merge(maneuv, sptraits, by='species')
for(i in 12:15){
	maneuv[, paste('rel', names(maneuv)[i], sep='.')] <- maneuv[ ,i] - maneuv[ ,i+24]
}

hist(maneuv$elevation) # birds were studied at <300m, and >3000m elev
maneuv$hi.elev <- maneuv$elevation > 3000

# standardize the behaviors and traits. STD = centered on 0 and scaled to have SD of 1
for(i in c(16:27,28:35,46:49)){
	maneuv[,paste(names(maneuv)[i],'STD',sep='.')] <- scale(maneuv[,i])
}
head(maneuv)

# get phylogeny:
tree <- read.nexus(file="hum294.tre") 
lookup <- read.csv('hum294_tiplabels.csv')
lookup <- subset(lookup, species!='' & species %in% maneuv$species)
lookup$tree <- factor(lookup$tree); lookup$species <- factor(lookup$species)
maneuvtree <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% lookup$tree])
plot(maneuvtree)
lookup <- lookup[match(maneuvtree$tip.label, lookup$tree), ] # reorder to match tree
maneuvtree$tip.label <- as.character(lookup$species)
plot(maneuvtree)
plottable <- unique(maneuv[,c('species','clade','color')])
plottable$pch <- c(15, 15, 16, 17, 15, 16, 15, 15, 15, 16, 17, 15, 16, 17, 16, 17, 18, 21, 22, 23, 15, 18, 21, 24, 25)
plottable # unique color-symbol combination for each species

# Figure 1 and Figure S1
# tree, individual traits, species averages.
# sample sizes

# tree

plotTree(maneuvtree, node.numbers=T) # painting branches using phytools package
COLmaneuvtree <- paintSubTree(maneuvtree, node=36, state="4") # 4.
COLmaneuvtree <- paintSubTree(COLmaneuvtree, node=46, state="3") # 3.
COLmaneuvtree <- paintSubTree(COLmaneuvtree, node=44, state="9") # 9
COLmaneuvtree <- paintSubTree(COLmaneuvtree, node=45, state="5") # 5.
COLmaneuvtree <- paintSubTree(COLmaneuvtree, node=41, state="8") # 8.
COLmaneuvtree <- paintSubTree(COLmaneuvtree, node=39, state="2") # 2.
COLmaneuvtree <- paintSubTree(COLmaneuvtree, node=27, state="6") # 6.
COLmaneuvtree <- paintSubTree(COLmaneuvtree, node=28, state="7") # 7.

# now let's plot using plotSimmap to ensure that the correct branches were painted
cols <- c("#000000", as.character(unique(plottable[,c('color')])))
names(cols) <- 1:9

sorttraits <- maneuv[,c('species','clade','mass','wing.area','wingload.gcm','AR')]
sorttraits$species <- factor(sorttraits$species, levels=maneuvtree$tip.label)
sorttraits$color <- plottable$color[match(sorttraits$clade, plottable$clade)]
sptraits$species <- factor(sptraits$species)
sptraits$species <- factor(sptraits$species, levels=maneuvtree$tip.label)

sprange <- summarize(LLsampleg, mass.lwr=min(mass, na.rm=T), mass.upr=max(mass, na.rm=T), winga.lwr=min(wing.area, na.rm=T), winga.upr=max(wing.area, na.rm=T), wingload.lwr=min(wing.load, na.rm=T), wingload.upr=max(wing.load, na.rm=T), AR.lwr=min(AR, na.rm=T), AR.upr=max(AR, na.rm=T), lift.mean=mean(lifted.beads, na.rm=T), lift.lwr=min(lifted.beads, na.rm=T), lift.upr=max(lifted.beads, na.rm=T))
sprange$species <- factor(sprange$species)
sprange$species <- factor(sprange$species, levels=maneuvtree$tip.label)
sprange <- sprange[with(sprange, order(species)),]

LLsampleg$species <- factor(LLsampleg$species)
LLsampleg$species <- factor(LLsampleg$species, levels=maneuvtree$tip.label)

elevationplot <- merge(unique(maneuv[,c('elevation','species')]), summarise(LLsampleg, mean.e=mean(elev, na.rm=T), min.e=min(elev, na.rm=T), max.e=max(elev, na.rm=T)), by='species')
elevationplot$species <- factor(elevationplot$species, levels=maneuvtree$tip.label)

# phylogeny, elevation, and trait data, for maneuvering and load lifting studies

dev.new(width=11,height=5)
layout(matrix(c(rep(1,6), rep(2,2), rep(3,2), rep(4,2), rep(5,2), rep(6,2), rep(7,2)), 2, 9, byrow=F))
plotSimmap(COLmaneuvtree, cols, pts=F, lwd=3, node.numbers=F, mar=c(4,0.1,0.1,0.1))
par(mar=c(4,2,0.25,0.25), las=1, bty='l', mgp=c(1.75, 0.5, 0), tck=-0.03)
plot(jitter(as.numeric(species), 0) ~ elevation, data=elevationplot, pch=plottable$pch, ylab='', xlab='elevation (m)', yaxt='n', col=as.character(plottable$color), xlim=c(-300, 3850), xaxt='n'); axis(2, at=1:25, labels=NA); axis(1, at=c(0,2000, 4000))
segments(y0=c(1:25)-0.45, y1=c(1:25)-0.25, x0=c(elevationplot$mean.e[match(levels(sorttraits$species),elevationplot$species)]), x1=c(elevationplot$mean.e[match(levels(sorttraits$species),elevationplot$species)]), col=as.character(plottable$color)[match(maneuvtree$tip.label, plottable$species)])
segments(x0=elevationplot$min.e[match(levels(sorttraits$species),elevationplot$species)], y0=c(1:25)-0.35, x1=elevationplot$max.e[match(levels(sorttraits$species),elevationplot$species)], y1=c(1:25)-0.35, col=as.character(plottable$color)[match(maneuvtree$tip.label, plottable$species)], lwd=0.5)
par(mar=c(4,2,0.25,0.25), las=1, bty='l', mgp=c(1.75, 0.5, 0), tck=-0.03)
plot(jitter(as.numeric(species), 0) ~ mass, sorttraits, pch=16, ylab='', xlab='body mass (g)', yaxt='n', col=(as.character(sorttraits$color))); axis(2, at=1:25, labels=NA)
segments(y0=c(1:25)-0.45, y1=c(1:25)-0.25, x0=c(sptraits$sp.mass[match(levels(sorttraits$species),sptraits$species)]), x1=c(sptraits$sp.mass[match(levels(sorttraits$species),sptraits$species)]), col=as.character(plottable$color)[match(maneuvtree$tip.label, plottable$species)])
segments(x0=c(sprange$mass.lwr), y0=c(1:25)-0.35, x1=c(sprange$mass.upr), y1=c(1:25)-0.35, col=as.character(plottable$color)[match(maneuvtree$tip.label, plottable$species)], lwd=0.5)
plot(jitter(as.numeric(species), 0) ~ c(wing.area/100), sorttraits, pch=16, ylab='', xlab='wing area (cm2)', yaxt='n', col=as.character(sorttraits$color)); axis(2, at=1:25, labels=NA)
segments(y=c(1:25)-0.45, y1=c(1:25)-0.25, x0=c(sptraits$sp.wing.area[match(levels(sorttraits$species),sptraits$species)]), x1=c(sptraits$sp.wing.area[match(levels(sorttraits$species),sptraits$species)]), col=as.character(plottable$color)[match(maneuvtree$tip.label, plottable$species)])
segments(x0=c(sprange$winga.lwr), y0=c(1:25)-0.35, x1=c(sprange$winga.upr), y1=c(1:25)-0.35, col=as.character(plottable$color)[match(maneuvtree$tip.label, plottable$species)], lwd=0.5)
plot(jitter(as.numeric(species), 0) ~ wingload.gcm, sorttraits, pch=16, ylab='', xlab='wing loading (g/cm2)', yaxt='n', col=as.character(sorttraits$color)); axis(2, at=1:25, labels=NA)
segments(y=c(1:25)-0.45, y1=c(1:25)-0.25, x0=c(sptraits$sp.wing.load[match(levels(sorttraits$species),sptraits$species)]), x1=c(sptraits$sp.wing.load[match(levels(sorttraits$species),sptraits$species)]), col=as.character(plottable$color)[match(maneuvtree$tip.label, plottable$species)])
segments(x0=c(sprange$wingload.lwr), y0=c(1:25)-0.35, x1=c(sprange$wingload.upr), y1=c(1:25)-0.35, col=as.character(plottable$color)[match(maneuvtree$tip.label, plottable$species)], lwd=0.5)
plot(jitter(as.numeric(species), 0) ~ AR, sorttraits, pch=16, ylab='', xlab='wing aspect ratio', yaxt='n', col=as.character(sorttraits$color), xlim=c(5.5,9.5)); axis(2, at=1:25, labels=NA)
segments(y0=c(1:25)-0.45, y1=c(1:25)-0.25, x0=c(sptraits$sp.AR[match(levels(sorttraits$species),sptraits$species)]), x1=c(sptraits$sp.AR[match(levels(sorttraits$species),sptraits$species)]), pch=2, cex=0.5, col=as.character(plottable$color)[match(maneuvtree$tip.label, plottable$species)])
segments(x0=c(sprange$AR.lwr), y0=c(1:25)-0.35, x1=c(sprange$AR.upr), y1=c(1:25)-0.35, col=as.character(plottable$color)[match(maneuvtree$tip.label, plottable$species)], lwd=0.5)
plot(jitter(as.numeric(species), 0) ~ lifted.beads, LLsampleg, type='n', pch=16, ylab='', xlab='load lifted (g)', yaxt='n', col=as.character(sorttraits$color)); axis(2, at=1:25, labels=NA)
segments(y=c(1:25)-0.45, y1=c(1:25)-0.25, x0=c(sprange$lift.mean[match(levels(sorttraits$species),sprange$species)]), x1=c(sprange$lift.mean[match(levels(sorttraits$species),sprange$species)]), pch=2, cex=0.5, col=as.character(plottable$color)[match(maneuvtree$tip.label, plottable$species)])
segments(x0=c(sprange$lift.lwr), y0=c(1:25)-0.35, x1=c(sprange$lift.upr), y1=c(1:25)-0.35, col=as.character(plottable$color)[match(maneuvtree$tip.label, plottable$species)], lwd=0.5)


# Figure 1 maneuver sample sizes
summary(maneuv[,28:35])
ss.maneuver <- data.frame(maneuver=c('V','HA','HD','PU','PD','Y','ARC','PRT'), mean=NA, lower=NA, upper=NA, SDlower=NA, SDupper=NA, sum=NA)
for(i in 1:8){
	ss.maneuver$mean[i] <- mean(maneuv[,27+i], na.rm=T)
	ss.maneuver$SDlower[i] <- mean(maneuv[,27+i], na.rm=T) - sd(maneuv[,27+i], na.rm=T)
	ss.maneuver$SDupper[i] <- mean(maneuv[,27+i], na.rm=T) + sd(maneuv[,27+i], na.rm=T)
	ss.maneuver$lower[i] <- min(maneuv[,27+i], na.rm=T)
	ss.maneuver$upper[i] <- max(maneuv[,27+i], na.rm=T)
	ss.maneuver$sum[i] <- sum(maneuv[,27+i], na.rm=T)
}

dev.new(width=3.75,height=3)
par(mfrow=c(1,2), mar=c(4,4,0.25,0.25), las=1, bty='l', mgp=c(1.75, 0.5, 0), tck=-0.03, cex=0.75)
plot(c(8:1) ~ sum, data=ss.maneuver, xlab='total # maneuvers', xlim=c(0,100000), yaxt='n', ylab='', xaxt='n')
axis(2, at=c(8:1), labels=c('Vel','AccHor','DecHor','PitchU','PitchD','Yaw','Arc','PRT'))
axis(1, at=c(0,50000,100000), labels=c('0','50,000','100,000'))
plot(c(8:1) ~ mean, data=ss.maneuver, xlab='# maneuvers/bird', xlim=c(0,1400), xaxt='n', yaxt='n', ylab='', type='n')
segments(y0=c(8:1)-0.2, y1=c(8:1)+0.2, x0=ss.maneuver$mean)
axis(2, at=c(8:1), labels=c('Vel','AccHor','DecHor','PitchU','PitchD','Yaw','Arc','PRT'))
axis(1, at=c(0, 1000, 2000), labels=c('0','1,000','2,000'))
for(i in 28:35){points(c(rep(-i+36,207) + runif(207,-0.1,0.1)) ~ maneuv[,i], pch=16, cex=0.4, col=rgb(0,0,0,0.2))}

# Figure S2 plot to illustrate the muscle capacity score
scores <- data.frame(score=ranef(mod1)[,1], species=rownames(ranef(mod1)))
agg <- aggregate(LLsample[,c('mass','lifted.beads','elev')], by=list(species=LLsample$species), FUN='mean', na.rm=T)
scores <- cbind(scores, agg[,2:4])
scores$col <- plottable$color; scores$pch <- plottable$pch
scores <- scores[with(scores, order(score)),]

dev.new(width=4,height=5)
par(mfrow=c(2,1), las=1, mar=c(5,3,0.25,0.25), mgp=c(1.5,0.5,0), cex=0.85, bty='l')
plot(scores$lifted.beads, pch=16, ylab='species average (g)', xaxt='n', xlab='', xaxt='n', yaxt='n', ylim=c(0,20), cex=0.5, col=as.character(scores$col))
axis(2, at=c(0,10,20), las=0)
points(scores$mass, cex=0.5, pch=1, col=as.character(scores$col))
text(x=1:25, y=scores$lifted.beads, round(scores$elev, -2), cex=0.35, srt=45, pos=3, offset=0.35)
legend('topleft', bty='n', pch=c(16,1), legend=c('beads lifted', 'body'), cex=0.5)
plot(scores$score, xlab='', xaxt='n', yaxt='n', ylab='muscle capacity score a.u.', cex=0.5, col=as.character(scores$col), pch=scores$pch)
axis(2, at=c(-4,0,4), las=0)
abline(h=0, lty=3)
for(i in 1:25) axis(1, at=i, labels=scores$species[i], cex.axis=0.5, las=2, col.axis=as.character(scores$col)[i])

# write.csv(maneuv, 'processed_maneuver_data_170724.csv')
# write.csv(plottable, 'plottable.csv')










