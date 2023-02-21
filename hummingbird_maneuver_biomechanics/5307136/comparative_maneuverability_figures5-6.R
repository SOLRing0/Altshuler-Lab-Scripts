
# Evolution reveals the biomechanical organization of maneuvering flight in hummingbirds
# Dakin, Segre, Straw and Altshuler
# Figure 5, Figure 6, and Figure S12

library(nlme)
library(visreg)
library(plotrix)

load('comparative_maneuverability_models.RData') # get model results & data

# re-fit models to raw (unstandardized) data for partial effect plots:
mod.listraw <- list()
mod.listraw[[1]] <- lme(hor_acel ~ hi.elev + sp.muscle + sp.mass + sp.wing.load + sp.AR + rel.mass + rel.wingload.gcm + rel.AR + ss.HA, random=~1|species, data=maneuv, na.action=na.omit, control=lmeControl(opt='optim'))
mod.listraw[[2]] <- lme(-hor_decel ~ hi.elev + sp.muscle + sp.mass + sp.wing.load + sp.AR + rel.mass + rel.wingload.gcm + rel.AR + ss.HD, random=~1|species, data=maneuv, na.action=na.omit, control=lmeControl(opt='optim'))
mod.listraw[[3]] <- lme(arc_force ~ hi.elev + sp.muscle + sp.mass + sp.wing.load + sp.AR + rel.mass + rel.wingload.gcm + rel.AR + ss.ARC, random=~1|species, data=maneuv, na.action=na.omit, control=lmeControl(opt='optim'))
mod.listraw[[4]] <- lme(-pitch_down ~ hi.elev + sp.muscle + sp.mass + sp.wing.load + sp.AR + rel.mass + rel.wingload.gcm + rel.AR + ss.PD, random=~1|species, data=maneuv, na.action=na.omit, control=lmeControl(opt='optim'))
mod.listraw[[5]] <- lme(yaw ~ hi.elev + sp.muscle + sp.mass + sp.wing.load + sp.AR + rel.mass + rel.wingload.gcm + rel.AR + ss.Y, random=~1|species, data=maneuv, na.action=na.omit, control=lmeControl(opt='optim'))
mod.listraw[[6]] <- lme(pPRT ~ hi.elev + sp.muscle + sp.mass + sp.wing.load + sp.AR + rel.mass + rel.wingload.gcm + rel.AR, random=~1|species, data=maneuv, na.action=na.omit, control=lmeControl(opt='optim'))
# order of decreasing effect sizes: yaw, pPRT, pitch down; accel, decel, arc; 

# Figure 5

# partial effect plots for species and individual wing load:
dev.new(width=6,height=8) 
par(las=1, mar=c(2,2,0.25,0.25), bty='l')
layout(matrix(c(1,3,5,2,4,6,7,9,11,8,10,12), nrow=4, ncol=3, byrow=T))

a <- visreg(mod.listraw[[5]], plot=F, type='conditional', cond=list('hi.elev'=T)) # yaw
b <- visreg(mod.listraw[[5]], plot=F, type='conditional', cond=list('hi.elev'=F))
temp2 <- aggregate(rbind(a[[4]]$res[,c('sp.wing.load','visregRes')], b[[4]]$res[,c('sp.wing.load','visregRes')]), by=list(species=rep(subset(maneuv,!is.na(mass)&!is.na(wingload.gcm))$species, 2)), 'mean')
temp2$color <- unique(maneuv[,c('species','color')])$color
temp2$n <- table(subset(maneuv,!is.na(mass)&!is.na(wingload.gcm))$species)
temp2 <- temp2[with(temp2, order(-n)),]
plot(visregRes ~ sp.wing.load, pch=16, data=temp2, xlim=c(0.1,0.4), ylim=c(1.3,2.3), ylab='residual', xlab='species wing loading', type='n', yaxt='n', xaxt='n')
axis(1, at=c(0.1,0.2,0.3,0.4)); axis(2, at=c(1.5,1.9,2.3))
for(i in 1:25){draw.circle(temp2$sp.wing.load[i], temp2$visregRes[i], radius=temp2$n[i]/650, col=as.character(temp2$color[i]))} # radius is based on x.
points(rowMeans(cbind(a[[4]]$fit$visregFit, b[[4]]$fit$visregFit)) ~ a[[4]]$fit$sp.wing.load, type='l', lwd=2.2)
plot(rowMeans(cbind(a[[7]]$res$visregRes, b[[7]]$res$visregRes)) ~ a[[7]]$res$rel.wingload.gcm, xlim=c(-0.1,0.2), ylim=c(1.25,2.25), ylab='residual', xlab='rel wing loading', type='p', col=as.character(temp2$color)[match(subset(maneuv,!is.na(mass)&!is.na(wingload.gcm))$species, temp2$species)], cex=0.5, pch=plottable$pch[match(subset(maneuv,!is.na(mass)&!is.na(wingload.gcm))$species, plottable$species)], xaxt='n', yaxt='n')
axis(1, at=c(-0.1,0,0.1,0.2)); axis(2, at=c(1.4,1.8,2.2)) 
points(rowMeans(cbind(a[[7]]$fit$visregFit, b[[7]]$fit$visregFit)) ~ a[[7]]$fit$rel.wingload.gcm, type='l', lty=3) 

a <- visreg(mod.listraw[[6]], plot=F, type='conditional', cond=list('hi.elev'=T)) # PRT%
b <- visreg(mod.listraw[[6]], plot=F, type='conditional', cond=list('hi.elev'=F)) 
temp2 <- aggregate(rbind(a[[4]]$res[,c('sp.wing.load','visregRes')], b[[4]]$res[,c('sp.wing.load','visregRes')]), by=list(species=rep(subset(maneuv,!is.na(mass)&!is.na(wingload.gcm))$species, 2)), 'mean')
temp2$color <- unique(maneuv[,c('species','color')])$color
temp2$n <- table(subset(maneuv,!is.na(mass)&!is.na(wingload.gcm))$species)
temp2 <- temp2[with(temp2, order(-n)),]
plot(visregRes ~ sp.wing.load, pch=16, data=temp2, xlim=c(0.1,0.4), ylim=c(0,1), ylab='residual', xlab='species wing loading', type='n', xaxt='n', yaxt='n')
axis(1, at=c(0.1,0.2,0.3,0.4)); axis(2, at=c(0,0.5,1), labels=c(0,50,100))
for(i in 1:25){draw.circle(temp2$sp.wing.load[i], temp2$visregRes[i], radius=temp2$n[i]/650, col=as.character(temp2$color[i]))} # radius is based on x.
points(rowMeans(cbind(a[[4]]$fit$visregFit, b[[4]]$fit$visregFit)) ~ a[[4]]$fit$sp.wing.load, type='l', lwd=2.2)
plot(rowMeans(cbind(a[[7]]$res$visregRes, b[[7]]$res$visregRes)) ~ a[[7]]$res$rel.wingload.gcm, xlim=c(-0.1,0.2), ylim=c(0,1), ylab='residual', xlab='rel wing loading', type='p', col=as.character(temp2$color)[match(subset(maneuv,!is.na(mass)&!is.na(wingload.gcm))$species, temp2$species)], cex=0.5, pch=plottable$pch[match(subset(maneuv,!is.na(mass)&!is.na(wingload.gcm))$species, plottable$species)], xaxt='n', yaxt='n')
axis(1, at=c(-0.1,0,0.1,0.2)); axis(2, at=c(0,0.5,1), labels=c(0,50,100)) 
points(rowMeans(cbind(a[[7]]$fit$visregFit, b[[7]]$fit$visregFit)) ~ a[[7]]$fit$rel.wingload.gcm, type='l', lty=3) 

a <- visreg(mod.listraw[[4]], plot=F, type='conditional', cond=list('hi.elev'=T)) # pitch down
b <- visreg(mod.listraw[[4]], plot=F, type='conditional', cond=list('hi.elev'=F)) 
temp2 <- aggregate(rbind(a[[4]]$res[,c('sp.wing.load','visregRes')], b[[4]]$res[,c('sp.wing.load','visregRes')]), by=list(species=rep(subset(maneuv,!is.na(mass)&!is.na(wingload.gcm))$species, 2)), 'mean')
temp2$color <- unique(maneuv[,c('species','color')])$color
temp2$n <- table(subset(maneuv,!is.na(mass)&!is.na(wingload.gcm))$species)
temp2 <- temp2[with(temp2, order(-n)),]
plot(visregRes ~ sp.wing.load, pch=16, data=temp2, xlim=c(0.1,0.4), ylim=c(0.6,1.6), ylab='residual', xlab='species wing loading', type='n', xaxt='n', yaxt='n')
axis(1, at=c(0.1,0.2,0.3,0.4)); axis(2, at=c(0.8,1.2,1.6))
for(i in 1:25){draw.circle(temp2$sp.wing.load[i], temp2$visregRes[i], radius=temp2$n[i]/650, col=as.character(temp2$color[i]))} # radius is based on x.
points(rowMeans(cbind(a[[4]]$fit$visregFit, b[[4]]$fit$visregFit)) ~ a[[4]]$fit$sp.wing.load, type='l', lwd=2.2)
plot(rowMeans(cbind(a[[7]]$res$visregRes, b[[7]]$res$visregRes)) ~ a[[7]]$res$rel.wingload.gcm, xlim=c(-0.1,0.2), ylim=c(0.6,1.6), ylab='residual', xlab='rel wing loading', type='p', col=as.character(temp2$color)[match(subset(maneuv,!is.na(mass)&!is.na(wingload.gcm))$species, temp2$species)], cex=0.5, pch=plottable$pch[match(subset(maneuv,!is.na(mass)&!is.na(wingload.gcm))$species, plottable$species)], xaxt='n', yaxt='n')
axis(1, at=c(-0.1,0,0.1,0.2)); axis(2, at=c(0.8,1.2,1.6)) 
points(rowMeans(cbind(a[[7]]$fit$visregFit, b[[7]]$fit$visregFit)) ~ a[[7]]$fit$rel.wingload.gcm, type='l', lty=3) 

a <- visreg(mod.listraw[[1]], plot=F, type='conditional', cond=list('hi.elev'=T)) # accel
b <- visreg(mod.listraw[[1]], plot=F, type='conditional', cond=list('hi.elev'=F)) 
temp2 <- aggregate(rbind(a[[4]]$res[,c('sp.wing.load','visregRes')], b[[4]]$res[,c('sp.wing.load','visregRes')]), by=list(species=rep(subset(maneuv,!is.na(mass)&!is.na(wingload.gcm))$species, 2)), 'mean')
temp2$color <- unique(maneuv[,c('species','color')])$color
temp2$n <- table(subset(maneuv,!is.na(mass)&!is.na(wingload.gcm))$species)
temp2 <- temp2[with(temp2, order(-n)),]
plot(visregRes ~ sp.wing.load, pch=16, data=temp2, xlim=c(0.1,0.4), ylim=c(2,8.1), ylab='residual', xlab='species wing loading', type='n', yaxt='n', xaxt='n')
axis(2, at=c(2,4,6,8)); axis(1, at=c(0.1,0.2,0.3,0.4))
for(i in 1:25){draw.circle(temp2$sp.wing.load[i], temp2$visregRes[i], radius=temp2$n[i]/650, col=as.character(temp2$color[i]))} # radius is based on x.
points(rowMeans(cbind(a[[4]]$fit$visregFit, b[[4]]$fit$visregFit)) ~ a[[4]]$fit$sp.wing.load, type='l', lty=3)
for(i in 1:3){draw.circle(x=0.25, y=7, radius=c(c(1,10,20)/650)[i], col=rgb(0,0,0,0.4), border=NULL)}
plot(rowMeans(cbind(a[[7]]$res$visregRes, b[[7]]$res$visregRes)) ~ a[[7]]$res$rel.wingload.gcm, xlim=c(-0.1,0.2), ylim=c(2,8.1), ylab='residual', xlab='rel wing loading', type='p', col=as.character(temp2$color)[match(subset(maneuv,!is.na(mass)&!is.na(wingload.gcm))$species, temp2$species)], cex=0.5, pch=plottable$pch[match(subset(maneuv,!is.na(mass)&!is.na(wingload.gcm))$species, plottable$species)], xaxt='n', yaxt='n')
axis(2, at=c(2,4,6,8)); axis(1, at=c(-0.1,0,0.1,0.2)) 
points(rowMeans(cbind(a[[7]]$fit$visregFit, b[[7]]$fit$visregFit)) ~ a[[7]]$fit$rel.wingload.gcm, type='l', lwd=2.2) 

a <- visreg(mod.listraw[[2]], plot=F, type='conditional', cond=list('hi.elev'=T)) # decel
b <- visreg(mod.listraw[[2]], plot=F, type='conditional', cond=list('hi.elev'=F))
temp2 <- aggregate(rbind(a[[4]]$res[,c('sp.wing.load','visregRes')], b[[4]]$res[,c('sp.wing.load','visregRes')]), by=list(species=rep(subset(maneuv,!is.na(mass)&!is.na(wingload.gcm))$species, 2)), 'mean')
temp2$color <- unique(maneuv[,c('species','color')])$color
temp2$n <- table(subset(maneuv,!is.na(mass)&!is.na(wingload.gcm))$species)
temp2 <- temp2[with(temp2, order(-n)),]
plot(visregRes ~ sp.wing.load, pch=16, data=temp2, xlim=c(0.1,0.4), ylim=c(2,8.1), ylab='residual', xlab='species wing loading', type='n', xaxt='n', yaxt='n')
axis(2, at=c(2,4,6,8)); axis(1, at=c(0.1,0.2,0.3,0.4))
for(i in 1:25){draw.circle(temp2$sp.wing.load[i], temp2$visregRes[i], radius=temp2$n[i]/650, col=as.character(temp2$color[i]))} # radius is based on x.
points(rowMeans(cbind(a[[4]]$fit$visregFit, b[[4]]$fit$visregFit)) ~ a[[4]]$fit$sp.wing.load, type='l', lty=3)
plot(rowMeans(cbind(a[[7]]$res$visregRes, b[[7]]$res$visregRes)) ~ a[[7]]$res$rel.wingload.gcm, xlim=c(-0.1,0.2), ylim=c(2,8.1), ylab='residual', xlab='rel wing loading', type='p', col=as.character(temp2$color)[match(subset(maneuv,!is.na(mass)&!is.na(wingload.gcm))$species, temp2$species)], cex=0.5, pch=plottable$pch[match(subset(maneuv,!is.na(mass)&!is.na(wingload.gcm))$species, plottable$species)], xaxt='n', yaxt='n') 
axis(2, at=c(2,4,6,8)); axis(1, at=c(-0.1,0,0.1,0.2)) 
points(rowMeans(cbind(a[[7]]$fit$visregFit, b[[7]]$fit$visregFit)) ~ a[[7]]$fit$rel.wingload.gcm, type='l', lwd=2.2) 

a <- visreg(mod.listraw[[3]], plot=F, type='conditional', cond=list('hi.elev'=T)) # centripetal
b <- visreg(mod.listraw[[3]], plot=F, type='conditional', cond=list('hi.elev'=F)) 
temp2 <- aggregate(rbind(a[[4]]$res[,c('sp.wing.load','visregRes')], b[[4]]$res[,c('sp.wing.load','visregRes')]), by=list(species=rep(subset(maneuv,!is.na(mass)&!is.na(wingload.gcm)&!is.na(arc_force))$species, 2)), 'mean')
temp2$color <- unique(maneuv[,c('species','color')])$color
temp2$n <- table(subset(maneuv,!is.na(mass)&!is.na(wingload.gcm)&!is.na(arc_force))$species)
temp2 <- temp2[with(temp2, order(-n)),]
plot(visregRes ~ sp.wing.load, pch=16, data=temp2, xlim=c(0.1,0.4), ylim=c(2,8.1), ylab='residual', xlab='species wing loading', type='n', xaxt='n', yaxt='n')
axis(2, at=c(2,4,6,8)); axis(1, at=c(0.1,0.2,0.3,0.4))
for(i in 1:25){draw.circle(temp2$sp.wing.load[i], temp2$visregRes[i], radius=temp2$n[i]/650, col=as.character(temp2$color[i]))} # radius is based on x.
points(rowMeans(cbind(a[[4]]$fit$visregFit, b[[4]]$fit$visregFit)) ~ a[[4]]$fit$sp.wing.load, type='l', lty=3)
plot(rowMeans(cbind(a[[7]]$res$visregRes, b[[7]]$res$visregRes)) ~ a[[7]]$res$rel.wingload.gcm, xlim=c(-0.1,0.2), ylim=c(2,8.1), ylab='residual', xlab='rel wing loading', type='p', col=as.character(temp2$color)[match(subset(maneuv,!is.na(mass)&!is.na(wingload.gcm)&!is.na(arc_force))$species, temp2$species)], cex=0.5, pch=plottable$pch[match(subset(maneuv,!is.na(mass)&!is.na(wingload.gcm)&!is.na(arc_force))$species, plottable$species)], xaxt='n', yaxt='n')
axis(2, at=c(2,4,6,8)); axis(1, at=c(-0.1,0,0.1,0.2)) 
points(rowMeans(cbind(a[[7]]$fit$visregFit, b[[7]]$fit$visregFit)) ~ a[[7]]$fit$rel.wingload.gcm, type='l', lty=3) 


# Figure 6

# partial effect plots for further models of PRT%:

tmod2plothi <- visreg(PRTmod1, cond=list(hi.elev=1), plot=F)
tmod2plotlo <- visreg(PRTmod1, cond=list(hi.elev=0), plot=F)

eleveffects <- eleveffects[with(eleveffects, order(abs(FULL.hi.elevTRUE))),]

dev.new(width=4.75,height=6)
par(mar=c(3,2,0.25,0.25), mgp=c(1.5,0.5,0), las=1, bty='l')
layout(matrix(c(rep(1,4), 2:5, 6,6,7,7), nrow=3, ncol=4, byrow=T))

plot(abs(eleveffects[,2]) ~ c(1:12), pch=16, ylim=c(0,1.75), ylab='|effect size|', yaxt='n', xaxt='n', xlab='')
segments(x0=c(1:12), y0=eleveffects$lower, y1=eleveffects$upper)
axis(2, at=c(0,0.5,1,1.5)); axis(1, at=1:12, labels=eleveffects$metric, cex.axis=0.5)

temp <- aggregate(tmod2plothi[[1]]$res[,c('hi.elev','visregRes')], by=list(species=subset(maneuv, !is.na(mass)&!is.na(wing.area)&!is.na(arc_force))$species), FUN='mean')
plot(visregRes ~ jitter(hi.elev+1), data=temp, ylim=c(0.2,1.1), xlim=c(0.5,2.5), ylab='', yaxt='n', xaxt='n', pch=plottable$pch[match(temp$species, plottable$species)], col=as.character(plottable$color[match(temp$species, plottable$species)]), xlab='elevation', cex=1.5, type='p') # 
axis(2, at=c(0.3,0.6,0.9), labels=c(30,60,90))
axis(1, at=c(1,2), labels=c('low', 'high'))
boxplot(temp$visregRes ~ temp$hi.elev, range=0, add=T, xaxt='n', yaxt='n')

temp1 <- aggregate(tmod2plothi[[2]]$res[,c('sp.mass','visregRes')], by=list(species=subset(maneuv, !is.na(mass)&!is.na(wing.area)&!is.na(arc_force))$species), FUN='mean')
temp2 <- aggregate(tmod2plotlo[[2]]$res[,c('sp.mass','visregRes')], by=list(species=subset(maneuv, !is.na(mass)&!is.na(wing.area)&!is.na(arc_force))$species), FUN='mean')
plot(rowMeans(cbind(temp1$visregRes, temp2$visregRes)) ~ temp1$sp.mass, ylim=c(0.2,1.1), ylab='', yaxt='n', xaxt='n', pch=plottable$pch[match(spturns$species, plottable$species)], col=as.character(plottable$color[match(spturns$species, plottable$species)]), xlab='sp mass (g)', cex=1.5)
axis(2, at=c(0.3,0.6,0.9), labels=c(30,60,90))
axis(1, at=c(3,6,9))
points(rowMeans(cbind(tmod2plotlo[[2]]$fit$visregFit, tmod2plothi[[2]]$fit$visregFit)) ~ tmod2plothi[[2]]$fit$sp.mass, type='l', lwd=2.2)

temp1 <- aggregate(tmod2plothi[[3]]$res[,c('sp.wing.load','visregRes')], by=list(species=subset(maneuv, !is.na(mass)&!is.na(wing.area)&!is.na(arc_force))$species), FUN='mean')
temp2 <- aggregate(tmod2plotlo[[3]]$res[,c('sp.wing.load','visregRes')], by=list(species=subset(maneuv, !is.na(mass)&!is.na(wing.area)&!is.na(arc_force))$species), FUN='mean')
plot(rowMeans(cbind(temp1$visregRes, temp2$visregRes)) ~ temp1$sp.wing.load, ylim=c(0.2,1.1), ylab='', yaxt='n', xaxt='n', pch=plottable$pch[match(spturns$species, plottable$species)], col=as.character(plottable$color[match(spturns$species, plottable$species)]), xlab='sp wing loading (g/cm2)', cex=1.5)
axis(2, at=c(0.3,0.6,0.9), labels=c(30,60,90))
axis(1, at=c(0.2,0.3))
points(rowMeans(cbind(tmod2plotlo[[3]]$fit$visregFit, tmod2plothi[[3]]$fit$visregFit)) ~ tmod2plothi[[3]]$fit$sp.wing.load, type='l', lwd=2.2)

temp1 <- aggregate(tmod2plothi[[5]]$res[,c('sp.arc_force','visregRes')], by=list(species=subset(maneuv, !is.na(mass)&!is.na(wing.area)&!is.na(arc_force))$species), FUN='mean')
temp2 <- aggregate(tmod2plotlo[[5]]$res[,c('sp.arc_force','visregRes')], by=list(species=subset(maneuv, !is.na(mass)&!is.na(wing.area)&!is.na(arc_force))$species), FUN='mean')
plot(rowMeans(cbind(temp1$visregRes, temp2$visregRes)) ~ temp1$sp.arc_force, ylim=c(0.2,1.1), xlim=c(2.5,8.5), ylab='', yaxt='n', xaxt='n', pch=plottable$pch[match(spturns$species, plottable$species)], col=as.character(plottable$color[match(spturns$species, plottable$species)]), xlab='sp arc cent (m/s2)', cex=1.5)
axis(2, at=c(0.3,0.6,0.9), labels=c(30,60,90))
axis(1, at=c(4,6,8))
points(rowMeans(cbind(tmod2plotlo[[5]]$fit$visregFit, tmod2plothi[[5]]$fit$visregFit)) ~ tmod2plotlo[[5]]$fit$sp.arc_force, type='l', lwd=2.2)

plot(rowMeans(cbind(tmod2plothi[[4]]$res$visregRes, tmod2plotlo[[4]]$res$visregRes)) ~ tmod2plothi[[4]]$res$rel.AR, ylab='', yaxt='n', xaxt='n', xlab='rel AR', xlim=c(-2,2), ylim=c(0.2,1.1), cex=0.75, col=as.character(maneuv$color[!is.na(maneuv$arc_force)&!is.na(maneuv$mass)&!is.na(maneuv$wing.area)]), pch=plottable$pch[match(maneuv$species[!is.na(maneuv$arc_force)&!is.na(maneuv$mass)&!is.na(maneuv$wing.area)], plottable$species)])
points(rowMeans(cbind(tmod2plotlo[[4]]$fit$visregFit, tmod2plothi[[4]]$fit$visregFit)) ~ tmod2plothi[[4]]$fit$rel.AR, type='l', lwd=2.2)
axis(2, at=c(0.3,0.6,0.9), labels=c(30,60,90))
axis(1, at=c(-2,0,2))

plot(rowMeans(cbind(tmod2plothi[[6]]$res$visregRes, tmod2plotlo[[6]]$res$visregRes)) ~ tmod2plothi[[6]]$res$rel.arc_force, ylab='', yaxt='n', xaxt='n', xlab='rel arc cent (m/s2)', ylim=c(0.2,1.1), xlim=c(-3,3), cex=0.75, col=as.character(maneuv$color[!is.na(maneuv$arc_force)&!is.na(maneuv$mass)&!is.na(maneuv$wing.area)]), pch=plottable$pch[match(maneuv$species[!is.na(maneuv$arc_force)&!is.na(maneuv$mass)&!is.na(maneuv$wing.area)], plottable$species)])
points(rowMeans(cbind(tmod2plotlo[[6]]$fit$visregFit, tmod2plothi[[6]]$fit$visregFit)) ~ tmod2plothi[[6]]$fit$rel.arc_force, type='l', lwd=2.2)
axis(2, at=c(0.3,0.6,0.9), labels=c(30,60,90))
axis(1, at=c(-2,0,2))




# Figure S12

# sample size results

# refit models on scale of original data:
mod.listraw[[7]] <- lme(PRT_degrees ~ hi.elev + sp.muscle + sp.mass + sp.wing.load + sp.AR + rel.mass + rel.wingload.gcm + rel.AR + ss.PRT, random=~1|species, data=maneuv, na.action=na.omit, control=lmeControl(opt='optim'))
mod.listraw[[8]] <- lme(pitch_up~ hi.elev + sp.muscle + sp.mass + sp.wing.load + sp.AR + rel.mass + rel.wingload.gcm + rel.AR + ss.PU, random=~1|species, data=maneuv, na.action=na.omit, control=lmeControl(opt='optim'))
mod.listraw[[9]] <- lme(total_vel ~ hi.elev + sp.muscle + sp.mass + sp.wing.load + sp.AR + rel.mass + rel.wingload.gcm + rel.AR + ss.V, random=~1|species, data=maneuv, na.action=na.omit, control=lmeControl(opt='optim'))
mod.listraw[[10]] <- lme(arc_radius~ hi.elev + sp.muscle + sp.mass + sp.wing.load + sp.AR + rel.mass + rel.wingload.gcm + rel.AR + ss.ARC, random=~1|species, data=maneuv, na.action=na.omit, control=lmeControl(opt='optim'))
mod.listraw[[11]] <- lme(-PRT_time ~ hi.elev + sp.muscle + sp.mass + sp.wing.load + sp.AR + rel.mass + rel.wingload.gcm + rel.AR + ss.PRT, random=~1|species, data=maneuv, na.action=na.omit, control=lmeControl(opt='optim'))

dev.new(width=4, height=6)
par(mar=c(2,2,2,0.25), las=1, bty='l', mgp=c(1.5,0.5,0), xpd=T)
layout(matrix(c(1,1,2,3,4,5,6,7,8,9), nrow=5, ncol=2, byrow=T))

plot(abs(FULL.ss.STD) ~ c(1:11), data=sseffects, xlim=c(1,12), ylim=c(-0.05,0.75), pch=16, ylab='|effect size|', xlab='', xaxt='n', col=c('#FA2800', '#00A0FA')[(sseffects$FULL.ss.STD > 0) + 1], main='# of maneuvers performed', cex.main=0.75); abline(h=0, lty=3)
axis(1, at=c(1:11), labels=F)
text(x=c(1:11), y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]), labels=sseffects$metric, srt=45, adj=1, xpd=TRUE, cex=0.5)
segments(x0=1:11, y0=sseffects$lower, y1=sseffects$upper, col=c('#FA2800', '#00A0FA')[(sseffects$FULL.ss.STD > 0) + 1])

temp1 <- visreg(mod.listraw[[3]], 'ss.ARC', cond=list(hi.elev=1), plot=F) # arc force 
temp2 <- visreg(mod.listraw[[3]], 'ss.ARC', cond=list(hi.elev=0), plot=F)
plot(rowMeans(cbind(temp1$res$visregRes, temp2$res$visregRes)) ~ temp1$res$ss.ARC, xlab='', ylab='', xlim=c(0,1300), xaxt='n', ylim=c(1,8.5), yaxt='n', main='arc cent,max (m/s2)', cex.main=0.5, cex=0.75, col=as.character(maneuv$color[!is.na(maneuv$mass)&!is.na(maneuv$wing.area)&!is.na(maneuv$arc_force)]), pch=plottable$pch[match(maneuv$species[!is.na(maneuv$mass)&!is.na(maneuv$wing.area)&!is.na(maneuv$arc_force)], plottable$species)])
axis(1, at=c(0,500,1000))
axis(2, at=c(2,5,8))
points(rowMeans(cbind(temp1$fit$visregFit, temp2$fit$visregFit)) ~ temp1$fit$ss.ARC, type='l', lwd=2)

temp1 <- visreg(mod.listraw[[7]], 'ss.PRT', cond=list(hi.elev=1), plot=F) # prt deg 
temp2 <- visreg(mod.listraw[[7]], 'ss.PRT', cond=list(hi.elev=0), plot=F) 
plot(rowMeans(cbind(temp1$res$visregRes, temp2$res$visregRes)) ~ temp1$res$ss.PRT, xlab='', ylab='', xlim=c(0,1300), xaxt='n', ylim=c(0,180), yaxt='n', main='PRT degrees (º)', cex.main=0.5, cex=0.75, col=as.character(maneuv$color[!is.na(maneuv$mass)&!is.na(maneuv$wing.area)]), pch=plottable$pch[match(maneuv$species[!is.na(maneuv$mass)&!is.na(maneuv$wing.area)], plottable$species)])
axis(1, at=c(0,500,1000))
axis(2, at=c(0,90,180))
points(rowMeans(cbind(temp1$fit$visregFit, temp2$fit$visregFit)) ~ temp1$fit$ss.PRT, type='l', lwd=2)

temp1 <- visreg(mod.listraw[[4]], 'ss.PD', cond=list(hi.elev=1), plot=F) # pitch down 
temp2 <- visreg(mod.listraw[[4]], 'ss.PD', cond=list(hi.elev=0), plot=F) 
plot(rowMeans(cbind(temp1$res$visregRes, temp2$res$visregRes)) ~ temp1$res$ss.PD, xlab='', ylab='', xlim=c(0,1300), xaxt='n', ylim=c(0.7,1.5), yaxt='n', main='pitch down (rev/s)', cex.main=0.5, cex=0.75, col=as.character(maneuv$color[!is.na(maneuv$mass)&!is.na(maneuv$wing.area)]), pch=plottable$pch[match(maneuv$species[!is.na(maneuv$mass)&!is.na(maneuv$wing.area)], plottable$species)])
axis(1, at=c(0,500,1000))
axis(2, at=c(0.7,1.1,1.5))
points(rowMeans(cbind(temp1$fit$visregFit, temp2$fit$visregFit)) ~ temp1$fit$ss.PD, type='l', lwd=2)

temp1 <- visreg(mod.listraw[[8]], 'ss.PU', cond=list(hi.elev=1), plot=F) # pitch up 
temp2 <- visreg(mod.listraw[[8]], 'ss.PU', cond=list(hi.elev=0), plot=F)
plot(rowMeans(cbind(temp1$res$visregRes, temp2$res$visregRes)) ~ temp1$res$ss.PU, xlab='', ylab='', xlim=c(0,1300), xaxt='n', ylim=c(0.9,1.7), yaxt='n', main='pitch up (rev/s)', cex.main=0.5, cex=0.75, col=as.character(maneuv$color[!is.na(maneuv$mass)&!is.na(maneuv$wing.area)]), pch=plottable$pch[match(maneuv$species[!is.na(maneuv$mass)&!is.na(maneuv$wing.area)], plottable$species)])
axis(1, at=c(0,500,1000))
axis(2, at=c(0.9,1.3,1.7))
points(rowMeans(cbind(temp1$fit$visregFit, temp2$fit$visregFit)) ~ temp1$fit$ss.PU, type='l', lwd=2)

temp1 <- visreg(mod.listraw[[9]], 'ss.V', cond=list(hi.elev=1), plot=F) # total vel 
temp2 <- visreg(mod.listraw[[9]], 'ss.V', cond=list(hi.elev=0), plot=F)
plot(rowMeans(cbind(temp1$res$visregRes, temp2$res$visregRes)) ~ temp1$res$ss.V, xlab='', ylab='', xlim=c(0,1300), xaxt='n', ylim=c(1,3), yaxt='n', main='vel (m/s)', cex.main=0.5, cex=0.75, col=as.character(maneuv$color[!is.na(maneuv$mass)&!is.na(maneuv$wing.area)]), pch=plottable$pch[match(maneuv$species[!is.na(maneuv$mass)&!is.na(maneuv$wing.area)], plottable$species)])
axis(1, at=c(0,500,1000))
axis(2, at=c(1,2,3))
points(rowMeans(cbind(temp1$fit$visregFit, temp2$fit$visregFit)) ~ temp1$fit$ss.V, type='l', lwd=2)

temp1 <- visreg(mod.listraw[[10]], 'ss.ARC', cond=list(hi.elev=1), plot=F) # arc radius 
temp2 <- visreg(mod.listraw[[10]], 'ss.ARC', cond=list(hi.elev=0), plot=F)
plot(rowMeans(cbind(temp1$res$visregRes, temp2$res$visregRes)) ~ temp1$res$ss.ARC, xlab='', ylab='', xlim=c(0,1300), xaxt='n', ylim=c(0,1), yaxt='n', main='arc radius (m)', cex.main=0.5, cex=0.75, col=as.character(maneuv$color[!is.na(maneuv$mass)&!is.na(maneuv$wing.area)&!is.na(maneuv$arc_force)]), pch=plottable$pch[match(maneuv$species[!is.na(maneuv$mass)&!is.na(maneuv$wing.area)&!is.na(maneuv$arc_force)], plottable$species)])
axis(1, at=c(0,500,1000))
axis(2, at=c(0,0.5,1))
points(rowMeans(cbind(temp1$fit$visregFit, temp2$fit$visregFit)) ~ temp1$fit$ss.ARC, type='l', lwd=2)

temp1 <- visreg(mod.listraw[[11]], 'ss.PRT', cond=list(hi.elev=1), plot=F) # prt time 
temp2 <- visreg(mod.listraw[[11]], 'ss.PRT', cond=list(hi.elev=0), plot=F)
plot(rowMeans(cbind(temp1$res$visregRes, temp2$res$visregRes)) ~ temp1$res$ss.PRT, xlab='', ylab='', xlim=c(0,1300), xaxt='n', ylim=c(-0.7,-0.3), yaxt='n', main='PRT time (s)', cex.main=0.5, cex=0.75, col=as.character(maneuv$color[!is.na(maneuv$mass)&!is.na(maneuv$wing.area)]), pch=plottable$pch[match(maneuv$species[!is.na(maneuv$mass)&!is.na(maneuv$wing.area)], plottable$species)])
axis(1, at=c(0,500,1000))
axis(2, at=c(-0.7,-0.5,-0.3), labels=c(0.7,0.5,0.3))
points(rowMeans(cbind(temp1$fit$visregFit, temp2$fit$visregFit)) ~ temp1$fit$ss.PRT, type='l', lwd=2)

temp1 <- visreg(mod.listraw[[1]], 'ss.HA', cond=list(hi.elev=1), plot=F) # hor acel 
temp2 <- visreg(mod.listraw[[1]], 'ss.HA', cond=list(hi.elev=0), plot=F)
plot(rowMeans(cbind(temp1$res$visregRes, temp2$res$visregRes)) ~ temp1$res$ss.HA, xlab='', ylab='', xlim=c(0,1300), xaxt='n', ylim=c(1,8.5), yaxt='n', main='acchor (m/s2)', cex.main=0.5, cex=0.75, col=as.character(maneuv$color[!is.na(maneuv$mass)&!is.na(maneuv$wing.area)]), pch=plottable$pch[match(maneuv$species[!is.na(maneuv$mass)&!is.na(maneuv$wing.area)], plottable$species)])
axis(1, at=c(0,500,1000))
axis(2, at=c(2,5,8))
points(rowMeans(cbind(temp1$fit$visregFit, temp2$fit$visregFit)) ~ temp1$fit$ss.HA, type='l', lwd=2)






