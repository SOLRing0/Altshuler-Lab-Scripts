
# Hummingbird evolution reveals the biomechanical organization of maneuverability
# Dakin, Segre, and Altshuler
# models, Figures S6, S7, S10, S11

library(ape)
library(MCMCglmm)
library(nlme)
library(lattice)
library(MuMIn)

# get maneuver data: 

maneuv <- read.csv('processed_maneuver_data_170724.csv')[,-1]

# get phylogeny:

tree <- read.nexus(file="hum294.tre") 
lookup <- read.csv('hum294_tiplabels.csv')
lookup <- subset(lookup, species!='' & species %in% maneuv$species)
lookup$tree <- factor(lookup$tree); lookup$species <- factor(lookup$species)
maneuvtree <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% lookup$tree])
lookup <- lookup[match(maneuvtree$tip.label, lookup$tree), ] # reorder to match tree
maneuvtree$tip.label <- as.character(lookup$species)
plot(maneuvtree)

# priors and phylo for MCMC models:

inv.phylo <- inverseA(maneuvtree, nodes="TIPS", scale=TRUE)
priorlist <- list(G=list(G1=list(V=1, nu=0.02), G2=list(V=1, nu=0.02)), R=list(V=1,nu=0.02))

# set up the lists to store models, scaling models (mass only), full models (+ other traits), x LME & MCMC
scalingmods <- list()
fullmods <- list()
nophyloscalingmods <- list()
nophylomods <- list()
burnin <- 3000
thin <- 500
nitt <- 1000000 + burnin

# function to run diagnostic checks for MCMC models:
mcmccheck = function(m1, m2){
	autocorr1 <- range(as.matrix(autocorr(m1$Sol))[seq(2, length(as.matrix(autocorr(m1$Sol))), by=5)])
	autocorr2 <- range(as.matrix(autocorr(m1$VCV))[seq(2, length(as.matrix(autocorr(m1$VCV))), by=5)])
	mbind <- lapply(list(m1, m2), function(m) m$Sol)
	mbind <- do.call(mcmc.list, mbind)
	geldiag <- gelman.diag(mbind)
	dev.new(width=nrow(summary(m1)$solutions), height=2.5)
	par(mfcol=c(2,nrow(summary(m1)$solutions)), mar=c(2,1,1,1), cex=0.25)
	mplot <- plot(mbind, ask=F, auto.layout=F)
	return(list(m1.autocorrelation = autocorr1, m2.autorcorrelation = autocorr2, Gelman.Rubin.convergence = geldiag, summary(m1), mplot))
}

# Vel,max
set.seed(1112)
scalingmods[[1]] <- MCMCglmm(total_vel.STD ~ hi.elev + sp.mass.STD + rel.mass.STD + ss.V.STD, random=~species+sp, family='gaussian', ginverse=list(species=inv.phylo$Ainv), prior=priorlist, data=subset(maneuv, !is.na(mass)), nitt=nitt, burnin=burnin, thin=thin); alarm()
scalingmods[[2]] <- MCMCglmm(total_vel.STD ~ hi.elev + sp.mass.STD + rel.mass.STD + ss.V.STD, random=~species+sp, family='gaussian', ginverse=list(species=inv.phylo$Ainv), prior=priorlist, data=subset(maneuv, !is.na(mass)), nitt=nitt, burnin=burnin, thin=thin); alarm()
mcmccheck(scalingmods[[1]], scalingmods[[2]])

set.seed(1212)
fullmods[[1]] <- MCMCglmm(total_vel.STD ~ hi.elev + sp.muscle.STD + sp.mass.STD + sp.wing.load.STD + sp.AR.STD + rel.mass.STD  + rel.wingload.gcm.STD + rel.AR.STD + ss.V.STD, random=~species+sp, family="gaussian", ginverse=list(species=inv.phylo$Ainv), prior=priorlist, data=subset(maneuv, !is.na(mass) & !is.na(wing.area)), nitt=nitt, burnin=burnin, thin=thin); alarm()
fullmods[[2]] <- MCMCglmm(total_vel.STD ~ hi.elev + sp.muscle.STD + sp.mass.STD + sp.wing.load.STD + sp.AR.STD + rel.mass.STD  + rel.wingload.gcm.STD + rel.AR.STD + ss.V.STD, random=~species+sp, family="gaussian", ginverse=list(species=inv.phylo$Ainv), prior=priorlist, data=subset(maneuv, !is.na(mass) & !is.na(wing.area)), nitt=nitt, burnin=burnin, thin=thin); alarm()
mcmccheck(fullmods[[1]], fullmods[[2]])

nophyloscalingmods[[1]] <- lme(total_vel.STD ~ hi.elev + sp.mass.STD + rel.mass.STD + ss.V.STD, random=~1|species, data=subset(maneuv, !is.na(mass) & !is.na(wing.area)), control=lmeControl(opt='optim'))
nophylomods[[1]] <- lme(total_vel.STD ~ hi.elev + sp.muscle.STD + sp.mass.STD + sp.wing.load.STD + sp.AR.STD + rel.mass.STD  + rel.wingload.gcm.STD + rel.AR.STD + ss.V.STD, random=~1|species, data=subset(maneuv, !is.na(mass) & !is.na(wing.area)), control=lmeControl(opt='optim'))

# AccHor,max
set.seed(2112)
scalingmods[[3]] <- MCMCglmm(hor_acel.STD ~ hi.elev + sp.mass.STD + rel.mass.STD + ss.HA.STD, random=~species+sp, family='gaussian', ginverse=list(species=inv.phylo$Ainv), prior=priorlist, data=subset(maneuv, !is.na(mass)), nitt=nitt, burnin=burnin, thin=thin); alarm()
scalingmods[[4]] <- MCMCglmm(hor_acel.STD ~ hi.elev + sp.mass.STD + rel.mass.STD + ss.HA.STD, random=~species+sp, family='gaussian', ginverse=list(species=inv.phylo$Ainv), prior=priorlist, data=subset(maneuv, !is.na(mass)), nitt=nitt, burnin=burnin, thin=thin); alarm()
mcmccheck(scalingmods[[3]], scalingmods[[4]])

set.seed(2212)
fullmods[[3]] <- MCMCglmm(hor_acel.STD ~ hi.elev + sp.muscle.STD + sp.mass.STD + sp.wing.load.STD + sp.AR.STD + rel.mass.STD  + rel.wingload.gcm.STD + rel.AR.STD + ss.HA.STD, random=~species+sp, family="gaussian", ginverse=list(species=inv.phylo$Ainv), prior=priorlist, data=subset(maneuv, !is.na(mass) & !is.na(wing.area)), nitt=nitt, burnin=burnin, thin=thin); alarm()
fullmods[[4]] <- MCMCglmm(hor_acel.STD ~ hi.elev + sp.muscle.STD + sp.mass.STD + sp.wing.load.STD + sp.AR.STD + rel.mass.STD  + rel.wingload.gcm.STD + rel.AR.STD + ss.HA.STD, random=~species+sp, family="gaussian", ginverse=list(species=inv.phylo$Ainv), prior=priorlist, data=subset(maneuv, !is.na(mass) & !is.na(wing.area)), nitt=nitt, burnin=burnin, thin=thin); alarm()
mcmccheck(fullmods[[3]], fullmods[[4]])

nophyloscalingmods[[2]] <- lme(hor_acel.STD ~ hi.elev + sp.mass.STD + rel.mass.STD + ss.HA.STD, random=~1|species, data=subset(maneuv, !is.na(mass) & !is.na(wing.area)), control=lmeControl(opt='optim'))
nophylomods[[2]] <- lme(hor_acel.STD ~ hi.elev + sp.muscle.STD + sp.mass.STD + sp.wing.load.STD + sp.AR.STD + rel.mass.STD  + rel.wingload.gcm.STD + rel.AR.STD + ss.HA.STD, random=~1|species, data=subset(maneuv, !is.na(mass) & !is.na(wing.area)), control=lmeControl(opt='optim'))

# DecHor,max
set.seed(3112)
scalingmods[[5]] <- MCMCglmm(-hor_decel.STD ~ hi.elev + sp.mass.STD + rel.mass.STD + ss.HD.STD, random=~species+sp, family='gaussian', ginverse=list(species=inv.phylo$Ainv), prior=priorlist, data=subset(maneuv, !is.na(mass)), nitt=nitt, burnin=burnin, thin=thin); alarm()
scalingmods[[6]] <- MCMCglmm(-hor_decel.STD ~ hi.elev + sp.mass.STD + rel.mass.STD + ss.HD.STD, random=~species+sp, family='gaussian', ginverse=list(species=inv.phylo$Ainv), prior=priorlist, data=subset(maneuv, !is.na(mass)), nitt=nitt, burnin=burnin, thin=thin); alarm()
mcmccheck(scalingmods[[5]], scalingmods[[6]])

set.seed(3212)
fullmods[[5]] <- MCMCglmm(-hor_decel.STD ~ hi.elev + sp.muscle.STD + sp.mass.STD + sp.wing.load.STD + sp.AR.STD + rel.mass.STD  + rel.wingload.gcm.STD + rel.AR.STD + ss.HD.STD, random=~species+sp, family="gaussian", ginverse=list(species=inv.phylo$Ainv), prior=priorlist, data=subset(maneuv, !is.na(mass) & !is.na(wing.area)), nitt=nitt, burnin=burnin, thin=thin); alarm()
fullmods[[6]] <- MCMCglmm(-hor_decel.STD ~ hi.elev + sp.muscle.STD + sp.mass.STD + sp.wing.load.STD + sp.AR.STD + rel.mass.STD  + rel.wingload.gcm.STD + rel.AR.STD + ss.HD.STD, random=~species+sp, family="gaussian", ginverse=list(species=inv.phylo$Ainv), prior=priorlist, data=subset(maneuv, !is.na(mass) & !is.na(wing.area)), nitt=nitt, burnin=burnin, thin=thin); alarm()
mcmccheck(fullmods[[5]], fullmods[[6]])

nophyloscalingmods[[3]] <- lme(-hor_decel.STD ~ hi.elev + sp.mass.STD + rel.mass.STD + ss.HD.STD, random=~1|species, data=subset(maneuv, !is.na(mass) & !is.na(wing.area)), control=lmeControl(opt='optim'))
nophylomods[[3]] <- lme(-hor_decel.STD ~ hi.elev + sp.muscle.STD + sp.mass.STD + sp.wing.load.STD + sp.AR.STD + rel.mass.STD  + rel.wingload.gcm.STD + rel.AR.STD + ss.HD.STD, random=~1|species, data=subset(maneuv, !is.na(mass) & !is.na(wing.area)), control=lmeControl(opt='optim'))

# PitchU,vel,avg
set.seed(4112)
scalingmods[[7]] <- MCMCglmm(pitch_up.STD ~ hi.elev + sp.mass.STD + rel.mass.STD + ss.PU.STD, random=~species+sp, family='gaussian', ginverse=list(species=inv.phylo$Ainv), prior=priorlist, data=subset(maneuv, !is.na(mass)), nitt=nitt, burnin=burnin, thin=thin); alarm()
scalingmods[[8]] <- MCMCglmm(pitch_up.STD ~ hi.elev + sp.mass.STD + rel.mass.STD + ss.PU.STD, random=~species+sp, family='gaussian', ginverse=list(species=inv.phylo$Ainv), prior=priorlist, data=subset(maneuv, !is.na(mass)), nitt=nitt, burnin=burnin, thin=thin); alarm()
mcmccheck(scalingmods[[7]], scalingmods[[8]])

set.seed(4212)
fullmods[[7]] <- MCMCglmm(pitch_up.STD ~ hi.elev + sp.muscle.STD + sp.mass.STD + sp.wing.load.STD + sp.AR.STD + rel.mass.STD  + rel.wingload.gcm.STD + rel.AR.STD + ss.PU.STD, random=~species+sp, family="gaussian", ginverse=list(species=inv.phylo$Ainv), prior=priorlist, data=subset(maneuv, !is.na(mass) & !is.na(wing.area)), nitt=nitt, burnin=burnin, thin=thin); alarm()
fullmods[[8]] <- MCMCglmm(pitch_up.STD ~ hi.elev + sp.muscle.STD + sp.mass.STD + sp.wing.load.STD + sp.AR.STD + rel.mass.STD  + rel.wingload.gcm.STD + rel.AR.STD + ss.PU.STD, random=~species+sp, family="gaussian", ginverse=list(species=inv.phylo$Ainv), prior=priorlist, data=subset(maneuv, !is.na(mass) & !is.na(wing.area)), nitt=nitt, burnin=burnin, thin=thin); alarm()
mcmccheck(fullmods[[7]], fullmods[[8]])

nophyloscalingmods[[4]] <- lme(pitch_up.STD ~ hi.elev + sp.mass.STD + rel.mass.STD + ss.PU.STD, random=~1|species, data=subset(maneuv, !is.na(mass) & !is.na(wing.area)), control=lmeControl(opt='optim'))
nophylomods[[4]] <- lme(pitch_up.STD ~ hi.elev + sp.muscle.STD + sp.mass.STD + sp.wing.load.STD + sp.AR.STD + rel.mass.STD  + rel.wingload.gcm.STD + rel.AR.STD + ss.PU.STD, random=~1|species, data=subset(maneuv, !is.na(mass) & !is.na(wing.area)), control=lmeControl(opt='optim'))

# PitchD,vel,avg
set.seed(5112)
scalingmods[[9]] <- MCMCglmm(-pitch_down.STD ~ hi.elev + sp.mass.STD + rel.mass.STD + ss.PD.STD, random=~species+sp, family='gaussian', ginverse=list(species=inv.phylo$Ainv), prior=priorlist, data=subset(maneuv, !is.na(mass)), nitt=nitt, burnin=burnin, thin=thin); alarm()
scalingmods[[10]] <- MCMCglmm(-pitch_down.STD ~ hi.elev + sp.mass.STD + rel.mass.STD + ss.PD.STD, random=~species+sp, family='gaussian', ginverse=list(species=inv.phylo$Ainv), prior=priorlist, data=subset(maneuv, !is.na(mass)), nitt=nitt, burnin=burnin, thin=thin); alarm()
mcmccheck(scalingmods[[9]], scalingmods[[10]])

set.seed(5212)
fullmods[[9]] <- MCMCglmm(-pitch_down.STD ~ hi.elev + sp.muscle.STD + sp.mass.STD + sp.wing.load.STD + sp.AR.STD + rel.mass.STD  + rel.wingload.gcm.STD + rel.AR.STD + ss.PD.STD, random=~species+sp, family="gaussian", ginverse=list(species=inv.phylo$Ainv), prior=priorlist, data=subset(maneuv, !is.na(mass) & !is.na(wing.area)), nitt=nitt, burnin=burnin, thin=thin); alarm()
fullmods[[10]] <- MCMCglmm(-pitch_down.STD ~ hi.elev + sp.muscle.STD + sp.mass.STD + sp.wing.load.STD + sp.AR.STD + rel.mass.STD  + rel.wingload.gcm.STD + rel.AR.STD + ss.PD.STD, random=~species+sp, family="gaussian", ginverse=list(species=inv.phylo$Ainv), prior=priorlist, data=subset(maneuv, !is.na(mass) & !is.na(wing.area)), nitt=nitt, burnin=burnin, thin=thin); alarm()
mcmccheck(fullmods[[9]], fullmods[[10]])

nophyloscalingmods[[5]] <- lme(-pitch_down.STD ~ hi.elev + sp.mass.STD + rel.mass.STD + ss.PD.STD, random=~1|species, data=subset(maneuv, !is.na(mass) & !is.na(wing.area)), control=lmeControl(opt='optim'))
nophylomods[[5]] <- lme(-pitch_down.STD ~ hi.elev + sp.muscle.STD + sp.mass.STD + sp.wing.load.STD + sp.AR.STD + rel.mass.STD  + rel.wingload.gcm.STD + rel.AR.STD + ss.PD.STD, random=~1|species, data=subset(maneuv, !is.na(mass) & !is.na(wing.area)), control=lmeControl(opt='optim'))

# Yaw,vel,avg
set.seed(6112)
scalingmods[[11]] <- MCMCglmm(yaw.STD ~ hi.elev + sp.mass.STD + rel.mass.STD + ss.Y.STD, random=~species+sp, family='gaussian', ginverse=list(species=inv.phylo$Ainv), prior=priorlist, data=subset(maneuv, !is.na(mass)), nitt=nitt, burnin=burnin, thin=thin); alarm()
scalingmods[[12]] <- MCMCglmm(yaw.STD ~ hi.elev + sp.mass.STD + rel.mass.STD + ss.Y.STD, random=~species+sp, family='gaussian', ginverse=list(species=inv.phylo$Ainv), prior=priorlist, data=subset(maneuv, !is.na(mass)), nitt=nitt, burnin=burnin, thin=thin); alarm()
mcmccheck(scalingmods[[11]], scalingmods[[12]])

set.seed(6212)
fullmods[[11]] <- MCMCglmm(yaw.STD ~ hi.elev + sp.muscle.STD + sp.mass.STD + sp.wing.load.STD + sp.AR.STD + rel.mass.STD  + rel.wingload.gcm.STD + rel.AR.STD + ss.Y.STD, random=~species+sp, family="gaussian", ginverse=list(species=inv.phylo$Ainv), prior=priorlist, data=subset(maneuv, !is.na(mass) & !is.na(wing.area)), nitt=nitt, burnin=burnin, thin=thin); alarm()
fullmods[[12]] <- MCMCglmm(yaw.STD ~ hi.elev + sp.muscle.STD + sp.mass.STD + sp.wing.load.STD + sp.AR.STD + rel.mass.STD  + rel.wingload.gcm.STD + rel.AR.STD + ss.Y.STD, random=~species+sp, family="gaussian", ginverse=list(species=inv.phylo$Ainv), prior=priorlist, data=subset(maneuv, !is.na(mass) & !is.na(wing.area)), nitt=nitt, burnin=burnin, thin=thin); alarm()
mcmccheck(fullmods[[11]], fullmods[[12]])

nophyloscalingmods[[6]] <- lme(yaw.STD ~ hi.elev + sp.mass.STD + rel.mass.STD + ss.Y.STD, random=~1|species, data=subset(maneuv, !is.na(mass) & !is.na(wing.area)), control=lmeControl(opt='optim'))
nophylomods[[6]] <- lme(yaw.STD ~ hi.elev + sp.muscle.STD + sp.mass.STD + sp.wing.load.STD + sp.AR.STD + rel.mass.STD  + rel.wingload.gcm.STD + rel.AR.STD + ss.Y.STD, random=~1|species, data=subset(maneuv, !is.na(mass) & !is.na(wing.area)), control=lmeControl(opt='optim'))

# Arc,radius
set.seed(7112)
scalingmods[[13]] <- MCMCglmm(arc_radius.STD ~ hi.elev + sp.mass.STD + rel.mass.STD + ss.ARC.STD, random=~species+sp, family='gaussian', ginverse=list(species=inv.phylo$Ainv), prior=priorlist, data=subset(maneuv, !is.na(mass) & !is.na(arc_radius)), nitt=nitt, burnin=burnin, thin=thin); alarm()
scalingmods[[14]] <- MCMCglmm(arc_radius.STD ~ hi.elev + sp.mass.STD + rel.mass.STD + ss.ARC.STD, random=~species+sp, family='gaussian', ginverse=list(species=inv.phylo$Ainv), prior=priorlist, data=subset(maneuv, !is.na(mass) & !is.na(arc_radius)), nitt=nitt, burnin=burnin, thin=thin); alarm()
mcmccheck(scalingmods[[13]], scalingmods[[14]])

set.seed(7212)
fullmods[[13]] <- MCMCglmm(arc_radius.STD ~ hi.elev + sp.muscle.STD + sp.mass.STD + sp.wing.load.STD + sp.AR.STD + rel.mass.STD  + rel.wingload.gcm.STD + rel.AR.STD + ss.ARC.STD, random=~species+sp, family="gaussian", ginverse=list(species=inv.phylo$Ainv), prior=priorlist, data=subset(maneuv, !is.na(mass) & !is.na(wing.area) & !is.na(arc_radius)), nitt=nitt, burnin=burnin, thin=thin); alarm()
fullmods[[14]] <- MCMCglmm(arc_radius.STD ~ hi.elev + sp.muscle.STD + sp.mass.STD + sp.wing.load.STD + sp.AR.STD + rel.mass.STD  + rel.wingload.gcm.STD + rel.AR.STD + ss.ARC.STD, random=~species+sp, family="gaussian", ginverse=list(species=inv.phylo$Ainv), prior=priorlist, data=subset(maneuv, !is.na(mass) & !is.na(wing.area) & !is.na(arc_radius)), nitt=nitt, burnin=burnin, thin=thin); alarm()
mcmccheck(fullmods[[13]], fullmods[[14]])

nophyloscalingmods[[7]] <- lme(arc_radius.STD ~ hi.elev + sp.mass.STD + rel.mass.STD + ss.ARC.STD, random=~1|species, data=subset(maneuv, !is.na(mass) & !is.na(wing.area) & !is.na(arc_radius)), control=lmeControl(opt='optim'))
nophylomods[[7]] <- lme(arc_radius.STD ~ hi.elev + sp.muscle.STD + sp.mass.STD + sp.wing.load.STD + sp.AR.STD + rel.mass.STD  + rel.wingload.gcm.STD + rel.AR.STD + ss.ARC.STD, random=~1|species, data=subset(maneuv, !is.na(mass) & !is.na(wing.area) & !is.na(arc_radius)), control=lmeControl(opt='optim'))

# Arc,vel,avg
set.seed(8112)
scalingmods[[15]] <- MCMCglmm(arc_avg_vel.STD ~ hi.elev + sp.mass.STD + rel.mass.STD + ss.ARC.STD, random=~species+sp, family='gaussian', ginverse=list(species=inv.phylo$Ainv), prior=priorlist, data=subset(maneuv, !is.na(mass) & !is.na(arc_radius)), nitt=nitt, burnin=burnin, thin=thin); alarm()
scalingmods[[16]] <- MCMCglmm(arc_avg_vel.STD ~ hi.elev + sp.mass.STD + rel.mass.STD + ss.ARC.STD, random=~species+sp, family='gaussian', ginverse=list(species=inv.phylo$Ainv), prior=priorlist, data=subset(maneuv, !is.na(mass) & !is.na(arc_radius)), nitt=nitt, burnin=burnin, thin=thin); alarm()
mcmccheck(scalingmods[[15]], scalingmods[[16]])

set.seed(8212)
fullmods[[15]] <- MCMCglmm(arc_avg_vel.STD ~ hi.elev + sp.muscle.STD + sp.mass.STD + sp.wing.load.STD + sp.AR.STD + rel.mass.STD  + rel.wingload.gcm.STD + rel.AR.STD + ss.ARC.STD, random=~species+sp, family="gaussian", ginverse=list(species=inv.phylo$Ainv), prior=priorlist, data=subset(maneuv, !is.na(mass) & !is.na(wing.area) & !is.na(arc_radius)), nitt=nitt, burnin=burnin, thin=thin); alarm()
fullmods[[16]] <- MCMCglmm(arc_avg_vel.STD ~ hi.elev + sp.muscle.STD + sp.mass.STD + sp.wing.load.STD + sp.AR.STD + rel.mass.STD  + rel.wingload.gcm.STD + rel.AR.STD + ss.ARC.STD, random=~species+sp, family="gaussian", ginverse=list(species=inv.phylo$Ainv), prior=priorlist, data=subset(maneuv, !is.na(mass) & !is.na(wing.area) & !is.na(arc_radius)), nitt=nitt, burnin=burnin, thin=thin); alarm()
mcmccheck(fullmods[[15]], fullmods[[16]])

nophyloscalingmods[[8]] <- lme(arc_avg_vel.STD ~ hi.elev + sp.mass.STD + rel.mass.STD + ss.ARC.STD, random=~1|species, data=subset(maneuv, !is.na(mass) & !is.na(wing.area) & !is.na(arc_radius)), control=lmeControl(opt='optim'))
nophylomods[[8]] <- lme(arc_avg_vel.STD ~ hi.elev + sp.muscle.STD + sp.mass.STD + sp.wing.load.STD + sp.AR.STD + rel.mass.STD  + rel.wingload.gcm.STD + rel.AR.STD + ss.ARC.STD, random=~1|species, data=subset(maneuv, !is.na(mass) & !is.na(wing.area) & !is.na(arc_radius)), control=lmeControl(opt='optim'))

# Arc,cent,max
set.seed(9112)
scalingmods[[17]] <- MCMCglmm(arc_force.STD ~ hi.elev + sp.mass.STD + rel.mass.STD + ss.ARC.STD, random=~species+sp, family='gaussian', ginverse=list(species=inv.phylo$Ainv), prior=priorlist, data=subset(maneuv, !is.na(mass) & !is.na(arc_radius)), nitt=nitt, burnin=burnin, thin=thin); alarm()
scalingmods[[18]] <- MCMCglmm(arc_force.STD ~ hi.elev + sp.mass.STD + rel.mass.STD + ss.ARC.STD, random=~species+sp, family='gaussian', ginverse=list(species=inv.phylo$Ainv), prior=priorlist, data=subset(maneuv, !is.na(mass) & !is.na(arc_radius)), nitt=nitt, burnin=burnin, thin=thin); alarm()
mcmccheck(scalingmods[[17]], scalingmods[[18]])

set.seed(9212)
fullmods[[17]] <- MCMCglmm(arc_force.STD ~ hi.elev + sp.muscle.STD + sp.mass.STD + sp.wing.load.STD + sp.AR.STD + rel.mass.STD  + rel.wingload.gcm.STD + rel.AR.STD + ss.ARC.STD, random=~species+sp, family="gaussian", ginverse=list(species=inv.phylo$Ainv), prior=priorlist, data=subset(maneuv, !is.na(mass) & !is.na(wing.area) & !is.na(arc_radius)), nitt=nitt, burnin=burnin, thin=thin); alarm()
fullmods[[18]] <- MCMCglmm(arc_force.STD ~ hi.elev + sp.muscle.STD + sp.mass.STD + sp.wing.load.STD + sp.AR.STD + rel.mass.STD  + rel.wingload.gcm.STD + rel.AR.STD + ss.ARC.STD, random=~species+sp, family="gaussian", ginverse=list(species=inv.phylo$Ainv), prior=priorlist, data=subset(maneuv, !is.na(mass) & !is.na(wing.area) & !is.na(arc_radius)), nitt=nitt, burnin=burnin, thin=thin); alarm()
mcmccheck(fullmods[[17]], fullmods[[18]])

nophyloscalingmods[[9]] <- lme(arc_force.STD ~ hi.elev + sp.mass.STD + rel.mass.STD + ss.ARC.STD, random=~1|species, data=subset(maneuv, !is.na(mass) & !is.na(wing.area) & !is.na(arc_radius)), control=lmeControl(opt='optim'))
nophylomods[[9]] <- lme(arc_force.STD ~ hi.elev + sp.muscle.STD + sp.mass.STD + sp.wing.load.STD + sp.AR.STD + rel.mass.STD  + rel.wingload.gcm.STD + rel.AR.STD + ss.ARC.STD, random=~1|species, data=subset(maneuv, !is.na(mass) & !is.na(wing.area) & !is.na(arc_radius)), control=lmeControl(opt='optim'))

# PRT,degrees
set.seed(10112)
scalingmods[[19]] <- MCMCglmm(PRT_degrees.STD ~ hi.elev + sp.mass.STD + rel.mass.STD + ss.PRT.STD, random=~species+sp, family='gaussian', ginverse=list(species=inv.phylo$Ainv), prior=priorlist, data=subset(maneuv, !is.na(mass)), nitt=nitt, burnin=burnin, thin=thin); alarm()
scalingmods[[20]] <- MCMCglmm(PRT_degrees.STD ~ hi.elev + sp.mass.STD + rel.mass.STD + ss.PRT.STD, random=~species+sp, family='gaussian', ginverse=list(species=inv.phylo$Ainv), prior=priorlist, data=subset(maneuv, !is.na(mass)), nitt=nitt, burnin=burnin, thin=thin); alarm()
mcmccheck(scalingmods[[19]], scalingmods[[20]])

set.seed(10212)
fullmods[[19]] <- MCMCglmm(PRT_degrees.STD ~ hi.elev + sp.muscle.STD + sp.mass.STD + sp.wing.load.STD + sp.AR.STD + rel.mass.STD  + rel.wingload.gcm.STD + rel.AR.STD + ss.PRT.STD, random=~species+sp, family="gaussian", ginverse=list(species=inv.phylo$Ainv), prior=priorlist, data=subset(maneuv, !is.na(mass) & !is.na(wing.area)), nitt=nitt, burnin=burnin, thin=thin); alarm()
fullmods[[20]] <- MCMCglmm(PRT_degrees.STD ~ hi.elev + sp.muscle.STD + sp.mass.STD + sp.wing.load.STD + sp.AR.STD + rel.mass.STD  + rel.wingload.gcm.STD + rel.AR.STD + ss.PRT.STD, random=~species+sp, family="gaussian", ginverse=list(species=inv.phylo$Ainv), prior=priorlist, data=subset(maneuv, !is.na(mass) & !is.na(wing.area)), nitt=nitt, burnin=burnin, thin=thin); alarm()
mcmccheck(fullmods[[19]], fullmods[[20]])

nophyloscalingmods[[10]] <- lme(PRT_degrees.STD ~ hi.elev + sp.mass.STD + rel.mass.STD + ss.PRT.STD, random=~1|species, data=subset(maneuv, !is.na(mass) & !is.na(wing.area)), control=lmeControl(opt='optim'))
nophylomods[[10]] <- lme(PRT_degrees.STD ~ hi.elev + sp.muscle.STD + sp.mass.STD + sp.wing.load.STD + sp.AR.STD + rel.mass.STD  + rel.wingload.gcm.STD + rel.AR.STD + ss.PRT.STD, random=~1|species, data=subset(maneuv, !is.na(mass) & !is.na(wing.area)), control=lmeControl(opt='optim'))

# PRT,time
set.seed(11112)
scalingmods[[21]] <- MCMCglmm(-PRT_time.STD ~ hi.elev + sp.mass.STD + rel.mass.STD + ss.PRT.STD, random=~species+sp, family='gaussian', ginverse=list(species=inv.phylo$Ainv), prior=priorlist, data=subset(maneuv, !is.na(mass)), nitt=nitt, burnin=burnin, thin=thin); alarm()
scalingmods[[22]] <- MCMCglmm(-PRT_time.STD ~ hi.elev + sp.mass.STD + rel.mass.STD + ss.PRT.STD, random=~species+sp, family='gaussian', ginverse=list(species=inv.phylo$Ainv), prior=priorlist, data=subset(maneuv, !is.na(mass)), nitt=nitt, burnin=burnin, thin=thin); alarm()
mcmccheck(scalingmods[[21]], scalingmods[[22]])

set.seed(11212)
fullmods[[21]] <- MCMCglmm(-PRT_time.STD ~ hi.elev + sp.muscle.STD + sp.mass.STD + sp.wing.load.STD + sp.AR.STD + rel.mass.STD  + rel.wingload.gcm.STD + rel.AR.STD + ss.PRT.STD, random=~species+sp, family="gaussian", ginverse=list(species=inv.phylo$Ainv), prior=priorlist, data=subset(maneuv, !is.na(mass) & !is.na(wing.area)), nitt=nitt, burnin=burnin, thin=thin); alarm()
fullmods[[22]] <- MCMCglmm(-PRT_time.STD ~ hi.elev + sp.muscle.STD + sp.mass.STD + sp.wing.load.STD + sp.AR.STD + rel.mass.STD  + rel.wingload.gcm.STD + rel.AR.STD + ss.PRT.STD, random=~species+sp, family="gaussian", ginverse=list(species=inv.phylo$Ainv), prior=priorlist, data=subset(maneuv, !is.na(mass) & !is.na(wing.area)), nitt=nitt, burnin=burnin, thin=thin); alarm()
mcmccheck(fullmods[[21]], fullmods[[22]])

nophyloscalingmods[[11]] <- lme(-PRT_time.STD ~ hi.elev + sp.mass.STD + rel.mass.STD + ss.PRT.STD, random=~1|species, data=subset(maneuv, !is.na(mass) & !is.na(wing.area)), control=lmeControl(opt='optim'))
nophylomods[[11]] <- lme(-PRT_time.STD ~ hi.elev + sp.muscle.STD + sp.mass.STD + sp.wing.load.STD + sp.AR.STD + rel.mass.STD  + rel.wingload.gcm.STD + rel.AR.STD + ss.PRT.STD, random=~1|species, data=subset(maneuv, !is.na(mass) & !is.na(wing.area)), control=lmeControl(opt='optim'))

# %PRT
set.seed(12112)
scalingmods[[23]] <- MCMCglmm(pPRT.STD ~ hi.elev + sp.mass.STD + rel.mass.STD, random=~species+sp, family='gaussian', ginverse=list(species=inv.phylo$Ainv), prior=priorlist, data=subset(maneuv, !is.na(mass)), nitt=nitt, burnin=burnin, thin=thin); alarm()
scalingmods[[24]] <- MCMCglmm(pPRT.STD ~ hi.elev + sp.mass.STD + rel.mass.STD, random=~species+sp, family='gaussian', ginverse=list(species=inv.phylo$Ainv), prior=priorlist, data=subset(maneuv, !is.na(mass)), nitt=nitt, burnin=burnin, thin=thin); alarm()
mcmccheck(scalingmods[[23]], scalingmods[[24]])

set.seed(12212)
fullmods[[23]] <- MCMCglmm(pPRT.STD ~ hi.elev + sp.muscle.STD + sp.mass.STD + sp.wing.load.STD + sp.AR.STD + rel.mass.STD  + rel.wingload.gcm.STD + rel.AR.STD, random=~species+sp, family="gaussian", ginverse=list(species=inv.phylo$Ainv), prior=priorlist, data=subset(maneuv, !is.na(mass) & !is.na(wing.area)), nitt=nitt, burnin=burnin, thin=thin); alarm()
fullmods[[24]] <- MCMCglmm(pPRT.STD ~ hi.elev + sp.muscle.STD + sp.mass.STD + sp.wing.load.STD + sp.AR.STD + rel.mass.STD  + rel.wingload.gcm.STD + rel.AR.STD, random=~species+sp, family="gaussian", ginverse=list(species=inv.phylo$Ainv), prior=priorlist, data=subset(maneuv, !is.na(mass) & !is.na(wing.area)), nitt=nitt, burnin=burnin, thin=thin); alarm()
mcmccheck(fullmods[[23]], fullmods[[24]])

nophyloscalingmods[[12]] <- lme(pPRT.STD ~ hi.elev + sp.mass.STD + rel.mass.STD, random=~1|species, data=subset(maneuv, !is.na(mass) & !is.na(wing.area)), control=lmeControl(opt='optim'))
nophylomods[[12]] <- lme(pPRT.STD ~ hi.elev + sp.muscle.STD + sp.mass.STD + sp.wing.load.STD + sp.AR.STD + rel.mass.STD  + rel.wingload.gcm.STD + rel.AR.STD, random=~1|species, data=subset(maneuv, !is.na(mass) & !is.na(wing.area)), control=lmeControl(opt='optim'))

# save.image('comparative_maneuverability_models.RData')

# compile the trait effect sizes

mcmc.posterior <- data.frame(matrix(NA, nrow=12, ncol=13))
names(mcmc.posterior) <- c(paste('SCALING',rownames(summary(scalingmods[[1]])$solutions[2:5,1:3]),sep='.'), paste('FULL',rownames(summary(fullmods[[1]])$solutions[2:10,1:3]),sep='.'))
names(mcmc.posterior)[4] <- 'SCALING.ss.STD'
names(mcmc.posterior)[13] <- 'FULL.ss.STD'
mcmc.posterior$metric <- names(maneuv)[c(16:27)]; mcmc.posterior <- mcmc.posterior[,c(14,1:13)]
mcmc.lower <- mcmc.posterior; mcmc.upper <- mcmc.posterior
nophylo <- mcmc.posterior; nophylo.lower <- mcmc.posterior; nophylo.upper <- mcmc.posterior

for(i in 1:12){
	if(i < 12){
		mcmc.posterior[i, 2:5] <- unname(summary(scalingmods[[i*2-1]])$solutions[2:5,1])
		mcmc.posterior[i, 6:14] <- unname(summary(fullmods[[i*2-1]])$solutions[2:10,1])
		mcmc.lower[i, 2:5] <- unname(summary(scalingmods[[i*2-1]])$solutions[2:5,2])
		mcmc.lower[i, 6:14] <- unname(summary(fullmods[[i*2-1]])$solutions[2:10,2])
		mcmc.upper[i, 2:5] <- unname(summary(scalingmods[[i*2-1]])$solutions[2:5,3])
		mcmc.upper[i, 6:14] <- unname(summary(fullmods[[i*2-1]])$solutions[2:10,3])		
		nophylo[i, 2:5] <- unname(summary(nophyloscalingmods[[i]])$tTable[2:5,'Value'])
		nophylo[i, 6:14] <- unname(summary(nophylomods[[i]])$tTable[2:10,'Value'])
		nophylo.lower[i, 2:5] <- unname(summary(nophyloscalingmods[[i]])$tTable[2:5,'Value']) - 1.96*unname(summary(nophyloscalingmods[[i]])$tTable[2:5,'Std.Error'])
		nophylo.lower[i, 6:14] <- unname(summary(nophylomods[[i]])$tTable[2:10,'Value']) - 1.96*unname(summary(nophylomods[[i]])$tTable[2:10,'Std.Error'])
		nophylo.upper[i, 2:5] <- unname(summary(nophyloscalingmods[[i]])$tTable[2:5,'Value']) + 1.96*unname(summary(nophyloscalingmods[[i]])$tTable[2:5,'Std.Error'])
		nophylo.upper[i, 6:14] <- unname(summary(nophylomods[[i]])$tTable[2:10,'Value']) + 1.96*unname(summary(nophylomods[[i]])$tTable[2:10,'Std.Error'])		
	} else {
		mcmc.posterior[i, 2:4] <- unname(summary(scalingmods[[i*2-1]])$solutions[2:4,1])
		mcmc.posterior[i, 6:13] <- unname(summary(fullmods[[i*2-1]])$solutions[2:9,1])
		mcmc.lower[i, 2:4] <- unname(summary(scalingmods[[i*2-1]])$solutions[2:4,2])
		mcmc.lower[i, 6:13] <- unname(summary(fullmods[[i*2-1]])$solutions[2:9,2])
		mcmc.upper[i, 2:4] <- unname(summary(scalingmods[[i*2-1]])$solutions[2:4,3])
		mcmc.upper[i, 6:13] <- unname(summary(fullmods[[i*2-1]])$solutions[2:9,3])
		nophylo[i, 2:4] <- unname(summary(nophyloscalingmods[[i]])$tTable[2:4,'Value'])
		nophylo[i, 6:13] <- unname(summary(nophylomods[[i]])$tTable[2:9,'Value'])
		nophylo.lower[i, 2:4] <- unname(summary(nophyloscalingmods[[i]])$tTable[2:4,'Value']) - 1.96*unname(summary(nophyloscalingmods[[i]])$tTable[2:4,'Std.Error'])
		nophylo.lower[i, 6:13] <- unname(summary(nophylomods[[i]])$tTable[2:9,'Value']) - 1.96*unname(summary(nophylomods[[i]])$tTable[2:9,'Std.Error'])
		nophylo.upper[i, 2:4] <- unname(summary(nophyloscalingmods[[i]])$tTable[2:4,'Value']) + 1.96*unname(summary(nophyloscalingmods[[i]])$tTable[2:4,'Std.Error'])		
		nophylo.upper[i, 6:13] <- unname(summary(nophylomods[[i]])$tTable[2:9,'Value']) + 1.96*unname(summary(nophylomods[[i]])$tTable[2:9,'Std.Error'])
	}
}

# Figure S7 plot all effect sizes:

dev.new(width=11, height=6)
par(mfcol=c(3,4), mar=c(4,4,2,0.25), mgp=c(2,0.5,0), las=1, cex=0.5, bty='l')
for(i in 1:12){
	plot(unlist(mcmc.posterior[i,3:4]) ~ c(3-0.3,4-0.3), xlim=c(0,9), ylim=c(-1,1), pch=16, col=1, ylab='effect size', xlab='', xaxt='n', main=mcmc.posterior$metric[i]); abline(h=0, lty=3)
	axis(1, at=c(1:9), labels=F)
	text(x=c(1:9), y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]), labels=gsub('.STD','',gsub('FULL.','',names(mcmc.posterior)[c(6,14,8,11,7,9,12,10,13)])), srt=45, adj=1, xpd=TRUE, cex=.75)
	points(unlist(mcmc.posterior[i,c(6,14,8,11,7,9,12,10,13)]) ~ c(1:9), pch=16, col=1)
	segments(x0=1:9, y0=unlist(mcmc.lower[i,c(6,14,8,11,7,9,12,10,13)]), y1=unlist(mcmc.upper[i,c(6,14,8,11,7,9,12,10,13)]), col=1)
	segments(x0=c(3-0.3,4-0.3), y0=unlist(mcmc.lower[i,3:4]), y1=unlist(mcmc.upper[i,3:4]), col=1)
	points(unlist(nophylo[i,7:10]) ~ c(c(5,3,6,8)+0.15), col='green', pch=16)
	segments(x0=c(5,3,6,8)+0.15, y0=unlist(nophylo.lower[i,7:10]), y1=unlist(nophylo.upper[i,7:10]), col='green')
	points(unlist(nophylo[i,3:4]) ~ c(c(3-0.3,4-0.3)+0.15), col='green', pch=16)
	segments(x0=c(3-0.3,4-0.3)+0.15, y0=unlist(nophylo.lower[i,3:4]), y1=unlist(nophylo.upper[i,3:4]), col='green')
}

for(i in 1:12){
	print(r.squaredGLMM(nophylomods[[i]])[1])
}

# compile elevation results:

eleveffects <- cbind(mcmc.posterior[1:12,c(1,6)],mcmc.lower[1:12,6],mcmc.upper[1:12,6])
eleveffects$lower <- ifelse(eleveffects$FULL.hi.elevTRUE < 0, -eleveffects[,3], eleveffects[,3])
eleveffects$upper <- ifelse(eleveffects$FULL.hi.elevTRUE < 0, -eleveffects[,4], eleveffects[,4])

# compile maneuver sample size results:

sseffects <- cbind(mcmc.posterior[1:11,c(1,14)],mcmc.lower[1:11,14],mcmc.upper[1:11,14])
sseffects <- sseffects[with(sseffects,order(abs(FULL.ss.STD))),]
sseffects$lower <- ifelse(sseffects$FULL.ss.STD < 0, -sseffects[,3], sseffects[,3])
sseffects$upper <- ifelse(sseffects$FULL.ss.STD < 0, -sseffects[,4], sseffects[,4])

# Figure S8 compare MCMC and LME results

dev.new(width=5, height=5)
par(mar=c(4,4,0.25,0.25), las=1)
plot(as.vector(as.matrix(mcmc.posterior[,c(11:13)])) ~ as.vector(as.matrix(nophylo[,c(11:13)])), pch=16, cex=0.5, xlab='effect size, LME model (no phylo)', ylab='effect size, MCMC model', ylim=c(-0.8,0.8), xlim=c(-0.8,0.8), asp=1)
points(as.vector(as.matrix(mcmc.posterior[,c(7:10)])) ~ as.vector(as.matrix(nophylo[,c(7:10)])), pch=16, cex=0.5, col='green')
points(as.vector(as.matrix(mcmc.posterior[,c(6)])) ~ as.vector(as.matrix(nophylo[,c(6)])), pch=16, cex=0.5, col='blue')
points(as.vector(as.matrix(mcmc.posterior[,c(14)])) ~ as.vector(as.matrix(nophylo[,c(14)])), pch=16, cex=0.5, col='red')
abline(a=0, b=1, lty=3); abline(h=0,v=0, lty=2)
legend('topleft', pch=16, cex=0.5, col=c('white','green','black','red','blue'), legend=c('legend','species average trait','individual trait, relative to the species average','number of maneuvers','elevation'), bty='n')

# save.image('comparative_maneuverability_models.RData')

# plot body mass results 

# refit LME models on the scale of the original data for plotting
nophyloscalingmods[[14]] <- lme(hor_acel ~ hi.elev + sp.mass + rel.mass + ss.HA, random=~1|species, data=subset(maneuv, !is.na(mass) & !is.na(wing.area)), control=lmeControl(opt='optim'))
nophyloscalingmods[[15]] <- lme(-hor_decel ~ hi.elev + sp.mass + rel.mass + ss.HD, random=~1|species, data=subset(maneuv, !is.na(mass) & !is.na(wing.area)), control=lmeControl(opt='optim'))
nophyloscalingmods[[21]] <- lme(arc_force ~ hi.elev + sp.mass + rel.mass + ss.ARC, random=~1|species, data=subset(maneuv, !is.na(mass) & !is.na(wing.area) & !is.na(arc_radius)), control=lmeControl(opt='optim'))
plottable <- read.csv('plottable.csv')

dev.new(width=7, height=2)
par(mfcol=c(1,3), mar=c(4,4,.25,.25), mgp=c(1.75,0.5,0), las=1)
plot(hor_acel ~ mass, data=maneuv, pch=plottable$pch[match(maneuv$species, plottable$species)], cex=0.45, col=as.character(color), xlim=c(2,10.5), ylim=c(2,9), type='p')
polygon(x=c(2,10.5,10.5,2), y=c(2,2,9,9), border=NA, col=rgb(1,1,1,0.6))
for(i in 1:17){
	mysamp <- subset(maneuv, !is.na(mass) & species==rownames(table(maneuv$species)[table(maneuv$species)>4])[i])
	mymod <- lm(hor_acel ~ mass, data=mysamp)
	segments(x0=c(min(mysamp$mass)), y0=predict(mymod, newdata=data.frame(mass=min(mysamp$mass))), x1=c(max(mysamp$mass)), y1=predict(mymod, newdata=data.frame(mass=max(mysamp$mass))), col=as.character(mysamp$color[1]), lwd=1.25, lty=3)
	points(mean(mysamp$hor_acel) ~ mean(mysamp$mass), pch=plottable$pch[match(mysamp$species[1],plottable$species)], col=as.character(mysamp$color[1]), cex=1.25)
}
points(predict(nophyloscalingmods[[14]], newdata=data.frame(hi.elev=c(F,F), sp.mass=range(maneuv$mass,na.rm=T), rel.mass=0, ss.HA=mean(maneuv$ss.HA)), level=0) ~ range(maneuv$mass, na.rm=T), type='l') 
plot(-hor_decel ~ mass, data=maneuv, pch=plottable$pch[match(maneuv$species, plottable$species)], cex=0.45, col=as.character(color), xlim=c(2,10.5), ylim=c(2,9), type='p')
polygon(x=c(2,10.5,10.5,2), y=c(2,2,9,9), border=NA, col=rgb(1,1,1,0.6))
for(i in 1:17){
	mysamp <- subset(maneuv, !is.na(mass) & species==rownames(table(maneuv$species)[table(maneuv$species)>4])[i])
	mymod <- lm(-hor_decel ~ mass, data=mysamp)
	segments(x0=c(min(mysamp$mass)), y0=predict(mymod, newdata=data.frame(mass=min(mysamp$mass))), x1=c(max(mysamp$mass)), y1=predict(mymod, newdata=data.frame(mass=max(mysamp$mass))), col=as.character(mysamp$color[1]), lwd=1.25, lty=3)
	points(mean(-mysamp$hor_decel) ~ mean(mysamp$mass), pch=plottable$pch[match(mysamp$species[1],plottable$species)], col=as.character(mysamp$color[1]), cex=1.25)
}
points(predict(nophyloscalingmods[[15]], newdata=data.frame(hi.elev=c(F,F), sp.mass=range(maneuv$mass,na.rm=T), rel.mass=0, ss.HD=mean(maneuv$ss.HD)), level=0) ~ range(maneuv$mass, na.rm=T), type='l') 
plot(arc_force ~ mass, data=maneuv, pch=plottable$pch[match(maneuv$species, plottable$species)], cex=0.45, col=as.character(color), xlim=c(2,10.5), ylim=c(2,9))
polygon(x=c(2,10.5,10.5,2), y=c(2,2,9,9), border=NA, col=rgb(1,1,1,0.6))
for(i in 1:17){
	mysamp <- subset(maneuv, !is.na(mass) & species==rownames(table(maneuv$species)[table(maneuv$species)>4])[i])
	mymod <- lm(arc_force ~ mass, data=mysamp, na.action=na.omit)
	segments(x0=c(min(mysamp$mass)), y0=predict(mymod, newdata=data.frame(mass=min(mysamp$mass))), x1=c(max(mysamp$mass)), y1=predict(mymod, newdata=data.frame(mass=max(mysamp$mass))), col=as.character(mysamp$color[1]), lwd=1.25, lty=3)
	points(mean(mysamp$arc_force) ~ mean(mysamp$mass), pch=plottable$pch[match(mysamp$species[1],plottable$species)], col=as.character(mysamp$color[1]), cex=1.25)
}
points(predict(nophyloscalingmods[[21]], newdata=data.frame(hi.elev=c(F,F), sp.mass=range(maneuv$mass,na.rm=T), rel.mass=0, ss.ARC=mean(maneuv$ss.ARC, na.rm=T)), level=0) ~ range(maneuv$mass, na.rm=T), type='l') 


# further analysis of the distribution of complex turns, %PRT

# at the species level, we consider elevation, body mass, wing loading, Arccent,max, and PRTtime
# at the individual level, we consider elevation, Arccent,max, and PRTtime

spturns <- aggregate(maneuv[,c('pPRT','sp.mass','sp.wing.load','arc_force','PRT_time')], by=list(hi.elev=maneuv$hi.elev, species=maneuv$species, sp=maneuv$sp), 'mean', na.rm=T)
names(spturns)[c(4,7,8)] <- paste('sp', names(spturns)[c(4,7,8)], sep='.')

maneuvT <- merge(maneuv, spturns[,c('species','sp.pPRT','sp.arc_force','sp.PRT_time')])
maneuvT$rel.pPRT <- maneuvT$pPRT - maneuvT$sp.pPRT
maneuvT$rel.arc_force <- maneuvT$arc_force - maneuvT$sp.arc_force
maneuvT$rel.PRT_time <- maneuvT$PRT_time - maneuvT$sp.PRT_time

PRTmod1 <- lme(pPRT ~ hi.elev + sp.mass + sp.wing.load + rel.AR + sp.arc_force + rel.arc_force, random=~1|species, data=subset(maneuvT, !is.na(mass) & !is.na(arc_force) &!is.na(wing.area)), control=lmeControl(opt='optim')) # no phylogeny
set.seed(13212)
PRTmod2 <- MCMCglmm(pPRT ~ hi.elev + sp.mass + sp.wing.load + rel.AR + sp.arc_force + rel.arc_force, random=~species+sp, family='gaussian', ginverse=list(species=inv.phylo$Ainv), prior=priorlist, data=subset(maneuvT, !is.na(mass) & !is.na(arc_force) &!is.na(wing.area)), nitt=nitt, burnin=burnin, thin=thin); alarm()
PRTmod3 <- MCMCglmm(pPRT ~ hi.elev + sp.mass + sp.wing.load + rel.AR + sp.arc_force + rel.arc_force, random=~species+sp, family='gaussian', ginverse=list(species=inv.phylo$Ainv), prior=priorlist, data=subset(maneuvT, !is.na(mass) & !is.na(arc_force) &!is.na(wing.area)), nitt=nitt, burnin=burnin, thin=thin); alarm()
mcmccheck(PRTmod2, PRTmod3)
summary(PRTmod1)
r.squaredGLMM(PRTmod1)

# In Fig 6B showing this result, there is a high leverage point for the species Eriocnemis luciani with n=1. Here we check that the results are robust to removing this species:

PRTmod1X <- lme(pPRT ~ hi.elev + sp.mass + sp.wing.load + rel.AR + sp.arc_force + rel.arc_force, random=~1|species, data=subset(maneuvT, !is.na(mass) & !is.na(arc_force) &!is.na(wing.area)&species!='Eriocnemis.luciani'), control=lmeControl(opt='optim')) # no phylogeny
set.seed(13212)
PRTmod2X <- MCMCglmm(pPRT ~ hi.elev + sp.mass + sp.wing.load + rel.AR + sp.arc_force + rel.arc_force, random=~species+sp, family='gaussian', ginverse=list(species=inv.phylo$Ainv), prior=priorlist, data=subset(maneuvT, !is.na(mass) & !is.na(arc_force) &!is.na(wing.area)&species!='Eriocnemis.luciani'), nitt=1000000 + burnin, burnin=burnin, thin=thin); alarm()
PRTmod3X <- MCMCglmm(pPRT ~ hi.elev + sp.mass + sp.wing.load + rel.AR + sp.arc_force + rel.arc_force, random=~species+sp, family='gaussian', ginverse=list(species=inv.phylo$Ainv), prior=priorlist, data=subset(maneuvT, !is.na(mass) & !is.na(arc_force) &!is.na(wing.area)&species!='Eriocnemis.luciani'), nitt=1000000 + burnin, burnin=burnin, thin=thin); alarm()
mcmccheck(PRTmod2X, PRTmod3X)
summary(PRTmod1X); plot(PRTmod1X$coef$fixed ~ PRTmod1$coef$fixed, pch=16, xlab='Results in text (all species)', ylab='Results with Eriocnemis luciani removed', las=1, bty='l'); abline(a=0,b=1,lty=3) # fixed effects are unchanged by removing this one species
cor.test(PRTmod1X$coef$fixed, PRTmod1$coef$fixed)
r.squaredGLMM(PRTmod1X) # 0.28, virtually same as before


# write.csv(nophylo, 'lme_model_results.csv')
# write.csv(mcmc.posterior, 'mcmc_model_results.csv')


# single-species models for effects of individual traits:

highsp  <- names(sort(table(subset(maneuv, !is.na(mass) & !is.na(arc_force) & !is.na(wing.area))$species), decreasing=T)[1:8]) # use 8 species with n = 10 or greater (n between 10 and 17)
ns <- as.vector(sort(table(subset(maneuv, !is.na(mass) & !is.na(arc_force) & !is.na(wing.area))$species), decreasing=T))[1:8]

(smallmods <- data.frame(metric=rep(names(maneuv)[c(52,53,58,59,60,62)][c(1:4, 1:2, 5:6)], each=length(highsp)), trait=rep(c(rep('rel.mass.STD',4), rep('rel.wingload.gcm.STD',2), rep('rel.AR.STD',2)), each=length(highsp)), species=highsp, n=ns, coef=NA, p=NA, lower=NA, upper=NA))
smallmods$metric <- ifelse(as.character(smallmods$metric)=='hor_decel.STD', '-hor_decel.STD', as.character(smallmods$metric))
for(i in 1:64){
	mod.single <- lm(as.formula(paste(smallmods$metric[i], '~', smallmods$trait[i])), data=subset(maneuv, species %in% smallmods$species[i]))
	smallmods$coef[i] <- summary(mod.single)$coef[2,1]
	smallmods$p[i] <- summary(mod.single)$coef[2,4]
	smallmods$lower[i] <- summary(mod.single)$coef[2,1] - 1.96*summary(mod.single)$coef[2,2]
	smallmods$upper[i] <- summary(mod.single)$coef[2,1] + 1.96*summary(mod.single)$coef[2,2]
}
smallmods <- merge(smallmods, plottable[,c(2,4,5)])
smallmods <- smallmods[with(smallmods, order(trait, metric, -p)),]

# Figure S10 single-species effect sizes:

dev.new(width=5, height=6)

par(mfrow=c(3,1), mar=c(3,3,0.5,0.25), mgp=c(1,0.5,0), bty='l', las=1)

plot(coef~c(1:8), pch=pch, col=as.character(color), data=subset(smallmods, trait=='rel.mass.STD' & metric=='hor_acel.STD'), ylim=c(-1.5,1.5), xlim=c(1,45), yaxt='n', xaxt='n', ylab='effect size', xlab='', main='body mass', cex.main=0.75, cex=0.5)
segments(x0=c(1:8), y0=subset(smallmods, trait=='rel.mass.STD' & metric=='hor_acel.STD')$lower, y1=subset(smallmods, trait=='rel.mass.STD' & metric=='hor_acel.STD')$upper, col=as.character(subset(smallmods, trait=='rel.mass.STD' & metric=='hor_acel.STD')$color))

points(coef~c(13:20), pch=pch, col=as.character(color), data=subset(smallmods, trait=='rel.mass.STD' & metric=='-hor_decel.STD'), cex=0.5)
segments(x0=c(13:20), y0=subset(smallmods, trait=='rel.mass.STD' & metric=='-hor_decel.STD')$lower, y1=subset(smallmods, trait=='rel.mass.STD' & metric=='-hor_decel.STD')$upper, col=as.character(subset(smallmods, trait=='rel.mass.STD' & metric=='-hor_decel.STD')$color))

points(coef~c(25:32), pch=pch, col=as.character(color), data=subset(smallmods, trait=='rel.mass.STD' & metric=='arc_avg_vel.STD'), cex=0.5)
segments(x0=c(25:32), y0=subset(smallmods, trait=='rel.mass.STD' & metric=='arc_avg_vel.STD')$lower, y1=subset(smallmods, trait=='rel.mass.STD' & metric=='arc_avg_vel.STD')$upper, col=as.character(subset(smallmods, trait=='rel.mass.STD' & metric=='arc_avg_vel.STD')$color))

points(coef~c(37:44), pch=pch, col=as.character(color), data=subset(smallmods, trait=='rel.mass.STD' & metric=='arc_force.STD'), cex=0.5)
segments(x0=c(37:44), y0=subset(smallmods, trait=='rel.mass.STD' & metric=='arc_force.STD')$lower, y1=subset(smallmods, trait=='rel.mass.STD' & metric=='arc_force.STD')$upper, col=as.character(subset(smallmods, trait=='rel.mass.STD' & metric=='arc_force.STD')$color))

abline(h=0, lty=3)
axis(2, at=c(-1,0,1))
points(nophylo$FULL.rel.mass.STD[c(2,3,8,9)] ~ c(9,21,33,45), pch=16)
segments(x0=c(9,21,33,45), y0=nophylo.lower$FULL.rel.mass.STD[c(2,3,8,9)], y1=nophylo.upper$FULL.rel.mass.STD[c(2,3,8,9)])
axis(1, at=c(9,21,33,45), labels=c('AccHor', 'DecHor', 'Arcavgvel', 'Arccentmax'), cex.axis=0.5, las=2)

plot(coef~c(1:8), pch=pch, col=as.character(color), data=subset(smallmods, trait=='rel.wingload.gcm.STD' & metric=='hor_acel.STD'), ylim=c(-1.5,1.5), xlim=c(1,45), yaxt='n', xaxt='n', ylab='effect size', xlab='', main='wing loading', cex.main=0.75, cex=0.5)
segments(x0=c(1:8), y0=subset(smallmods, trait=='rel.wingload.gcm.STD' & metric=='hor_acel.STD')$lower, y1=subset(smallmods, trait=='rel.wingload.gcm.STD' & metric=='hor_acel.STD')$upper, col=as.character(subset(smallmods, trait=='rel.wingload.gcm.STD' & metric=='hor_acel.STD')$color))

points(coef~c(13:20), pch=pch, col=as.character(color), data=subset(smallmods, trait=='rel.wingload.gcm.STD' & metric=='-hor_decel.STD'), cex=0.5)
segments(x0=c(13:20), y0=subset(smallmods, trait=='rel.wingload.gcm.STD'& metric=='-hor_decel.STD')$lower, y1=subset(smallmods, trait=='rel.wingload.gcm.STD'& metric=='-hor_decel.STD')$upper, col=as.character(subset(smallmods, trait=='rel.wingload.gcm.STD'& metric=='-hor_decel.STD')$color))

abline(h=0, lty=3)
axis(2, at=c(-1,0,1))
points(nophylo$FULL.rel.wingload.gcm.STD[c(2,3)] ~ c(9,21), pch=16)
segments(x0=c(9,21), y0=nophylo.lower$FULL.rel.wingload.gcm.STD[c(2,3)], y1=nophylo.upper$FULL.rel.wingload.gcm.STD[c(2,3)])
axis(1, at=c(9,21), labels=c('AccHor', 'DecHor'), cex.axis=0.5, las=2)

plot(coef~c(1:8), pch=pch, col=as.character(color), data=subset(smallmods, trait=='rel.AR.STD' & metric=='PRT_degrees.STD'), ylim=c(-1.5,1.5), xlim=c(1,45), yaxt='n', xaxt='n', ylab='effect size', xlab='', main='AR', cex.main=0.75, cex=0.5)
segments(x0=c(1:8), y0=subset(smallmods, trait=='rel.AR.STD'& metric=='PRT_degrees.STD')$lower, y1=subset(smallmods, trait=='rel.AR.STD'& metric=='PRT_degrees.STD')$upper, col=as.character(subset(smallmods, trait=='rel.AR.STD'& metric=='PRT_degrees.STD')$color))

points(coef~c(13:20), pch=pch, col=as.character(color), data=subset(smallmods, trait=='rel.AR.STD' & metric=='pPRT.STD'), cex=0.5)
segments(x0=c(13:20), y0=subset(smallmods, trait=='rel.AR.STD'& metric=='pPRT.STD')$lower, y1=subset(smallmods, trait=='rel.AR.STD'& metric=='pPRT.STD')$upper, col=as.character(subset(smallmods, trait=='rel.AR.STD'& metric=='pPRT.STD')$color))

abline(h=0, lty=3)
axis(2, at=c(-1,0,1))
points(nophylo$FULL.rel.AR.STD[c(10,12)] ~ c(9,21), pch=16)
segments(x0=c(9,21), y0=nophylo.lower$FULL.rel.AR.STD[c(10,12)], y1=nophylo.upper$FULL.rel.AR.STD[c(10,12)])
axis(1, at=c(9,21), labels=c('PRT_degrees', 'PRT%'), cex.axis=0.5, las=2)


# further analysis of the effect of individual wing shape (aspect ratio):

summary(nophylomods[[12]]) # is the effect of individual (relative) wing shape species-specific? refit model with a random slope for each species for wing shape
ARslopemod <- update(nophylomods[[12]], random=~rel.AR.STD|species)
spslopes1 <- ranef(ARslopemod)$rel.AR.STD # a slope estimate for each species
nitt <- 10000
permslopes1 <- data.frame(matrix(NA, nrow=nitt, ncol=25))
set.seed(1234)
for(i in 1:nitt){
	print(i)
	shuffdat <- maneuv[complete.cases(maneuv[,c('mass','wing.area')]),]
	shuffindex <- sample(1:length(shuffdat[,1]))
	shuffdat$rel.AR.STD <- shuffdat$rel.AR.STD[shuffindex]
	shuffmod <- try(  update(ARslopemod, data=shuffdat, na.action=na.omit)  )
	if(class(shuffmod) == 'try-error'){
		next
	} else {
		permslopes1[i,] <- ranef(shuffmod)$rel.AR.STD
	}	
}
permslopes1$mad <- apply(permslopes1[,1:25], 1, FUN='mad', constant=1)

summary(nophylomods[[2]]) # rel. wingload and accel, has comparable effect size to rel. AR and pPRT
WLslopemod <- update(nophylomods[[2]], random=~rel.wingload.gcm.STD|species)
spslopes2 <- ranef(WLslopemod)$rel.wingload.gcm.STD # a slope estimate for each species
permslopes2 <- data.frame(matrix(NA, nrow=nitt, ncol=25))
set.seed(5678)
for(i in 1:nitt){
	print(i)
	shuffdat <- maneuv[complete.cases(maneuv[,c('mass','wing.area')]),]
	shuffindex <- sample(1:length(shuffdat[,1]))
	shuffdat$rel.wingload.gcm.STD <- shuffdat$rel.wingload.gcm.STD[shuffindex]
	shuffmod <- try(  update(WLslopemod, data=shuffdat, na.action=na.omit)  )
	if(class(shuffmod) == 'try-error'){
		next
	} else {
		permslopes2[i,] <- ranef(shuffmod)$rel.wingload.gcm.STD
	}	
}
permslopes2$mad <- apply(permslopes2[,1:25], 1, FUN='mad', constant=1)

# Figure S11

dev.new(width=6,height=6)
layout(matrix(c(1,1,2,3,3,4), nrow=2, ncol=3, byrow=T))
par(las=1, mar=c(3,3,2,0.5), mgp=c(1.5,0.5,0))
hist(spslopes1, main='PRT% ~ AR', xlab='species-specific slopes', ylim=c(0,10), ylab='', xlim=c(-0.5,0.5), yaxt='n'); axis(2, at=c(0,5,10)); abline(v=0, lty=3)
abline(v=mad(spslopes1, constant=1), col='red')
abline(v=median(spslopes1), col='green')
hist(permslopes1$mad, main='', xlab='median abs. deviation', ylab='', xlim=c(0,0.5), ylim=c(0,10000), yaxt='n'); axis(2, at=c(0,5000,10000))
abline(v=mad(spslopes1, constant=1), col='red')

hist(spslopes2, main='AccHor ~ WL', xlab='species-specific slopes', ylim=c(0,10), ylab='', xlim=c(-0.5,0.5), yaxt='n'); axis(2, at=c(0,5,10)); abline(v=0, lty=3)
abline(v=mad(spslopes2, constant=1), col='red')
abline(v=median(spslopes2), col='green')
hist(permslopes2$mad, main='', xlab='median abs. deviation', ylab='', xlim=c(0,0.5), ylim=c(0,10000), yaxt='n'); axis(2, at=c(0,5000,10000))
abline(v=mad(spslopes2, constant=1), col='red')

# within different species, wing shape has a positive or negative association with PRT%:
plot1 <- xyplot(pPRT ~ rel.AR|sp, type=c('p','r'), pch=16, data=subset(maneuvT, species %in% names(sort(table(maneuvT$species), T))[1:15]))
dev.new(width=4,height=3); par(mar=c(1,1,0.25,0.25)); update(plot1, index.cond = function(x, y) coef(lm(y ~ x))[2], layout=c(5,3), strip=F, xlab='Relative AR', ylab='PRT%', as.table=T, ylim=c(0,1.1), xlim=c(-1.4,1.4), scales=list(alternating=F, cex=0.75, y=list(at=c(0.3,0.6,0.9)), x=list(at=c(-1,0,1)))) 

plot2 <- xyplot(hor_acel ~ rel.wingload.gcm|sp, type=c('p','r'), pch=16, data=subset(maneuvT, species %in% names(sort(table(maneuvT$species), T))[1:15]))
dev.new(width=4,height=3); update(plot2, index.cond = function(x, y) coef(lm(y ~ x))[2], layout=c(5,3), strip=F, xlab='Relative wing load', ylab='AccHormax', as.table=T, ylim=c(2,8), xlim=c(-0.13,0.19), scales=list(alternating=F, cex=0.75, y=list(at=c(3,5,7)), x=list(at=c(-0.1,0,0.1)))) 

# the negative correlation between PRT% and Arc centripetal force is robust across species:
plot3 <- xyplot(pPRT ~ arc_force|sp, type=c('p', 'r'), pch=16, data=subset(maneuvT, species %in% names(sort(table(maneuvT$species), T))[1:15]))
dev.new(width=4,height=3); update(plot3, index.cond = function(x, y) coef(lm(y ~ x))[2], layout=c(5,3), strip=F, xlab='Arccent,max (m/s2)', ylab='PRT%', as.table=T, ylim=c(0,1.1), scales=list(alternating=F, cex=0.75, y=list(at=c(0.3,0.6,0.9)), x=list(at=c(2,4,6)))) 

# save.image('comparative_maneuverability_models.RData')

