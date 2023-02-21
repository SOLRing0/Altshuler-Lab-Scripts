
# Evolution reveals the biomechanical organization of maneuvering flight in hummingbirds
# Dakin, Segre, Straw and Altshuler
# Figure 3 and Figure S5

library(adegenet)
library(phylobase)
library(phylosignal)
library(lattice)
library(MASS)
library(ape)

# get maneuver data: 

maneuv <- read.csv('processed_maneuver_data_170724.csv')[,-1]
plottable <- read.csv('plottable.csv')[,-1]

# set up for (DAPC) LDA analyses and cross-validation

maneuv.LDA <- subset(maneuv, !is.na(arc_avg_vel)&!is.na(wingload.gcm))[,c(1,12:15,16:35)] # use individuals with complete data only
dim(maneuv.LDA) # n = 180 individuals
rownames(maneuv.LDA) <- 1:dim(maneuv.LDA)[1]
maneuv.LDA$row <- rownames(maneuv.LDA)
for(i in 18:25){
	maneuv.LDA[,i] <- log(maneuv.LDA[,i]) # transform maneuver sample sizes to improve normality
	names(maneuv.LDA)[i] <- paste('log',names(maneuv.LDA)[i],sep='.')
}

nreps = 10000
results.morph.LDA <- data.frame(percent.correct=rep(NA, nreps), clade.correct=NA)
results.morph.rand <- data.frame(percent.correct=rep(NA, nreps), clade.correct=NA)
results.behav.LDA <- data.frame(percent.correct=rep(NA, nreps), clade.correct=NA)
results.behav.rand <- data.frame(percent.correct=rep(NA, nreps), clade.correct=NA)

# morphology LDA, 3 traits (mass, wing area, wing AR):

model1.morph <- dapc(x=maneuv.LDA[,2:4], grp=maneuv.LDA$species, n.pca=3, n.da=3) # model using all complete records

# cross-validation for morphology:
set.seed(101)
for(j in 1:nreps){
	thing <- aggregate(maneuv.LDA$row, by=list(maneuv.LDA$species), sample)$x # take a sample of each multi-individual species, for cross validation
	index <- NULL
	for(i in 1:25){
		index <- c(index, as.numeric(as.character(thing[[i]][1:ceiling((length(thing[[i]])/1.5))])))
	}
	index <- sort(index)
	build.dat <- maneuv.LDA[index, 1:4] # n = 129
	test.dat <- maneuv.LDA[-index, 1:4] # n = 51
	mymodel <- dapc(x=build.dat[,2:4], grp=build.dat$species, n.pca=3, n.da=3)
	pred <- predict(mymodel, newdata=test.dat[,c(2:4)])
	results.morph.LDA$percent.correct[j] <- sum(pred$assign == test.dat$species)/length(pred$assign)
	results.morph.LDA$clade.correct[j] <- sum(plottable$clade[match(pred$assign, plottable$species)] == plottable$clade[match(test.dat$species, plottable$species)])/length(pred$assign) 
	print(j)
}
summary(results.morph.LDA) # 65% correctly classified

# randomization test for morphology, randomizing species labels:
set.seed(102)
for(j in 1:nreps){
	shuffle <- sample(1:dim(maneuv.LDA)[1])
	maneuv.LDA.rand <- maneuv.LDA
	maneuv.LDA.rand$species2 <- maneuv.LDA.rand$species[shuffle] # randomize the species labels
	thing <- aggregate(maneuv.LDA.rand$row, by=list(maneuv.LDA.rand$species2), sample)$x
	index <- NULL
	for(i in 1:25){
		index <- c(index, as.numeric(as.character(thing[[i]][1:ceiling((length(thing[[i]])/1.5))])))
	}
	index <- sort(index)
	build.dat <- maneuv.LDA.rand[index, c(27,2:4)] 
	test.dat <- maneuv.LDA.rand[-index, c(27,2:4)] 
	mymodel <- dapc(x=build.dat[,2:4], grp=build.dat$species2, n.pca=3, n.da=3)
	pred <- predict(mymodel, newdata=test.dat[,c(2:4)])
	results.morph.rand$percent.correct[j] <- sum(pred$assign == test.dat$species2)/length(pred$assign)
	results.morph.rand$clade.correct[j] <- sum(plottable$clade[match(pred$assign, plottable$species)] == plottable$clade[match(test.dat$species2, plottable$species)])/length(pred$assign) 
	print(j)
}

# maneuver LDA, 20 variables (12 metrics + 8 sample sizes):

model2.behav <- dapc(x=maneuv.LDA[,6:25], grp=maneuv.LDA$species, n.pca=20, n.da=20) # model using all data

# cross-validation for maneuvers:
set.seed(103)
for(j in 1:nreps){
	thing <- aggregate(maneuv.LDA$row, by=list(maneuv.LDA$species), sample)$x
	index <- NULL
	for(i in 1:25){
		index <- c(index, as.numeric(as.character(thing[[i]][1:ceiling((length(thing[[i]])/1.5))])))
	}
	index <- sort(index)
	build.dat <- maneuv.LDA[index, c(1,6:25)] 
	test.dat <- maneuv.LDA[-index, c(1,6:25)] 
	mymodel <- dapc(x=build.dat[,2:21], grp=build.dat$species, n.pca=20, n.da=20)
	pred <- predict(mymodel, newdata=test.dat[,c(2:21)])
	results.behav.LDA$percent.correct[j] <- sum(pred$assign == test.dat$species)/length(pred$assign)
	results.behav.LDA$clade.correct[j] <- sum(plottable$clade[match(pred$assign, plottable$species)] == plottable$clade[match(test.dat$species, plottable$species)])/length(pred$assign)
	print(j)
}
summary(results.behav.LDA) # 34% correctly classified

# randomization test for maneuvers:
set.seed(104)
for(j in 1:nreps){
	shuffle <- sample(1:dim(maneuv.LDA)[1])
	maneuv.LDA.rand <- maneuv.LDA
	maneuv.LDA.rand$species2 <- maneuv.LDA.rand$species[shuffle]
	thing <- aggregate(maneuv.LDA.rand$row, by=list(maneuv.LDA.rand$species2), sample)$x
	index <- NULL
	for(i in 1:25){
		index <- c(index, as.numeric(as.character(thing[[i]][1:ceiling((length(thing[[i]])/1.5))])))
	}
	index <- sort(index)
	build.dat <- maneuv.LDA.rand[index, c(27,6:25)] 
	test.dat <- maneuv.LDA.rand[-index, c(27,6:25)] 
	mymodel <- dapc(x=build.dat[,2:21], grp=build.dat$species2, n.pca=20, n.da=20)
	pred <- predict(mymodel, newdata=test.dat[,c(2:21)])
	results.behav.rand$percent.correct[j] <- sum(pred$assign == test.dat$species2)/length(pred$assign)
	results.behav.rand$clade.correct[j] <- sum(plottable$clade[match(pred$assign, plottable$species)] == plottable$clade[match(test.dat$species2, plottable$species)])/length(pred$assign)
	print(j)
}

# Cohen's kappa

# morphology
(mean(results.morph.LDA$percent.correct) - mean(results.morph.rand$percent.correct))/(1 - mean(results.morph.rand$percent.correct)) # 0.61

# maneuvers
(mean(results.behav.LDA$percent.correct) - mean(results.behav.rand$percent.correct))/(1 - mean(results.behav.rand$percent.correct)) # 0.30

# Figure 3 scatterplots with insets:

# inset panels for % correct in figure 3
myinset1 <- function(){
	par(las=1, bty='n', cex.lab=0.5, cex.axis=0.5, mgp=c(0.75, 0.05, 0), tck=-0.03)
	hist(results.morph.rand$percent.correct, xlim=c(0,1), xlab='percent correct', main='', ylab='', yaxt='n', xaxt='n', ylim=c(0,4000))
	axis(1, at=c(0, 0.5, 1), labels=c(0, 50, 100))
	points(c(500) ~ quantile(results.morph.LDA$percent.correct, 0.5))
	segments(x0=quantile(results.morph.LDA$percent.correct, 0.025), x1=quantile(results.morph.LDA$percent.correct, 0.975), y0=500, y1=500)
}
myinset2 <- function(){
	par(las=1, bty='n', cex.lab=0.5, cex.axis=0.5, mgp=c(0.75, 0.05, 0), tck=-0.03)
	hist(results.behav.rand$percent.correct, xlim=c(0,1), xlab='percent correct', main='', ylab='', yaxt='n', xaxt='n', ylim=c(0,4000))
	axis(1, at=c(0, 0.5, 1), labels=c(0, 50, 100))
	points(c(500) ~ quantile(results.behav.LDA$percent.correct, 0.5))
	segments(x0=quantile(results.behav.LDA$percent.correct, 0.025), x1=quantile(results.behav.LDA$percent.correct, 0.975), y0=500, y1=500)
}

dev.new(width=4,height=6)
par(mfrow=c(2,1))
scatter.dapc(model2.behav, label=NULL, scree.da=F, leg=F, col=plottable$color, pch=plottable$pch, cell=0, cstar=0, cex=0.75, ylim=c(5,-5))
add.scatter(myinset2(), ratio=0.25, inset=c(-0.15,-0.22), bg.col=rgb(1,1,1,0))
scatter.dapc(model1.morph, label=NULL, scree.da=F, leg=F, col=plottable$color, pch=plottable$pch, cell=0, cstar=0, cex=0.75)
add.scatter(myinset1(), ratio=0.25, inset=c(-0.15,-0.22), bg.col=rgb(1,1,1,0))
# axes represent DF1 (x) & DF2 (y) (a.u., arbitrary units)

# top loadings for DF1s:

dev.new(width=2,height=3)
par(las=1, mfrow=c(2,1), mar=c(4,0.25,0.25,4))
plot(sort(round(model2.behav$var.contr[,'LD1'], 10), decreasing=T)[1:4], ylab='', xlab='', xaxt='n', yaxt='n', ylim=c(0,1), xlim=c(0.5,4.5), bty='l', pch=16, cex=0.5) # > only 0.8
axis(1, at=1:4, labels=names(sort(round(model2.behav$var.contr[,'LD1'], 5), decreasing=T))[1:4], las=2, cex.axis=0.5)
axis(4, at=c(0,1), las=3, cex.axis=0.5)
segments(x0=1:4, x1=1:4, y1=0, y0=sort(round(model2.behav$var.contr[,'LD1'], 5), decreasing=T)[1:4])
plot(sort(round(model1.morph$var.contr[,'LD1'], 10), decreasing=T), ylab='', xlab='', xaxt='n', yaxt='n', ylim=c(0,1), xlim=c(0.5,4.5), bty='l', pch=16, cex=0.5)
axis(1, at=1:3, labels=names(sort(round(model1.morph$var.contr[,'LD1'], 5), decreasing=T)), las=2, cex.axis=0.5)
axis(4, at=c(0,1), las=3, cex.axis=0.5)
segments(x0=1:3, x1=1:3, y1=0, y0=sort(round(model1.morph$var.contr[,'LD1'], 5), decreasing=T))

# Figure S5 all loadings plot:

# heatmap of loadings with cell widths scaled to the eigenvalues of each DF

model2.behav$var.contr # loadings
model2.behav$eig/sum(model2.behav$eig) # column widths
model2.behav$eig/sum(model2.behav$eig)*20

heatmatrix <- t(model2.behav$var.contr)[,20:1]
heatmatrix2 <- t(model1.morph$var.contr)[,3:1]
heatmatrix <- cbind(heatmatrix, maxmin=c(rep(0,10),rep(1,10)))
heatmatrix2 <- cbind(heatmatrix2, maxmin=c(rep(0,1),rep(1,2)))

# using levelplot from the lattice package. because of the way row.values works, have to manually solve for the break points:

(morph.ends <- cumsum(model1.morph$eig/sum(model1.morph$eig))*3)
(morph.breaks <- c(0, cumsum(model1.morph$eig/sum(model1.morph$eig))*3)[-4]+diff(c(0, cumsum(model1.morph$eig/sum(model1.morph$eig))*3))/2)
morph.breaks[1] <- morph.ends[1] - (morph.breaks[2] - morph.ends[1])
morph.breaks[3] <- morph.ends[2] + (morph.ends[2] - morph.breaks[2])
(behav.ends <- cumsum(model2.behav$eig/sum(model2.behav$eig))*20)
(behav.breaks <- c(0, cumsum(model2.behav$eig/sum(model2.behav$eig))*20)[-21]+diff(c(0, cumsum(model2.behav$eig/sum(model2.behav$eig))*20))/2)
behav.breaks[1] <- behav.ends[1] - (behav.breaks[2] - behav.ends[1])
behav.breaks[20] <- behav.ends[19] + (behav.ends[19] - behav.breaks[19])

# the following plots need to be manually tidied, to bring the start of the first column to 0 and end of the last column to the max (i.e., 3 for morphology, 20 for maneuvers)
# the first row in each plot (maxmin) is there as a reference only and is deleted from the final version
dev.new(width=5, height=5)
levelplot(x=heatmatrix, row.values=behav.breaks, xlab='', ylab='', col.regions=colorRampPalette(c('white','lemonchiffon2','goldenrod','red','purple','darkblue'))(50), region=T, cuts=50, xlim=c(0-0.1*20, 20+0.1*20))

dev.new(width=5, height=5)
levelplot(x=heatmatrix2, row.values=morph.breaks, xlab='', ylab='', col.regions=colorRampPalette(c('white','lemonchiffon2','goldenrod','red','purple','darkblue'))(50), region=T, cuts=50, xlim=c(0-0.1*3, 3+0.1*3))
# row.values works by drawing n-1 lines halfway between the specified values. then adds the end lines. you will have to adjust the end lines to be the correct size

# phylogenetic lambda for 1st DF axis as a trait (using sp. average coordinate values):

tree <- read.nexus(file="hum294.tre") 
lookup <- read.csv('hum294_tiplabels.csv')
lookup <- subset(lookup, species!='' & species %in% maneuv$species)
lookup$tree <- factor(lookup$tree); lookup$species <- factor(lookup$species)
maneuvtree <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% lookup$tree])
lookup <- lookup[match(maneuvtree$tip.label, lookup$tree), ] # reorder to match tree
maneuvtree$tip.label <- as.character(lookup$species)
plot(maneuvtree)

# species average coordinate values:
morph.coords <- data.frame(model1.morph$grp.coord)
morph.coords$species <- factor(rownames(morph.coords), levels=maneuvtree$tip.label)
morph.coords <- morph.coords[with(morph.coords, order(species)),]
behav.coords <- data.frame(model2.behav$grp.coord)
behav.coords$species <- factor(rownames(behav.coords), levels=maneuvtree$tip.label)
behav.coords <- behav.coords[with(behav.coords, order(species)),]

phylo.morph <- phylo4d(maneuvtree, tip.data=morph.coords$LD1)
phyloSignal(phylo.morph) # lambda = 1.1, p = 0.03
phylo.maneuv <- phylo4d(maneuvtree, tip.data=behav.coords$LD1)
phyloSignal(phylo.maneuv) # lambda = 0.8, p = 0.003

# check, how well can the species be classified on the basis of just 1 trait/performance metric at a time?

set.seed(123)
one.trait <- data.frame(matrix(NA, nrow=nreps, ncol=23))
names(one.trait) <- names(maneuv.LDA)[c(2:4,6:25)]
for(k in 1:3){
	for(j in 1:nreps){
		thing <- aggregate(maneuv.LDA$row, by=list(maneuv.LDA$species), sample)$x
		index <- NULL
		for(i in 1:25){
			index <- c(index, as.numeric(as.character(thing[[i]][1:ceiling((length(thing[[i]])/1.5))])))
		}
		index <- sort(index)
		build.dat <- maneuv.LDA[index, c('species',names(one.trait)[k])] 
		test.dat <- maneuv.LDA[-index, c('species',names(one.trait)[k])] 
		mymodel <- lda(as.formula(paste('species ~ ',names(one.trait)[k],sep='')), data=build.dat)
		pred <- predict(mymodel, newdata=test.dat)$class
		one.trait[j, k] <- sum(pred == test.dat$species)/length(pred)
		print(c(k, j))
	}
}
for(k in 4:23){
	for(j in 1:nreps){
		thing <- aggregate(maneuv.LDA$row, by=list(maneuv.LDA$species), sample)$x
		index <- NULL
		for(i in 1:25){
			index <- c(index, as.numeric(as.character(thing[[i]][1:ceiling((length(thing[[i]])/1.5))])))
		}
		index <- sort(index)
		build.dat <- maneuv.LDA[index, c('species',names(one.trait)[k])] 
		test.dat <- maneuv.LDA[-index, c('species',names(one.trait)[k])] 
		mymodel <- lda(as.formula(paste('species ~ ',names(one.trait)[k],sep='')), data=build.dat)
		pred <- predict(mymodel, newdata=test.dat)$class
		one.trait[j, k] <- sum(pred == test.dat$species)/length(pred)
		print(c(k, j))
	}
}
range(colMeans(one.trait[,1:3])) # for morphology, 20-47% correct using single traits
range(colMeans(one.trait[,4:23])) # for behaviour, 8-19% correct using single metrics

# save.image('comparative_maneuverability_figure3.RData')









