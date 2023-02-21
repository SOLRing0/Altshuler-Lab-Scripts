
# Evolution reveals the biomechanical organization of maneuvering flight in hummingbirds
# Dakin, Segre, Straw and Altshuler
# Figure 2, 4, Figure S5, and Figure S9

library(plotrix)
library(qgraph)
library(dplyr)
library(cluster)
library(pvclust)
library(nlme)
library(ape)

# get model results:
mcmcres <- read.csv('mcmc_model_results.csv')[,-1] # mcmc results
lmeres <- read.csv('lme_model_results.csv')[,-1] # lme results
rownames(mcmcres) <- mcmcres$metric; mcmcres <- mcmcres[,c(3,4,8,11,7,9:10,12:13)]
rownames(lmeres) <- lmeres$metric; lmeres <- lmeres[,c(3,4,8,11,7,9:10,12:13)]

# matrix of posterior coefficients:
formodplot <- t(mcmcres)
formodplot <- formodplot[c(5,4,8,3,6,7,9,1,2),c(1:3,9,4,5:6,12,7,8,10,11)]
formodplot <- t(formodplot)
lmeforplot <- t(lmeres)
lmeforplot <- lmeforplot[c(5,4,8,3,6,7,9,1,2),c(1:3,9,4,5:6,12,7,8,10,11)]
lmeforplot <- t(lmeforplot)

# hierarchical clustering analysis:
clustmat <- t(lmeres[,c(3:9)]) # lme reuslts
clustres <- hclust(as.dist(1-cor(clustmat, method='pearson')), method='average') # distances are in terms of correlations
plot(clustres)
coef.hclust(clustres) # agglomerative coefficient 0.84
bootclust <- pvclust(clustmat, method.hclust='average', method.dist='correlation', nboot=5000, iseed=101) # get bootstrap support for cluster nodes. takes a minute or two.
plot(bootclust); bootclust
bootclust$edges

# null distribution of agglomerative coefficient, from randomly permuted matrices:
nreps <- 10000
coef.null <- rep(NA, nreps)
set.seed(101)
for(i in 1:nreps){
	shuffmat <- clustmat
	shuffmat[] <- sample(clustmat)
	shuffclust <- hclust(as.dist(1-cor(shuffmat, method='pearson')), method='average')
	coef.null[i] <- coef.hclust(shuffclust)
	print(i)
}
quantile(coef.null, c(0.025,0.975)) # 95% central range for randomized data = 0.54-0.76

# average effect size magnitudes, per trait:
rowMeans(abs(clustmat))

# to get null distributions for average effect sizes per trait, refit LME models to randomly permuted data:
maneuv <- read.csv('processed_maneuver_data_170724.csv')[,-1] # get original data
formulalist <- list() # list of model formulas for each maneuver
for(i in 1:12){
	formulalist[[i]] <- paste(names(maneuv)[50+i] , ' ~ ', 'hi.elev + sp.mass.STD + rel.mass.STD + sp.muscle.STD + sp.wing.load.STD + sp.AR.STD + rel.wingload.gcm.STD + rel.AR.STD + ', c(names(maneuv)[63:68], rep(names(maneuv)[69],3), rep(names(maneuv)[70],2), 'hi')[i], sep='')
}
formulalist[[12]] <- paste(names(maneuv)[50+12] , ' ~ ', 'hi.elev + sp.mass.STD + rel.mass.STD + sp.muscle.STD + sp.wing.load.STD + sp.AR.STD + rel.wingload.gcm.STD + rel.AR.STD', sep='')
formulalist[[3]] <- paste('-',formulalist[[3]],sep='')
formulalist[[5]] <- paste('-',formulalist[[5]],sep='')
formulalist[[11]] <- paste('-',formulalist[[11]],sep='')

# randomize data and refit model results:
nreps <- 10000
nulltrait <- data.frame(matrix(NA, nrow=nreps, ncol=7))
names(nulltrait) <- gsub('FULL.', '', rownames(clustmat))
nulltrait$coef <- NA
set.seed(101)
for(i in 1:nreps){
	shuffmat <- clustmat
	shuffdat <- cbind(maneuv[,c(1,50:62)], maneuv[,c(41,71,45,44,43,74,73)], maneuv[,c(63:70)])
	shuffdat[,2:29] <- apply(shuffdat[,2:29], 2, sample) # permute each column, separately
	for(k in 1:12){	
		shuffmat[,k] <- abs(summary(lme(as.formula(formulalist[[k]]), random=~1|species, data=shuffdat, na.action=na.omit))$tTable[3:9,'Value'])
	}
	nulltrait[i,1:7] <- rowMeans(abs(shuffmat))
	shuffclust <- hclust(as.dist(1-cor(shuffmat, method='pearson')), method='average') # added
	nulltrait$coef[i] <- coef.hclust(shuffclust) # added
	print(i)
}

# Figure 4 effect size matrix, with dendrogram and trait avg. results:

dendorder <- rownames(lmeres)[c(9,1:3,4,6,5,12,11,10,8,7)]
formodplot <- formodplot[match(dendorder,rownames(formodplot)), c(1,4,5,6,2,3,7,8,9)]
lmeforplot <- lmeforplot[match(dendorder,rownames(lmeforplot)), c(1,4,5,6,2,3,7,8,9)]
nulltrait <- nulltrait[ ,match(gsub('FULL.','',colnames(lmeforplot)[1:7]),  names(nulltrait))]
traitavgs <- rowMeans(abs(clustmat))[match(colnames(lmeforplot)[1:7], rownames(clustmat))]

dev.new(width=9, height=6)
layout(matrix(c(rep(c(rep(1,3),rep(2,11)),3), rep(3,3), 4:14), 4, 14, byrow = T))
plot(as.phylo(clustres), cex=0.9, label.offset=0.1)
# nodelabels() # identify the node labels
nodelabels(text=(round(bootclust$edges$au[1:10],2))*100, node=c(21,23,22,20,17,19,18,14,16,15), frame='none', cex=0.75, adj=c(1.1,-0.3))
par(las=1)
plot(c(1:12), type='n', xaxt='n', yaxt='n', ylim=c(12.5,0.5), xlim=c(0.5,11.5), xlab='', ylab='',  mar=c(0.25,8,7,0.25)) # axis is reversed
axis(3, at=1:11, labels=(gsub('FULL.','',gsub('.STD','',colnames(formodplot)[c(1:7,8,2,9,5)]))), cex.axis=0.5)
axis(2, at=1:12, labels=rownames(formodplot), cex.axis=0.5)
range(round(formodplot + 0.4314587, 2)*100+1) 
range(round(lmeforplot + 0.4314587, 2)*100+1) # therefore, use a scale from 1-82
myramp <- colorRampPalette(colors=c('#FA2800','white','#00A0FA'))

for(i in 1:7){
	for(j in 1:12){
		draw.circle(i, j, radius=abs(lmeforplot[j,i])*0.6, col=myramp(82)[round(lmeforplot[j,i] + 0.4314587, 2)*100+1], border=0)
		draw.circle(i, j, radius=abs(formodplot[j,i])*0.6)
	}
}
for(i in 8:11){
	index <- c(rep(NA, 7), 8,2,9,5)
	for(j in 1:12){
		draw.circle(i, j, radius=abs(lmeforplot[j,index[i]])*0.6, col=myramp(82)[round(lmeforplot[j,index[i]] + 0.4314587, 2)*100+1], border=0)
		draw.circle(i, j, radius=abs(formodplot[j,index[i]])*0.6)
	}
}
for(i in 1:3) draw.circle(11.5, 12, radius=c(0.1,0.2,0.4)*0.6[i])
abline(h=c(0.5,12.5), lty=4); abline(v=c(0.5,7.5,11.5), lty=4)
par(mar=c(7,1,1,0.25), las=1, cex.axis=0.5, mgp=c(1.5,0.1,0), tck=-0.05, las=0, cex.main=0.5)
plot.new()
hist(coef.null, main='coef', xlab='', ylab='', xlim=c(0.5,1), ylim=c(0,4000), xaxt='n', yaxt='n'); abline(v=coef.hclust(clustres), col='red')
polygon(x=c(quantile(coef.null, c(0.025,0.975)), rev(quantile(coef.null, c(0.025,0.975)))), y=c(0,0,4000,4000), border=NA, col=rgb(0,0,0,0.1))
axis(1, at=c(0.5,1)); axis(2, at=c(0,4000))
for(i in 1:7){
	hist(nulltrait[,i], main=gsub('.STD','',names(nulltrait)[i]), xlim=c(0,0.22), ylim=c(0,4000), xlab='', ylab='', xaxt='n', yaxt='n')
	axis(1, at=c(0,0.2))
	axis(2, at=c(0,4000))
	abline(v=traitavgs[i], col='red')
	polygon(x=c(quantile(nulltrait[,i], c(0.025,0.975)), rev(quantile(nulltrait[,i], c(0.025,0.975)))), y=c(0,0,4000,4000), border=NA, col=rgb(0,0,0,0.1))
}

color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
	scale = (length(lut)-1)/(max-min)
    dev.new(width=1.75, height=5)
    plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
    axis(2, ticks, las=1)
    for (i in 1:(length(lut)-1)) {
     y = (i-1)/scale + min
     rect(0,y,10,y+1/scale, col=lut[i], border=NA)
    }
}
color.bar(myramp(82), min(lmeforplot), max=max(formodplot), ticks=seq(-0.4,0.4,by=0.1))


# Figure S9 trait-maneuver bipartite network:

fullgraph <- matrix(0, nrow=19, ncol=19, dimnames=list(c(colnames(formodplot)[1:7],rownames(formodplot)), c(colnames(formodplot)[1:7],rownames(formodplot))))
fullgraph[1:7,8:19] <- t(as.matrix(formodplot[,1:7]))
fullgraph[8:19,1:7] <- t(fullgraph[1:7,8:19])
fullgraph[abs(fullgraph)<0.15] <- 0
shapes <- c('triangle','square','circle','circle')
allg <- list(traits=c(1:7), trans=c(9:11), rot=c(12:14), turns=c(8,15:19))

dev.new(width=10,height=5)
par(mfrow=c(2,2))
spgraph <- qgraph(fullgraph*40, directed=F, groups=allg, label.cex=1.5, labels=gsub('.STD','',gsub('FULL.','',rownames(fullgraph))), label.scale.equal=T, maximum=0.42, layout='spring', mode='direct', shape=shapes[c(rep(1,7),4,rep(2,3),rep(3,3),rep(4,5))], posCol='#00A0FA', negCol='#FA2800', trans=T, fade=T)
# generate a key for the coefficients, from -0.5 to 0.5:
# http://sachaepskamp.com/qgraph/examples#EXAMPLES_FOR_EDGES_UNDER_DIFFERENT_ARGUMENTS
max(abs(fullgraph))
(dat.3 <- data.frame(matrix(NA, ncol=42, nrow=42)))
seq(-0.5,0.5,by=0.05)
for(i in 1:21){
	dat.3[2*i-1, 2*i] <- seq(-0.5,0.5,by=0.05)[i]
}
L.3 <- matrix(c(rep(1:21, each=2), rep(c(1,20), 21)), ncol=2, nrow = 42, byrow=F) # create grid layout
leggraph <- qgraph(dat.3*40, directed=F, edge.labels=T, maximum=0.42, layout=L.3, cut=NULL, mode='direct', posCol='#00A0FA', negCol='#FA2800', trans=T, fade=T)


# maneuver correlation network:

maneuver <- matrix(NA, nrow=12, ncol=12, dimnames=list(names(maneuv)[16:27], names(maneuv)[16:27]))
for(i in 1:12){
	for(j in 1:12){
		maneuver[i,j] <- cor(abs(maneuv[!is.na(maneuv$arc_radius),i+15]), abs(maneuv[!is.na(maneuv$arc_radius),j+15]))
	}
}
maneuver[,'PRT_time'] <- -maneuver[,'PRT_time']
maneuver['PRT_time',] <- -maneuver['PRT_time',]
diag(maneuver) <- 0
maneuver[lower.tri(maneuver)] <- 0
shapes <- c('triangle','square','circle','circle')
mang <- list(traits=NA, trans=c(1:3), rot=c(4:6), turns=c(7:12))

# average correlation coef (3 types), centrality, shortest path length (3 categories). plot and then randomize.
maneuver[lower.tri(maneuver)] <- t(maneuver)[lower.tri(t(maneuver))]
isSymmetric(maneuver)
m <- data.frame(matrix(as.matrix(maneuver), dimnames=list(t(outer(colnames(maneuver), rownames(maneuver), FUN=paste)), NULL)))
names(m) <- 'coef'
m$type1 <- c(rep('trans',36), rep('rot',36), rep('turn',72))
m$type2 <- c(rep('trans',3),rep('rot',3),rep('turn',6))
m$coef[m$coef==0] <- NA
m$pairtype <- paste(m$type1, m$type2, sep='-')
m$pairtype <- ifelse(m$pairtype=='rot-trans','trans-rot',ifelse(m$pairtype=='turn-trans','trans-turn',ifelse(m$pairtype=='turn-rot','rot-turn',m$pairtype)))
mstore <- m[with(m, order(abs(coef), decreasing=T)),]
m <- group_by(m, pairtype)
m$coef <- ifelse(duplicated(m$coef), NA, m$coef)
expected.r <- mean(maneuver[upper.tri(maneuver)])
m <- summarize(m[!is.na(m$coef),], mean.r=mean(abs(coef)), sd.r=sd(abs(coef)), n=n(), lower.r=mean.r-1.96*sd.r/sqrt(n), upper.r=mean.r+1.96*sd.r/sqrt(n))
m <- m[with(m, order(mean.r, decreasing=T)),]

mangraph <- qgraph(maneuver*20, directed=F, groups=mang, label.cex=1.5, labels=rownames(maneuver), label.scale.equal=T, maximum=1, layout='spring', mode='direct', shape=shapes[c(rep(2,3),rep(3,3),rep(4,6))], posCol='#00A0FA', negCol='#FA2800', trans=T, fade=T)

# within-individual correlation analysis to plot in Figure 2:

library(stringr)
library(dplyr)

mdata4 <- read.csv('mdata4.csv')[,-1]
head(mdata4)

mdata4$id <- paste(mdata4$species, ifelse(nchar(mdata4$id)==1, paste('00', mdata4$id,sep=''), paste('0', mdata4$id,sep='')), sep='')
mdata4$mid.frame <- rowMeans(mdata4[,c('framee','frames')])
mdata4$id <- factor(mdata4$id)
mdata4$video <- factor(mdata4$video)
mdata4$cat <- factor(mdata4$cat)

vid.table <- table(factor(mdata4$id), mdata4$video)!=0
rownames(vid.table)[rowSums(vid.table)>1] # these 7 individuals have multi-part videos

timediff <- 1/60 # find pairs of maneuvers within 1 second of each other. units in minutes
# list of 55 possible pairwise correlations
mypairs <- subset(expand.grid(levels(mdata4$cat), levels(mdata4$cat)), Var1!=Var2)
mypairs <- subset(mypairs, Var1!=Var2)
mypairs <- mypairs[c(1:10,12:20,23:30,34:40,45:50,56:60,67:70,78:80,89:90,100),]
rownames(mypairs) <- 1:dim(mypairs)[1]

within.cor <- expand.grid(levels(factor(mdata4$id)), paste(mypairs$Var1, mypairs$Var2, sep='-'))
within.cor$id <- within.cor$Var1
within.cor$pair <- within.cor$Var2
within.cor$Var1 <- str_split(within.cor$pair, '-', simplify=T)[,1]
within.cor$Var2 <- str_split(within.cor$pair, '-', simplify=T)[,2]
dim(within.cor) # 11,385 correlations to obtain. 55 combinations. 207 birds but 7 have complication of >1 video
within.cor$r <- 0

# change decel, pitch down, and PRT time to negatives for interpretation:
mdata4$val <- ifelse(mdata4$cat=='hor_decel' | mdata4$cat=='pitch_down' | mdata4$cat=='PRT_time', -mdata4$val, mdata4$val)

for(i in 1:207){ # get 55 within-individual correlations for each of 207 individuals
	temp <- subset(mdata4, id==levels(factor(mdata4$id))[i])
	temp <- temp[with(temp, order(video, mid.frame)), ]
	if(length(levels(factor(temp$video))) > 1){ # renumber frames for multi-video
		corre <- aggregate(temp$frames, by=list(video=temp$video), 'max')
		corre$startS <- as.numeric(substr(corre$video, 10,11)) + 60*as.numeric(substr(corre$video, 8,9)) + 60*60*as.numeric(substr(corre$video, 6,7))
		corre$frame.after <- 200*(corre$startS - min(corre$startS))
		corre$frame.after <- ifelse(corre$frame.after>0, corre$frame.after+temp$frames[1], corre$frame.after)
		temp <- merge(temp, corre[,c('video','frame.after')])
		temp$mid.frame <- temp$mid.frame + temp$frame.after
	}

	temp$cut <- cut(temp$mid.frame/200/60, seq(0,max(temp$mid.frame/200/60)+1, by=timediff))
	for(k in 1:55){
		tempA <- subset(temp, cat %in% mypairs[k,1])
		tempB <- subset(temp, cat %in% mypairs[k,2])
		tempM <- merge(tempA[,c('cut','cat','val','mid.frame')], tempB[,c('cut','cat','val','mid.frame')], by='cut')
		tempM <- tempM[!duplicated(tempM$mid.frame.x),] # use each maneuver once
		tempM <- tempM[!duplicated(tempM$mid.frame.y),]
		if(dim(tempM)[1] > 10){ # set minimum sample size as 10 time-paired maneuvers
			within.cor[within.cor$id==levels(factor(temp$id)) & within.cor$Var1==tempM$cat.x[1] & within.cor$Var2==tempM$cat.y[1],]$r <- cor(tempM$val.x, tempM$val.y)
		} else {
			within.cor[within.cor$id==levels(factor(temp$id)) & within.cor$Var1==mypairs[k,1] & within.cor$Var2==mypairs[k,2],]$r <- NA
		}		
	}
	print(i)
}
hist(within.cor$r)

# summarize the average within-individual correlations for maneuver pairs within 1s of each other. we assume at least 10 such pairs are needed to calculate a correlation. 95% CIs that exlcude 0 are noted as "sig"
cor.summary <- summarize(group_by(within.cor, pair), mean.r=mean(r,na.rm=T), sd.r=sd(r,na.rm=T), n=sum(!is.na(r)))
cor.summary$lower <- cor.summary$mean.r - 1.96* cor.summary$sd.r/sqrt(cor.summary$n)
cor.summary$upper <- cor.summary$mean.r + 1.96* cor.summary$sd.r/sqrt(cor.summary$n)
cor.summary <- cor.summary[with(cor.summary, order(abs(mean.r), decreasing=T)),]
cor.summary$exclude <- ifelse(sign(cor.summary$lower)==sign(cor.summary$upper), 'sig', NA)
data.frame(cor.summary)

mstore$pair <- gsub(' ', '-', rownames(mstore))
mstore <- subset(mstore, pair %in% cor.summary$pair) # 55 pairwise corr
rownames(mstore) <- 1:55
cor.compare <- merge(cor.summary, mstore)
names(cor.compare)[c(2,8)] <- c('within.individual.r','among.individual.r')
cor.compare <- cor.compare[with(cor.compare, order(abs(among.individual.r), decreasing=T)),]
# eliminate the data from the exact same maneuver:
cor.compare <- subset(cor.compare, !pair %in% c('PRT_time-PRT_degrees', 'arc_radius-arc_force', 'arc_radius-arc_avg_vel', 'arc_force-arc_avg_vel'))


# Figure S5 showing correlation network :

dev.new(width=10,height=5)
layout(matrix(data=c(1,1,2,2,1,1,2,2,3,3,3,3,3,3,3,3), nrow=4, ncol=4, byrow=T))
mangraph <- qgraph(maneuver*20, directed=F, groups=mang, label.cex=1.5, labels=rownames(maneuver), label.scale.equal=T, maximum=1, layout='spring', mode='direct', shape=shapes[c(rep(2,3),rep(3,3),rep(4,6))], posCol='#00A0FA', negCol='#FA2800', trans=T, fade=T)

# key for the correlation coefs, from -1 to 1:
(dat.3 <- data.frame(matrix(NA, ncol=42, nrow=42)))
seq(-1,1,by=0.1)
for(i in 1:21){
	dat.3[2*i-1, 2*i] <- seq(-1,1,by=0.1)[i]
}
L.3 <- matrix(c(rep(1:21, each=2), rep(c(1,20), 21)), ncol=2, nrow = 42, byrow=F) # create grid layout
leggraph <- qgraph(dat.3*20, directed=F, edge.labels=T, maximum=1, layout=L.3, cut=NULL, mode='direct', posCol='#00A0FA', negCol='#FA2800', trans=T, fade=T)
par(las=1, bty='l')

plot((cor.compare$among.individual.r), pch=1, cex=0.35, ylim=c(-1,1), ylab='r', xlab='', xaxt='n')
abline(h=0, lty=3)
axis(1, at=1:51, labels=cor.compare$pair, cex.axis=0.5, las=2)
points((cor.compare$within.individual.r), pch=16, cex=0.5)
segments(x0=1:51, y0=cor.compare$lower, y1=cor.compare$upper)
legend('topright', pch=c(1,16), legend=c('between','within'), bty='n', cex=0.75)


# Figure 2 showing correlation analysis:

dev.new(width=5,height=2.5)
par(mfrow=c(1,2), las=1, mar=c(3,3,0.25,0.25), bty='l')
plot(mean.r ~ c(6:1), m, xlim=c(6,1), xaxt='n', ylab='|r|', xlab='', ylim=c(-0.05,1), pch=16)
abline(h=expected.r)
segments(y0=c(m$lower.r), y1=c(m$upper.r), x0=c(6:1))
axis(1, at=6:1, labels=m$pairtype, cex.axis=0.26)

plot(abs(cor.compare$among.individual.r), pch=1, cex=0.5, ylim=c(-0.05,1), ylab='abs r', xlab='', xaxt='n')
abline(h=0, lty=3)
points(abs(cor.compare$within.individual.r), pch=16, cex=0.5)
segments(x0=1:51, y0=ifelse(sign(cor.compare$lower)!=sign(cor.compare$upper)&cor.compare$within.individual.r<0, -cor.compare$lower, ifelse(cor.compare$within.individual.r>0, cor.compare$lower, abs(cor.compare$lower))), y1=ifelse(sign(cor.compare$lower)!=sign(cor.compare$upper)&cor.compare$within.individual.r<0, -cor.compare$upper, ifelse(cor.compare$within.individual.r>0, cor.compare$upper, abs(cor.compare$upper))))
legend('topright', pch=c(1,16), legend=c('between','within'), bty='n', cex=0.75)


# further checks of the influence of chamber size 

# check whether the within-individual correlations depend on body mass:

mass.corr <- merge(within.cor, maneuv[,c('id','mass','species')])
head(mass.corr)
mass.corr.results <- data.frame(p=rep(NA, 55), combination=levels(factor(mass.corr$pair)), intercept=NA, slope=NA)
for(i in 1:55){
	temptest <- 99
	temp <- subset(mass.corr, pair==levels(factor(mass.corr$pair))[i])
	try(temptest <- lme(r ~ mass, random=~1|species, data=temp, na.action=na.omit))
	if(is.list(temptest)){
		mass.corr.results$p[i] <- summary(temptest)$tTable[2,5]
		mass.corr.results$intercept[i] <- summary(temptest)$tTable[1,1]
		mass.corr.results$slope[i] <- summary(temptest)$tTable[2,1]
		print(i)
	} else {
		next
	}
}
head(mass.corr.results)
# remove the results for non-distinct maneuvers:
mass.corr.results <- mass.corr.results[-c(1,2,11,50),]
rownames(mass.corr.results) <- 1:51
mass.corr.results[mass.corr.results$p<0.05&!is.na(mass.corr.results$p),] # the only two with sig. slope have a weak negative slope

# check if the max and/or span of arcing turn radii is related to body mass:

plottable <- read.csv('plottable.csv')[,-1]
dim(subset(mdata4, cat=='arc_radius')) # 9,302
range(subset(mdata4, cat=='arc_radius')$val) # 8 cm to 2.0 m
archeck <- merge(subset(mdata4, cat=='arc_radius')[,c('id','val')], maneuv[,c('id','mass','species')], by='id')
dim(archeck) # 9,302
archeck <- merge(archeck, plottable)
archeck <- summarize(group_by(archeck, species, id, color, pch), mean=mean(val), min=min(val), max=max(val), range=max-min, mass=mean(mass))

summary(lme(range ~ mass, random=~1|species, data=archeck, na.action=na.omit)) 
summary(lme(max ~ mass, random=~1|species, data=archeck, na.action=na.omit)) 

# plot results for Figure S5:

dev.new(width=6, height=2)
par(mfrow=c(1,3), mar=c(4,4,0.25,0.25), mgp=c(2,0.5,0), las=1)
plot(1:10, type='n', bty='l', xlim=c(2.2,10.5), ylim=c(-1,1), xlab='Body mass (g)', ylab='Correlation within flight bouts')
for(i in as.numeric(rownames(mass.corr.results)[!is.na(mass.corr.results$intercept)])){
	abline(a=mass.corr.results$intercept[i], b=mass.corr.results$slope[i], col=rgb(0,0,0,0.5))
}
plot(max ~ mass, archeck, bty='l', las=1, xlab='Body mass (g)', ylab='Individual max Arc radius (m)', ylim=c(0,2), pch=pch, col=as.character(color))
plot(range ~ mass, archeck, bty='l', las=1, xlab='Body mass (g)', ylab='Individual range of Arc radii (m)', ylim=c(0,2), pch=pch, col=as.character(color))


# save.image('comparative_maneuverability_figure4.RData')











