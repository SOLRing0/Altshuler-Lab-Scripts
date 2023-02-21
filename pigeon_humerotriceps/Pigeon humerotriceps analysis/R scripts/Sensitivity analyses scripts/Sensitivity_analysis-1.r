### FIRST SENSITIVITY ANALYSIS--RAISED CUT-OFF ###
### STATS FOR GAM W/ CYCLIC COMPONENT ###
# load packages
library(mgcv)
library(plotrix)
library(visreg)

getwd()
setwd(choose.dir()) # select 'Experiments' folder
data <- read.csv("pigeon HT power by freq and phase.csv", stringsAsFactors=FALSE)

# create a column with power (mW) normalized to muscle mass (g) and converted to
# standard units of measure: W/kg
data$norm.Power_W..kg<-((data$Power_3workLoops_mW/1000)/(data$HuTri_mass_g/1000))

# FIT GENERALIZED ADDITIVE MODELS
#
# We fit two models to each data subset, by stimulus duration (normalized to cycle frequency).
#
# RESPONSE: the raw power estimated from 3 successive lengthening cycles
# PREDICTORS: 
#     (1) smoothed cyclic function of the stimulus phase
#     (2) individual as a factor
#     (3) linear function of the time since the start of the experiment
#     (4) linear functions of body mass, and humerotriceps length (no apparent effect on fitting)

out_50pc_b6 <- gam( norm.Power_W..kg ~ 
                  s(Phase_pc, bs="cc") + #factor(Bird_code) +
                  Time_Since_Start_s + Mass_g + HuTri_length_mm,
             data=subset(data, c(Freq_Hz==8.6 & Stim_dur_pc_cyc==50 & Bird_code=="COLLI-041")) )

out_69pc <- gam( norm.Power_W..kg ~ 
                  s(Phase_pc, bs="cc") + factor(Bird_code) +
                  Time_Since_Start_s + Mass_g + HuTri_length_mm,
             data=subset(data, c(Freq_Hz==8.6 & Stim_dur_pc_cyc==69)) )

# check the reasonableness of the fits
gam.check(out_50pc_b6)
gam.check(out_69pc)

plot(out_50pc_b6,pages=1)
plot(out_69pc,pages=1)

### COMPARE THE EFFECTS OF STIMULUS ONSET PHASE
### QUESTION 1 ###
#### 69% stimulus duration at 8.6 Hz ####
# check model with and without autocorrelation
phasesreduced <- lme(norm.Power_W..kg ~ factor(Phase_pc), 
		random=~1|Bird_code, 
data=subset(data, Freq_Hz==8.6 & Stim_dur_pc_cyc==69), method="ML")
# model with autocorrelation
phases <- lme(norm.Power_W..kg ~ factor(Phase_pc), 
		random=~1|Bird_code, 
		correlation=corAR1(form=~Time_Since_Start_s),
	data=subset(data, Freq_Hz==8.6 & Stim_dur_pc_cyc==69), method="ML")
# compare the two
AIC(phases, phasesreduced) # not much difference, w/o is minutely better.
# keeping autocorrelation to account for force degradation since we are not
# applying a correction factor like many others do with this type of data.

hist(residuals(phases)) # good
qqnorm(residuals(phases))
shapiro.test(residuals(phases)) # not normally distributed, however anova is robust against violating normality
summary(phases)
anova(phases)

# create a null model to do likelihood ratio testing
phases.null <- lme(norm.Power_W..kg ~ 1, 
		random=~1|Bird_code, 
		correlation=corAR1(form=~Time_Since_Start_s),
	data=subset(data, Freq_Hz==8.6 & Stim_dur_pc_cyc==69), method="ML")
hist(residuals(phases.null))
qqnorm(residuals(phases.null))
shapiro.test(residuals(phases.null)) # not normally distributed, however anova is robust against violating normality
# Likelihood ratio test
out.stim.phases69<-anova(phases.null, phases)

### 50% stimulus duration at 8.6 Hz ###
 phases2 <- lme(norm.Power_W..kg ~ factor(Phase_pc), 
		random=~1|Bird_code, 
		correlation=corAR1(form=~Time_Since_Start_s),
	data=subset(data, Freq_Hz==8.6 & Stim_dur_pc_cyc==50 & Bird_code=="COLLI-041"), method="ML")
 hist(residuals(phases2)) # good
 shapiro.test(residuals(phases2)) # not normally distributed, however anova is robust against violating normality
 summary(phases2)
 anova(phases2)

# null model
phases2.null <- lme(norm.Power_W..kg ~ 1,  
		random=~1|Bird_code, 
		correlation=corAR1(form=~Time_Since_Start_s),
	data=subset(data, Freq_Hz==8.6 & Stim_dur_pc_cyc==50 & Bird_code=="COLLI-041"), method="ML")
shapiro.test(residuals(phases2.null)) # good

# Likelihood ratio test
out.stim.phases50_b6<-anova(phases2.null, phases2)

##### calculate the summary statistics from the fit for each individual and then run an anova--lm() on summary stats~treatment

#### calculate the summary statistics for 8.6 Hz, 50% stimulus duration--Bird 6 ####
head(out_50pc_b6)
rm(data50_b6)
data50_b6<- subset(data, c(Freq_Hz==8.6 & Stim_dur_pc_cyc==50 & Bird_code=="COLLI-041"))
data50_b6$fitted<-fitted(out_50pc_b6)
# structure of fitted data prevents grouping by Bird_code therefore specify 
# using positions that correspond to each individual
# Check which fitted values belong to which subject
data50_b6$Bird_code[1:22]
# data50_b6$fitted[1:22] is COLLI-041 aka Bird 6

summaryPower50_b6 <- data.frame( stim = 50, Subject = "Bird 6",
	mean = mean(data50_b6$fitted[1:22]),
	min = range(data50_b6$fitted[1:22])[1],
	max = range(data50_b6$fitted[1:22])[2],
	range = diff(range(data50_b6$fitted[1:22])),
	Act = (length(data50_b6$Phase_pc[which(data50_b6$fitted[1:22] >= 0)])/(length(unique(data50_b6$Phase_pc))))*100 )

#### calculate the summary statistics for 8.6 Hz, 69% stimulus duration ####
head(out_69pc)

data69 <- subset(data, c(Freq_Hz==8.6 & Stim_dur_pc_cyc==69))
data69$fitted<-fitted(out_69pc)
# structure of fitted data prevents grouping by Bird_code therefore specify 
# using positions that correspond to each individual
# Check which fitted values belong to which subject
data69$Bird_code[1:22]
data69$Bird_code[23:44]
data69$Bird_code[45:62]
data69$Bird_code[63:84]
# data69$fitted[1:22] is COLLI-031 aka Bird 1
# data69$fitted[23:44] is COLLI-035 aka Bird 3
# data69$fitted[45:62] is COLLI-036 aka Bird 4
# data69$fitted[63:84] is COLLI-033 aka Bird 2

summaryPower69 <- data.frame( stim = 69, Subject = c("Bird 1", "Bird 3", 
		"Bird 4", "Bird 2"),
	mean = c(mean(data69$fitted[1:22]), mean(data69$fitted[23:44]),
		 mean(data69$fitted[45:62]),mean(data69$fitted[63:84])),
	min = c(range(data69$fitted[1:22])[1], range(data69$fitted[23:44])[1],
		 range(data69$fitted[45:62])[1],range(data69$fitted[63:84])[1]),
	max = c(range(data69$fitted[1:22])[2], range(data69$fitted[23:44])[2],
		 range(data69$fitted[45:62])[2],range(data69$fitted[63:84])[2]),		 
	range = c(diff(range(data69$fitted[1:22])),diff(range(data69$fitted[23:44])),
		  diff(range(data69$fitted[45:62])),diff(range(data69$fitted[63:84]))),
	Act = c((length(data69$Phase_pc[which(data69$fitted[1:22] >= 0)])/(length(unique(data69$Phase_pc))))*100,
		 (length(data69$Phase_pc[which(data69$fitted[23:44] >= 0)])/(length(unique(data69$Phase_pc))))*100,
		 (length(data69$Phase_pc[which(data69$fitted[45:62] >= 0)])/(length(unique(data69$Phase_pc))))*100,
		 (length(data69$Phase_pc[which(data69$fitted[63:84] >= 0)])/(length(unique(data69$Phase_pc))))*100) )

## Calculate means and standard errors of each summary statistic for reporting 
## in results

stim.summstatstable <- data.frame( stim = c(50,69),
avg_min = c(mean(summaryPower50_b6$min), mean(summaryPower69$min)),
SE_min = c(std.error(summaryPower50_b6$min), std.error(summaryPower69$min)),
avg_max = c(mean(summaryPower50_b6$max), mean(summaryPower69$max)),
SE_max = c(std.error(summaryPower50_b6$max), std.error(summaryPower69$max)),
avg_mean = c(mean(summaryPower50_b6$mean), mean(summaryPower69$mean)),
SE_mean = c(std.error(summaryPower50_b6$mean), std.error(summaryPower69$mean)),
avg_act = c(mean(summaryPower50_b6$Act), mean(summaryPower69$Act)),
SE_act = c(std.error(summaryPower50_b6$Act), std.error(summaryPower69$Act)) )

## create a single dataframe containing all of the summary stats of interest
Pow.summ.stats <-rbind(summaryPower50_b6, summaryPower69)

### COMPARE SUMMARY STATS ###
### means ###
means <-lm(mean~stim, data = Pow.summ.stats)
hist(residuals(means))
qqnorm(residuals(means))
shapiro.test(residuals(means)) # good
summary(means) # look at coefficients
out.stim.means<-anova(means) # still insignificant difference

### max ###
max.z <-lm(max~stim, data = Pow.summ.stats)
hist(residuals(max.z))
qqnorm(residuals(max.z))
shapiro.test(residuals(max.z)) # good
summary(max.z) # look at coefficients
out.stim.max<-anova(max.z) # slightly higher p-value: 0.1916 vs 0.076 but still insignificant

### min ###
min.z <-lm(min~stim, data = Pow.summ.stats)
hist(residuals(min.z))
qqnorm(residuals(min.z)) # good
shapiro.test(residuals(min.z)) # good
summary(min.z) # look at coefficients
out.stim.min<-anova(min.z) # still significant difference

### Act ###
Acts <-lm(Act~stim, data = Pow.summ.stats)
hist(residuals(Acts))
qqnorm(residuals(Acts))
shapiro.test(residuals(Acts)) # not normally distributed, however anova is robust against violating normality
summary(Acts) # look at coefficients
out.stim.acts<-anova(Acts) # still insignificant difference

# write a row to append to the summary stats table, containing the p-values 
# from the comparison of that summary stat across the two stimulus durations
summ.p_values<-cbind(out.stim.min$Pr[1], out.stim.max$Pr[1], out.stim.means$Pr[1], out.stim.acts$Pr[1]) 
# append to summary stats table
stim.summstatstable[3,c(2,4,6,8)]<-summ.p_values
# add row name
stim.summstatstable[3,1]<-"p-values-summary stats comparisons"
# add column with phase effects p-values
stim.summstatstable$p_values.phase_effects<-c(out.stim.phases50_b6$'p-value'[2],out.stim.phases69$'p-value'[2],NA)

# write csv of summary stats and p-values
write.csv(stim.summstatstable,"summary stats and p-values--first sensitivity analysis.csv",row.names=FALSE)

# SET UP THE PLOTTING WINDOW
plot.mat <- matrix( 1, ncol=2, nrow=6 )

plot.mat[1] <-0
plot.mat[2:3] <- 1
plot.mat[4:5] <- 2
plot.mat[6] <- 0
plot.mat[7:8] <- 3
plot.mat[9:10] <- 4
plot.mat[11:12] <- 5


# some global options
ylims=c(-800,200)
xlims=c(-50,50)
presid.cex = 0


# PLOTTING
layout(plot.mat)
par( mar=c(2,5,1,1) )

# visreg used for its plotting format and for visualising GAM fit (if desired)--
# to see model fit, change visreg line parameters 
## If model fit line plotting turned off, plotting error returned which can be ignored.
# Plot 8.6 Hz at 69% norm.stim duration 
visreg( out_69pc, "Phase_pc",
      frame=FALSE,
      xaxt="n", yaxt="n",
      ylab=" ",
      xlab="Stimulus onset (% cycle)",
      ylim=ylims, xlim=xlims,
      points=list(cex=0),
      band=F,
      line=list(col="purple", lwd=0, lty=" ")) # remove lwd= and lty= to see model fit
points( norm.Power_W..kg ~ Phase_pc,
        col="purple", pch=(as.numeric(factor(Bird_code))+14), cex=0.75,
        data=subset(data, Freq_Hz==8.6 & Stim_dur_pc_cyc==69) )
#Add axes manually
axis(2, at=seq(ylims[1],ylims[2], by=ylims[2]), labels=TRUE, tcl=0.5)
abline(h=0, lty="dashed")

# Add legend
legend(-48,-300, legend = c("Bird 1","Bird 2","Bird 3","Bird 4"), 
        pch = c(15,18,16,17), col = "purple", bty="n")
mtext(side = 2, at=-200, text = "8.6 Hz, 69% cycle", line =3.5)
mtext(side = 2, at=-750, text = "Power (W/kg)", line =2.2, cex=0.75)

# Plot 8.6 Hz at 50% norm.stim duration 
visreg( out_50pc_b6, "Phase_pc",
      frame=FALSE, 
      xaxt="n", yaxt="n",
      xlab="Stimulus onset (% cycle)",
      ylab=" ",
      ylim=ylims, xlim=xlims,
      points=list(cex=0),
      band=F,
      line=list(col="red", lwd=0, lty=" ")) # remove lwd= and lty= to see model fit
#Add axes manually
axis(1, at=seq(-50,50,by=25), labels=TRUE, tcl=0.5)
points( norm.Power_W..kg ~ Phase_pc,
        col="red", pch=(as.numeric(factor(Bird_code))+1), cex=0.75,
        data=subset(data, Freq_Hz==8.6 & Stim_dur_pc_cyc==50 & Bird_code=="COLLI-041") )
#Add axes manually
axis(2, at=seq(ylims[1],ylims[2], by=ylims[2]), labels=TRUE, tcl=0.5)
abline(h=0, lty="dashed")

# Add legend
legend(-48,-400, legend = "Bird 6", 
        pch = 2, col = "red", bty="n")
mtext(side = 2, at=-200, text = "8.6 Hz, 50% cycle", line =3.5)

# plot mean vs stimulus intensity
plot( mean~stim, data=Pow.summ.stats,
      cex=1,
      col=c("red","purple")[as.factor(Pow.summ.stats$stim)],
      pch=c(2,15:18)[as.factor(Pow.summ.stats$Subject)], bty="n", tcl=0.5,
      ylim=c(-250,5), xlim=c(40,80),
      xaxt="n",
      ylab="Mean power (W/kg)" )
abline(h=0, lty="dashed")
axis(1, at=c(40,50,69,80), labels=FALSE, tcl=0.5)

# plot Max & Min power vs stimulus intensity
plot( max~stim, data=Pow.summ.stats,
      cex=1,
      col=c("red","purple")[as.factor(Pow.summ.stats$stim)],
      pch=c(2,15:18)[as.factor(Pow.summ.stats$Subject)], bty="n", tcl=0.5,
      ylim=c(-600,200), xlim=c(40,80),
      xaxt="n",
      ylab="Power (W/kg)" )
points( min~stim, data=Pow.summ.stats,
	cex=1,
	col=c("red","purple")[as.factor(Pow.summ.stats$stim)],
	pch=c(2,15:18)[as.factor(Pow.summ.stats$Subject)])
abline(h=0, lty="dashed")
axis(1, at=c(40,50,69,80), labels=FALSE, tcl=0.5)

# plot actuator range vs stimulus intensity, coloured by stimulus intensity 
plot( Act~stim, data=Pow.summ.stats,
      cex=1,
      col=c("red","purple")[as.factor(Pow.summ.stats$stim)],
      pch=c(2,15:18)[as.factor(Pow.summ.stats$Subject)], bty="n", tcl=0.5,
      xlim=c(40,80), ylim=c(0,50),      
      xaxt="n",
      xlab="Stimulus intensity (%)", ylab="Actuator range (% cycle)")
axis(1, at=c(40,50,69,80), labels=TRUE, tcl=0.5)


