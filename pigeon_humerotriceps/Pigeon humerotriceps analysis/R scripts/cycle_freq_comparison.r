##### COMPARE THE EFFECTS OF CYCLE FREQUENCY ##### 
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
# We fit three models to each data subset, by frequency.
#
# RESPONSE: the raw power estimated from 3 successive lengthening cycles
# PREDICTORS: 
#     (1) smoothed cyclic function of the stimulus phase
#     (2) individual as a factor
#     (3) linear function of the time since the start of the experiment
#     (4) linear functions of body mass, and humerotriceps length (no apparent effect on fitting)

out_6.1_50 <- gam( I(norm.Power_W..kg) ~ 
                  s(Phase_pc, bs="cc") + factor(Bird_code) +
                  Time_Since_Start_s + Mass_g + HuTri_length_mm,
             data=subset(data, c(Freq_Hz==6.1 & Stim_dur_pc_cyc==49)) )

out_8.6_50 <- gam( I(norm.Power_W..kg) ~ 
                  s(Phase_pc, bs="cc") + factor(Bird_code) +
                  Time_Since_Start_s + Mass_g + HuTri_length_mm,
             data=subset(data, c(Freq_Hz==8.6 & Stim_dur_pc_cyc==50)) )

out_10.1_50 <- gam( I(norm.Power_W..kg) ~ 
                   s(Phase_pc, bs="cc") + factor(Bird_code) +
                   Time_Since_Start_s + Mass_g + HuTri_length_mm,
             data=subset(data, c(Freq_Hz==10.1 & Stim_dur_pc_cyc==50)) )

# check the reasonableness of the fits
gam.check(out_6.1_50)
gam.check(out_8.6_50)
gam.check(out_10.1_50)

plot(out_6.1_50,pages=1)
plot(out_8.6_50,pages=1)
plot(out_10.1_50,pages=1)

### COMPARE THE EFFECTS OF STIMULUS ONSET PHASE
### QUESTION 1 ###
#### 6.1 Hz ####
 phases3 <- lme(norm.Power_W..kg ~ factor(Phase_pc), 
		random=~1|Bird_code, 
		correlation=corAR1(form=~Time_Since_Start_s),
	data=subset(data, Freq_Hz==6.1 & Stim_dur_pc_cyc==49), method="ML")
 hist(residuals(phases3)) # good
 shapiro.test(residuals(phases3)) # good
 summary(phases3)
 anova(phases3)

# null model
phases3.null <- lme(norm.Power_W..kg ~ 1, 
		random=~1|Bird_code, 
		correlation=corAR1(form=~Time_Since_Start_s),
	data=subset(data, Freq_Hz==6.1 & Stim_dur_pc_cyc==49), method="ML") 
shapiro.test(residuals(phases3.null)) # borderline non-normal, however anova is robust against violating normality
# Likelihood ratio test
out.freq.phases6<-anova(phases3.null, phases3)

#### 8.6 Hz ####
 phases4 <- lme(norm.Power_W..kg ~ factor(Phase_pc), 
		random=~1|Bird_code, 
		correlation=corAR1(form=~Time_Since_Start_s),
	data=subset(data, Freq_Hz==8.6 & Stim_dur_pc_cyc==50), method="ML")
 hist(residuals(phases4)) # good
 shapiro.test(residuals(phases4)) # good
 summary(phases4)
 anova(phases4)

# null model
phases4.null <- lme(norm.Power_W..kg ~ 1, 
		random=~1|Bird_code, 
		correlation=corAR1(form=~Time_Since_Start_s),
	data=subset(data, Freq_Hz==8.6 & Stim_dur_pc_cyc==50), method="ML")
shapiro.test(residuals(phases4.null)) # good
# Likelihood ratio test
out.freq.phases8<-anova(phases4.null, phases4)

#### 10.1 Hz ####
 phases5 <- lme(norm.Power_W..kg ~ factor(Phase_pc), 
		random=~1|Bird_code, 
		correlation=corAR1(form=~Time_Since_Start_s),
	data=subset(data, Freq_Hz==10.1 & Stim_dur_pc_cyc==50), method="ML")
 hist(residuals(phases5)) # good
 shapiro.test(residuals(phases5)) # good
 summary(phases5)
 anova(phases5)

# null model
phases5.null <- lme(norm.Power_W..kg ~ 1, 
		random=~1|Bird_code, 
		correlation=corAR1(form=~Time_Since_Start_s),
	data=subset(data, Freq_Hz==10.1 & Stim_dur_pc_cyc==50), method="ML")
shapiro.test(residuals(phases5.null)) # not normally distributed, however anova is robust against violating normality
# Likelihood ratio test
out.freq.phases10<-anova(phases5.null, phases5)

##### calculate the summary statistics from the fit for each individual and then run an anova--lm() on summary stats~treatment

#### calculate the summary statistics for 6.1 Hz ####
head(out_6.1_50)

data6 <- subset(data, c(Freq_Hz==6.1 & Stim_dur_pc_cyc==49))
data6$fitted<-fitted(out_6.1_50)
# structure of fitted data prevents grouping by Bird_code therefore specify 
# using positions that correspond to each individual 
# Check which fitted values belong to which subject 
data6$Bird_code[1:22]
data6$Bird_code[23:44]
# data6$fitted[1:22] is COLLI-031 aka Bird 1
# data6$fitted[23:44] is COLLI-036 aka Bird 4

summaryPower6 <- data.frame( freq = 6.1, Subject=c("Bird 1", "Bird 4"),
	mean = c(mean(data6$fitted[1:22]), mean(data6$fitted[23:44])),
	min = c(range(data6$fitted[1:22])[1], range(data6$fitted[23:44])[1]),
	max = c(range(data6$fitted[1:22])[2], range(data6$fitted[23:44])[2]),
	range = c(diff(range(data6$fitted[1:22])),diff(range(data6$fitted[23:44]))),
	Act = c((length(data6$Phase_pc[which(data6$fitted[1:22] >= 0)])/(length(unique(data6$Phase_pc))))*100,
		 (length(data6$Phase_pc[which(data6$fitted[23:44] >= 0)])/(length(unique(data6$Phase_pc))))*100)  )

#### calculate the summary statistics for 8.6 Hz ####
head(out_8.6_50)

data8 <- subset(data, c(Freq_Hz==8.6 & Stim_dur_pc_cyc==50))
data8$fitted<-fitted(out_8.6_50)
# structure of fitted data prevents grouping by Bird_code therefore specify 
# using positions that correspond to each individual
# Check which fitted values belong to which subject
data8$Bird_code[1:22]
data8$Bird_code[23:44]
# data8$fitted[1:22] is COLLI-039 aka Bird 5
# data8$fitted[23:44] is COLLI-041 aka Bird 6

summaryPower8 <- data.frame( freq = 8.6, Subject=c("Bird 5", "Bird 6"),
	mean = c(mean(data8$fitted[1:22]), mean(data8$fitted[23:44])),
	min = c(range(data8$fitted[1:22])[1], range(data8$fitted[23:44])[1]),
	max = c(range(data8$fitted[1:22])[2], range(data8$fitted[23:44])[2]),
	range = c(diff(range(data8$fitted[1:22])),diff(range(data8$fitted[23:44]))),
	Act = c((length(data8$Phase_pc[which(data8$fitted[1:22] >= 0)])/(length(unique(data8$Phase_pc))))*100,
		 (length(data8$Phase_pc[which(data8$fitted[23:44] >= 0)])/(length(unique(data8$Phase_pc))))*100) )

#### calculate the summary statistics for 10.1 Hz ####
head(out_10.1_50)

data10 <- subset(data, c(Freq_Hz==10.1 & Stim_dur_pc_cyc==50))
data10$fitted<-fitted(out_10.1_50)
# structure of fitted data prevents grouping by Bird_code therefore specify 
# using positions that correspond to each individual
# Check which fitted values belong to which subject
data10$Bird_code[1:22]
data10$Bird_code[23:44]
data10$Bird_code[45:66]
# data10$fitted[1:22] is COLLI-042 aka Bird 7
# data10$fitted[23:44] is COLLI-045 aka Bird 8
# data10$fittef[45:66] is COLLI-046 aka Bird 9

summaryPower10 <- data.frame( freq = 10.1, Subject=c("Bird 7","Bird 8","Bird 9"),
	mean = c(mean(data10$fitted[1:22]), mean(data10$fitted[23:44]),
		 mean(data10$fitted[45:66])),
	min = c(range(data10$fitted[1:22])[1], range(data10$fitted[23:44])[1],
		 range(data10$fitted[45:66])[1]),
	max = c(range(data10$fitted[1:22])[2], range(data10$fitted[23:44])[2],
		 range(data10$fitted[45:66])[2]),		 
	range = c(diff(range(data10$fitted[1:22])),diff(range(data10$fitted[23:44])),
		  diff(range(data10$fitted[45:66]))),
	Act = c((length(data10$Phase_pc[which(data10$fitted[1:22] >= 0)])/(length(unique(data10$Phase_pc))))*100,
		 (length(data10$Phase_pc[which(data10$fitted[23:44] >= 0)])/(length(unique(data10$Phase_pc))))*100,
		 (length(data10$Phase_pc[which(data10$fitted[45:66] >= 0)])/(length(unique(data10$Phase_pc))))*100) )

## Calculate means and standard errors of each summary statistic for reporting 
## in results

freq.summstatstable <- data.frame( freq = c(6.1, 8.6, 10.1),
avg_min = c(mean(summaryPower6$min), mean(summaryPower8$min),
	     mean(summaryPower10$min)),
SE_min = c(std.error(summaryPower6$min), std.error(summaryPower8$min),
	    std.error(summaryPower10$min)),
avg_max = c(mean(summaryPower6$max), mean(summaryPower8$max),
	     mean(summaryPower10$max)),
SE_max = c(std.error(summaryPower6$max), std.error(summaryPower8$max),
	    std.error(summaryPower10$max)),
avg_mean = c(mean(summaryPower6$mean), mean(summaryPower8$mean),
	     mean(summaryPower10$mean)),
SE_mean = c(std.error(summaryPower6$mean), std.error(summaryPower8$mean),
	    std.error(summaryPower10$mean)),
avg_act = c(mean(summaryPower6$Act), mean(summaryPower8$Act),
	     mean(summaryPower10$Act)),
SE_act = c(std.error(summaryPower6$Act), std.error(summaryPower8$Act),
	    std.error(summaryPower10$Act)) )

## create a single dataframe containing all of the summary stats of interest
Pow.summ.stats2 <-rbind(summaryPower6, summaryPower8, summaryPower10)

### COMPARE SUMMARY STATS ###
### means ###
means <-lm(mean~freq, data = Pow.summ.stats2)
hist(residuals(means))
shapiro.test(residuals(means)) # good
summary(means) # look at coefficients
out.freq.means<-anova(means)

### max ###
max.z <-lm(max~freq, data = Pow.summ.stats2)
hist(residuals(max.z))
shapiro.test(residuals(max.z)) # good
summary(max.z) # look at coefficients
out.freq.max<-anova(max.z)

### min ###
min.z <-lm(min~freq, data = Pow.summ.stats2)
hist(residuals(min.z))
shapiro.test(residuals(min.z)) # good
summary(min.z) # look at coefficients
out.freq.min<-anova(min.z)

### Act ###
Acts <-lm(Act~freq, data = Pow.summ.stats2)
hist(residuals(Acts))
shapiro.test(residuals(Acts)) # good
summary(Acts) # look at coefficients
out.freq.acts<-anova(Acts)

# write a row to append to the summary stats table, containing the p-values 
# from the comparison of that summary stat across the two stimulus durations
summ.p_values<-cbind(out.freq.min$Pr[1], out.freq.max$Pr[1], out.freq.means$Pr[1], out.freq.acts$Pr[1]) 
# append to summary stats table
freq.summstatstable[4,c(2,4,6,8)]<-summ.p_values
# add row name
freq.summstatstable[4,1]<-"p-values-summary stats comparisons"
# add column with phase effects p-values
freq.summstatstable$p_values.phase_effects<-c(out.freq.phases6$'p-value'[2],out.freq.phases8$'p-value'[2],out.freq.phases10$'p-value'[2],NA)

# write csv of summary stats and p-values
write.csv(freq.summstatstable,"summary stats and p-values--frequencies.csv",row.names=FALSE)

################ PLOTTING #################

# SET UP THE PLOTTING WINDOW
plot.mat <- matrix( 1, ncol=2, nrow=6 )

plot.mat[3:4] <- 2
plot.mat[5:6] <- 3
plot.mat[7:8] <- 4
plot.mat[9:10] <- 5
plot.mat[11:12] <- 6


# some global options
ylims=c(-800,200)
xlims=c(-50,50)
presid.cex = 0


#### PLOTTING FOR FIGURE 6 ####
layout(plot.mat)
par( mar=c(2,5,1,1) )

# visreg used for its plotting format and for visualising GAM fit (if desired)--
# to see model fit, change visreg line parameters
## If model fit line plotting turned off, plotting error returned which can be ignored.
# Plot 6.1 Hz, 50% cyc 
visreg( out_6.1_50, "Phase_pc",
      frame=FALSE,
      xaxt="n", yaxt="n",
      ylab="", xlab="",
      ylim=ylims, xlim=xlims,
      points=list(cex=0),
      band=F,
      line=list(col="darkgreen", lwd=0, lty=" ")) # remove lwd= and lty= to see fit
# check that symbols are correct (pch 15 --bird 1 and 17 -- bird 4
points( norm.Power_W..kg ~ Phase_pc,
        col="darkgreen", pch=c((as.numeric(factor(Bird_code)[1:22]))+14,
	 (as.numeric(factor(Bird_code)[23:44]))+15), cex=0.75,
        data=subset(data, Freq_Hz==6.1 & Stim_dur_pc_cyc==49) )
#Add legend
legend(-48,-550, legend = c("Bird 1","Bird 4"), 
        pch = c(15,17), col = "darkgreen", bty="n")
axis(2, at=seq(ylims[1],ylims[2], by=ylims[2]), labels=TRUE, tcl=0.5)
abline(h=0, lty="dashed")
mtext(side = 2, at=-200, text = "6.1 Hz, 50% cycle", col="black", line =3.5)

# plot 8.6 Hz, 50% cyc
visreg( out_8.6_50, "Phase_pc",
        frame=FALSE,
        xaxt="n", yaxt="n",
        ylab="", xlab="",
        ylim=ylims, xlim=xlims,
        points=list(cex=0),
        band=F,
        line=list(col="red",lwd=0, lty=" ")) # remove lwd= and lty= to see fit
points( norm.Power_W..kg ~ Phase_pc,
        col="red", pch=(as.numeric(factor(Bird_code))), cex=0.75,
        data=subset(data, Freq_Hz==8.6 & Stim_dur_pc_cyc==50) )
#Add legend
legend(-48,-550, legend = c("Bird 5","Bird 6"), 
        pch = c(1,2), col = "red", bty="n")
axis(2, at=seq(ylims[1],ylims[2], by=ylims[2]), labels=TRUE, tcl=0.5)
abline(h=0, lty="dashed")
mtext(side = 2, at=-200, text = "8.6 Hz, 50% cycle", line =3.5)
mtext(side = 2, at=-200, text = "Power (W/kg)", line =2.2, cex=0.75)

# plot 10.1 Hz
visreg( out_10.1_50, "Phase_pc",
        frame=FALSE,
        xaxt="n", yaxt="n",
        ylab="", xlab="",
        ylim=ylims, xlim=xlims,
        points=list(cex=0),
        band=F,
	 line=list(col="blue", lwd=0, lty=" ")) # remove lwd= and lty= to see fit
points( norm.Power_W..kg ~ Phase_pc,
        col="blue", pch=c((as.numeric(factor(Bird_code)[1:22]))-1,
	 (as.numeric(factor(Bird_code)[23:44]))+3, 
	 (as.numeric(factor(Bird_code)[45:66]))+3), cex=0.75,
        data=subset(data, Freq_Hz==10.1 & Stim_dur_pc_cyc==50) )
#Add legend
legend(-48,-450, legend = c("Bird 7","Bird 8","Bird 9"), 
        pch = c(0,5,6), col = "blue", bty="n")
axis(2, at=seq(ylims[1],ylims[2], by=ylims[2]), labels=TRUE, tcl=0.5)
abline(h=0, lty="dashed")
mtext(side = 2, at=-200, text = "10.1 Hz, 50% cycle", line =3.5)

axis(1, at=seq(-50,50,by=25), labels=TRUE, tcl=0.5)

# plot mean vs frequency
plot( mean~freq, data=Pow.summ.stats2,
      cex=1,
      col=c("darkgreen","red","blue")[as.factor(Pow.summ.stats2$freq)],
      pch=c(15,17,1,2,0,5,6)[as.factor(Pow.summ.stats2$Subject)], bty="n", tcl=0.5,
      ylim=c(-250,5), xlim=c(4,12),
      xaxt="n",
      ylab="Mean power (W/kg)" )
abline(h=0, lty="dashed")
axis(1, at=c(4,6.1,8.6,10.1,12), labels=FALSE, tcl=0.5)

# plot Max & Min power vs frequency
plot( max~freq, data=Pow.summ.stats2,
      cex=1,
      col=c("darkgreen","red","blue")[as.factor(Pow.summ.stats2$freq)],
      pch=c(15,17,1,2,0,5,6)[as.factor(Pow.summ.stats2$Subject)], bty="n", tcl=0.5,
      ylim=c(-600,200), xlim=c(4,12),
      xaxt="n",
      ylab="Power (W/kg)" )
points( min~freq, data=Pow.summ.stats2,
	cex=1,
	col=c("darkgreen","red","blue")[as.factor(Pow.summ.stats2$freq)],
	pch=c(15,17,1,2,0,5,6)[as.factor(Pow.summ.stats2$Subject)])
abline(h=0, lty="dashed")
axis(1, at=c(4,6.1,8.6,10.1,12), labels=FALSE, tcl=0.5)

# plot actuator range vs frequency, coloured by frequency 
plot( Act~freq, data=Pow.summ.stats2,
      cex=1,
      col=c("darkgreen","red","blue")[as.factor(Pow.summ.stats2$freq)],
      pch=c(15,17,1,2,0,5,6)[as.factor(Pow.summ.stats2$Subject)], bty="n", tcl=0.5,
      xlim=c(4,12), ylim=c(0,50),      
      xaxt="n",
      xlab="Frequency (Hz)", ylab="Actuator range (% cycle)")
axis(1, at=c(4,6.1,8.6,10.1,12), labels=TRUE, tcl=0.5)




