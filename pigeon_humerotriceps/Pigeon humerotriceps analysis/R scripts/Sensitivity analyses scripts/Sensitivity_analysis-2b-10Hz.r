### SECOND SENSITIVITY ANALYSIS--COMPARE INDIVIDUALS AT 10.1 Hz, 50% duty cycle ### 
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

out_10.1_50_b7 <- gam( I(norm.Power_W..kg) ~ 
                   s(Phase_pc, bs="cc") + #factor(Bird_code) +
                   Time_Since_Start_s + Mass_g + HuTri_length_mm,
             data=subset(data, c(Freq_Hz==10.1 & Stim_dur_pc_cyc==50 & Bird_code=="COLLI-042")) )

out_10.1_50_b8 <- gam( I(norm.Power_W..kg) ~ 
                   s(Phase_pc, bs="cc") + #factor(Bird_code) +
                   Time_Since_Start_s + Mass_g + HuTri_length_mm,
             data=subset(data, c(Freq_Hz==10.1 & Stim_dur_pc_cyc==50 & Bird_code=="COLLI-045")) )

out_10.1_50_b9 <- gam( I(norm.Power_W..kg) ~ 
                   s(Phase_pc, bs="cc") + #factor(Bird_code) +
                   Time_Since_Start_s + Mass_g + HuTri_length_mm,
             data=subset(data, c(Freq_Hz==10.1 & Stim_dur_pc_cyc==50 & Bird_code=="COLLI-046")) )

# check the reasonableness of the fits
gam.check(out_10.1_50_b7)
gam.check(out_10.1_50_b8)
gam.check(out_10.1_50_b9)

plot(out_10.1_50_b7,pages=1)
plot(out_10.1_50_b8,pages=1)
plot(out_10.1_50_b9,pages=1)

################ PLOTTING #################

# SET UP THE PLOTTING WINDOW
plot.mat <- matrix( 1, ncol=1, nrow=4 )

plot.mat[1] <- 0
plot.mat[2:3] <- 1
plot.mat[4] <- 0

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
# Plot 10.1 Hz, 50% cyc for bird 7 
visreg( out_10.1_50_b7, "Phase_pc",
      frame=FALSE,
      xaxt="n", yaxt="n",
      ylab="", xlab="",
      ylim=ylims, xlim=xlims,
      points=list(cex=0),
      band=T,
      fill.par=list(col="grey50"),
      line=list(col="black", lwd=0.75)) # remove lwd= and lty= to see fit

points( norm.Power_W..kg ~ Phase_pc,
        col="black", pch=0, cex=0.75,
        data=subset(data, Freq_Hz==10.1 & Stim_dur_pc_cyc==50 & Bird_code=="COLLI-042") )
#Add legend
legend(-48,-400, legend = "Bird 7", 
        pch = 0, col = "black", bty="n")
axis(2, at=seq(ylims[1],ylims[2], by=ylims[2]), labels=TRUE, tcl=0.5)
abline(h=0, lty="dashed")
mtext(side = 2, at=-200, text = "10.1 Hz, 50% cycle", col="black", line =3.5)

par(new=T)
# Plot 10.1 Hz, 50% cyc for bird 8 
visreg( out_10.1_50_b8, "Phase_pc",
      frame=FALSE,
      xaxt="n", yaxt="n",
      ylab="", xlab="",
      ylim=ylims, xlim=xlims,
      points=list(cex=0),
      band=T,
      fill.par=list(col="aquamarine"),
      line=list(col="deepskyblue3", lwd=0.75, lty=2)) # remove lwd= and lty= to see fit

points( norm.Power_W..kg ~ Phase_pc,
        col="deepskyblue3", pch=5, cex=0.75,
        data=subset(data, Freq_Hz==10.1 & Stim_dur_pc_cyc==50 & Bird_code=="COLLI-045") )
#Add legend
legend(-48,-450, legend = "Bird 8", 
        pch = 5, col = "deepskyblue3", bty="n")
axis(2, at=seq(ylims[1],ylims[2], by=ylims[2]), labels=TRUE, tcl=0.5)
abline(h=0, lty="dashed")
mtext(side = 2, at=-200, text = "10.1 Hz, 50% cycle", col="deepskyblue3", line =3.5)

par(new=T)
# plot 10.1 Hz for bird 9
visreg( out_10.1_50_b9, "Phase_pc",
        frame=FALSE,
        xaxt="n", yaxt="n",
        ylab="", xlab="",
        ylim=ylims, xlim=xlims,
        points=list(cex=0),
        band=T,
        fill.par=list(col="steelblue1"),
	 line=list(col="blue", lwd=0.75, lty=3)) # remove lwd= and lty= to see fit
points( norm.Power_W..kg ~ Phase_pc,
        col="blue", pch=6, cex=0.75,
        data=subset(data, Freq_Hz==10.1 & Stim_dur_pc_cyc==50 & Bird_code=="COLLI-046") )
#Add legend
legend(-48,-500, legend = "Bird 9", 
        pch = 6, col = "blue", bty="n")
axis(2, at=seq(ylims[1],ylims[2], by=ylims[2]), labels=TRUE, tcl=0.5)
abline(h=0, lty="dashed")
mtext(side = 2, at=-200, text = "10.1 Hz, 50% cycle", line =3.5)

axis(1, at=seq(-50,50,by=10), labels=TRUE, tcl=0.5)

