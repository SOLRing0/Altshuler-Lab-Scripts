### SECOND SENSITIVITY ANALYSIS--COMPARE INDIVIDUALS AT 8.6 Hz, 50% duty cycle ###
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

out_50pc_b5 <- gam( norm.Power_W..kg ~ 
                  s(Phase_pc, bs="cc") + #factor(Bird_code) +
                  Time_Since_Start_s + Mass_g + HuTri_length_mm,
             data=subset(data, c(Freq_Hz==8.6 & Stim_dur_pc_cyc==50 & Bird_code=="COLLI-039")) )

# check the reasonableness of the fits
gam.check(out_50pc_b6)
gam.check(out_50pc_b5)

plot(out_50pc_b6,pages=1)
plot(out_50pc_b5,pages=1)

# SET UP THE PLOTTING WINDOW
plot.mat <- matrix( 1, ncol=1, nrow=4 )

plot.mat[1] <-0
plot.mat[2:3] <- 1
plot.mat[4] <- 0

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

# Plot 8.6 Hz at 50% norm.stim duration for bird 6
visreg( out_50pc_b6, "Phase_pc",
      frame=FALSE,
      xaxt="n", yaxt="n",
      xlab="Stimulus onset (% cycle)",
      ylab=" ",
      ylim=ylims, xlim=xlims,
      points=list(cex=0),
      band=T, 
      fill.par=list(col="hotpink1"),
      line=list(col="red", lwd=0.75)) # remove lwd= and lty= to see model fit
#Add axes manually
axis(1, at=seq(-50,50,by=10), labels=TRUE, tcl=0.5)
points( norm.Power_W..kg ~ Phase_pc,
        col="red", pch=2, cex=0.75,
        data=subset(data, Freq_Hz==8.6 & Stim_dur_pc_cyc==50 & Bird_code=="COLLI-041") )
#Add axes manually
axis(2, at=seq(ylims[1],ylims[2], by=ylims[2]), labels=TRUE, tcl=0.5)
abline(h=0, lty="dashed")

# Add legend
legend(-48,-450, legend = "Bird 6", 
        pch = 2, col = "red", bty="n")
mtext(side = 2, at=-200, text = "8.6 Hz, 50% cycle", line =3.5)

# Plot 8.6 Hz at 50% norm.stim duration for bird 5
par(new=T)
visreg( out_50pc_b5, "Phase_pc",
      frame=FALSE, 
      xaxt="n", yaxt="n",
      xlab="Stimulus onset (% cycle)",
      ylab=" ",
      ylim=ylims, xlim=xlims,
      points=list(cex=0),
      band=T,
      fill.par=list(col="yellow"),
      line=list(col="darkorange2", lwd=0.75, lty=2)) # remove lwd= and lty= to see model fit
#Add axes manually
axis(1, at=seq(-50,50,by=10), labels=TRUE, tcl=0.5)
points( norm.Power_W..kg ~ Phase_pc,
        col="darkorange2", pch=1, cex=0.75,
        data=subset(data, Freq_Hz==8.6 & Stim_dur_pc_cyc==50 & Bird_code=="COLLI-039") )
#Add axes manually
axis(2, at=seq(ylims[1],ylims[2], by=ylims[2]), labels=TRUE, tcl=0.5)
abline(h=0, lty="dashed")

# Add legend
legend(-48,-400, legend = "Bird 5", 
        pch = 1, col = "darkorange2", bty="n")
mtext(side = 2, at=-200, text = "8.6 Hz, 50% cycle", line =3.5)

