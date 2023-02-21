### Writes a function to read in twitch0XX.ddf files and obtain all the necessary parameters from the metadata

read.ddf <- 
  function( filename, 
          renameColumns=list( c(2,3), c("Position","Force")),
          deleteColumns=4:11)
{
  
  # pull useful metadata from the header
  header <- readLines(filename, n=25)
  
  Sample_Frequency <- as.numeric(sub( "Sample Frequency \\(Hz\\): ", "", header[2] ))
  Reference_Area <- as.numeric( sub( " sq. mm","", sub( "Reference Area: ", "", header[3] ) ))
  Reference_Force <- as.numeric( sub( " mN","", sub( "Reference Force: ", "", header[4] ) ))
  Reference_Length <- as.numeric( sub( " mm","", sub( "Reference Length: ", "", header[5] ) ))
  
  calibration <- strsplit(header[7],"\t",fixed=TRUE)[[1]][-1]
  units <- strsplit(header[8],"\t",fixed=TRUE)[[1]][-1]
  scale <- as.numeric(strsplit(header[9],"\t",fixed=TRUE)[[1]][-1]) # units per volt
  offset <- as.numeric(strsplit(header[10],"\t",fixed=TRUE)[[1]][-1]) # volts
  tads <- strsplit(header[11],"\t",fixed=TRUE)[[1]][-1] # volts
  
  protocol <- sapply(header[17:18], function(x) paste0(x,"\n"))
  
  # last line of header is the column names of the data
  colNames <- strsplit(header[23],"\t",fixed=TRUE)[[1]]
  
  # read the body and convert to matrix
  body <- scan( filename, skip=23, what="numeric", sep="\t", quiet=TRUE )
  bodymat <- matrix( as.numeric(body), ncol=length(colNames), byrow=TRUE )
  
  colnames(bodymat) <- colNames
  
  # rename columns, if desired
  if(!is.null(renameColumns))
    colNames[renameColumns[[1]]] <- renameColumns[[2]]
  
  # rescale the data from volts to units
  bodyUnits <- sapply( 1:length(calibration), function(i)
    ( bodymat[,calibration[i]] + offset[i] ) * scale[i] )
  
  colnames(bodyUnits) <- colNames[-c(1,length(colNames))]
  
  bodydf <- data.frame( Time = (1:nrow(bodyUnits))/Sample_Frequency,
                        bodyUnits,
                        Stim = bodymat[,"Stim"] )
  
  bodydf <- bodydf[,-deleteColumns]
  units <- c("s", units[-deleteColumns+1], "TTL") # +1 b/c no Sample column
  
  attr(bodydf,"Sample_Frequency") <- Sample_Frequency
  attr(bodydf,"Reference_Area") <- Reference_Area
  attr(bodydf,"Reference_Force") <- Reference_Force
  attr(bodydf,"Reference_Length") <- Reference_Length
  attr(bodydf,"units") <- units
  
  bodydf
}

# assign matrix to an object that can then be analysed 
## Time is in seconds, Postion is in millimeters, Force is in milliNewtons and Stim is 
## all (1) or none (0)

#### PLOT SAMPLE ISOMETRIC CONTRACTION TRACE ####
# get working directory
getwd()
# change working directory
setwd(choose.dir()) # select 'Experiments' folder
# Use '2017-06-27/twitch003.ddf' for figure 2A
twitch<-read.ddf('2017-06-27/twitch003.ddf')
head(twitch)
# adjust forces for top hole of motor arm and convert from mN to N
twitch$Force <- ((twitch$Force)*4)/1000
# convert time to ms
twitch$Time <- twitch$Time*1000
plot(twitch$Force ~ twitch$Time, xlim=c(0,300), ylab="Force (N)", xlab="Time (ms)", 
	bty="n", type="l", col="gray35", lwd=2.5)
# add lines to mark 90% and 10% peak force
abline(h=(0.9*(max(twitch$Force))), lty="dashed")
abline(h=(0.1*(max(twitch$Force))), lty="dashed")

########### Summarize data from isometric contractions (twitch) ###############
getwd() # make sure returned to 'Experiments' folder
# get a list of data folders from relevant dates
folders <- list.dirs(path = dir(pattern = "^2017-0", full.names = T, all.files=F, 
	recursive=FALSE, ignore.case=T, include.dirs=F))

# get a list of relevant subfolders within each of those data folders and remove unwanted subfolders
subfolders <- subset(folders, !grepl("length-tension|length-tension2", folders))

# get a list of all the 'twitch' files within each of those data folders
allfiles <- list.files(subfolders, pattern ="^twitch" , 
	full.names=T, recursive=FALSE, ignore.case=T, include.dirs=F)
# get a list of only the 'twitch0XX.ddf' files to be analysed
files <- subset(allfiles, grepl(".ddf", allfiles))

# Write a function to summarize the data in the *.ddf files

twitchkin <- function(filename) {
  # read in *.ddf file using function defined above.
  dat <- read.ddf(filename)
  # Adjust Forces for top hole on motor arm
  dat$Force <- dat$Force*4
  # convert time to ms
  dat$Time <- dat$Time*1000
  # Find the stimuli times and their associated forces 
  stim.period <- subset(dat, Stim == 1, select = c("Time", "Force"))
  # Find the maximum force output and time that it occurs (in ms)
  PF <- subset(dat, Force == max(dat$Force) & 
	Time == min(dat$Time[which(dat$Force==max(dat$Force))]), select = c("Time", "Force"))
  # split twitch into two categories
  Forcerise <- subset(dat, Time <= PF$Time, select = c("Time", "Force"))
  Forcefall <- subset(dat, Time >= PF$Time, select = c("Time", "Force"))
  # Find 10% and 90% PF
  # sort forces 'non-decreasingly'
  sortFR <- Forcerise[order(Forcerise$Force),]
  sortFF <- Forcefall[order(Forcefall$Force),]
  # Developed force calculations
  PF10 <- sortFR[findInterval(PF$Force*0.1, sortFR$Force),]
  PF90 <- sortFR[findInterval(PF$Force*0.9, sortFR$Force),]
  # Calculate the change in force from resting to peak force (developed force) 
  DF <- PF90$Force-PF10$Force
  # and the time it takes to reach 90% of peak force from rest 
  DFT <- PF90$Time-PF10$Time
  # Relaxation calculations
  Rx10 <- sortFF[findInterval(PF$Force*0.1, sortFF$Force),]
  Rx50 <- sortFF[findInterval(PF$Force*0.5, sortFF$Force),]
  Rx90 <- sortFF[findInterval(PF$Force*0.9, sortFF$Force),]
  # calculate the time it takes to relax following peak force  
  RD10 <- Rx10$Time-Rx90$Time
  RD50 <- Rx50$Time-Rx90$Time
  # Create a matrix and populate it with the results of the above calculations of interest for each
  # *.ddf file.
  matrix(c(filename, PF$Force, DF, DFT, RD10, RD50), nrow = length(filename), ncol=6 , byrow= F)
  }

# loop through all the *.ddf files and get the desired summary data 
out<-lapply(files, function(filename) {
	x <- twitchkin(filename)})
# create a data frame containing the summary data from all *.ddf files of interest
data<-data.frame(matrix(unlist(out), nrow= length(files), ncol= 6, byrow= T), stringsAsFactors = F)
# rename the columns to properly identify the metrics
colnames(data) <- c("Date.and.Treatment", "Peak.force", "Developed.force", "Time.to.DF", "Time.to.90pc.Relax", "Time.to.50pc.Relax")

#check data structure--make sure data are numeric since converting from list to dataframe
str(data) # numbers are in character format--convert data to numeric
data$Peak.force <-as.numeric(data$Peak.force)
data$Developed.force <-as.numeric(data$Developed.force)
data$Time.to.DF <-as.numeric(data$Time.to.DF)
data$Time.to.90pc.Relax <-as.numeric(data$Time.to.90pc.Relax)
data$Time.to.50pc.Relax <-as.numeric(data$Time.to.50pc.Relax)

# check it
str(data)

# reassign to new name for summary, plotting and analysis
tkin <- data

head(tkin)
# Replace the file names with the subject IDs 
levels(tkin$Date.and.Treatment) <- c(levels(tkin$Date.and.Treatment), "Bird 1", "Bird 2",
	"Bird 3", "Bird 4", "Bird 5", "Bird 6", "Bird 7", "Bird 8", "Bird 9")
tkin$Date.and.Treatment[grepl("2017-04-05", tkin$Date.and.Treatment)]<-"Bird 1"
tkin$Date.and.Treatment[grepl("2017-04-10", tkin$Date.and.Treatment)]<-"Bird 3"
tkin$Date.and.Treatment[grepl("2017-04-12", tkin$Date.and.Treatment)]<-"Bird 4"
tkin$Date.and.Treatment[grepl("2017-04-19", tkin$Date.and.Treatment)]<-"Bird 2"
tkin$Date.and.Treatment[grepl("2017-06-06", tkin$Date.and.Treatment)]<-"Bird 5"
tkin$Date.and.Treatment[grepl("2017-06-11", tkin$Date.and.Treatment)]<-"Bird 6"
tkin$Date.and.Treatment[grepl("2017-06-14", tkin$Date.and.Treatment)]<-"Bird 7"
tkin$Date.and.Treatment[grepl("2017-06-27", tkin$Date.and.Treatment)]<-"Bird 8"
tkin$Date.and.Treatment[grepl("2017-06-28", tkin$Date.and.Treatment)]<-"Bird 9"

# check it
head(tkin)
tail(tkin)

# Rename the column
colnames(tkin)[1] <- "Subject_ID"
head(tkin)
# convert Subject ID to factor
str(tkin)
tkin$Subject_ID<-as.factor(tkin$Subject_ID)
# check it
str(tkin)

#####################################################################################
####### Summary and analysis of pigeon humerotriceps twitch kinetics data ###########
#####################################################################################
# estimate full developed force duration 
tkin$Full_DF.time <- tkin$Time.to.DF * 1.2
# estimate full relaxation duration
tkin$Full_Relax.time <- tkin$Time.to.90pc.Relax * 1.2
# estimate total duration for both contraction and relaxation
tkin$Full_CnR.time <- tkin$Full_DF.time + tkin$Full_Relax.time
# calculate the maximum contractile frequency 
tkin$max_freq <- 1/(tkin$Full_CnR.time/1000)
# calculate the mean max frequency by individual 
tapply(tkin$max_freq, INDEX=tkin$Subject_ID, FUN=mean)

# Calculate mean and SE for each metric by individual
library(plotrix)

## Mean time to DF by individual
T2PF <- with(tkin, aggregate(Time.to.DF ~ Subject_ID, FUN = function(x) 
	 (AVG = mean(x)) ))
## SE time to DF by individual
T2PF.SE <- with(tkin, aggregate(Time.to.DF ~ Subject_ID, FUN = function(x) 
	 (SE = std.error(x)) ))
## Mean 90% relaxation duration by individual
T2R.90 <- with(tkin, aggregate(Time.to.90pc.Relax ~ Subject_ID, FUN = function(x) 
	(AVG = mean(x)) ))
## SE 90% relaxation duration by individual
T2R.90.SE <- with(tkin, aggregate(Time.to.90pc.Relax ~ Subject_ID, FUN = function(x) 
	(SE = std.error(x)) ))

# Combine two summary data frames into one
tk <- cbind(T2PF, T2PF.SE$Time.to.DF, T2R.90$Time.to.90pc.Relax, T2R.90.SE$Time.to.90pc.Relax)
# Rename the columns
colnames(tk) <- c("Subject_ID","m.t2pf","se.t2pf","m.t2r90","se.t2r90") 
# Check it
head(tk)

# Check t-test assumptions
hist(tk$m.t2pf)
shapiro.test(tk$m.t2pf) # just barely good
hist(tk$m.t2r90)
shapiro.test(tk$m.t2r90) # good

# Run paired t-test on means
out.tk <- t.test(tk$m.t2r90, tk$m.t2pf, paired = T) 
# write a matrix with the means of time to develop 90% peak force, time to 
# relaxation, and p-value from t-test
t.out<-data.frame(matrix(nrow=1,ncol=3))
t.out$X1<-mean(tk$m.t2pf) 
t.out$X2<-mean(tk$m.t2r90) 
t.out$X3<-out.tk$p.value
#rename columns to write csv of t-test output for twitch kinetics
colnames(t.out)<-c("mean_time_to_develop_90pc_peak_force","mean_time_to_relax_following_90pc_peak_force","p-value")

## write text file of twitch kinetics t-test results/p-value
write.table(t.out, "twitch kinetics means and t-test p-value.txt", row.names=FALSE)

## calculate means and sems for each parameter
# summarize parameters by individual

## mean full developed force duration
z1 <- with(tkin, aggregate(Full_DF.time ~ Subject_ID, FUN = function(x) 
	 (AVG = mean(x)) ))
## SE full developed force duration
z2 <- with(tkin, aggregate(Full_DF.time ~ Subject_ID, FUN = function(x) 
	 (SE = std.error(x)) ))
## mean full relaxation duration
z3 <- with(tkin, aggregate(Full_Relax.time ~ Subject_ID, FUN = function(x) 
	 (AVG = mean(x)) ))
## SE full relaxation duration
z4 <- with(tkin, aggregate(Full_Relax.time ~ Subject_ID, FUN = function(x) 
	 (SE = std.error(x)) ))
## mean total duration for both contraction and relaxation
z5 <- with(tkin, aggregate(Full_CnR.time ~ Subject_ID, FUN = function(x) 
	 (AVG = mean(x)) ))
## SE total duration for both contraction and relaxation
z6 <- with(tkin, aggregate(Full_CnR.time ~ Subject_ID, FUN = function(x) 
	 (SE = std.error(x)) ))
## calculate mean the maximum contractile frequency 
z7 <- with(tkin, aggregate(max_freq ~ Subject_ID, FUN = function(x) 
	 (AVG = mean(x)) ))
## calculate SE the maximum contractile frequency 
z8 <- with(tkin, aggregate(max_freq ~ Subject_ID, FUN = function(x) 
	 (SE = std.error(x)) ))

tk$m.full_DF.time <-z1$Full_DF.time
tk$SE.full_DF.time <-z2$Full_DF.time
tk$m.full_Relax.time <- z3$Full_Relax.time
tk$SE.full_Relax.time <- z4$Full_Relax.time
tk$m.full_CnR.time <- z5$Full_CnR.time
tk$SE.full_CnR.time <- z6$Full_CnR.time
tk$m.max_freq <- z7$max_freq
tk$SE.max_freq <- z8$max_freq

## calculate mean of the means of each parameters and their SEs
## put into one matrix to be appended as the last row of tk
Mean.of.means<-cbind((MoM.t2pf<-mean(tk$m.t2pf)),
(SEM.t2pf<-std.error(tk$m.t2pf)), 
(MoM.t2r90<-mean(tk$m.t2r90)),
(SEM.t2r90<-std.error(tk$m.t2r90)),
(MoM.full_DF.time <-mean(tk$m.full_DF.time)),
(SEM.full_DF.time <-std.error(tk$m.full_DF.time)),
(MoM.full_Relax.time<-mean(tk$m.full_Relax.time)),
(SEM.full_Relax.time <- std.error(tk$m.full_Relax.time)),
(MoM.full_CnR.time<-mean(tk$m.full_CnR.time)),
(SEM.full_CnR.time<-std.error(tk$m.full_CnR.time)),
(MoM.max_freq<-mean(tk$m.max_freq)), 
(SEM.max_freq<-std.error(tk$m.max_freq)) )

# check the new matrix
head(Mean.of.means)
str(Mean.of.means)

## append row name and matrix as last row of tk
tk[10,c(2:13)] <-rbind(Mean.of.means)
# add factor level to give last row its name
levels(tk$Subject_ID)<-c(levels(tk$Subject_ID),"Mean of means")
# add row name
tk[10,1]<-"Mean of means"

## check it
head(tk)
tail(tk)
str(tk)

# Rename the columns for writing csv
colnames(tk) <- c("Subject_ID","mean_time_to_develop_90pc_peak_force.ms","SE_time_to_90pc_pf.ms",
			"mean_time_to_relax_from_90pc_peak_force.ms","SE_time_to_90pc_relax.ms",
			"mean_full_time_to_develop_force.ms","SE_full_time_to_DF.ms",
			"mean_full_time_to_relax.ms","SE_full_time_to_relax.ms",
			"mean_full_contract&relax_time.ms","SE_full_contract&relax_time.ms",
			"mean_max_freq.Hz","SE_max_freq.Hz") 

# check it
head(tk)

#### PLOT Fig.2B -- twitch kinetics by bird ####
# find range of twitch kinetic values to set y-axis limits
range(tkin$Time.to.90pc.Relax)
range(tkin$Time.to.DF)

# initiate plotting for stripchart of twitch kinetics--set y limits from 20-90 to capture all data points
stripchart(1, vertical=TRUE, frame.plot=FALSE, xlab="Subject (Bird #)", ylab="Time (ms)",
	ylim=c(20,90), xlim=c(1,9), xaxt="n")
# add developed force duration points
points(tkin$Time.to.DF~tkin$Subject_ID, pch=c(15,18,16,17,1,2,0,5,6)[as.factor(tkin$Subject_ID)],
	col="gray25")
# add relaxation points
points(tkin$Time.to.90pc.Relax~tkin$Subject_ID, pch=c(15,18,16,17,1,2,0,5,6)[as.factor(tkin$Subject_ID)],
	col="gray55")
# add x-axis to identify individuals by subject ID number (shorter than bird code)
axis(1, at=seq(1,9, by=1))
legend("topleft", legend = c("Developed force", "Relaxation"), pch=15, 
	col= c("gray25", "gray55"), bty="n")


### tk is the means of each parameter by individual--to be used in 
### conjunction with 'tkin' to generate twitch kinetics summary table
### write both data frames to csvs

### Used for Table 2 ###
write.csv(tk, "summarized twitch kinetics and estimates.csv", row.names=FALSE)

## all twitch kinetic data per individual
write.csv(tkin, "twitch kinetics.csv", row.names=FALSE) 

