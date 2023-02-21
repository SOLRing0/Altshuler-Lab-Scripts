## Script to get time data for each of the trials within an experiment, using their file timestamps

############################## Get all relevant files ##############################################

# get working directory
getwd()
# change working directory
setwd(choose.dir())
# get a list of data folders from relevant dates
folders <- list.dirs(path = dir(pattern = "^2017-0", full.names = T, all.files=F, 
	recursive=FALSE, ignore.case=T, include.dirs=F))

# get a list of relevant subfolders within each of those data folders and remove unwanted subfolders
subfolders <- subset(folders, !grepl("length-tension|length-tension2", folders))
# list only treatment folders, not experiment folders that contain them
subfolders1 <- list.files(subfolders, pattern ="cycle",  
	full.names=T, recursive=FALSE, ignore.case=T, include.dirs=F)

# get a list of all the *.ddf files to be analysed
ddf_files <- list.files(subfolders, pattern ="*.ddf", 
	full.names= T, recursive= F, ignore.case= T, include.dirs= F)
# ignore WL, twitch and tetanus files -- only want trial ddf files
files <- subset(ddf_files, !grepl("WL|twitch|tetanus", ddf_files))

############################## Loop through and get timestamps ######################################
# Write a function to summarize the data in the *.ddf files
timedata <- function(filename) {
  Time <-file.mtime(filename)
  data <- matrix(c(filename, Time), nrow=length(filename), ncol=2, byrow=F)
  }

# loop through all the *.ddf files and pull the time stamps 
out<-lapply(files, function(filename) {
	x <- timedata(filename)})
# create a data frame containing the summary data from all *.ddf files of interest
Timestamps<-data.frame(matrix(unlist(out), nrow= length(files), ncol= 2, byrow= T), stringsAsFactors = F)# rename the columns to properly identify the metrics
colnames(Timestamps) <- c("Filename", "Time")
head(Timestamps)

# Timestamps are output in seconds a default origin  
# check file creation dates and times
dts = Timestamps$Time
mydates = structure(dts,class=c('POSIXt','POSIXct'))
mydates

###### Processing file paths and names to sort data by phase (same as power summary data) #######
# Split file pathnames 
z<-data.frame(matrix(unlist(strsplit(Timestamps$Filename, "[/]")), nrow= length(files), ncol= 4, byrow= T), stringsAsFactors = F)
timedat<-subset(z, select = c(X2, X3, X4))
colnames(timedat) <- c("Date", "Treatment", "File")

# get everything in the filename before the '%' (end of phase in filename)
z1<-data.frame(matrix(unlist(strsplit(timedat$File, "%")), nrow=length(files), ncol=2, byrow=T),stringsAsFactors = F)
timedat$File<-z1$X1

# get everything in the filename after the second '-' (indicates start of phase)
## save as numeric in new vector of trial phase 
timedat$t_phase<-as.numeric(gsub('^(?:[^-]*-){2}','',timedat$File))

# add timestamp column to timedat dataframe before sorting by ascending phase
timedat$Time<-Timestamps$Time

# Sort timedat by Date, then by Treatment, and then by t_phase
s_timedat<-timedat[with(timedat, order(Date, Treatment, t_phase)),]
# check structure of new dataframe
str(s_timedat)
# Time is class 'chr'--convert to numeric before passing to power csv
s_timedat$Time<-as.numeric(s_timedat$Time)
# check
str(s_timedat)

########## Transform timestamps to reflect times relative to start of each experiment (Date) ########
library(plyr)
# group data by experiment date and then mutate the Time to determine when relative to the start of 
# each experiment each trial was run (in seconds).

x <- ddply(s_timedat, .(Date), mutate,
	time_since_start.s = Time - min(Time)
	)

# read in power by frequency and phase data csv to add the fully processed timestamp data.
data <- read.csv("pigeon HT power by freq and phase.csv", stringsAsFactors=FALSE)
head(data)

# assign the vector containing the time calculated times to a new vector within 
# the power by phase and frequency data frame.
data$Time_Since_Start_s <- x$time_since_start.s
# check it
head(data)

# write to csv - overwrites old power by frequency and phase data csv
write.csv(data,"pigeon HT power by freq and phase.csv", row.names=FALSE)




