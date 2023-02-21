### Write a function to read in work loop *.ddf files and obtain all the necessary parameters from the metadata

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
  
  protocol <- sapply(header[17:20], function(x) paste0(x,"\n"))
  
  # last line of header is the column names of the data
  colNames <- strsplit(header[25],"\t",fixed=TRUE)[[1]]
  
  # read the body and convert to matrix
  body <- scan( filename, skip=25, what="numeric", sep="\t", quiet=TRUE )
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

# get working directory
getwd()
# change working directory
setwd(choose.dir()) # select 'Manuscript experiments' folder
# get a list of data folders from relevant dates
folders <- list.dirs(path = dir(pattern = "^2017-0", full.names = T, all.files=F, 
	recursive=FALSE, ignore.case=T, include.dirs=F))

# get a list of relevant subfolders within each of those data folders and remove unwanted subfolders
subfolders1 <- list.files(subfolders, pattern ="cycle",  
	full.names=T, recursive=FALSE, ignore.case=T, include.dirs=F)

### 8.6 Hz ###
# get a list of 8.6 Hz subfolders within each of those data folders
subfolders2 <- subset(subfolders1, grepl("8.6 Hz", subfolders1))

# get a list of all the 8.6 Hz *.ddf files to be analysed
files8 <- list.files(subfolders2, pattern = "*.ddf", 
	full.names=T, recursive=FALSE, ignore.case=T, include.dirs=F)

### 6.1 Hz ###
# get a list of 6.1 Hz subfolders within each of those data folders
subfolders3 <- subset(subfolders1, grepl("6.1 Hz", subfolders1))

# get a list of all 6.1 Hz the *.ddf files to be analysed
files6 <- list.files(subfolders3, pattern = "*.ddf", 
	full.names=T, recursive=FALSE, ignore.case=T, include.dirs=F)

### 10.1 Hz ###
# get a list of 10.1 Hz subfolders within each of those data folders
subfolders4 <- subset(subfolders1, grepl("10.1 Hz", subfolders1))

# get a list of all the 10.1 Hz *.ddf files to be analysed
files10 <- list.files(subfolders4, pattern = "*.ddf", 
	full.names=T, recursive=FALSE, ignore.case=T, include.dirs=F)

# Write a function to summarize the 8.6 Hz data in the *.ddf files

strain8 <- function(filename) {
  # read in *.ddf file using function defined above.
  dat <- read.ddf(filename)
  WLstr <-subset(dat, select = c("Time", "Position"))
  WL.p<-subset(WLstr, Time>=0.03 & Time <= 0.731) # clips the recordings to remove extra recording time on either side of a trial (buffer to ensure full recording obtained).
  WL.p$Position <- (WL.p$Position/4) # calculates the true position--location of the hole on the servo motor arm (1/4 of arm length) causes position values to be quadrupled.
# Create a matrix and populate it with the strain data of each *.ddf file.
  z<-matrix(nrow=length(WL.p$Time), ncol=3, byrow= F)
  z[,1]=rep(filename, length(WL.p$Time))
  z[,2]=as.numeric(WL.p$Time)
  z[,3]=as.numeric(WL.p$Position)
  return(z)
  }

# Write a function to summarize the 6.1 Hz data in the *.ddf files

strain6 <- function(filename) {
  # read in *.ddf file using function defined above.
  dat <- read.ddf(filename)
  WLstr <-subset(dat, select = c("Time", "Position")) 
  WL.p<-subset(WLstr, Time>=0.0455 & Time <= 1.0375) # clips the recordings to remove extra recording time on either side of a trial (buffer to ensure full recording obtained).
  WL.p$Position <- (WL.p$Position/4)
# Create a matrix and populate it with the strain data of each *.ddf file.
  z<-matrix(nrow=length(WL.p$Time), ncol=3, byrow= F)
  z[,1]=rep(filename, length(WL.p$Time))
  z[,2]=as.numeric(WL.p$Time)
  z[,3]=as.numeric(WL.p$Position)
  return(z)
  }

# Write a function to summarize the 10.1 Hz data in the *.ddf files

strain10 <- function(filename) {
  # read in *.ddf file using function defined above.
  dat <- read.ddf(filename)
  WLstr <-subset(dat, select = c("Time", "Position"))
  WL.p<-subset(WLstr, Time>=0.0308 & Time <= 0.6252) # clips the recordings to remove extra recording time on either side of a trial (buffer to ensure full recording obtained).
  WL.p$Position <- (WL.p$Position/4) # calculates the true position--location of the hole on the servo motor arm (1/4 of arm length) causes position values to be quadrupled.
# Create a matrix and populate it with the strain data of each *.ddf file.
  z<-matrix(nrow=length(WL.p$Time), ncol=3, byrow= F)
  z[,1]=rep(filename, length(WL.p$Time))
  z[,2]=as.numeric(WL.p$Time)
  z[,3]=as.numeric(WL.p$Position)
  return(z)
  }

### 8.6 Hz ###
# loop through all the *.ddf files and get the desired summary data 
out8<-lapply(files8, function(filename) {
	x <- strain8(filename)})

### 6.1 Hz ###
# loop through all the *.ddf files and get the desired summary data 
out6<-lapply(files6, function(filename) {
	x <- strain6(filename)})

### 10.1 Hz ###
# loop through all the *.ddf files and get the desired summary data 
out10<-lapply(files10, function(filename) {
	x <- strain10(filename)})

# create a data frame containing the strain data from all *.ddf files of interest
### 8.6 Hz ###
p8<-NULL
for (i in 1:length(out8)){
	w8 <- data.frame(out8[[i]], stringsAsFactors=F)
	p8 <- rbind(p8, w8)
	}
	
# rename the columns to properly identify the metrics
colnames(p8) <- c("Rep_num", "Time", "Position")

# Check it
head(p8)
tail(p8)
str(p8)
p8$Time <- as.numeric(p8$Time)
p8$Position <- as.numeric(p8$Position)
p8$Rep_num <- as.numeric(as.factor(p8$Rep_num))

### 6.1 Hz ###
p6<-NULL
for (i in 1:length(out6)){
	w6 <- data.frame(out6[[i]], stringsAsFactors=F)
	p6 <- rbind(p6, w6)
	}
	
# rename the columns to properly identify the metrics
colnames(p6) <- c("Rep_num", "Time", "Position")

# Check it
head(p6)
tail(p6)
str(p6)
p6$Time <- as.numeric(p6$Time)
p6$Position <- as.numeric(p6$Position)
p6$Rep_num <- as.numeric(as.factor(p6$Rep_num))

### 10.1 Hz ###
p10<-NULL
for (i in 1:length(out10)){
	w10 <- data.frame(out10[[i]], stringsAsFactors=F)
	p10 <- rbind(p10, w10)
	}
	
# rename the columns to properly identify the metrics
colnames(p10) <- c("Rep_num", "Time", "Position")

# Check it
head(p10)
tail(p10)
str(p10)
p10$Time <- as.numeric(p10$Time)
p10$Position <- as.numeric(p10$Position)
p10$Rep_num <- as.numeric(as.factor(p10$Rep_num))

# write csvs of stacked in vitro strains to be used in FFT analysis in matlab
write.csv(p8, file = "all-8.6-Hz_in situ_strains.csv", row.names=F)

write.csv(p6, file = "all-6.1-Hz_in situ_strains.csv", row.names=F)

write.csv(p10, file = "all-10.1-Hz_in situ_strains.csv", row.names=F)
