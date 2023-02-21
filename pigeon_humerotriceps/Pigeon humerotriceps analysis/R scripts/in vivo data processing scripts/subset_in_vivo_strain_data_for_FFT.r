# Script to subset Robertson and Biewener in vivo strain data for FFT analysis
# get working directory
getwd()
# change working directory
## select subfolder within all in vivo strain data folder containing all in 
## vivo strain data that has been processed for analysis ('in vivo strain 
## data for analysis' folder)
setwd(choose.dir())

# read in csvs for Bird Www

www_L1<-read.csv("BirdWww1 in vivo_HT strain.csv")
www_T1<-read.csv("BirdWww2 in vivo_HT strain.csv")
www_T2<-read.csv("BirdWww3 in vivo_HT strain.csv")
www_T3<-read.csv("BirdWww4 in vivo_HT strain.csv")
www_T4<-read.csv("BirdWww5 in vivo_HT strain.csv")

# read in csvs for Bird Xxx

xxx_L1<-read.csv("BirdXxx1 in vivo_HT strain.csv")
xxx_T1<-read.csv("BirdXxx2 in vivo_HT strain.csv")
xxx_T2<-read.csv("BirdXxx3 in vivo_HT strain.csv")
xxx_T3<-read.csv("BirdXxx4 in vivo_HT strain.csv")
xxx_T4<-read.csv("BirdXxx5 in vivo_HT strain.csv")

# read in csvs for Bird Uuu

uuu_T1<-read.csv("BirdUuu1 in vivo_HT strain.csv")
uuu_T2<-read.csv("BirdUuu2 in vivo_HT strain.csv")
uuu_L1<-read.csv("BirdUuu3 in vivo_HT strain.csv")
uuu_L2<-read.csv("BirdUuu4 in vivo_HT strain.csv")
uuu_L3<-read.csv("BirdUuu5 in vivo_HT strain.csv")

# read in csvs for Bird Tie

tie_T1<-read.csv("BirdTie1 in vivo_HT strain.csv")
tie_L1<-read.csv("BirdTie2 in vivo_HT strain.csv")
tie_L2<-read.csv("BirdTie3 in vivo_HT strain.csv")
tie_L3<-read.csv("BirdTie4 in vivo_HT strain.csv")
tie_L4<-read.csv("BirdTie5 in vivo_HT strain.csv")

# read in csvs for Bird Yyy

yyy_L1<-read.csv("BirdYyy1 in vivo_HT strain.csv")
yyy_T1<-read.csv("BirdYyy2 in vivo_HT strain.csv")
yyy_T2<-read.csv("BirdYyy3 in vivo_HT strain.csv")
yyy_T3<-read.csv("BirdYyy4 in vivo_HT strain.csv")
yyy_T4<-read.csv("BirdYyy5 in vivo_HT strain.csv")

#### Bird Www ####
## trial 1
# subset trial to only have time and HT sono data
Birdwww1 <- subset(www_L1, select=c(Time, HT.sono))
# rename columns
colnames(Birdwww1) <- c("time","humtri_SONO")
# check it
head(Birdwww1)

## trial 2
# subset trial to only have time and HT sono data
Birdwww2 <- subset(www_T1, select=c(Time, HT.sono))
colnames(Birdwww2) <- c("time","humtri_SONO")
# check it
head(Birdwww2)

## trial 3
# subset trial to only have time and HT sono data
Birdwww3 <- subset(www_T2, select=c(Time, HT.sono))
colnames(Birdwww3) <- c("time","humtri_SONO")
# check it
head(Birdwww3)

## trial 4
# subset trial to only have time and HT sono data
Birdwww4 <- subset(www_T3, select=c(Time, HT.sono))
colnames(Birdwww4) <- c("time","humtri_SONO")
# check it
head(Birdwww4)

## trial 5
# subset trial to only have time and HT sono data
Birdwww5 <- subset(www_T4, select=c(Time, HT.sono))
# rename columns
colnames(Birdwww5) <- c("time","humtri_SONO")
#check it
head(Birdwww5)

#### Bird Xxx ####
## trial 1
# subset trial to only have time and HT sono data
Birdxxx1 <- subset(xxx_L1, select=c(Time, HT.sono))
colnames(Birdxxx1) <- c("time","humtri_SONO")
# check it
head(Birdxxx1)

## trial 2
# subset trial to only have time and HT sono data
Birdxxx2 <- subset(xxx_T1, select=c(Time, HT.sono))
colnames(Birdxxx2) <- c("time","humtri_SONO")
# check it
head(Birdxxx2)

## trial 3
# subset trial to only have time and HT sono data
Birdxxx3 <- subset(xxx_T2, select=c(Time, HT.sono))
colnames(Birdxxx3) <- c("time","humtri_SONO")
# check it
head(Birdxxx3)

## trial 4
# subset trial to only have time and HT sono data
Birdxxx4 <- subset(xxx_T3, select=c(Time, HT.sono))
colnames(Birdxxx4) <- c("time","humtri_SONO")
# check it
head(Birdxxx4)

## trial 5
# subset trial to only have time and HT sono data
Birdxxx5 <- subset(xxx_T4, select=c(Time, HT.sono))
colnames(Birdxxx5) <- c("time","humtri_SONO")
# check it
head(Birdxxx5)

#### Bird Uuu ####
## trial 1
# subset trial to only have time and HT sono data
Birduuu1 <- subset(uuu_T1, select=c(Time, HT.sono))
colnames(Birduuu1) <- c("time","humtri_SONO")
# check it
head(Birduuu1)

## trial 2
# subset trial to only have time and HT sono data
Birduuu2 <- subset(uuu_T2, select=c(Time, HT.sono))
colnames(Birduuu2) <- c("time","humtri_SONO")
# check it
head(Birduuu2)

## trial 3
# subset trial to only have time and HT sono data
Birduuu3 <- subset(uuu_L1, select=c(Time, HT.sono))
colnames(Birduuu3) <- c("time","humtri_SONO")
# check it
head(Birduuu3)

## trial 4
# subset trial to only have time and HT sono data
Birduuu4 <- subset(uuu_L2, select=c(Time, HT.sono))
colnames(Birduuu4) <- c("time","humtri_SONO")
# check it
head(Birduuu4)

## trial 5
# subset trial to only have time and HT sono data
Birduuu5 <- subset(uuu_L3, select=c(Time, HT.sono))
colnames(Birduuu5) <- c("time","humtri_SONO")
# check it
head(Birduuu5)

#### Bird Tie ####
## trial 1
# subset trial to only have time and HT sono data
Birdtie1 <- subset(tie_T1, select=c(Time, HT.sono))
colnames(Birdtie1) <- c("time","humtri_SONO")
# check it
head(Birdtie1)

## trial 2
# subset trial to only have time and HT sono data
Birdtie2 <- subset(tie_L1, select=c(Time, HT.sono))
colnames(Birdtie2) <- c("time","humtri_SONO")
# check it
head(Birdtie2)

## trial 3
# subset trial to only have time and HT sono data
Birdtie3 <- subset(tie_L2, select=c(Time, HT.sono))
colnames(Birdtie3) <- c("time","humtri_SONO")
# check it
head(Birdtie3)

## trial 4
# subset trial to only have time and HT sono data
Birdtie4 <- subset(tie_L3, select=c(Time, HT.sono))
colnames(Birdtie4) <- c("time","humtri_SONO")
# check it
head(Birdtie4)

## trial 5
# subset trial to only have time and HT sono data
Birdtie5 <- subset(tie_L4, select=c(Time, HT.sono))
colnames(Birdtie5) <- c("time","humtri_SONO")
# check it
head(Birdtie5)

#### Bird Yyy ####
## trial 1
# subset trial to only have time and HT sono data
Birdyyy1 <- subset(yyy_L1, select=c(Time, HT.sono))
colnames(Birdyyy1) <- c("time","humtri_SONO")
# check it
head(Birdyyy1)

## trial 2
# subset trial to only have time and HT sono data
Birdyyy2 <- subset(yyy_T1, select=c(Time, HT.sono))
colnames(Birdyyy2) <- c("time","humtri_SONO")
# check it
head(Birdyyy2)

## trial 3
# subset trial to only have time and HT sono data
Birdyyy3 <- subset(yyy_T2, select=c(Time, HT.sono))
colnames(Birdyyy3) <- c("time","humtri_SONO")
# check it
head(Birdyyy3)

## trial 4
# subset trial to only have time and HT sono data
Birdyyy4 <- subset(yyy_T3, select=c(Time, HT.sono))
colnames(Birdyyy4) <- c("time","humtri_SONO")
# check it
head(Birdyyy4)

## trial 5
# subset trial to only have time and HT sono data
Birdyyy5 <- subset(yyy_T4, select=c(Time, HT.sono))
colnames(Birdyyy5) <- c("time","humtri_SONO")
# check it
head(Birdyyy5)

# change working directory
setwd(choose.dir()) # select 'Manuscript experiments' folder 

# write csv for strain data from Bird Www
write.csv(Birdwww1, "BirdWww1 in vivo_HT strain.csv", row.names=FALSE)
write.csv(Birdwww2, "BirdWww2 in vivo_HT strain.csv", row.names=FALSE)
write.csv(Birdwww3, "BirdWww3 in vivo_HT strain.csv", row.names=FALSE)
write.csv(Birdwww4, "BirdWww4 in vivo_HT strain.csv", row.names=FALSE)
write.csv(Birdwww5, "BirdWww5 in vivo_HT strain.csv", row.names=FALSE)

# write csv for strain data from Bird Xxx
write.csv(Birdxxx1, "BirdXxx1 in vivo_HT strain.csv", row.names=FALSE)
write.csv(Birdxxx2, "BirdXxx2 in vivo_HT strain.csv", row.names=FALSE)
write.csv(Birdxxx3, "BirdXxx3 in vivo_HT strain.csv", row.names=FALSE)
write.csv(Birdxxx4, "BirdXxx4 in vivo_HT strain.csv", row.names=FALSE)
write.csv(Birdxxx5, "BirdXxx5 in vivo_HT strain.csv", row.names=FALSE)

# write csv for strain data from Bird Uuu
write.csv(Birduuu1, "BirdUuu1 in vivo_HT strain.csv", row.names=FALSE)
write.csv(Birduuu2, "BirdUuu2 in vivo_HT strain.csv", row.names=FALSE)
write.csv(Birduuu3, "BirdUuu3 in vivo_HT strain.csv", row.names=FALSE)
write.csv(Birduuu4, "BirdUuu4 in vivo_HT strain.csv", row.names=FALSE)
write.csv(Birduuu5, "BirdUuu5 in vivo_HT strain.csv", row.names=FALSE)

# write csv for strain data from Bird Tie
write.csv(Birdtie1, "BirdTie1 in vivo_HT strain.csv", row.names=FALSE)
write.csv(Birdtie2, "BirdTie2 in vivo_HT strain.csv", row.names=FALSE)
write.csv(Birdtie3, "BirdTie3 in vivo_HT strain.csv", row.names=FALSE)
write.csv(Birdtie4, "BirdTie4 in vivo_HT strain.csv", row.names=FALSE)
write.csv(Birdtie5, "BirdTie5 in vivo_HT strain.csv", row.names=FALSE)

# write csv for strain data from Bird Yyy
write.csv(Birdyyy1, "BirdYyy1 in vivo_HT strain.csv", row.names=FALSE)
write.csv(Birdyyy2, "BirdYyy2 in vivo_HT strain.csv", row.names=FALSE)
write.csv(Birdyyy3, "BirdYyy3 in vivo_HT strain.csv", row.names=FALSE)
write.csv(Birdyyy4, "BirdYyy4 in vivo_HT strain.csv", row.names=FALSE)
write.csv(Birdyyy5, "BirdYyy5 in vivo_HT strain.csv", row.names=FALSE)

