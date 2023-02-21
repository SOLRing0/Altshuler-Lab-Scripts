# Script to calculate lengthening to shortening ratios for all 
# Robertson and Biewener in vivo strain data
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

# plot sono data by time for five flights from each bird to check it
# Bird Www
plot(www_L1$HT.sono~www_L1$Time, type="l")
plot(www_T1$HT.sono~www_T1$Time, type="l")
plot(www_T2$HT.sono~www_T2$Time, type="l")
plot(www_T3$HT.sono~www_T3$Time, type="l")
plot(www_T4$HT.sono~www_T4$Time, type="l")

# Bird Xxx
plot(xxx_L1$HT.sono~xxx_L1$Time, type="l") 
plot(xxx_T1$HT.sono~xxx_T1$Time, type="l") 
plot(xxx_T2$HT.sono~xxx_T2$Time, type="l") 
plot(xxx_T3$HT.sono~xxx_T3$Time, type="l") 
plot(xxx_T4$HT.sono~xxx_T4$Time, type="l") 

# Bird Uuu
plot(uuu_T1$HT.sono~uuu_T1$Time, type="l")
plot(uuu_T2$HT.sono~uuu_T2$Time, type="l")
plot(uuu_L1$HT.sono~uuu_L1$Time, type="l")
plot(uuu_L2$HT.sono~uuu_L2$Time, type="l")
plot(uuu_L3$HT.sono~uuu_L3$Time, type="l")

# Bird Tie
plot(tie_T1$HT.sono~tie_T1$Time, type="l")
plot(tie_L1$HT.sono~tie_L1$Time, type="l")
plot(tie_L2$HT.sono~tie_L2$Time, type="l")
plot(tie_L3$HT.sono~tie_L3$Time, type="l")
plot(tie_L4$HT.sono~tie_L4$Time, type="l")

# Bird Yyy
plot(yyy_L1$HT.sono~yyy_L1$Time, type="l")
plot(yyy_T1$HT.sono~yyy_T1$Time, type="l")
plot(yyy_T2$HT.sono~yyy_T2$Time, type="l")
plot(yyy_T3$HT.sono~yyy_T3$Time, type="l")
plot(yyy_T4$HT.sono~yyy_T4$Time, type="l")

# load pracma package for findpeaks function
library(pracma)

# Define find peaks function 
find_peaks <- function (x, m = 3){
     shape <- diff(sign(diff(x, na.pad = FALSE)))
     pks <- sapply(which(shape < 0), FUN = function(i){
        z <- i - m + 1
        z <- ifelse(z > 0, z, 1)
        w <- i + m + 1
        w <- ifelse(w < length(x), w, length(x))
        if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
     })
      pks <- unlist(pks)
      pks
 }

### Find peaks for Bird Www's five humerotriceps SONO traces ###

### Flight 1 ###
www1=www_L1$HT.sono 
plot(www1, type="l")

## Bird Www - flight 1 HT peaks
pw1=find_peaks(www1,m=200)
points(pw1, www1[pw1], col="blue", pch=19) # plots peaks in blue onto strain profile plot
pw1
length(pw1) # counts peaks

## Bird Www - flight 1 HT troughs
tw1=find_peaks(-www1,m=200)
points(tw1, www1[tw1], col="blue") # plots open blue points for trough locations
tw1
length(tw1) # counts troughs

# HT length change cycle duration for Bird Www flight 1
HTpeaksw1 <- www_L1$Time[pw1] 
HTtroughsw1 <- www_L1$Time[tw1] 
HTcycDw1 <- diff(HTpeaksw1)

#  Last marked point is a trough, for simplicity, 
## calculate Bird Www shortening ratio first 
shortw1 <- HTtroughsw1[-1] - HTpeaksw1 # first trough removed--no peak prior to it
# remove last shortw1 -- no cycle duration 
shortw1 <- shortw1[-(length(shortw1))]
shortw1_pc <- (shortw1/HTcycDw1)*100
avg_Sw1_pc <- mean(shortw1_pc) 

# Calculate Bird Www lengthening ratio using cycle duration
longw1 <- HTcycDw1-shortw1
longw1_pc <- (longw1/HTcycDw1)*100
avg_Lw1_pc <- mean(longw1_pc) 

# Check the math
shortw1_pc + longw1_pc # should all sum to 100% -- good
avg_Sw1_pc + avg_Lw1_pc # should also sum to ~100% -- good

### Flight 2 ###
www2=www_T1$HT.sono
plot(www2, type="l")

## Bird Www - flight 2 HT peaks
pw2=find_peaks(www2, m = 200)
points(pw2, www2[pw2], col="blue", pch=19) # plots peaks in blue onto strain profile plot
pw2
length(pw2) # counts peaks

# too many peaks -- remove duplicates
# since values are both time and position duplicates, take first instance
# of its occurence if sono values are the same
pw2<-pw2[!duplicated(www_T1$HT.sono[pw2])]
# check it
plot(www2, type="l")
points(pw2, www2[pw2], col="blue", pch=19) # plots peaks in blue onto strain profile plot
length(pw2)

## Bird Www - flight 2 HT troughs
tw2=find_peaks(-www2, m = 200)
points(tw2, www2[tw2], col="blue") # plots open blue points for trough locations
tw2
length(tw2) # counts troughs

# too many troughs -- remove duplicates
# since values are both time and position duplicates, take first instance
# of its occurence if sono values are the same
tw2<-tw2[!duplicated(www_T1$HT.sono[tw2])]
# check it
plot(www2, type="l")
points(pw2, www2[pw2], col="blue", pch=19) # plots peaks in blue onto strain profile plot
points(tw2, www2[tw2], col="blue") # plots open blue points for trough locations
length(tw2)

# HT length change cycle duration for Bird Www flight 2
HTpeaksw2 <- www_T1$Time[pw2] 
HTtroughsw2 <- www_T1$Time[tw2] 
HTcycDw2 <- diff(HTpeaksw2)

#  Last marked point is a trough, for simplicity, 
## calculate Bird Www shortening ratio first 
shortw2 <- HTtroughsw2[-1] - HTpeaksw2 # first trough removed--no peak prior to it
# remove last shortw2 -- no cycle duration 
shortw2 <- shortw2[-(length(shortw2))]
shortw2_pc <- (shortw2/HTcycDw2)*100
avg_Sw2_pc <- mean(shortw2_pc) 

# Calculate Bird Www lengthening ratio using cycle duration
longw2 <- HTcycDw2-shortw2
longw2_pc <- (longw2/HTcycDw2)*100
avg_Lw2_pc <- mean(longw2_pc) 

# Check the math
shortw2_pc + longw2_pc # should all sum to 100% -- good
avg_Sw2_pc + avg_Lw2_pc # should also sum to ~100% -- good

### Flight 3 ###
www3=www_T2$HT.sono
plot(www3, type="l")

## Bird Www - flight 3 HT peaks
pw3=find_peaks(www3, m = 200)
points(pw3, www3[pw3], col="blue", pch=19) # plots peaks in blue onto strain profile plot
pw3
length(pw3) # counts peaks

# too many peaks -- remove duplicates
# since values are both time and position duplicates, take first instance
# of its occurence if sono values are the same
pw3<-pw3[!duplicated(www_T2$HT.sono[pw3])]
# check it
plot(www3, type="l")
points(pw3, www3[pw3], col="blue", pch=19) # plots peaks in blue onto strain profile plot
length(pw3)

## Bird Www - flight 3 HT troughs
tw3=find_peaks(-www3, m = 200)
points(tw3, www3[tw3], col="blue") # plots open blue points for trough locations
tw3
length(tw3) # counts troughs

# too many troughs -- remove duplicates
# since values are both time and position duplicates, take first instance
# of its occurence if sono values are the same
tw3<-tw3[!duplicated(www_T2$HT.sono[tw3])]
# check it
plot(www3, type="l")
points(pw3, www3[pw3], col="blue", pch=19) # plots peaks in blue onto strain profile plot
points(tw3, www3[tw3], col="blue") # plots open blue points for trough locations
length(tw3)

# HT length change cycle duration for Bird Www flight 3
HTpeaksw3 <- www_T2$Time[pw3] 
HTtroughsw3 <- www_T2$Time[tw3] 
HTcycDw3 <- diff(HTpeaksw3)

#  Last marked point is a peak, for simplicity, 
## calculate Bird Www shortening ratio first 
shortw3 <- HTtroughsw3[-1] - HTpeaksw3[-15] # first trough removed--no peak prior to it
shortw3_pc <- (shortw3/HTcycDw3)*100
avg_Sw3_pc <- mean(shortw3_pc) 

# Calculate Bird Www lengthening ratio using cycle duration
longw3 <- HTcycDw3-shortw3
longw3_pc <- (longw3/HTcycDw3)*100
avg_Lw3_pc <- mean(longw3_pc) 

# Check the math
shortw3_pc + longw3_pc # should all sum to 100% -- good
avg_Sw3_pc + avg_Lw3_pc # should also sum to ~100% -- good

### Flight 4 ###
www4=www_T3$HT.sono
plot(www4, type="l")

## Bird Www - flight 4 HT peaks
pw4=find_peaks(www4, m = 200)
points(pw4, www4[pw4], col="blue", pch=19) # plots peaks in blue onto strain profile plot
pw4
length(pw4) # counts peaks

# too many peaks -- remove duplicates
# since values are both time and position duplicates, take first instance
# of its occurence if sono values are the same
pw4<-pw4[!duplicated(www_T3$HT.sono[pw4])]
# check it
plot(www4, type="l")
points(pw4, www4[pw4], col="blue", pch=19) # plots peaks in blue onto strain profile plot
length(pw4)

## Bird Www - flight 4 HT troughs
tw4=find_peaks(-www4, m = 200)
points(tw4, www4[tw4], col="blue") # plots open blue points for trough locations
tw4
length(tw4) # counts troughs

# too many troughs -- remove duplicates
# since values are both time and position duplicates, take first instance
# of its occurence if sono values are the same
tw4<-tw4[!duplicated(www_T3$HT.sono[tw4])]
# missing a midpoint trough -- add it back manually
tw4<-c(tw4[1:7],4352,tw4[8:16])
# check it
plot(www4, type="l")
points(pw4, www4[pw4], col="blue", pch=19) # plots peaks in blue onto strain profile plot
points(tw4, www4[tw4], col="blue") # plots open blue points for trough locations
length(tw4)

# HT length change cycle duration for Bird Www flight 4
HTpeaksw4 <- www_T3$Time[pw4] 
HTtroughsw4 <- www_T3$Time[tw4] 
HTcycDw4 <- diff(HTpeaksw4)

#  Last marked point is a trough, for simplicity, 
## calculate Bird Www shortening ratio first 
shortw4 <- HTtroughsw4[-1] - HTpeaksw4 # first trough removed--no peak prior to it
# remove last shortw4 -- no cycle duration 
shortw4 <- shortw4[-(length(shortw4))]
shortw4_pc <- (shortw4/HTcycDw4)*100
avg_Sw4_pc <- mean(shortw4_pc) 

# Calculate Bird Www lengthening ratio using cycle duration
longw4 <- HTcycDw4-shortw4
longw4_pc <- (longw4/HTcycDw4)*100
avg_Lw4_pc <- mean(longw4_pc) 

# Check the math
shortw4_pc + longw4_pc # should all sum to 100% -- good
avg_Sw4_pc + avg_Lw4_pc # should also sum to ~100% -- good

### Flight 5 ###
www5=www_T4$HT.sono
plot(www5, type="l")

## Bird Www - flight 5 HT peaks
pw5=find_peaks(www5, m = 200)
points(pw5, www5[pw5], col="blue", pch=19) # plots peaks in blue onto strain profile plot
pw5
length(pw5) # counts peaks

# too many peaks -- remove duplicates
# since values are both time and position duplicates, take first instance
# of its occurence if sono values are the same
pw5<-pw5[!duplicated(www_T4$HT.sono[pw5])]
# check it
plot(www5, type="l")
points(pw5, www5[pw5], col="blue", pch=19) # plots peaks in blue onto strain profile plot
length(pw5)

## Bird Www - flight 5 HT troughs
tw5=find_peaks(-www5, m = 200)
points(tw5, www5[tw5], col="blue") # plots open blue points for trough locations
tw5
length(tw5) # counts troughs

# too many troughs -- remove duplicates
# since values are both time and position duplicates, take first instance
# of its occurence if sono values are the same
tw5<-tw5[!duplicated(www_T4$HT.sono[tw5])]
# missing two troughs, add them back manually
tw5 <- c(tw5[1:5],2969,tw5[6:13],8641,tw5[14])
# check it
plot(www5, type="l")
points(pw5, www5[pw5], col="blue", pch=19) # plots peaks in blue onto strain profile plot
points(tw5, www5[tw5], col="blue") # plots open blue points for trough locations
length(tw5)

# HT length change cycle duration for Bird Www flight 5
HTpeaksw5 <- www_T4$Time[pw5] 
HTtroughsw5 <- www_T4$Time[tw5] 
HTcycDw5 <- diff(HTpeaksw5)

#  Last marked point is a peak, for simplicity, 
## calculate Bird Www shortening ratio first and remove last peak
shortw5 <- HTtroughsw5[-1] - HTpeaksw5[-(length(HTpeaksw5))] # first trough removed--no peak prior to it
shortw5_pc <- (shortw5/HTcycDw5)*100
avg_Sw5_pc <- mean(shortw5_pc) 

# Calculate Bird Www lengthening ratio using cycle duration
longw5 <- HTcycDw5-shortw5
longw5_pc <- (longw5/HTcycDw5)*100
avg_Lw5_pc <- mean(longw5_pc) 

# Check the math
shortw5_pc + longw5_pc # should all sum to 100% -- good
avg_Sw5_pc + avg_Lw5_pc # should also sum to ~100% -- good

### Find peaks for Bird Xxx's five humerotriceps SONO traces ###

### Flight 1 ###
xxx1=xxx_L1$HT.sono
plot(xxx1, type="l")

## Bird Xxx - flight 1 HT peaks
# use pracma findpeaks function for this trace as there are too many duplicates to do manually
# also already tested when comparing filtered to unfiltered data
px1=sort(findpeaks(xxx1,minpeakdistance=500)[,2])
points(px1, xxx1[px1], col="blue", pch=19) # plots peaks in blue onto strain profile plot
px1
length(px1) # counts peaks

## Bird Xxx - flight 1 HT troughs
tx1=sort(findpeaks(-xxx1,minpeakdistance=500)[,2])
points(tx1, xxx1[tx1], col="blue") # plots open blue points for trough locations
tx1
length(tx1) # counts troughs

# HT length change cycle duration for Bird Xxx flight 1
HTpeaksx1 <- xxx_L1$Time[px1] 
HTtroughsx1 <- xxx_L1$Time[tx1] 
HTcycDx1 <- diff(HTpeaksx1)

#  Last marked point is a trough, for simplicity, 
## calculate Bird Xxx shortening ratio first 
shortx1 <- HTtroughsx1[-1] - HTpeaksx1 # first trough removed--no peak prior to it
# remove last shortx1 -- no cycle duration 
shortx1 <- shortx1[-(length(shortx1))]
shortx1_pc <- (shortx1/HTcycDx1)*100
avg_Sx1_pc <- mean(shortx1_pc) 

# Calculate Bird Xxx lengthening ratio using cycle duration
longx1 <- HTcycDx1-shortx1
longx1_pc <- (longx1/HTcycDx1)*100
avg_Lx1_pc <- mean(longx1_pc) 

# Check the math
shortx1_pc + longx1_pc # should all sum to 100% -- good
avg_Sx1_pc + avg_Lx1_pc # should also sum to ~100% -- good

### Flight 2 ###
xxx2=xxx_T1$HT.sono
plot(xxx2, type="l")

## Bird Xxx - flight 2 HT peaks
px2=sort(findpeaks(xxx2,minpeakdistance=500)[,2])
points(px2, xxx2[px2], col="blue", pch=19) # plots peaks in blue onto strain profile plot
px2
length(px2) # counts peaks

## Bird Xxx - flight 2 HT troughs
tx2=sort(findpeaks(-xxx2,minpeakdistance=500)[,2])
points(tx2, xxx2[tx2], col="blue") # plots open blue points for trough locations
tx2
length(tx2) # counts troughs

# HT length change cycle duration for Bird Xxx flight 2
HTpeaksx2 <- xxx_T1$Time[px2] 
HTtroughsx2 <- xxx_T1$Time[tx2] 
HTcycDx2 <- diff(HTpeaksx2)

#  Last marked point is a trough, for simplicity, 
## calculate Bird Xxx shortening ratio first 
shortx2 <- HTtroughsx2[-1] - HTpeaksx2 # first trough removed--no peak prior to it
# remove last shortx2 -- no cycle duration 
shortx2 <- shortx2[-(length(shortx2))]
shortx2_pc <- (shortx2/HTcycDx2)*100
avg_Sx2_pc <- mean(shortx2_pc) 

# Calculate Bird Xxx lengthening ratio using cycle duration
longx2 <- HTcycDx2-shortx2
longx2_pc <- (longx2/HTcycDx2)*100
avg_Lx2_pc <- mean(longx2_pc) 

# Check the math
shortx2_pc + longx2_pc # should all sum to 100% -- good
avg_Sx2_pc + avg_Lx2_pc # should also sum to ~100% -- good

### Flight 3 ###
xxx3=xxx_T2$HT.sono
plot(xxx3, type="l")

## Bird Xxx - flight 3 HT peaks
px3=sort(findpeaks(xxx3,minpeakdistance=500)[,2])
points(px3, xxx3[px3], col="blue", pch=19) # plots peaks in blue onto strain profile plot
px3
length(px3) # counts peaks

## Bird Xxx - flight 3 HT troughs
tx3=sort(findpeaks(-xxx3,minpeakdistance=500)[,2])
points(tx3, xxx3[tx3], col="blue") # plots open blue points for trough locations
# adjust i=17--only one that looks off therefore do manually rather than change method
tx3[17]<-9900 # slightly too far ahead therefore move back 44 indices
tx3
length(tx3) # counts troughs

# HT length change cycle duration for Bird Xxx flight 3
HTpeaksx3 <- xxx_T2$Time[px3] 
HTtroughsx3 <- xxx_T2$Time[tx3] 
HTcycDx3 <- diff(HTpeaksx3)

#  Last marked point is a trough, for simplicity, 
## calculate Bird Xxx shortening ratio first 
shortx3 <- HTtroughsx3[-1] - HTpeaksx3 # first trough removed--no peak prior to it
# remove last shortx3 -- no cycle duration 
shortx3 <- shortx3[-(length(shortx3))]
shortx3_pc <- (shortx3/HTcycDx3)*100
avg_Sx3_pc <- mean(shortx3_pc) 

# Calculate Bird Xxx lengthening ratio using cycle duration
longx3 <- HTcycDx3-shortx3
longx3_pc <- (longx3/HTcycDx3)*100
avg_Lx3_pc <- mean(longx3_pc) 

# Check the math
shortx3_pc + longx3_pc # should all sum to 100% -- good
avg_Sx3_pc + avg_Lx3_pc # should also sum to ~100% -- good

### Flight 4 ###
xxx4=xxx_T3$HT.sono
plot(xxx4, type="l")

## Bird Xxx - flight 4 HT peaks
px4=sort(findpeaks(xxx4,minpeakdistance=500)[,2])
points(px4, xxx4[px4], col="blue", pch=19) # plots peaks in blue onto strain profile plot
px4
length(px4) # counts peaks

## Bird Xxx - flight 4 HT troughs
tx4=sort(findpeaks(-xxx4,minpeakdistance=500)[,2])
points(tx4, xxx4[tx4], col="blue") # plots open blue points for trough locations
tx4
length(tx4) # counts troughs

# HT length change cycle duration for Bird Xxx flight 4
HTpeaksx4 <- xxx_T3$Time[px4] 
HTtroughsx4 <- xxx_T3$Time[tx4] 
HTcycDx4 <- diff(HTpeaksx4)

#  Last marked point is a trough, for simplicity, 
## calculate Bird Xxx shortening ratio first 
shortx4 <- HTtroughsx4[-1] - HTpeaksx4 # first trough removed--no peak prior to it
# remove last shortx4 -- no cycle duration 
shortx4 <- shortx4[-(length(shortx4))]
shortx4_pc <- (shortx4/HTcycDx4)*100
avg_Sx4_pc <- mean(shortx4_pc) 

# Calculate Bird Xxx lengthening ratio using cycle duration
longx4 <- HTcycDx4-shortx4
longx4_pc <- (longx4/HTcycDx4)*100
avg_Lx4_pc <- mean(longx4_pc) 

# Check the math
shortx4_pc + longx4_pc # should all sum to 100% -- good
avg_Sx4_pc + avg_Lx4_pc # should also sum to ~100% -- good

### Flight 5 ###
xxx5=xxx_T4$HT.sono
plot(xxx5, type="l")

## Bird Xxx - flight 5 HT peaks
px5=sort(findpeaks(xxx5,minpeakdistance=500)[,2])
points(px5, xxx5[px5], col="blue", pch=19) # plots peaks in blue onto strain profile plot
px5
length(px5) # counts peaks

## Bird Xxx - flight 5 HT troughs
tx5=sort(findpeaks(-xxx5,minpeakdistance=500)[,2])
points(tx5, xxx5[tx5], col="blue") # plots open blue points for trough locations
tx5
length(tx5) # counts troughs

# HT length change cycle duration for Bird Xxx flight 5
HTpeaksx5 <- xxx_T4$Time[px5] 
HTtroughsx5 <- xxx_T4$Time[tx5] 
HTcycDx5 <- diff(HTpeaksx5)

#  Last marked point is a trough, for simplicity, 
## calculate Bird Xxx shortening ratio first 
shortx5 <- HTtroughsx5[-1] - HTpeaksx5 # first trough removed--no peak prior to it
# remove last shortx5 -- no cycle duration 
shortx5 <- shortx5[-(length(shortx5))]
shortx5_pc <- (shortx5/HTcycDx5)*100
avg_Sx5_pc <- mean(shortx5_pc) 

# Calculate Bird Xxx lengthening ratio using cycle duration
longx5 <- HTcycDx5-shortx5
longx5_pc <- (longx5/HTcycDx5)*100
avg_Lx5_pc <- mean(longx5_pc) 

# Check the math
shortx5_pc + longx5_pc # should all sum to 100% -- good
avg_Sx5_pc + avg_Lx5_pc # should also sum to ~100% -- good

### Find peaks for Bird Uuu's five humerotriceps SONO traces ###

### Flight 1 ###
uuu1=uuu_T1$HT.sono
plot(uuu1, type="l")

## Bird Uuu - flight 1 HT peaks
pu1=find_peaks(uuu1, m = 200)
points(pu1, uuu1[pu1], col="blue", pch=19) # plots peaks in blue onto strain profile plot
pu1
length(pu1) # counts peaks

## Bird Uuu - flight 1 HT troughs
tu1=find_peaks(-uuu1, m = 200)
points(tu1, uuu1[tu1], col="blue") # plots open blue points for trough locations
tu1
length(tu1) # counts troughs

# HT length change cycle duration for Bird Uuu flight 1
HTpeaksu1 <- uuu_T1$Time[pu1] 
HTtroughsu1 <- uuu_T1$Time[tu1] 
HTcycDu1 <- diff(HTpeaksu1[-1]) # exclude first peak--not true peak

#  Last marked point is a trough, for simplicity, 
## calculate Bird Uuu shortening ratio first 
shortu1 <- HTtroughsu1[-1] - HTpeaksu1[-1] # first trough removed--no true peak prior to it and remove first 'peak'-- not true peak
# remove last shortu1 -- no cycle duration 
shortu1 <- shortu1[-(length(shortu1))]
shortu1_pc <- (shortu1/HTcycDu1)*100
avg_Su1_pc <- mean(shortu1_pc) 

# Calculate Bird Uuu lengthening ratio using cycle duration
longu1 <- HTcycDu1-shortu1
longu1_pc <- (longu1/HTcycDu1)*100
avg_Lu1_pc <- mean(longu1_pc) 

# Check the math
shortu1_pc + longu1_pc # should all sum to 100% -- good
avg_Su1_pc + avg_Lu1_pc # should also sum to ~100% -- good

### Flight 2 ###
uuu2=uuu_T2$HT.sono
plot(uuu2, type="l")

## Bird Uuu - flight 2 HT peaks
pu2=find_peaks(uuu2, m = 200)
points(pu2, uuu2[pu2], col="blue", pch=19) # plots peaks in blue onto strain profile plot
pu2
length(pu2) # counts peaks

## Bird Uuu - flight 2 HT troughs
tu2=find_peaks(-uuu2, m = 200)
points(tu2, uuu2[tu2], col="blue") # plots open blue points for trough locations
tu2
length(tu2) # counts troughs

# HT length change cycle duration for Bird Uuu flight 2
HTpeaksu2 <- uuu_T2$Time[pu2] 
HTtroughsu2 <- uuu_T2$Time[tu2] 
HTcycDu2 <- diff(HTpeaksu2) 

#  Last marked point is a trough, for simplicity, 
## calculate Bird Uuu shortening ratio first 
shortu2 <- HTtroughsu2[-1] - HTpeaksu2 # first trough removed--no peak prior to it 
# remove last shortu2 -- no cycle duration 
shortu2 <- shortu2[-(length(shortu2))]
shortu2_pc <- (shortu2/HTcycDu2)*100
avg_Su2_pc <- mean(shortu2_pc) 

# Calculate Bird Uuu lengthening ratio using cycle duration
longu2 <- HTcycDu2-shortu2
longu2_pc <- (longu2/HTcycDu2)*100
avg_Lu2_pc <- mean(longu2_pc) 

# Check the math
shortu2_pc + longu2_pc # should all sum to 100% -- good
avg_Su2_pc + avg_Lu2_pc # should also sum to ~100% -- good

### Flight 3 ###
uuu3=uuu_L1$HT.sono
plot(uuu3, type="l")

## Bird Uuu - flight 3 HT peaks
pu3=find_peaks(uuu3, m = 200)
points(pu3, uuu3[pu3], col="blue", pch=19) # plots peaks in blue onto strain profile plot
pu3
length(pu3) # counts peaks

## Bird Uuu - flight 3 HT troughs
tu3=find_peaks(-uuu3, m = 200)
points(tu3, uuu3[tu3], col="blue") # plots open blue points for trough locations
tu3
length(tu3) # counts troughs

# HT length change cycle duration for Bird Uuu flight 3
HTpeaksu3 <- uuu_L1$Time[pu3] 
HTtroughsu3 <- uuu_L1$Time[tu3] 
HTcycDu3 <- diff(HTpeaksu3) 

#  Last marked point is a peak, for simplicity, 
## calculate Bird Uuu shortening ratio first and remove last peak
shortu3 <- HTtroughsu3[-1] - HTpeaksu3[-(length(HTpeaksu3))] # first trough removed--no peak prior to it 
shortu3_pc <- (shortu3/HTcycDu3)*100
avg_Su3_pc <- mean(shortu3_pc) 

# Calculate Bird Uuu lengthening ratio using cycle duration
longu3 <- HTcycDu3-shortu3
longu3_pc <- (longu3/HTcycDu3)*100
avg_Lu3_pc <- mean(longu3_pc) 

# Check the math
shortu3_pc + longu3_pc # should all sum to 100% -- good
avg_Su3_pc + avg_Lu3_pc # should also sum to ~100% -- good

### Flight 4 ###
uuu4=uuu_L2$HT.sono
plot(uuu4, type="l")

## Bird Uuu - flight 4 HT peaks
pu4=find_peaks(uuu4, m = 200)
points(pu4, uuu4[pu4], col="blue", pch=19) # plots peaks in blue onto strain profile plot
pu4
length(pu4) # counts peaks

## Bird Uuu - flight 4 HT troughs
tu4=find_peaks(-uuu4, m = 200)
points(tu4, uuu4[tu4], col="blue") # plots open blue points for trough locations
tu4
length(tu4) # counts troughs

# HT length change cycle duration for Bird Uuu flight 4
HTpeaksu4 <- uuu_L2$Time[pu4] 
HTtroughsu4 <- uuu_L2$Time[tu4] 
HTcycDu4 <- diff(HTpeaksu4[-1]) # exclude first peak--not true peak

#  Last marked point is a peak, for simplicity, 
## calculate Bird Uuu shortening ratio first and remove last peak
shortu4 <- HTtroughsu4[-1] - HTpeaksu4[-c(1,(length(HTpeaksu4)))] # first trough removed--no true peak prior to it and remove first 'peak'-- not true peak
shortu4_pc <- (shortu4/HTcycDu4)*100
avg_Su4_pc <- mean(shortu4_pc) 

# Calculate Bird Uuu lengthening ratio using cycle duration
longu4 <- HTcycDu4-shortu4
longu4_pc <- (longu4/HTcycDu4)*100
avg_Lu4_pc <- mean(longu4_pc) 

# Check the math
shortu4_pc + longu4_pc # should all sum to 100% -- good
avg_Su4_pc + avg_Lu4_pc # should also sum to ~100% -- good

### Flight 5 ###
uuu5=uuu_L3$HT.sono
plot(uuu5, type="l")

## Bird Uuu - flight 5 HT peaks
pu5=find_peaks(uuu5, m = 200)
points(pu5, uuu5[pu5], col="blue", pch=19) # plots peaks in blue onto strain profile plot
pu5
length(pu5) # counts peaks

## Bird Uuu - flight 5 HT troughs
tu5=find_peaks(-uuu5, m = 200)
points(tu5, uuu5[tu5], col="blue") # plots open blue points for trough locations
tu5
length(tu5) # counts troughs

# HT length change cycle duration for Bird Uuu flight 5
HTpeaksu5 <- uuu_L3$Time[pu5] 
HTtroughsu5 <- uuu_L3$Time[tu5] 
HTcycDu5 <- diff(HTpeaksu5) 

#  Last marked point is a trough, for simplicity, 
## calculate Bird Uuu shortening ratio first
shortu5 <- HTtroughsu5[-1] - HTpeaksu5 # first trough removed--no peak prior to it 
# remove last shortu5 -- no cycle duration 
shortu5 <- shortu5[-(length(shortu5))]
shortu5_pc <- (shortu5/HTcycDu5)*100
avg_Su5_pc <- mean(shortu5_pc) 

# Calculate Bird Uuu lengthening ratio using cycle duration
longu5 <- HTcycDu5-shortu5
longu5_pc <- (longu5/HTcycDu5)*100
avg_Lu5_pc <- mean(longu5_pc) 

# Check the math
shortu5_pc + longu5_pc # should all sum to 100% -- good
avg_Su5_pc + avg_Lu5_pc # should also sum to ~100% -- good

### Find peaks for Bird Tie's five humerotriceps SONO traces ###

### Flight 1 ###
tie1=tie_T1$HT.sono
plot(tie1, type="l")

## Bird Tie - flight 1 HT peaks
pt1=find_peaks(tie1, m = 200)
points(pt1, tie1[pt1], col="blue", pch=19) # plots peaks in blue onto strain profile plot
pt1
length(pt1) # counts peaks

## Bird Tie - flight 1 HT troughs
tt1=find_peaks(-tie1, m = 200)
points(tt1, tie1[tt1], col="blue") # plots open blue points for trough locations
tt1
length(tt1) # counts troughs

# HT length change cycle duration for Bird Tie flight 1
HTpeakst1 <- tie_T1$Time[pt1] 
HTtroughst1 <- tie_T1$Time[tt1] 
HTcycDt1 <- diff(HTpeakst1) 

#  Last marked point is a trough, for simplicity, 
## calculate Bird Tie shortening ratio first
shortt1 <- HTtroughst1 - HTpeakst1 
# remove last shortt1 -- no cycle duration 
shortt1 <- shortt1[-(length(shortt1))]
shortt1_pc <- (shortt1/HTcycDt1)*100
avg_St1_pc <- mean(shortt1_pc) 

# Calculate Bird Tie lengthening ratio using cycle duration
longt1 <- HTcycDt1-shortt1
longt1_pc <- (longt1/HTcycDt1)*100
avg_Lt1_pc <- mean(longt1_pc) 

# Check the math
shortt1_pc + longt1_pc # should all sum to 100% -- good
avg_St1_pc + avg_Lt1_pc # should also sum to ~100% -- good

### Flight 2 ###
tie2=tie_L1$HT.sono
plot(tie2, type="l")

## Bird Tie - flight 2 HT peaks
pt2=find_peaks(tie2, m = 200)
points(pt2, tie2[pt2], col="blue", pch=19) # plots peaks in blue onto strain profile plot
pt2
length(pt2) # counts peaks

## Bird Tie - flight 2 HT troughs
tt2=find_peaks(-tie2, m = 200)
points(tt2, tie2[tt2], col="blue") # plots open blue points for trough locations
tt2
length(tt2) # counts troughs

# HT length change cycle duration for Bird Tie flight 2
HTpeakst2 <- tie_L1$Time[pt2] 
HTtroughst2 <- tie_L1$Time[tt2] 
HTcycDt2 <- diff(HTpeakst2) 

#  Last marked point is a trough, for simplicity, 
## calculate Bird Tie shortening ratio first and remove last peak
shortt2 <- HTtroughst2[-1] - HTpeakst2[-(length(HTpeakst2))] # first trough removed--no peak prior to it
shortt2_pc <- (shortt2/HTcycDt2)*100
avg_St2_pc <- mean(shortt2_pc) 

# Calculate Bird Tie lengthening ratio using cycle duration
longt2 <- HTcycDt2-shortt2
longt2_pc <- (longt2/HTcycDt2)*100
avg_Lt2_pc <- mean(longt2_pc) 

# Check the math
shortt2_pc + longt2_pc # should all sum to 100% -- good
avg_St2_pc + avg_Lt2_pc # should also sum to ~100% -- good

### Flight 3 ###
tie3=tie_L2$HT.sono
plot(tie3, type="l")

## Bird Tie - flight 3 HT peaks
pt3=find_peaks(tie3, m = 200)
points(pt3, tie3[pt3], col="blue", pch=19) # plots peaks in blue onto strain profile plot
pt3
length(pt3) # counts peaks

## Bird Tie - flight 3 HT troughs
tt3=find_peaks(-tie3, m = 200)
points(tt3, tie3[tt3], col="blue") # plots open blue points for trough locations
tt3
length(tt3) # counts troughs

# HT length change cycle duration for Bird Tie flight 3
HTpeakst3 <- tie_L2$Time[pt3] 
HTtroughst3 <- tie_L2$Time[tt3] 
HTcycDt3 <- diff(HTpeakst3) 

#  Last marked point is a trough, for simplicity, 
## calculate Bird Tie shortening ratio first
shortt3 <- HTtroughst3[-1] - HTpeakst3 # first trough removed--no peak prior to it
# remove last shortt3 -- no cycle duration 
shortt3 <- shortt3[-(length(shortt3))]
shortt3_pc <- (shortt3/HTcycDt3)*100
avg_St3_pc <- mean(shortt3_pc) 

# Calculate Bird Tie lengthening ratio using cycle duration
longt3 <- HTcycDt3-shortt3
longt3_pc <- (longt3/HTcycDt3)*100
avg_Lt3_pc <- mean(longt3_pc) 

# Check the math
shortt3_pc + longt3_pc # should all sum to 100% -- good
avg_St3_pc + avg_Lt3_pc # should also sum to ~100% -- good

### Flight 4 ###
tie4=tie_L3$HT.sono
plot(tie4, type="l")

## Bird Tie - flight 4 HT peaks
pt4=find_peaks(tie4, m = 200)
points(pt4, tie4[pt4], col="blue", pch=19) # plots peaks in blue onto strain profile plot
pt4
length(pt4) # counts peaks

## Bird Tie - flight 4 HT troughs
tt4=find_peaks(-tie4, m = 200)
points(tt4, tie4[tt4], col="blue") # plots open blue points for trough locations
tt4
length(tt4) # counts troughs

# HT length change cycle duration for Bird Tie flight 4
HTpeakst4 <- tie_L3$Time[pt4] 
HTtroughst4 <- tie_L3$Time[tt4] 
HTcycDt4 <- diff(HTpeakst4) 

#  Last marked point is a peak, for simplicity, 
## calculate Bird Tie shortening ratio first and remove last peak
shortt4 <- HTtroughst4[-1] - HTpeakst4[-(length(HTpeakst4))] # first trough removed--no peak prior to it
shortt4_pc <- (shortt4/HTcycDt4)*100
avg_St4_pc <- mean(shortt4_pc) 

# Calculate Bird Tie lengthening ratio using cycle duration
longt4 <- HTcycDt4-shortt4
longt4_pc <- (longt4/HTcycDt4)*100
avg_Lt4_pc <- mean(longt4_pc) 

# Check the math
shortt4_pc + longt4_pc # should all sum to 100% -- good
avg_St4_pc + avg_Lt4_pc # should also sum to ~100% -- good

### Flight 5 ###
tie5=tie_L4$HT.sono
plot(tie5, type="l")

## Bird Tie - flight 5 HT peaks
pt5=find_peaks(tie5, m = 200)
points(pt5, tie5[pt5], col="blue", pch=19) # plots peaks in blue onto strain profile plot
pt5
length(pt5) # counts peaks

## Bird Tie - flight 5 HT troughs
tt5=find_peaks(-tie5, m = 200)
points(tt5, tie5[tt5], col="blue") # plots open blue points for trough locations
tt5
length(tt5) # counts troughs

# HT length change cycle duration for Bird Tie flight 5
HTpeakst5 <- tie_L4$Time[pt5] 
HTtroughst5 <- tie_L4$Time[tt5] 
HTcycDt5 <- diff(HTpeakst5) 

#  Last marked point is a peak, for simplicity, 
## calculate Bird Tie shortening ratio first 
shortt5 <- HTtroughst5[-1] - HTpeakst5 # first trough removed--no peak prior to it
# remove last shortt5 -- no cycle duration 
shortt5 <- shortt5[-(length(shortt5))]
shortt5_pc <- (shortt5/HTcycDt5)*100
avg_St5_pc <- mean(shortt5_pc) 

# Calculate Bird Tie lengthening ratio using cycle duration
longt5 <- HTcycDt5-shortt5
longt5_pc <- (longt5/HTcycDt5)*100
avg_Lt5_pc <- mean(longt5_pc) 

# Check the math
shortt5_pc + longt5_pc # should all sum to 100% -- good
avg_St5_pc + avg_Lt5_pc # should also sum to ~100% -- good

### Find peaks for Bird Yyy's five humerotriceps SONO traces ###

### Flight 1 ###
yyy1=yyy_L1$HT.sono
plot(yyy1, type="l")

## Bird Yyy - flight 1 HT peaks
py1=find_peaks(yyy1, m = 200)
points(py1, yyy1[py1], col="blue", pch=19) # plots peaks in blue onto strain profile plot
py1
length(py1) # counts peaks
# last point is a duplicate--remove
py1<-py1[-(length(py1))]
# check it
plot(yyy1, type="l")
points(py1, yyy1[py1], col="blue", pch=19) # plots peaks in blue onto strain profile plot
py1
length(py1) # counts peaks

## Bird Yyy - flight 1 HT troughs
ty1=find_peaks(-yyy1, m = 200)
points(ty1, yyy1[ty1], col="blue") # plots open blue points for trough locations
ty1
length(ty1) # counts troughs

# too many troughs -- remove duplicates
# since values are both time and position duplicates, take first instance
# of its occurence if sono values are the same
ty1<-ty1[!duplicated(yyy_L1$HT.sono[ty1])]

# check it
ty1
length(ty1) # counts troughs
plot(yyy1, type="l")
points(py1, yyy1[py1], col="blue", pch=19) # plots peaks in blue onto strain profile plot
points(ty1, yyy1[ty1], col="blue") # plots open blue points for trough locations

# HT length change cycle duration for Bird Yyy flight 1
HTpeaksy1 <- yyy_L1$Time[py1] 
HTtroughsy1 <- yyy_L1$Time[ty1] 
HTcycDy1 <- diff(HTpeaksy1)

#  Last marked point is a trough, for simplicity, 
## calculate Bird Yyy shortening ratio first 
shorty1 <- HTtroughsy1[-1] - HTpeaksy1 # first trough removed--no peak prior to it
# remove last shorty1 -- no cycle duration 
shorty1 <- shorty1[-(length(shorty1))]
shorty1_pc <- (shorty1/HTcycDy1)*100
avg_Sy1_pc <- mean(shorty1_pc) 

# Calculate Bird Yyy lengthening ratio using cycle duration
longy1 <- HTcycDy1-shorty1
longy1_pc <- (longy1/HTcycDy1)*100
avg_Ly1_pc <- mean(longy1_pc) 

# Check the math
shorty1_pc + longy1_pc # should all sum to 100% -- good
avg_Sy1_pc + avg_Ly1_pc # should also sum to ~100% -- good

### Flight 2 ###
yyy2=yyy_T1$HT.sono
plot(yyy2, type="l")

## Bird Yyy - flight 2 HT peaks
py2=find_peaks(yyy2, m = 200)
points(py2, yyy2[py2], col="blue", pch=19) # plots peaks in blue onto strain profile plot
py2
length(py2) # counts peaks
# too many peaks -- remove duplicates
# since values are both time and position duplicates, take first instance
# of its occurence if sono values are the same
py2<-py2[!duplicated(yyy_T1$HT.sono[py2])]
# check it
py2
length(py2) # counts peaks
plot(yyy2, type="l")
points(py2, yyy2[py2], col="blue", pch=19) # plots peaks in blue onto strain profile plot

## Bird Yyy - flight 2 HT troughs
ty2=find_peaks(-yyy2, m = 200)
points(ty2, yyy2[ty2], col="blue") # plots open blue points for trough locations
ty2
length(ty2) # counts troughs
# too many troughs -- remove duplicates
# since values are both time and position duplicates, take first instance
# of its occurence if sono values are the same
ty2<-ty2[!duplicated(yyy_T1$HT.sono[ty2])]
# check it
ty2
length(ty2) # counts troughs
# missing 3rd last trough, add back manually
ty2<-c(ty2[1:14],8259,ty2[15:16])
# check it
ty2
length(ty2) # counts peaks
plot(yyy2, type="l")
points(py2, yyy2[py2], col="blue", pch=19) # plots peaks in blue onto strain profile plot
points(ty2, yyy2[ty2], col="blue") # plots open blue points for trough locations

# HT length change cycle duration for Bird Yyy flight 2
HTpeaksy2<- yyy_T1$Time[py2] 
HTtroughsy2 <- yyy_T1$Time[ty2] 
HTcycDy2 <- diff(HTpeaksy2)

#  Last marked point is a peak, for simplicity, 
## calculate Bird Yyy shortening ratio first and remove last peak
shorty2 <- HTtroughsy2[-1] - HTpeaksy2[-(length(HTpeaksy2))] # first trough removed--no peak prior to it
shorty2_pc <- (shorty2/HTcycDy2)*100
avg_Sy2_pc <- mean(shorty2_pc) 

# Calculate Bird Yyy lengthening ratio using cycle duration
longy2 <- HTcycDy2-shorty2
longy2_pc <- (longy2/HTcycDy2)*100
avg_Ly2_pc <- mean(longy2_pc) 

# Check the math
shorty2_pc + longy2_pc # should all sum to 100% -- good
avg_Sy2_pc + avg_Ly2_pc # should also sum to 100% -- good

### Flight 3 ###
yyy3=yyy_T2$HT.sono
plot(yyy3, type="l")

## Bird Yyy - flight 3 HT peaks
py3=find_peaks(yyy3, m = 200)
points(py3, yyy3[py3], col="blue", pch=19) # plots peaks in blue onto strain profile plot
py3
length(py3) # counts peaks

## Bird Yyy - flight 3 HT troughs
ty3=find_peaks(-yyy3, m = 200)
points(ty3, yyy3[ty3], col="blue") # plots open blue points for trough locations
ty3
length(ty3) # counts troughs

# HT length change cycle duration for Bird Yyy flight 3
HTpeaksy3 <- yyy_T2$Time[py3] 
HTtroughsy3 <- yyy_T2$Time[ty3] 
HTcycDy3 <- diff(HTpeaksy3[-1]) # first peak not true peak--remove

#  Last marked point is a trough, for simplicity, 
## calculate Bird Yyy shortening ratio first 
shorty3 <- HTtroughsy3[-1] - HTpeaksy3[-1] # first trough removed--no peak prior to it and remove first peak--not true peak
# remove last shorty3 -- no cycle duration 
shorty3 <- shorty3[-(length(shorty3))]
shorty3_pc <- (shorty3/HTcycDy3)*100
avg_Sy3_pc <- mean(shorty3_pc) 

# Calculate Bird Yyy lengthening ratio using cycle duration
longy3 <- HTcycDy3-shorty3
longy3_pc <- (longy3/HTcycDy3)*100
avg_Ly3_pc <- mean(longy3_pc) 

# Check the math
shorty3_pc + longy3_pc # should all sum to 100% -- good
avg_Sy3_pc + avg_Ly3_pc # should also sum to ~100% -- good

### Flight 4 ###
yyy4=yyy_T3$HT.sono
plot(yyy4, type="l")

## Bird Yyy - flight 4 HT peaks
py4=find_peaks(yyy4, m = 200)
points(py4, yyy4[py4], col="blue", pch=19) # plots peaks in blue onto strain profile plot
py4
length(py4) # counts peaks

## Bird Yyy - flight 4 HT troughs
ty4=find_peaks(-yyy4, m = 200)
points(ty4, yyy4[ty4], col="blue") # plots open blue points for trough locations
ty4
length(ty4) # counts troughs

# HT length change cycle duration for Bird Yyy flight 4
HTpeaksy4 <- yyy_T3$Time[py4] 
HTtroughsy4 <- yyy_T3$Time[ty4] 
HTcycDy4 <- diff(HTpeaksy4[-1]) # first peak not true peak--remove

#  Last marked point is a trough, for simplicity, 
## calculate Bird Yyy shortening ratio first 
shorty4 <- HTtroughsy4[-1] - HTpeaksy4[-1] # first trough removed--no peak prior to it and remove first peak--not true peak
# remove last shorty4 -- no cycle duration 
shorty4 <- shorty4[-(length(shorty4))]
shorty4_pc <- (shorty4/HTcycDy4)*100
avg_Sy4_pc <- mean(shorty4_pc) 

# Calculate Bird Yyy lengthening ratio using cycle duration
longy4 <- HTcycDy4-shorty4
longy4_pc <- (longy4/HTcycDy4)*100
avg_Ly4_pc <- mean(longy4_pc) 

# Check the math
shorty4_pc + longy4_pc # should all sum to 100% -- good
avg_Sy4_pc + avg_Ly4_pc # should also sum to ~100% -- good

### Flight 5 ###
yyy5=yyy_T4$HT.sono
plot(yyy5, type="l")

## Bird Yyy - flight 5 HT peaks
py5=find_peaks(yyy5, m = 200)
points(py5, yyy5[py5], col="blue", pch=19) # plots peaks in blue onto strain profile plot
py5
length(py5) # counts peaks
# one duplicate--remove manually
py5<-py5[-18]
# check it
py5
length(py5) # counts peaks
plot(yyy5, type="l")
points(py5, yyy5[py5], col="blue", pch=19) # plots peaks in blue onto strain profile plot

## Bird Yyy - flight 5 HT troughs
ty5=find_peaks(-yyy5, m = 200)
points(ty5, yyy5[ty5], col="blue") # plots open blue points for trough locations
ty5
length(ty5) # counts troughs

# HT length change cycle duration for Bird Yyy flight 5
HTpeaksy5 <- yyy_T4$Time[py5] 
HTtroughsy5 <- yyy_T4$Time[ty5] 
HTcycDy5 <- diff(HTpeaksy5)

#  Last marked point is a peak, for simplicity, 
## calculate Bird Yyy shortening ratio first and remove last peak
shorty5 <- HTtroughsy5[-1] - HTpeaksy5[-(length(HTpeaksy5))] # first trough removed--no peak prior to it
shorty5_pc <- (shorty5/HTcycDy5)*100
avg_Sy5_pc <- mean(shorty5_pc) 

# Calculate Bird Yyy lengthening ratio using cycle duration
longy5 <- HTcycDy5-shorty5
longy5_pc <- (longy5/HTcycDy5)*100
avg_Ly5_pc <- mean(longy5_pc) 

# Check the math
shorty5_pc + longy5_pc # should all sum to 100% -- good
avg_Sy5_pc + avg_Ly5_pc # should also sum to ~100% -- good

#################################################################################
#####################   Write results to text file   ############################
#################################################################################


# create a matrix with as many rows as trials and populate it
# with the average lengthening and shortening percentages of per trial, per bird and add a 3rd 
# column with the ID of the individual that the strain data were taken from.
x <- matrix(nrow=26, ncol=3, byrow=F)
x[,1] <- c("Www1","Www2","Www3","Www4","Www5","Xxx1","Xxx2","Xxx3","Xxx4","Xxx5",
		"Uuu1","Uuu2","Uuu3","Uuu4","Uuu5","Tie1","Tie2","Tie3","Tie4","Tie5",
		"Yyy1","Yyy2","Yyy3","Yyy4","Yyy5", "Grand means")
x[,2] <- c(avg_Sw1_pc,avg_Sw2_pc,avg_Sw3_pc,avg_Sw4_pc,avg_Sw5_pc,
		avg_Sx1_pc,avg_Sx2_pc,avg_Sx3_pc,avg_Sx4_pc,avg_Sx5_pc,
		avg_Su1_pc,avg_Su2_pc,avg_Su3_pc,avg_Su4_pc,avg_Su5_pc,
		avg_St1_pc,avg_St2_pc,avg_St3_pc,avg_St4_pc,avg_St5_pc,
		avg_Sy1_pc,avg_Sy2_pc,avg_Sy3_pc,avg_Sy4_pc,avg_Sy5_pc,
		mean(avg_Sw1_pc,avg_Sw2_pc,avg_Sw3_pc,avg_Sw4_pc,avg_Sw5_pc,
		avg_Sx1_pc,avg_Sx2_pc,avg_Sx3_pc,avg_Sx4_pc,avg_Sx5_pc,
		avg_Su1_pc,avg_Su2_pc,avg_Su3_pc,avg_Su4_pc,avg_Su5_pc,
		avg_St1_pc,avg_St2_pc,avg_St3_pc,avg_St4_pc,avg_St5_pc,
		avg_Sy1_pc,avg_Sy2_pc,avg_Sy3_pc,avg_Sy4_pc,avg_Sy5_pc))
x[,3] <- c(avg_Lw1_pc,avg_Lw2_pc,avg_Lw3_pc,avg_Lw4_pc,avg_Lw5_pc,
		avg_Lx1_pc,avg_Lx2_pc,avg_Lx3_pc,avg_Lx4_pc,avg_Lx5_pc,
		avg_Lu1_pc,avg_Lu2_pc,avg_Lu3_pc,avg_Lu4_pc,avg_Lu5_pc,
		avg_Lt1_pc,avg_Lt2_pc,avg_Lt3_pc,avg_Lt4_pc,avg_Lt5_pc,
		avg_Ly1_pc,avg_Ly2_pc,avg_Ly3_pc,avg_Ly4_pc,avg_Ly5_pc,
		mean(avg_Lw1_pc,avg_Lw2_pc,avg_Lw3_pc,avg_Lw4_pc,avg_Lw5_pc,
		avg_Lx1_pc,avg_Lx2_pc,avg_Lx3_pc,avg_Lx4_pc,avg_Lx5_pc,
		avg_Lu1_pc,avg_Lu2_pc,avg_Lu3_pc,avg_Lu4_pc,avg_Lu5_pc,
		avg_Lt1_pc,avg_Lt2_pc,avg_Lt3_pc,avg_Lt4_pc,avg_Lt5_pc,
		avg_Ly1_pc,avg_Ly2_pc,avg_Ly3_pc,avg_Ly4_pc,avg_Ly5_pc))

# create a data frame from the above generated matrix (converts values to numeric).
# can also be exported as csv if need be.
avg.L_Sratios <- data.frame(x, stringsAsFactors = F)
# rename the columns to properly identify the metrics
colnames(avg.L_Sratios) <- c("Bird/trial ID", "Mean shortening (%)", "Mean lengthening (%)")

# move up one folder level (out of 'in vivo strain data for analysis' folder) before saving text file
setwd(choose.dir())

## write .txt file of in vivo summary data for lengthening to shortening ratios 
write.table(avg.L_Sratios, "in vivo strain ratio data.txt", row.names=FALSE)