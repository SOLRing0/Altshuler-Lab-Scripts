# Script to determine the offset of the timing for peak length of the 
# humerotriceps relative to peak length of the pectoralis from all of 
# the Robertson and Biewener in vivo strain data
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

## Bird Www - flight 1 PC peaks
wpc1=www_L1$PC.sono
par(new=T)
plot(wpc1, type="l", col="red")
## Bird Www - flight 1 PC peaks
pcw1=find_peaks(wpc1,m=200)
points(pcw1, wpc1[pcw1], col="black", pch=19) # plots peaks in blue onto strain profile plot
pcw1
length(pcw1) # counts peaks

# PC length change cycle duration for Bird Www flight 1
HTpeaksw1 <- www_L1$Time[pw1] 
PCpeaksw1 <- www_L1$Time[pcw1] 
PCcycDw1 <- diff(PCpeaksw1)

# calculate the timing offset between HT peak length (PL) and pect PL
PLoffsetw1 <- HTpeaksw1 - PCpeaksw1[-(length(PCpeaksw1))] # remove last PC peak

# determine what that offset is as a percentage of the length change cycle
HT_pkoffsetw1 <- (PLoffsetw1/PCcycDw1)*100
avg_takeoffw1<-mean(HT_pkoffsetw1[1:5])
avg_midflightw1<-mean(HT_pkoffsetw1[6:((length(HT_pkoffsetw1))-5)])
avg_landingw1<-mean(HT_pkoffsetw1[((length(HT_pkoffsetw1))-4):(length(HT_pkoffsetw1))])
avg_HT_pkoffsetw1 <-mean(HT_pkoffsetw1)

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

## Bird Www - flight 2 PC peaks
wpc2=www_T1$PC.sono
par(new=T)
plot(wpc2, type="l", col="red")
## Bird Www - flight 2 PC peaks
pcw2=find_peaks(wpc2,m=200)
points(pcw2, wpc2[pcw2], col="black", pch=19) # plots peaks in blue onto strain profile plot
pcw2
length(pcw2) # counts peaks

# too many peaks -- remove duplicates
# since values are both time and position duplicates, take first instance
# of its occurence if sono values are the same
pcw2<-pcw2[!duplicated(www_T1$PC.sono[pcw2])]
# 7th, 10th and 2nd last peaks missing--add back manually
pcw2<-c(pcw2[1:6],3577,pcw2[7:9],5709,pcw2[10:13],9499,pcw2[14])
# check it
pcw2
length(pcw2) # counts peaks
plot(www2, type="l")
points(pw2, www2[pw2], col="blue", pch=19) # plots peaks in blue onto strain profile plot
par(new=T)
plot(wpc2, type="l", col="red")
points(pcw2, wpc2[pcw2], col="black", pch=19) # plots peaks in blue onto strain profile plot

# PC length change cycle duration for Bird Www flight 2
HTpeaksw2 <- www_T1$Time[pw2] 
PCpeaksw2 <- www_T1$Time[pcw2] 
PCcycDw2 <- diff(PCpeaksw2)

# calculate the timing offset between HT peak length (PL) and pect PL
PLoffsetw2 <- HTpeaksw2 - PCpeaksw2[-(length(PCpeaksw2))] # remove last PC peak

# determine what that offset is as a percentage of the length change cycle
HT_pkoffsetw2 <- (PLoffsetw2/PCcycDw2)*100
avg_takeoffw2<-mean(HT_pkoffsetw2[1:5])
avg_midflightw2<-mean(HT_pkoffsetw2[6:((length(HT_pkoffsetw2))-5)])
avg_landingw2<-mean(HT_pkoffsetw2[((length(HT_pkoffsetw2))-4):(length(HT_pkoffsetw2))])
avg_HT_pkoffsetw2 <-mean(HT_pkoffsetw2)

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

## Bird Www - flight 3 PC peaks
wpc3=www_T2$PC.sono
par(new=T)
plot(wpc3, type="l", col="red")
## Bird Www - flight 3 PC peaks
pcw3=find_peaks(wpc3,m=200)
points(pcw3, wpc3[pcw3], col="black", pch=19) # plots peaks in blue onto strain profile plot
pcw3
length(pcw3) # counts peaks

# too many peaks -- remove duplicates
# since values are both time and position duplicates, take first instance
# of its occurence if sono values are the same
pcw3<-pcw3[!duplicated(www_T2$PC.sono[pcw3])]
# check it
pcw3
length(pcw3) # counts peaks
plot(www3, type="l")
points(pw3, www3[pw3], col="blue", pch=19) # plots peaks in blue onto strain profile plot
par(new=T)
plot(wpc3, type="l", col="red")
points(pcw3, wpc3[pcw3], col="black", pch=19) # plots peaks in blue onto strain profile plot

# PC length change cycle duration for Bird Www flight 3
HTpeaksw3 <- www_T2$Time[pw3] 
PCpeaksw3 <- www_T2$Time[pcw3] 
PCcycDw3 <- diff(PCpeaksw3)

# calculate the timing offset between HT peak length (PL) and pect PL
PLoffsetw3 <- HTpeaksw3 - PCpeaksw3[-(length(PCpeaksw3))] # remove last PC peak
 
# determine what that offset is as a percentage of the length change cycle
HT_pkoffsetw3 <- (PLoffsetw3/PCcycDw3)*100
avg_takeoffw3<-mean(HT_pkoffsetw3[1:5])
avg_midflightw3<-mean(HT_pkoffsetw3[6:((length(HT_pkoffsetw3))-5)])
avg_landingw3<-mean(HT_pkoffsetw3[((length(HT_pkoffsetw3))-4):(length(HT_pkoffsetw3))])
avg_HT_pkoffsetw3 <-mean(HT_pkoffsetw3)

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

## Bird Www - flight 4 PC peaks
wpc4=www_T3$PC.sono
par(new=T)
plot(wpc4, type="l", col="red")
## Bird Www - flight 1 PC peaks
pcw4=find_peaks(wpc4,m=200)
points(pcw4, wpc4[pcw4], col="black", pch=19) # plots peaks in blue onto strain profile plot
pcw4
length(pcw4) # counts peaks

# too many peaks -- remove duplicates
# since values are both time and position duplicates, take first instance
# of its occurence if sono values are the same
pcw4<-pcw4[!duplicated(www_T3$PC.sono[pcw4])]
# 9th and 2nd last peaks missing--add back manually
pcw4<-c(pcw4[1:8],4810,pcw4[9:14],9588,pcw4[15])
# check it
pcw4
length(pcw4) # counts peaks
plot(www4, type="l")
points(pw4, www4[pw4], col="blue", pch=19) # plots peaks in blue onto strain profile plot
par(new=T)
plot(wpc4, type="l", col="red")
points(pcw4, wpc4[pcw4], col="black", pch=19) # plots peaks in blue onto strain profile plot

# PC length change cycle duration for Bird Www flight 4
HTpeaksw4 <- www_T3$Time[pw4] 
PCpeaksw4 <- www_T3$Time[pcw4] 
PCcycDw4 <- diff(PCpeaksw4)

# calculate the timing offset between HT peak length (PL) and pect PL
PLoffsetw4 <- HTpeaksw4 - PCpeaksw4[-(length(PCpeaksw4))] # remove last PC peak
 
# determine what that offset is as a percentage of the length change cycle
HT_pkoffsetw4 <- (PLoffsetw4/PCcycDw4)*100
avg_takeoffw4<-mean(HT_pkoffsetw4[1:5])
avg_midflightw4<-mean(HT_pkoffsetw4[6:((length(HT_pkoffsetw4))-5)])
avg_landingw4<-mean(HT_pkoffsetw4[((length(HT_pkoffsetw4))-4):(length(HT_pkoffsetw4))])
avg_HT_pkoffsetw4 <-mean(HT_pkoffsetw4)

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

## Bird Www - flight 5 PC peaks
wpc5=www_T4$PC.sono
par(new=T)
plot(wpc5, type="l", col="red")
## Bird Www - flight 5 PC peaks
pcw5=find_peaks(wpc5,m=200)
points(pcw5, wpc5[pcw5], col="black", pch=19) # plots peaks in blue onto strain profile plot
pcw5
length(pcw5) # counts peaks
# too many peaks -- remove duplicates
# since values are both time and position duplicates, take first instance
# of its occurence if sono values are the same
pcw5<-pcw5[!duplicated(www_T4$PC.sono[pcw5])]
# 11th peak missing--add back manually
pcw5<-c(pcw5[1:10],6071,pcw5[11:15])
# check it
pcw5
length(pcw5) # counts peaks
plot(www5, type="l")
points(pw5, www5[pw5], col="blue", pch=19) # plots peaks in blue onto strain profile plot
par(new=T)
plot(wpc5, type="l", col="red")
points(pcw5, wpc5[pcw5], col="black", pch=19) # plots peaks in blue onto strain profile plot

# PC length change cycle duration for Bird Www flight 5
HTpeaksw5 <- www_T4$Time[pw5] 
PCpeaksw5 <- www_T4$Time[pcw5] 
PCcycDw5 <- diff(PCpeaksw5)

# calculate the timing offset between HT peak length (PL) and pect PL
PLoffsetw5 <- HTpeaksw5 - PCpeaksw5 

# determine what that offset is as a percentage of the length change cycle
HT_pkoffsetw5 <- (PLoffsetw5[-(length(PLoffsetw5))]/PCcycDw5)*100 # remove last offset--no PCcyc duration
avg_takeoffw5<-mean(HT_pkoffsetw5[1:5])
avg_midflightw5<-mean(HT_pkoffsetw5[6:((length(HT_pkoffsetw5))-5)])
avg_landingw5<-mean(HT_pkoffsetw5[((length(HT_pkoffsetw5))-4):(length(HT_pkoffsetw5))])
avg_HT_pkoffsetw5 <-mean(HT_pkoffsetw5)


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

## Bird Xxx - flight 1 PC peaks
xpc1=xxx_L1$PC.sono
par(new=T)
plot(xpc1, type="l", col="red")
## Bird Xxx - flight 1 PC peaks
pcx1=find_peaks(xpc1,m=200)
points(pcx1, xpc1[pcx1], col="black", pch=19) # plots peaks in blue onto strain profile plot
pcx1
length(pcx1) # counts peaks
# too many peaks -- remove duplicates
# since values are both time and position duplicates, take first instance
# of its occurence if sono values are the same
pcx1<-pcx1[!duplicated(xxx_L1$PC.sono[pcx1])]
# 7th and 16th peaks missing--add back manually
pcx1<-c(pcx1[1:6],3604,pcx1[7:16],10504,pcx1[17])
# check it
pcx1
length(pcx1) # counts peaks
plot(xxx1, type="l")
points(px1, xxx1[px1], col="blue", pch=19) # plots peaks in blue onto strain profile plot
par(new=T)
plot(xpc1, type="l", col="red")
points(pcx1, xpc1[pcx1], col="black", pch=19) # plots peaks in blue onto strain profile plot

# PC length change cycle duration for Bird Xxx flight 1
HTpeaksx1 <- xxx_L1$Time[px1] 
PCpeaksx1 <- xxx_L1$Time[pcx1] 
PCcycDx1 <- diff(PCpeaksx1)

# calculate the timing offset between HT peak length (PL) and pect PL
PLoffsetx1 <- HTpeaksx1 - PCpeaksx1[-(length(PCpeaksx1))] 

# determine what that offset is as a percentage of the length change cycle
HT_pkoffsetx1 <- (PLoffsetx1/PCcycDx1)*100
avg_takeoffx1<-mean(HT_pkoffsetx1[1:5])
avg_midflightx1<-mean(HT_pkoffsetx1[6:((length(HT_pkoffsetx1))-5)])
avg_landingx1<-mean(HT_pkoffsetx1[((length(HT_pkoffsetx1))-4):(length(HT_pkoffsetx1))])
avg_HT_pkoffsetx1 <-mean(HT_pkoffsetx1)

### Flight 2 ###
xxx2=xxx_T1$HT.sono
plot(xxx2, type="l")

## Bird Xxx - flight 2 HT peaks
px2=sort(findpeaks(xxx2,minpeakdistance=500)[,2])
points(px2, xxx2[px2], col="blue", pch=19) # plots peaks in blue onto strain profile plot
px2
length(px2) # counts peaks

## Bird Xxx - flight 2 PC peaks
xpc2=xxx_T1$PC.sono
par(new=T)
plot(xpc2, type="l", col="red")
## Bird Xxx - flight 2 PC peaks
pcx2=find_peaks(xpc2,m=200)
points(pcx2, xpc2[pcx2], col="black", pch=19) # plots peaks in blue onto strain profile plot
pcx2
length(pcx2) # counts peaks
# too many peaks -- remove duplicates
# since values are both time and position duplicates, take first instance
# of its occurence if sono values are the same
pcx2<-pcx2[!duplicated(xxx_T1$PC.sono[pcx2])]
# check it
pcx2
length(pcx2) # counts peaks
plot(xxx2, type="l")
points(px2, xxx2[px2], col="blue", pch=19) # plots peaks in blue onto strain profile plot
par(new=T)
plot(xpc2, type="l", col="red")
points(pcx2, xpc2[pcx2], col="black", pch=19) # plots peaks in blue onto strain profile plot

# PC length change cycle duration for Bird Xxx flight 2
HTpeaksx2 <- xxx_T1$Time[px2] 
PCpeaksx2 <- xxx_T1$Time[pcx2] 
PCcycDx2 <- diff(PCpeaksx2)

# calculate the timing offset between HT peak length (PL) and pect PL
PLoffsetx2 <- HTpeaksx2 - PCpeaksx2[-(length(PCpeaksx2))] 

# determine what that offset is as a percentage of the length change cycle
HT_pkoffsetx2 <- (PLoffsetx2/PCcycDx2)*100
avg_takeoffx2<-mean(HT_pkoffsetx2[1:5])
avg_midflightx2<-mean(HT_pkoffsetx2[6:((length(HT_pkoffsetx2))-5)])
avg_landingx2<-mean(HT_pkoffsetx2[((length(HT_pkoffsetx2))-4):(length(HT_pkoffsetx2))])
avg_HT_pkoffsetx2 <-mean(HT_pkoffsetx2)

### Flight 3 ###
xxx3=xxx_T2$HT.sono
plot(xxx3, type="l")

## Bird Xxx - flight 3 HT peaks
px3=sort(findpeaks(xxx3,minpeakdistance=500)[,2])
points(px3, xxx3[px3], col="blue", pch=19) # plots peaks in blue onto strain profile plot
px3
length(px3) # counts peaks

## Bird Xxx - flight 3 PC peaks
xpc3=xxx_T2$PC.sono
par(new=T)
plot(xpc3, type="l", col="red")
## Bird Xxx - flight 3 PC peaks
pcx3=find_peaks(xpc3,m=200)
points(pcx3, xpc3[pcx3], col="black", pch=19) # plots peaks in blue onto strain profile plot
pcx3
length(pcx3) # counts peaks
# too many peaks -- remove duplicates
# since values are both time and position duplicates, take first instance
# of its occurence if sono values are the same
pcx3<-pcx3[!duplicated(xxx_T2$PC.sono[pcx3])]
# 11th and 17th-20th peaks missing--add back manually
pcx3<-c(pcx3[1:10],6051,pcx3[11:15],9783,10392,11018,11632,pcx3[16])
# check it
pcx3
length(pcx3) # counts peaks
plot(xxx3, type="l")
points(px3, xxx3[px3], col="blue", pch=19) # plots peaks in blue onto strain profile plot
par(new=T)
plot(xpc3, type="l", col="red")
points(pcx3, xpc3[pcx3], col="black", pch=19) # plots peaks in blue onto strain profile plot

# PC length change cycle duration for Bird Xxx flight 3
HTpeaksx3 <- xxx_T2$Time[px3] 
PCpeaksx3 <- xxx_T2$Time[pcx3] 
PCcycDx3 <- diff(PCpeaksx3)

# calculate the timing offset between HT peak length (PL) and pect PL
PLoffsetx3 <- HTpeaksx3 - PCpeaksx3[-(length(PCpeaksx3))] 

# determine what that offset is as a percentage of the length change cycle
HT_pkoffsetx3 <- (PLoffsetx3/PCcycDx3)*100
avg_takeoffx3<-mean(HT_pkoffsetx3[1:5])
avg_midflightx3<-mean(HT_pkoffsetx3[6:((length(HT_pkoffsetx3))-5)])
avg_landingx3<-mean(HT_pkoffsetx3[((length(HT_pkoffsetx3))-4):(length(HT_pkoffsetx3))])
avg_HT_pkoffsetx3 <-mean(HT_pkoffsetx3)

### Flight 4 ###
xxx4=xxx_T3$HT.sono
plot(xxx4, type="l")

## Bird Xxx - flight 4 HT peaks
px4=sort(findpeaks(xxx4,minpeakdistance=500)[,2])
points(px4, xxx4[px4], col="blue", pch=19) # plots peaks in blue onto strain profile plot
px4
length(px4) # counts peaks

## Bird Xxx - flight 4 PC peaks
xpc4=xxx_T3$PC.sono
par(new=T)
plot(xpc4, type="l", col="red")
## Bird Xxx - flight 4 PC peaks
pcx4=find_peaks(xpc4,m=200)
points(pcx4, xpc4[pcx4], col="black", pch=19) # plots peaks in blue onto strain profile plot
pcx4
length(pcx4) # counts peaks
# too many peaks -- remove duplicates
# since values are both time and position duplicates, take first instance
# of its occurence if sono values are the same
pcx4<-pcx4[!duplicated(xxx_T3$PC.sono[pcx4])]
# check it
pcx4
length(pcx4) # counts peaks
plot(xxx4, type="l")
points(px4, xxx4[px4], col="blue", pch=19) # plots peaks in blue onto strain profile plot
par(new=T)
plot(xpc4, type="l", col="red")
points(pcx4, xpc4[pcx4], col="black", pch=19) # plots peaks in blue onto strain profile plot

# PC length change cycle duration for Bird Xxx flight 4
HTpeaksx4 <- xxx_T3$Time[px4] 
PCpeaksx4 <- xxx_T3$Time[pcx4] 
PCcycDx4 <- diff(PCpeaksx4)

# calculate the timing offset between HT peak length (PL) and pect PL
PLoffsetx4 <- HTpeaksx4 - PCpeaksx4[-(length(PCpeaksx4))] 

# determine what that offset is as a percentage of the length change cycle
HT_pkoffsetx4 <- (PLoffsetx4/PCcycDx4)*100
avg_takeoffx4<-mean(HT_pkoffsetx4[1:5])
avg_midflightx4<-mean(HT_pkoffsetx4[6:((length(HT_pkoffsetx4))-5)])
avg_landingx4<-mean(HT_pkoffsetx4[((length(HT_pkoffsetx4))-4):(length(HT_pkoffsetx4))])
avg_HT_pkoffsetx4 <-mean(HT_pkoffsetx4)

### Flight 5 ###
xxx5=xxx_T4$HT.sono
plot(xxx5, type="l")

## Bird Xxx - flight 5 HT peaks
px5=sort(findpeaks(xxx5,minpeakdistance=500)[,2])
points(px5, xxx5[px5], col="blue", pch=19) # plots peaks in blue onto strain profile plot
px5
length(px5) # counts peaks

## Bird Xxx - flight 5 PC peaks
xpc5=xxx_T4$PC.sono
par(new=T)
plot(xpc5, type="l", col="red")
## Bird Xxx - flight 5 PC peaks
pcx5=find_peaks(xpc5,m=200)
points(pcx5, xpc5[pcx5], col="black", pch=19) # plots peaks in blue onto strain profile plot
pcx5
length(pcx5) # counts peaks
# too many peaks -- remove duplicates
# since values are both time and position duplicates, take first instance
# of its occurence if sono values are the same
pcx5<-pcx5[!duplicated(xxx_T4$PC.sono[pcx5])]
# 3rd and 19th peaks missing--add back manually
pcx5<-c(pcx5[1:2],1102,pcx5[3:17],11025,pcx5[18:19])
# check it
pcx5
length(pcx5) # counts peaks
plot(xxx5, type="l")
points(px5, xxx5[px5], col="blue", pch=19) # plots peaks in blue onto strain profile plot
par(new=T)
plot(xpc5, type="l", col="red")
points(pcx5, xpc5[pcx5], col="black", pch=19) # plots peaks in blue onto strain profile plot

# PC length change cycle duration for Bird Xxx flight 5
HTpeaksx5 <- xxx_T4$Time[px5] 
PCpeaksx5 <- xxx_T4$Time[pcx5] 
PCcycDx5 <- diff(PCpeaksx5)

# calculate the timing offset between HT peak length (PL) and pect PL
PLoffsetx5 <- HTpeaksx5 - PCpeaksx5[-(length(PCpeaksx5))] 

# determine what that offset is as a percentage of the length change cycle
HT_pkoffsetx5 <- (PLoffsetx5/PCcycDx5)*100
avg_takeoffx5<-mean(HT_pkoffsetx5[1:5])
avg_midflightx5<-mean(HT_pkoffsetx5[6:((length(HT_pkoffsetx5))-5)])
avg_landingx5<-mean(HT_pkoffsetx5[((length(HT_pkoffsetx5))-4):(length(HT_pkoffsetx5))])
avg_HT_pkoffsetx5 <-mean(HT_pkoffsetx5)

### Find peaks for Bird Uuu's five humerotriceps SONO traces ###

### Flight 1 ###
uuu1=uuu_T1$HT.sono
plot(uuu1, type="l")

## Bird Uuu - flight 1 HT peaks
pu1=find_peaks(uuu1, m = 200)
points(pu1, uuu1[pu1], col="blue", pch=19) # plots peaks in blue onto strain profile plot
pu1
length(pu1) # counts peaks
# remove first HT peak, not true peak
pu1<-pu1[-1]
length(pu1) # counts peaks
pu1

## Bird Uuu - flight 1 PC peaks
upc1=uuu_T1$PC.sono
par(new=T)
plot(upc1, type="l", col="red")
## Bird Uuu - flight 1 PC peaks
pcu1=find_peaks(upc1,m=200)
points(pcu1, upc1[pcu1], col="black", pch=19) # plots peaks in blue onto strain profile plot
pcu1
length(pcu1) # counts peaks

# PC length change cycle duration for Bird Uuu flight 1
HTpeaksu1 <- uuu_T1$Time[pu1] 
PCpeaksu1 <- uuu_T1$Time[pcu1] 
PCcycDu1 <- diff(PCpeaksu1)

# calculate the timing offset between HT peak length (PL) and pect PL
PLoffsetu1 <- HTpeaksu1 - PCpeaksu1[-(length(PCpeaksu1))] 

# determine what that offset is as a percentage of the length change cycle
HT_pkoffsetu1 <- (PLoffsetu1/PCcycDu1)*100
avg_takeoffu1<-mean(HT_pkoffsetu1[1:5])
avg_midflightu1<-mean(HT_pkoffsetu1[6:((length(HT_pkoffsetu1))-5)])
avg_landingu1<-mean(HT_pkoffsetu1[((length(HT_pkoffsetu1))-4):(length(HT_pkoffsetu1))])
avg_HT_pkoffsetu1 <-mean(HT_pkoffsetu1)

### Flight 2 ###
uuu2=uuu_T2$HT.sono
plot(uuu2, type="l")

## Bird Uuu - flight 2 HT peaks
pu2=find_peaks(uuu2, m = 200)
points(pu2, uuu2[pu2], col="blue", pch=19) # plots peaks in blue onto strain profile plot
pu2
length(pu2) # counts peaks

## Bird Uuu - flight 2 PC peaks
upc2=uuu_T2$PC.sono
par(new=T)
plot(upc2, type="l", col="red")
## Bird Uuu - flight 2 PC peaks
pcu2=find_peaks(upc2,m=200)
points(pcu2, upc2[pcu2], col="black", pch=19) # plots peaks in blue onto strain profile plot
pcu2
length(pcu2) # counts peaks

# PC length change cycle duration for Bird Uuu flight 2
HTpeaksu2 <- uuu_T2$Time[pu2] 
PCpeaksu2 <- uuu_T2$Time[pcu2] 
PCcycDu2 <- diff(PCpeaksu2)

# calculate the timing offset between HT peak length (PL) and pect PL
PLoffsetu2 <- HTpeaksu2 - PCpeaksu2[-(length(PCpeaksu2))] 

# determine what that offset is as a percentage of the length change cycle
HT_pkoffsetu2 <- (PLoffsetu2/PCcycDu2)*100
avg_takeoffu2<-mean(HT_pkoffsetu2[1:5])
avg_midflightu2<-mean(HT_pkoffsetu2[6:((length(HT_pkoffsetu2))-5)])
avg_landingu2<-mean(HT_pkoffsetu2[((length(HT_pkoffsetu2))-4):(length(HT_pkoffsetu2))])
avg_HT_pkoffsetu2 <-mean(HT_pkoffsetu2)

### Flight 3 ###
uuu3=uuu_L1$HT.sono
plot(uuu3, type="l")

## Bird Uuu - flight 3 HT peaks
pu3=find_peaks(uuu3, m = 200)
points(pu3, uuu3[pu3], col="blue", pch=19) # plots peaks in blue onto strain profile plot
pu3
length(pu3) # counts peaks

## Bird Uuu - flight 3 PC peaks
upc3=uuu_L1$PC.sono
par(new=T)
plot(upc3, type="l", col="red")
## Bird Uuu - flight 3 PC peaks
pcu3=find_peaks(upc3,m=200)
points(pcu3, upc3[pcu3], col="black", pch=19) # plots peaks in blue onto strain profile plot
pcu3
length(pcu3) # counts peaks

# PC length change cycle duration for Bird Uuu flight 3
HTpeaksu3 <- uuu_L1$Time[pu3] 
PCpeaksu3 <- uuu_L1$Time[pcu3] 
PCcycDu3 <- diff(PCpeaksu3)

# calculate the timing offset between HT peak length (PL) and pect PL
PLoffsetu3 <- HTpeaksu3 - PCpeaksu3[-(length(PCpeaksu3))] 

# determine what that offset is as a percentage of the length change cycle
HT_pkoffsetu3 <- (PLoffsetu3/PCcycDu3)*100
avg_takeoffu3<-mean(HT_pkoffsetu3[1:5])
avg_midflightu3<-mean(HT_pkoffsetu3[6:((length(HT_pkoffsetu3))-5)])
avg_landingu3<-mean(HT_pkoffsetu3[((length(HT_pkoffsetu3))-4):(length(HT_pkoffsetu3))])
avg_HT_pkoffsetu3 <-mean(HT_pkoffsetu3)

### Flight 4 ###
uuu4=uuu_L2$HT.sono
plot(uuu4, type="l")

## Bird Uuu - flight 4 HT peaks
pu4=find_peaks(uuu4, m = 200)
points(pu4, uuu4[pu4], col="blue", pch=19) # plots peaks in blue onto strain profile plot
pu4
length(pu4) # counts peaks
# remove first HT peak -- not true peak
pu4<-pu4[-1]
pu4
length(pu4)
# re-plot 
plot(uuu4, type="l")
points(pu4, uuu4[pu4], col="blue", pch=19) # plots peaks in blue onto strain profile plot

## Bird Uuu - flight 4 PC peaks
upc4=uuu_L2$PC.sono
par(new=T)
plot(upc4, type="l", col="red")
## Bird Uuu - flight 4 PC peaks
pcu4=find_peaks(upc4,m=200)
points(pcu4, upc4[pcu4], col="black", pch=19) # plots peaks in blue onto strain profile plot
pcu4
length(pcu4) # counts peaks

# PC length change cycle duration for Bird Uuu flight 4
HTpeaksu4 <- uuu_L2$Time[pu4] 
PCpeaksu4 <- uuu_L2$Time[pcu4] 
PCcycDu4 <- diff(PCpeaksu4)

# calculate the timing offset between HT peak length (PL) and pect PL
PLoffsetu4 <- HTpeaksu4 - PCpeaksu4[-(length(PCpeaksu4))] 

# determine what that offset is as a percentage of the length change cycle
HT_pkoffsetu4 <- (PLoffsetu4/PCcycDu4)*100
avg_takeoffu4<-mean(HT_pkoffsetu4[1:5])
avg_midflightu4<-mean(HT_pkoffsetu4[6:((length(HT_pkoffsetu4))-5)])
avg_landingu4<-mean(HT_pkoffsetu4[((length(HT_pkoffsetu4))-4):(length(HT_pkoffsetu4))])
avg_HT_pkoffsetu4 <-mean(HT_pkoffsetu4)

### Flight 5 ###
uuu5=uuu_L3$HT.sono
plot(uuu5, type="l")

## Bird Uuu - flight 5 HT peaks
pu5=find_peaks(uuu5, m = 200)
points(pu5, uuu5[pu5], col="blue", pch=19) # plots peaks in blue onto strain profile plot
pu5
length(pu5) # counts peaks

## Bird Uuu - flight 5 PC peaks
upc5=uuu_L3$PC.sono
par(new=T)
plot(upc5, type="l", col="red")
## Bird Uuu - flight 5 PC peaks
pcu5=find_peaks(upc5,m=200)
points(pcu5, upc5[pcu5], col="black", pch=19) # plots peaks in blue onto strain profile plot
pcu5
length(pcu5) # counts peaks

# PC length change cycle duration for Bird Uuu flight 5
HTpeaksu5 <- uuu_L3$Time[pu5] 
PCpeaksu5 <- uuu_L3$Time[pcu5] 
PCcycDu5 <- diff(PCpeaksu5)

# calculate the timing offset between HT peak length (PL) and pect PL
PLoffsetu5 <- HTpeaksu5 - PCpeaksu5[-(length(PCpeaksu5))] 

# determine what that offset is as a percentage of the length change cycle
HT_pkoffsetu5 <- (PLoffsetu5/PCcycDu5)*100
avg_takeoffu5<-mean(HT_pkoffsetu5[1:5])
avg_midflightu5<-mean(HT_pkoffsetu5[6:((length(HT_pkoffsetu5))-5)])
avg_landingu5<-mean(HT_pkoffsetu5[((length(HT_pkoffsetu5))-4):(length(HT_pkoffsetu5))])
avg_HT_pkoffsetu5 <-mean(HT_pkoffsetu5)

### Find peaks for Bird Tie's five humerotriceps SONO traces ###

### Flight 1 ###
tie1=tie_T1$HT.sono
plot(tie1, type="l")

## Bird Tie - flight 1 HT peaks
pt1=find_peaks(tie1, m = 200)
points(pt1, tie1[pt1], col="blue", pch=19) # plots peaks in blue onto strain profile plot
pt1
length(pt1) # counts peaks

## Bird Tie - flight 1 PC peaks
tpc1=tie_T1$PC.sono
par(new=T)
plot(tpc1, type="l", col="red")
## Bird Tie - flight 1 PC peaks
pct1=find_peaks(tpc1,m=200)
points(pct1, tpc1[pct1], col="black", pch=19) # plots peaks in blue onto strain profile plot
pct1
length(pct1) # counts peaks

# PC length change cycle duration for Bird Tie flight 1
HTpeakst1 <- tie_T1$Time[pt1] 
PCpeakst1 <- tie_T1$Time[pct1] 
PCcycDt1 <- diff(PCpeakst1)

# calculate the timing offset between HT peak length (PL) and pect PL
PLoffsett1 <- HTpeakst1 - PCpeakst1 

# determine what that offset is as a percentage of the length change cycle
HT_pkoffsett1 <- (PLoffsett1[-(length(PLoffsett1))]/PCcycDt1)*100
avg_takeofft1<-mean(HT_pkoffsett1[1:5])
avg_midflightt1<-mean(HT_pkoffsett1[6:((length(HT_pkoffsett1))-5)])
avg_landingt1<-mean(HT_pkoffsett1[((length(HT_pkoffsett1))-4):(length(HT_pkoffsett1))])
avg_HT_pkoffsett1 <-mean(HT_pkoffsett1)

### Flight 2 ###
tie2=tie_L1$HT.sono
plot(tie2, type="l")

## Bird Tie - flight 2 HT peaks
pt2=find_peaks(tie2, m = 200)
points(pt2, tie2[pt2], col="blue", pch=19) # plots peaks in blue onto strain profile plot
pt2
length(pt2) # counts peaks

## Bird Tie - flight 2 PC peaks
tpc2=tie_L1$PC.sono
par(new=T)
plot(tpc2, type="l", col="red")
## Bird Tie - flight 2 PC peaks
pct2=find_peaks(tpc2,m=200)
points(pct2, tpc2[pct2], col="black", pch=19) # plots peaks in blue onto strain profile plot
pct2
length(pct2) # counts peaks

# PC length change cycle duration for Bird Tie flight 2
HTpeakst2 <- tie_L1$Time[pt2] 
PCpeakst2 <- tie_L1$Time[pct2] 
PCcycDt2 <- diff(PCpeakst2)

# calculate the timing offset between HT peak length (PL) and pect PL
PLoffsett2 <- HTpeakst2 - PCpeakst2 

# determine what that offset is as a percentage of the length change cycle
HT_pkoffsett2 <- (PLoffsett2[-(length(PLoffsett2))]/PCcycDt2)*100
avg_takeofft2<-mean(HT_pkoffsett2[1:5])
avg_midflightt2<-mean(HT_pkoffsett2[6:((length(HT_pkoffsett2))-5)])
avg_landingt2<-mean(HT_pkoffsett2[((length(HT_pkoffsett2))-4):(length(HT_pkoffsett2))])
avg_HT_pkoffsett2 <-mean(HT_pkoffsett2)

### Flight 3 ###
tie3=tie_L2$HT.sono
plot(tie3, type="l")

## Bird Tie - flight 3 HT peaks
pt3=find_peaks(tie3, m = 200)
points(pt3, tie3[pt3], col="blue", pch=19) # plots peaks in blue onto strain profile plot
pt3
length(pt3) # counts peaks

## Bird Tie - flight 3 PC peaks
tpc3=tie_L2$PC.sono
par(new=T)
plot(tpc3, type="l", col="red")
## Bird Tie - flight 3 PC peaks
pct3=find_peaks(tpc3,m=200)
points(pct3, tpc3[pct3], col="black", pch=19) # plots peaks in blue onto strain profile plot
pct3
length(pct3) # counts peaks

# PC length change cycle duration for Bird Tie flight 3
HTpeakst3 <- tie_L2$Time[pt3] 
PCpeakst3 <- tie_L2$Time[pct3] 
PCcycDt3 <- diff(PCpeakst3)

# calculate the timing offset between HT peak length (PL) and pect PL
PLoffsett3 <- HTpeakst3 - PCpeakst3 

# determine what that offset is as a percentage of the length change cycle
HT_pkoffsett3 <- (PLoffsett3[-(length(PLoffsett3))]/PCcycDt3)*100
avg_takeofft3<-mean(HT_pkoffsett3[1:5])
avg_midflightt3<-mean(HT_pkoffsett3[6:((length(HT_pkoffsett3))-5)])
avg_landingt3<-mean(HT_pkoffsett3[((length(HT_pkoffsett3))-4):(length(HT_pkoffsett3))])
avg_HT_pkoffsett3 <-mean(HT_pkoffsett3)

### Flight 4 ###
tie4=tie_L3$HT.sono
plot(tie4, type="l")

## Bird Tie - flight 4 HT peaks
pt4=find_peaks(tie4, m = 200)
points(pt4, tie4[pt4], col="blue", pch=19) # plots peaks in blue onto strain profile plot
pt4
length(pt4) # counts peaks

## Bird Tie - flight 4 PC peaks
tpc4=tie_L3$PC.sono
par(new=T)
plot(tpc4, type="l", col="red")
## Bird Tie - flight 4 PC peaks
pct4=find_peaks(tpc4,m=200)
points(pct4, tpc4[pct4], col="black", pch=19) # plots peaks in blue onto strain profile plot
pct4
length(pct4) # counts peaks

# PC length change cycle duration for Bird Tie flight 4
HTpeakst4 <- tie_L3$Time[pt4] 
PCpeakst4 <- tie_L3$Time[pct4] 
PCcycDt4 <- diff(PCpeakst4)

# calculate the timing offset between HT peak length (PL) and pect PL
PLoffsett4 <- HTpeakst4 - PCpeakst4 

# determine what that offset is as a percentage of the length change cycle
HT_pkoffsett4 <- (PLoffsett4[-(length(PLoffsett4))]/PCcycDt4)*100
avg_takeofft4<-mean(HT_pkoffsett4[1:5])
avg_midflightt4<-mean(HT_pkoffsett4[6:((length(HT_pkoffsett4))-5)])
avg_landingt4<-mean(HT_pkoffsett4[((length(HT_pkoffsett4))-4):(length(HT_pkoffsett4))])
avg_HT_pkoffsett4 <-mean(HT_pkoffsett4)

### Flight 5 ###
tie5=tie_L4$HT.sono
plot(tie5, type="l")

## Bird Tie - flight 5 HT peaks
pt5=find_peaks(tie5, m = 200)
points(pt5, tie5[pt5], col="blue", pch=19) # plots peaks in blue onto strain profile plot
pt5
length(pt5) # counts peaks

## Bird Tie - flight 5 PC peaks
tpc5=tie_L4$PC.sono
par(new=T)
plot(tpc5, type="l", col="red")
## Bird Tie - flight 5 PC peaks
pct5=find_peaks(tpc5,m=200)
points(pct5, tpc5[pct5], col="black", pch=19) # plots peaks in blue onto strain profile plot
pct5
length(pct5) # counts peaks

# PC length change cycle duration for Bird Tie flight 5
HTpeakst5 <- tie_L4$Time[pt5] 
PCpeakst5 <- tie_L4$Time[pct5] 
PCcycDt5 <- diff(PCpeakst5)

# calculate the timing offset between HT peak length (PL) and pect PL
PLoffsett5 <- HTpeakst5 - PCpeakst5 

# determine what that offset is as a percentage of the length change cycle
HT_pkoffsett5 <- (PLoffsett5[-(length(PLoffsett5))]/PCcycDt5)*100
avg_takeofft5<-mean(HT_pkoffsett5[1:5])
avg_midflightt5<-mean(HT_pkoffsett5[6:((length(HT_pkoffsett5))-5)])
avg_landingt5<-mean(HT_pkoffsett5[((length(HT_pkoffsett5))-4):(length(HT_pkoffsett5))])
avg_HT_pkoffsett5 <-mean(HT_pkoffsett5)

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

## Bird Yyy - flight 1 PC peaks
ypc1=yyy_L1$PC.sono
par(new=T)
plot(ypc1, type="l", col="red")
## Bird Yyy - flight 1 PC peaks
pcy1=find_peaks(ypc1,m=200)
points(pcy1, ypc1[pcy1], col="black", pch=19) # plots peaks in blue onto strain profile plot
pcy1
length(pcy1) # counts peaks
# too many peaks -- remove duplicates
# since values are both time and position duplicates, take first instance
# of its occurence if sono values are the same
pcy1<-pcy1[!duplicated(yyy_L1$PC.sono[pcy1])]
# 6th, 11th, 13th and 15th peaks missing--add back manually
pcy1<-c(pcy1[1:5],2760,pcy1[6:9],5768,pcy1[10],7056,pcy1[11],8334,pcy1[12:15])
# check it
pcy1
length(pcy1) # counts peaks
plot(yyy1, type="l")
points(py1, yyy1[py1], col="blue", pch=19) # plots peaks in blue onto strain profile plot
par(new=T)
plot(ypc1, type="l", col="red")
points(pcy1, ypc1[pcy1], col="black", pch=19) # plots peaks in blue onto strain profile plot

# PC length change cycle duration for Bird Yyy flight 1
HTpeaksy1 <- yyy_L1$Time[py1] 
PCpeaksy1 <- yyy_L1$Time[pcy1] 
PCcycDy1 <- diff(PCpeaksy1)

# calculate the timing offset between HT peak length (PL) and pect PL
PLoffsety1 <- HTpeaksy1 - PCpeaksy1[-(length(PCpeaksy1))] 

# determine what that offset is as a percentage of the length change cycle
HT_pkoffsety1 <- (PLoffsety1/PCcycDy1)*100
avg_takeoffy1<-mean(HT_pkoffsety1[1:5])
avg_midflighty1<-mean(HT_pkoffsety1[6:((length(HT_pkoffsety1))-5)])
avg_landingy1<-mean(HT_pkoffsety1[((length(HT_pkoffsety1))-4):(length(HT_pkoffsety1))])
avg_HT_pkoffsety1 <-mean(HT_pkoffsety1)

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

## Bird Yyy - flight 2 PC peaks
ypc2=yyy_T1$PC.sono
par(new=T)
plot(ypc2, type="l", col="red")
## Bird Yyy - flight 2 PC peaks
pcy2=find_peaks(ypc2,m=200)
points(pcy2, ypc2[pcy2], col="black", pch=19) # plots peaks in blue onto strain profile plot
pcy2
length(pcy2) # counts peaks
# too many peaks -- remove duplicates
# since values are both time and position duplicates, take first instance
# of its occurence if sono values are the same
pcy2<-pcy2[!duplicated(yyy_T1$PC.sono[pcy2])]
# 5th, 7th, 12th and 13th peaks missing--add back manually
pcy2<-c(pcy2[1:4],2212,pcy2[5],3325,pcy2[6:9],6363,6984,pcy2[10:13])
# check it
pcy2
length(pcy2) # counts peaks
plot(yyy2, type="l")
points(py2, yyy2[py2], col="blue", pch=19) # plots peaks in blue onto strain profile plot
par(new=T)
plot(ypc2, type="l", col="red")
points(pcy2, ypc2[pcy2], col="black", pch=19) # plots peaks in blue onto strain profile plot

# PC length change cycle duration for Bird Yyy flight 2
HTpeaksy2 <- yyy_T1$Time[py2] 
PCpeaksy2 <- yyy_T1$Time[pcy2] 
PCcycDy2 <- diff(PCpeaksy2)

# calculate the timing offset between HT peak length (PL) and pect PL
PLoffsety2 <- HTpeaksy2 - PCpeaksy2 

# determine what that offset is as a percentage of the length change cycle
HT_pkoffsety2 <- (PLoffsety2[-(length(PLoffsety2))]/PCcycDy2)*100
avg_takeoffy2<-mean(HT_pkoffsety2[1:5])
avg_midflighty2<-mean(HT_pkoffsety2[6:((length(HT_pkoffsety2))-5)])
avg_landingy2<-mean(HT_pkoffsety2[((length(HT_pkoffsety2))-4):(length(HT_pkoffsety2))])
avg_HT_pkoffsety2 <-mean(HT_pkoffsety2)

### Flight 3 ###
yyy3=yyy_T2$HT.sono
plot(yyy3, type="l")

## Bird Yyy - flight 3 HT peaks
py3=find_peaks(yyy3, m = 200)
points(py3, yyy3[py3], col="blue", pch=19) # plots peaks in blue onto strain profile plot
py3
length(py3) # counts peaks
# remove first HT peak --  not true peak
py3<-py3[-1]
py3
length(py3)
# re-plot
plot(yyy3, type="l")
points(py3, yyy3[py3], col="blue", pch=19) # plots peaks in blue onto strain profile plot

## Bird Yyy - flight 3 PC peaks
ypc3=yyy_T2$PC.sono
par(new=T)
plot(ypc3, type="l", col="red")
## Bird Yyy - flight 3 PC peaks
pcy3=find_peaks(ypc3,m=200)
points(pcy3, ypc3[pcy3], col="black", pch=19) # plots peaks in blue onto strain profile plot
pcy3
length(pcy3) # counts peaks

# PC length change cycle duration for Bird Yyy flight 3
HTpeaksy3 <- yyy_T2$Time[py3] 
PCpeaksy3 <- yyy_T2$Time[pcy3] 
PCcycDy3 <- diff(PCpeaksy3)

# calculate the timing offset between HT peak length (PL) and pect PL
PLoffsety3 <- HTpeaksy3 - PCpeaksy3[-(length(PCpeaksy3))] 

# determine what that offset is as a percentage of the length change cycle
HT_pkoffsety3 <- (PLoffsety3/PCcycDy3)*100
avg_takeoffy3<-mean(HT_pkoffsety3[1:5])
avg_midflighty3<-mean(HT_pkoffsety3[6:((length(HT_pkoffsety3))-5)])
avg_landingy3<-mean(HT_pkoffsety3[((length(HT_pkoffsety3))-4):(length(HT_pkoffsety3))])
avg_HT_pkoffsety3 <-mean(HT_pkoffsety3)

### Flight 4 ###
yyy4=yyy_T3$HT.sono
plot(yyy4, type="l")

## Bird Yyy - flight 4 HT peaks
py4=find_peaks(yyy4, m = 200)
points(py4, yyy4[py4], col="blue", pch=19) # plots peaks in blue onto strain profile plot
py4
length(py4) # counts peaks
# remove first HT peak -- not true peak
py4<-py4[-1]
py4
length(py4)
# re-plot
plot(yyy4, type="l")
points(py4, yyy4[py4], col="blue", pch=19) # plots peaks in blue onto strain profile plot

## Bird Yyy - flight 4 PC peaks
ypc4=yyy_T3$PC.sono
par(new=T)
plot(ypc4, type="l", col="red")
## Bird Yyy - flight 4 PC peaks
pcy4=find_peaks(ypc4,m=200)
points(pcy4, ypc4[pcy4], col="black", pch=19) # plots peaks in blue onto strain profile plot
pcy4
length(pcy4) # counts peaks

# PC length change cycle duration for Bird Yyy flight 4
HTpeaksy4 <- yyy_T3$Time[py4] 
PCpeaksy4 <- yyy_T3$Time[pcy4] 
PCcycDy4 <- diff(PCpeaksy4)

# calculate the timing offset between HT peak length (PL) and pect PL
PLoffsety4 <- HTpeaksy4 - PCpeaksy4[-(length(PCpeaksy4))] 

# determine what that offset is as a percentage of the length change cycle
HT_pkoffsety4 <- (PLoffsety4/PCcycDy4)*100
avg_takeoffy4<-mean(HT_pkoffsety4[1:5])
avg_midflighty4<-mean(HT_pkoffsety4[6:((length(HT_pkoffsety4))-5)])
avg_landingy4<-mean(HT_pkoffsety4[((length(HT_pkoffsety4))-4):(length(HT_pkoffsety4))])
avg_HT_pkoffsety4 <-mean(HT_pkoffsety4)

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

## Bird Yyy - flight 5 PC peaks
ypc5=yyy_T4$PC.sono
par(new=T)
plot(ypc5, type="l", col="red")
## Bird Yyy - flight 5 PC peaks
pcy5=find_peaks(ypc5,m=200)
points(pcy5, ypc5[pcy5], col="black", pch=19) # plots peaks in blue onto strain profile plot
pcy5
length(pcy5) # counts peaks
# too many peaks -- remove duplicates
# since values are both time and position duplicates, take first instance
# of its occurence if sono values are the same
pcy5<-pcy5[!duplicated(yyy_T4$PC.sono[pcy5])]
# 9th and 12th peaks missing--add back manually
pcy5<-c(pcy5[1:8],4680,pcy5[9:10],6580,pcy5[11:17])
# check it
pcy5
length(pcy5) # counts peaks
plot(yyy5, type="l")
points(py5, yyy5[py5], col="blue", pch=19) # plots peaks in blue onto strain profile plot
par(new=T)
plot(ypc5, type="l", col="red")
points(pcy5, ypc5[pcy5], col="black", pch=19) # plots peaks in blue onto strain profile plot

# PC length change cycle duration for Bird Yyy flight 5
HTpeaksy5 <- yyy_T4$Time[py5] 
PCpeaksy5 <- yyy_T4$Time[pcy5] 
PCcycDy5 <- diff(PCpeaksy5)

# calculate the timing offset between HT peak length (PL) and pect PL
PLoffsety5 <- HTpeaksy5 - PCpeaksy5[-(length(PCpeaksy5))] 

# determine what that offset is as a percentage of the length change cycle
HT_pkoffsety5 <- (PLoffsety5/PCcycDy5)*100
avg_takeoffy5<-mean(HT_pkoffsety5[1:5])
avg_midflighty5<-mean(HT_pkoffsety5[6:((length(HT_pkoffsety5))-5)])
avg_landingy5<-mean(HT_pkoffsety5[((length(HT_pkoffsety5))-4):(length(HT_pkoffsety5))])
avg_HT_pkoffsety5 <-mean(HT_pkoffsety5)

#################################################################################
#####################   Write results to text file   ############################
#################################################################################


# create a matrix with as many rows as trials and populate it
# with the average lengthening and shortening percentages of per trial, per bird and add a 3rd 
# column with the ID of the individual that the strain data were taken from.
x <- matrix(nrow=26, ncol=5, byrow=F)
x[,1] <- c("Www1","Www2","Www3","Www4","Www5","Xxx1","Xxx2","Xxx3","Xxx4","Xxx5",
		"Uuu1","Uuu2","Uuu3","Uuu4","Uuu5","Tie1","Tie2","Tie3","Tie4","Tie5",
		"Yyy1","Yyy2","Yyy3","Yyy4","Yyy5", "Grand means")
x[,2] <- c(avg_HT_pkoffsetw1,avg_HT_pkoffsetw2,avg_HT_pkoffsetw3,avg_HT_pkoffsetw4,avg_HT_pkoffsetw5,
		avg_HT_pkoffsetx1,avg_HT_pkoffsetx2,avg_HT_pkoffsetx3,avg_HT_pkoffsetx4,avg_HT_pkoffsetx5,
		avg_HT_pkoffsetu1,avg_HT_pkoffsetu2,avg_HT_pkoffsetu3,avg_HT_pkoffsetu4,avg_HT_pkoffsetu5,
		avg_HT_pkoffsett1,avg_HT_pkoffsett2,avg_HT_pkoffsett3,avg_HT_pkoffsett4,avg_HT_pkoffsett5,
		avg_HT_pkoffsety1,avg_HT_pkoffsety2,avg_HT_pkoffsety3,avg_HT_pkoffsety4,avg_HT_pkoffsety5,
		mean(avg_HT_pkoffsetw1,avg_HT_pkoffsetw2,avg_HT_pkoffsetw3,avg_HT_pkoffsetw4,avg_HT_pkoffsetw5,
		avg_HT_pkoffsetx1,avg_HT_pkoffsetx2,avg_HT_pkoffsetx3,avg_HT_pkoffsetx4,avg_HT_pkoffsetx5,
		avg_HT_pkoffsetu1,avg_HT_pkoffsetu2,avg_HT_pkoffsetu3,avg_HT_pkoffsetu4,avg_HT_pkoffsetu5,
		avg_HT_pkoffsett1,avg_HT_pkoffsett2,avg_HT_pkoffsett3,avg_HT_pkoffsett4,avg_HT_pkoffsett5,
		avg_HT_pkoffsety1,avg_HT_pkoffsety2,avg_HT_pkoffsety3,avg_HT_pkoffsety4,avg_HT_pkoffsety5))

x[,3] <- c(avg_takeoffw1,avg_takeoffw2,avg_takeoffw3,avg_takeoffw4,avg_takeoffw5,
		avg_takeoffx1,avg_takeoffx2,avg_takeoffx3,avg_takeoffx4,avg_takeoffx5,
		avg_takeoffu1,avg_takeoffu2,avg_takeoffu3,avg_takeoffu4,avg_takeoffu5,
		avg_takeofft1,avg_takeofft2,avg_takeofft3,avg_takeofft4,avg_takeofft5,
		avg_takeoffy1,avg_takeoffy2,avg_takeoffy3,avg_takeoffy4,avg_takeoffy5,
		mean(avg_takeoffw1,avg_takeoffw2,avg_takeoffw3,avg_takeoffw4,avg_takeoffw5,
		avg_takeoffx1,avg_takeoffx2,avg_takeoffx3,avg_takeoffx4,avg_takeoffx5,
		avg_takeoffu1,avg_takeoffu2,avg_takeoffu3,avg_takeoffu4,avg_takeoffu5,
		avg_takeofft1,avg_takeofft2,avg_takeofft3,avg_takeofft4,avg_takeofft5,
		avg_takeoffy1,avg_takeoffy2,avg_takeoffy3,avg_takeoffy4,avg_takeoffy5))

x[,4] <- c(avg_midflightw1,avg_midflightw2,avg_midflightw3,avg_midflightw4,avg_midflightw5,
		avg_midflightx1,avg_midflightx2,avg_midflightx3,avg_midflightx4,avg_midflightx5,
		avg_midflightu1,avg_midflightu2,avg_midflightu3,avg_midflightu4,avg_midflightu5,
		avg_midflightt1,avg_midflightt2,avg_midflightt3,avg_midflightt4,avg_midflightt5,
		avg_midflighty1,avg_midflighty2,avg_midflighty3,avg_midflighty4,avg_midflighty5,
		mean(avg_midflightw1,avg_midflightw2,avg_midflightw3,avg_midflightw4,avg_midflightw5,
		avg_midflightx1,avg_midflightx2,avg_midflightx3,avg_midflightx4,avg_midflightx5,
		avg_midflightu1,avg_midflightu2,avg_midflightu3,avg_midflightu4,avg_midflightu5,
		avg_midflightt1,avg_midflightt2,avg_midflightt3,avg_midflightt4,avg_midflightt5,
		avg_midflighty1,avg_midflighty2,avg_midflighty3,avg_midflighty4,avg_midflighty5))

x[,5] <- c(avg_landingw1,avg_landingw2,avg_landingw3,avg_landingw4,avg_landingw5,
		avg_landingx1,avg_landingx2,avg_landingx3,avg_landingx4,avg_landingx5,
		avg_landingu1,avg_landingu2,avg_landingu3,avg_landingu4,avg_landingu5,
		avg_landingt1,avg_landingt2,avg_landingt3,avg_landingt4,avg_landingt5,
		avg_landingy1,avg_landingy2,avg_landingy3,avg_landingy4,avg_landingy5,
		mean(avg_landingw1,avg_landingw2,avg_landingw3,avg_landingw4,avg_landingw5,
		avg_landingx1,avg_landingx2,avg_landingx3,avg_landingx4,avg_landingx5,
		avg_landingu1,avg_landingu2,avg_landingu3,avg_landingu4,avg_landingu5,
		avg_landingt1,avg_landingt2,avg_landingt3,avg_landingt4,avg_landingt5,
		avg_landingy1,avg_landingy2,avg_landingy3,avg_landingy4,avg_landingy5))

# create a data frame from the above generated matrix (converts values to numeric).
# can also be exported as csv if need be.
avg.HT_PL_Offsets <- data.frame(x, stringsAsFactors = F)
# rename the columns to properly identify the metrics
colnames(avg.HT_PL_Offsets) <- c("Bird/trial ID", "Mean HT peak length offset from PC (% of PC cyc)", 
					"Mean takeoff offset (%)", "Mean midflight offset (%)", 
					"Mean landing offset (%)")

# move up one folder level (out of processed in vivo strain data folder) before saving text file
setwd(choose.dir())

## write .txt file of in vivo summary data for lengthening to shortening ratios 
write.table(avg.HT_PL_Offsets, "in vivo HT peak length offset data.txt", row.names=FALSE)