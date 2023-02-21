# Script to crop Robertson and Biewener strain data from HT and Pec for FFT, 
# onset phase analyses and lengthening to shortening phase ratios

# get and set working directory 
getwd()
## select folder containing all in vivo strain data for processing 
setwd(choose.dir())

# read in csvs for Bird Www

www_L1<-read.csv("Bird Www/WwwMMr15L_EMGfiltdSONOFiltd_N.csv")
www_T1<-read.csv("Bird Www/WwwMMr07T_EMGfiltdSONOfixedFiltd_N.csv")
www_T2<-read.csv("Bird Www/WwwMMr08T_EMGfiltdSONOfixedFiltd_N.csv")
www_T3<-read.csv("Bird Www/WwwMMr09T_EMGfiltdSONOfixedFiltd_N.csv")
www_T4<-read.csv("Bird Www/WwwMMr10T_EMGfiltdSONOfixedFiltd_N.csv")

# read in csvs for Bird Xxx

xxx_L1<-read.csv("Bird Xxx/XxxMMr36L_EMGfiltdSONOfiltd_N.csv")
xxx_T1<-read.csv("Bird Xxx/XxxMMr30T_EMGfiltdSONOfixedFiltd_N.csv")
xxx_T2<-read.csv("Bird Xxx/XxxMMr35T_EMGfiltdSONOfixedFiltd_N.csv")
xxx_T3<-read.csv("Bird Xxx/XxxMMr37T_EMGfiltdSONOfixedFiltd_N.csv")
xxx_T4<-read.csv("Bird Xxx/XxxMMr39T_EMGfiltdSONOfiltd_N.csv")

# read in csvs for Bird Uuu

uuu_T1<-read.csv("Bird Uuu/TuuuMM06axo_EMGfiltdSONOfixedFiltd_N.csv")
uuu_T2<-read.csv("Bird Uuu/TuuuMM07axo_EMGfiltdSONOfixedFiltd_N.csv")
uuu_L1<-read.csv("Bird Uuu/LuuuMM06axo_EMGfiltdSONOfixedFiltd_N.csv")
uuu_L2<-read.csv("Bird Uuu/LuuuMM07axo_EMGfiltdSONOfixedFiltd_N.csv")
uuu_L3<-read.csv("Bird Uuu/LuuuMM08axo_EMGfiltdSONOfixedFiltd_N.csv")

# read in csvs for Bird Tie

tie_T1<-read.csv("Bird Tie/TtieMM07axo_EMGfiltdSONOFiltd_N.csv")
tie_L1<-read.csv("Bird Tie/LtieMM04axo_EMGfiltdSONOfixedFiltd_N.csv")
tie_L2<-read.csv("Bird Tie/LtieMM07axo_EMGfiltdSONOfixedFiltd_N.csv")
tie_L3<-read.csv("Bird Tie/LtieMM10axo_EMGfiltdSONOFiltd_N.csv")
tie_L4<-read.csv("Bird Tie/LtieMM12axo_EMGfiltdSONOFiltd_N.csv")

# read in csvs for Bird Yyy

yyy_L1<-read.csv("Bird Yyy/YyyMMr13L_EMGfiltdSONOfixedFiltd_N.csv")
yyy_T1<-read.csv("Bird Yyy/YyyMMr08T_EMGfiltdSONOfixedFiltd_N.csv")
yyy_T2<-read.csv("Bird Yyy/YyyMMr10T_EMGfiltdSONOfixedFiltd_N.csv")
yyy_T3<-read.csv("Bird Yyy/YyyMMr16T_EMGfiltdSONOfixedFiltd_N.csv")
yyy_T4<-read.csv("Bird Yyy/YyyMMr18T_EMGfiltdSONOfixedFiltd_N.csv")

######################## Find points to subset with! #######################

# start all strains at first trough, 
# if first wingbeat is visibly very different than the rest then start at
# second trough, and end them at second last or last trough (visually inspect 
# last wingbeat for obvious differences as well)  
## INSPECT BOTH HT AND PC STRAINS!! ##

##### PLOTTING #####
### HT ###

# Bird Www Landing 
head(www_L1)
plot(www_L1$HT.sono~www_L1$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points for subsetting in vivo HT strain data
www_i1<-identify(www_L1$Time, www_L1$HT.sono, labels= seq_along(www_L1$HT.sono))
# subset to crop noise and first and last wingbeats from strain
c.www_L1<-subset(www_L1, as.numeric(rownames(www_L1)) >=5265 
	& as.numeric(rownames(www_L1))<=15390, select=c(Time, PC.sono, PC.EMG, 
	HT.sono, HT.EMG))
# check it
plot(c.www_L1$HT.sono~c.www_L1$Time, type="l")
par(new=T)
plot(c.www_L1$PC.sono~c.www_L1$Time, col="blue", type="l")

# Bird Www Take-off (x4) all on same panel
head(www_T1)
plot(www_T1$HT.sono~www_T1$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points for subsetting in vivo HT strain data
www_i2<-identify(www_T1$Time, www_T1$HT.sono, labels= seq_along(www_T1$HT.sono))
# subset to crop noise and first and last wingbeats from strain
c.www_T1<-subset(www_T1, as.numeric(rownames(www_T1)) >=2050 
	& as.numeric(rownames(www_T1))<=12360, select=c(Time, PC.sono, PC.EMG, 
	HT.sono, HT.EMG))
# check it
plot(c.www_T1$HT.sono~c.www_T1$Time, type="l")
par(new=T)
plot(c.www_T1$PC.sono~c.www_T1$Time, col="blue", type="l")

head(www_T2)
plot(www_T2$HT.sono~www_T2$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points for subsetting in vivo HT strain data
www_i3<-identify(www_T2$Time, www_T2$HT.sono, labels= seq_along(www_T2$HT.sono))
# missing first PC peak, therefore plot PC strain to get starting point
plot(www_T2$PC.sono~www_T2$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
www_in3<-identify(www_T2$Time, www_T2$PC.sono, labels= seq_along(www_T2$PC.sono))
# subset to crop noise and first and last wingbeats from strain
c.www_T2<-subset(www_T2, as.numeric(rownames(www_T2)) >=4900 
	& as.numeric(rownames(www_T2))<=14420, select=c(Time, PC.sono, PC.EMG, 
	HT.sono, HT.EMG))
# check it
plot(c.www_T2$HT.sono~c.www_T2$Time, type="l")
par(new=T)
plot(c.www_T2$PC.sono~c.www_T2$Time, col="blue", type="l")

head(www_T3)
plot(www_T3$HT.sono~www_T3$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points for subsetting in vivo HT strain data
www_i4<-identify(www_T3$Time, www_T3$HT.sono, labels= seq_along(www_T3$HT.sono))
# subset to crop noise and first and last wingbeats from strain
c.www_T3<-subset(www_T3, as.numeric(rownames(www_T3)) >=3100 
	& as.numeric(rownames(www_T3))<=13560, select=c(Time, PC.sono, PC.EMG, 
	HT.sono, HT.EMG))
# check it
plot(c.www_T3$HT.sono~c.www_T3$Time, type="l")
par(new=T)
plot(c.www_T3$PC.sono~c.www_T3$Time, col="blue", type="l")

head(www_T4)
plot(www_T4$HT.sono~www_T4$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points for subsetting in vivo HT strain data
www_i5<-identify(www_T4$Time, www_T4$HT.sono, labels= seq_along(www_T4$HT.sono))
# subset to crop noise and first and last wingbeats from strain
c.www_T4<-subset(www_T4, as.numeric(rownames(www_T4)) >=4110 
	& as.numeric(rownames(www_T4))<=13890, select=c(Time, PC.sono, PC.EMG, 
	HT.sono, HT.EMG))
# check it
plot(c.www_T4$HT.sono~c.www_T4$Time, type="l")
par(new=T)
plot(c.www_T4$PC.sono~c.www_T4$Time, col="blue", type="l")

# Bird Xxx Landing 
head(xxx_L1)
plot(xxx_L1$HT.sono~xxx_L1$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points for subsetting in vivo HT strain data
xxx_i1<-identify(xxx_L1$Time, xxx_L1$HT.sono, labels= seq_along(xxx_L1$HT.sono))
# missing first PC peak, plot PC strain to get starting point
plot(xxx_L1$PC.sono~xxx_L1$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
xxx_in1<-identify(xxx_L1$Time, xxx_L1$PC.sono, labels= seq_along(xxx_L1$PC.sono))
# subset to crop noise and first and last wingbeats from strain
c.xxx_L1<-subset(xxx_L1, as.numeric(rownames(xxx_L1)) >=3820 
	& as.numeric(rownames(xxx_L1))<=15200, select=c(Time, PC.sono, PC.EMG, 
	HT.sono, HT.EMG))
# check it
plot(c.xxx_L1$HT.sono~c.xxx_L1$Time, type="l")
par(new=T)
plot(c.xxx_L1$PC.sono~c.xxx_L1$Time, col="blue", type="l")

# Bird Xxx Take-off (x4) all on same panel
head(xxx_T1)
plot(xxx_T1$HT.sono~xxx_T1$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points for subsetting in vivo HT strain data
xxx_i2<-identify(xxx_T1$Time, xxx_T1$HT.sono, labels= seq_along(xxx_T1$HT.sono))
# subset to crop noise and first and last wingbeats from strain
c.xxx_T1<-subset(xxx_T1, as.numeric(rownames(xxx_T1)) >=5080 
	& as.numeric(rownames(xxx_T1))<=17500, select=c(Time, PC.sono, PC.EMG, 
	HT.sono, HT.EMG))
# check it
plot(c.xxx_T1$HT.sono~c.xxx_T1$Time, type="l")
par(new=T)
plot(c.xxx_T1$PC.sono~c.xxx_T1$Time, col="blue", type="l")

head(xxx_T2)
plot(xxx_T2$HT.sono~xxx_T2$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points for subsetting in vivo HT strain data
xxx_i3<-identify(xxx_T2$Time, xxx_T2$HT.sono, labels= seq_along(xxx_T2$HT.sono))
# subset to crop noise and first and last wingbeats from strain
c.xxx_T2<-subset(xxx_T2, as.numeric(rownames(xxx_T2)) >=3850 
	& as.numeric(rownames(xxx_T2))<=16300, select=c(Time, PC.sono, PC.EMG, 
	HT.sono, HT.EMG))
# check it
plot(c.xxx_T2$HT.sono~c.xxx_T2$Time, type="l")
par(new=T)
plot(c.xxx_T2$PC.sono~c.xxx_T2$Time, col="blue", type="l")

head(xxx_T3)
plot(xxx_T3$HT.sono~xxx_T3$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points for subsetting in vivo HT strain data
xxx_i4<-identify(xxx_T3$Time, xxx_T3$HT.sono, labels= seq_along(xxx_T3$HT.sono))
# subset to crop noise and first and last wingbeats from strain
c.xxx_T3<-subset(xxx_T3, as.numeric(rownames(xxx_T3)) >=7520 
	& as.numeric(rownames(xxx_T3))<=20140, select=c(Time, PC.sono, PC.EMG, 
	HT.sono, HT.EMG))
# check it
plot(c.xxx_T3$HT.sono~c.xxx_T3$Time, type="l")
par(new=T)
plot(c.xxx_T3$PC.sono~c.xxx_T3$Time, col="blue", type="l")

head(xxx_T4)
plot(xxx_T4$HT.sono~xxx_T4$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points for subsetting in vivo HT strain data
xxx_i5<-identify(xxx_T4$Time, xxx_T4$HT.sono, labels= seq_along(xxx_T4$HT.sono))
# subset to crop noise and first and last wingbeats from strain
c.xxx_T4<-subset(xxx_T4, as.numeric(rownames(xxx_T4)) >=3520 
	& as.numeric(rownames(xxx_T4))<=16050, select=c(Time, PC.sono, PC.EMG, 
	HT.sono, HT.EMG))
# check it
plot(c.xxx_T4$HT.sono~c.xxx_T4$Time, type="l")
par(new=T)
plot(c.xxx_T4$PC.sono~c.xxx_T4$Time, col="blue", type="l")

# Bird Uuu Take-off (x2)
head(uuu_T1)
plot(uuu_T1$HT.sono~uuu_T1$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points for subsetting in vivo HT strain data
uuu_i1<-identify(uuu_T1$Time, uuu_T1$HT.sono, labels= seq_along(uuu_T1$HT.sono))
# subset to crop noise and first and last wingbeats from strain
c.uuu_T1<-subset(uuu_T1, as.numeric(rownames(uuu_T1)) >=2560 
	& as.numeric(rownames(uuu_T1))<=12050, select=c(Time, PC.sono, PC.EMG, 
	HT.sono, HT.EMG))
# check it
plot(c.uuu_T1$HT.sono~c.uuu_T1$Time, type="l")
par(new=T)
plot(c.uuu_T1$PC.sono~c.uuu_T1$Time, col="blue", type="l")

head(uuu_T2)
plot(uuu_T2$HT.sono~uuu_T2$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points for subsetting in vivo HT strain data
uuu_i2<-identify(uuu_T2$Time, uuu_T2$HT.sono, labels= seq_along(uuu_T2$HT.sono))
# subset to crop noise and first and last wingbeats from strain
c.uuu_T2<-subset(uuu_T2, as.numeric(rownames(uuu_T2)) >=2720 
	& as.numeric(rownames(uuu_T2))<=12795, select=c(Time, PC.sono, PC.EMG, 
	HT.sono, HT.EMG))
# check it
plot(c.uuu_T2$HT.sono~c.uuu_T2$Time, type="l")
par(new=T)
plot(c.uuu_T2$PC.sono~c.uuu_T2$Time, col="blue", type="l")

# Bird Uuu Landing (x3) all on same panel
head(uuu_L1)
plot(uuu_L1$HT.sono~uuu_L1$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points for subsetting in vivo HT strain data
uuu_i3<-identify(uuu_L1$Time, uuu_L1$HT.sono, labels= seq_along(uuu_L1$HT.sono))
# missing first PC peak, plot PC strain to find starting point
plot(uuu_L1$PC.sono~uuu_L1$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
uuu_in3<-identify(uuu_L1$Time, uuu_L1$PC.sono, labels= seq_along(uuu_L1$PC.sono))
# skip first PC peak because includes too much noise from start of HT strain
# subset to crop noise and first and last wingbeats from strain
c.uuu_L1<-subset(uuu_L1, as.numeric(rownames(uuu_L1)) >=5330 
	& as.numeric(rownames(uuu_L1))<=14720, select=c(Time, PC.sono, PC.EMG, 
	HT.sono, HT.EMG))
# check it
plot(c.uuu_L1$HT.sono~c.uuu_L1$Time, type="l")
par(new=T)
plot(c.uuu_L1$PC.sono~c.uuu_L1$Time, col="blue", type="l")

head(uuu_L2)
plot(uuu_L2$HT.sono~uuu_L2$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points for subsetting in vivo HT strain data
uuu_i4<-identify(uuu_L2$Time, uuu_L2$HT.sono, labels= seq_along(uuu_L2$HT.sono))
# subset to crop noise and first and last wingbeats from strain
c.uuu_L2<-subset(uuu_L2, as.numeric(rownames(uuu_L2)) >=4130 
	& as.numeric(rownames(uuu_L2))<=13650, select=c(Time, PC.sono, PC.EMG, 
	HT.sono, HT.EMG))
# check it
plot(c.uuu_L2$HT.sono~c.uuu_L2$Time, type="l")
par(new=T)
plot(c.uuu_L2$PC.sono~c.uuu_L2$Time, col="blue", type="l")

head(uuu_L3)
plot(uuu_L3$HT.sono~uuu_L3$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points for subsetting in vivo HT strain data
uuu_i5<-identify(uuu_L3$Time, uuu_L3$HT.sono, labels= seq_along(uuu_L3$HT.sono))

# subset to crop noise and first and last wingbeats from strain
## skip first peaks -- includes too much noise from HT strain start
c.uuu_L3<-subset(uuu_L3, as.numeric(rownames(uuu_L3)) >=2520 
	& as.numeric(rownames(uuu_L3))<=11625, select=c(Time, PC.sono, PC.EMG, 
	HT.sono, HT.EMG))
# check it
plot(c.uuu_L3$HT.sono~c.uuu_L3$Time, type="l")
par(new=T)
plot(c.uuu_L3$PC.sono~c.uuu_L3$Time, col="blue", type="l")

# Bird Tie Take-off 
head(tie_T1)
plot(tie_T1$HT.sono~tie_T1$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points for subsetting in vivo HT strain data
tie_i1<-identify(tie_T1$Time, tie_T1$HT.sono, labels= seq_along(tie_T1$HT.sono))
# small offset b/w PC and HT--check full PC strain
plot(tie_T1$PC.sono~tie_T1$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
tie_in1<-identify(tie_T1$Time, tie_T1$PC.sono, labels= seq_along(tie_T1$PC.sono))
# subset to crop noise and first and last wingbeats from strain
c.tie_T1<-subset(tie_T1, as.numeric(rownames(tie_T1)) >=3950 
	& as.numeric(rownames(tie_T1))<=11240, select=c(Time, PC.sono, PC.EMG, 
	HT.sono, HT.EMG))
# check it
plot(c.tie_T1$HT.sono~c.tie_T1$Time, type="l")
par(new=T)
plot(c.tie_T1$PC.sono~c.tie_T1$Time, col="blue", type="l")

# Bird Tie Landing (x4) all on same panel
head(tie_L1)
plot(tie_L1$HT.sono~tie_L1$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points for subsetting in vivo HT strain data
tie_i2<-identify(tie_L1$Time, tie_L1$HT.sono, labels= seq_along(tie_L1$HT.sono))
# small offset b/w PC and HT--check full PC strain
plot(tie_L1$PC.sono~tie_L1$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
tie_in2<-identify(tie_L1$Time, tie_L1$PC.sono, labels= seq_along(tie_L1$PC.sono))
# subset to crop noise and first and last wingbeats from strain
c.tie_L1<-subset(tie_L1, as.numeric(rownames(tie_L1)) >=6840 
	& as.numeric(rownames(tie_L1))<=15820, select=c(Time, PC.sono, PC.EMG, 
	HT.sono, HT.EMG))
# check it
plot(c.tie_L1$HT.sono~c.tie_L1$Time, type="l")
par(new=T)
plot(c.tie_L1$PC.sono~c.tie_L1$Time, col="blue", type="l")

head(tie_L2)
plot(tie_L2$HT.sono~tie_L2$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points for subsetting in vivo HT strain data
tie_i3<-identify(tie_L2$Time, tie_L2$HT.sono, labels= seq_along(tie_L2$HT.sono))
# subset to crop noise and first and last wingbeats from strain
c.tie_L2<-subset(tie_L2, as.numeric(rownames(tie_L2)) >=8030 
	& as.numeric(rownames(tie_L2))<=16910, select=c(Time, PC.sono, PC.EMG, 
	HT.sono, HT.EMG))
# check it
plot(c.tie_L2$HT.sono~c.tie_L2$Time, type="l")
par(new=T)
plot(c.tie_L2$PC.sono~c.tie_L2$Time, col="blue", type="l")

head(tie_L3)
plot(tie_L3$HT.sono~tie_L3$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points for subsetting in vivo HT strain data
tie_i4<-identify(tie_L3$Time, tie_L3$HT.sono, labels= seq_along(tie_L3$HT.sono))
# subset to crop noise and first and last wingbeats from strain
c.tie_L3<-subset(tie_L3, as.numeric(rownames(tie_L3)) >=9760 
	& as.numeric(rownames(tie_L3))<=17230, select=c(Time, PC.sono, PC.EMG, 
	HT.sono, HT.EMG))
# check it
plot(c.tie_L3$HT.sono~c.tie_L3$Time, type="l")
par(new=T)
plot(c.tie_L3$PC.sono~c.tie_L3$Time, col="blue", type="l")

head(tie_L4)
plot(tie_L4$HT.sono~tie_L4$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points for subsetting in vivo HT strain data
tie_i5<-identify(tie_L4$Time, tie_L4$HT.sono, labels= seq_along(tie_L4$HT.sono))
# subset to crop noise and first and last wingbeats from strain
c.tie_L4<-subset(tie_L4, as.numeric(rownames(tie_L4)) >=5370 
	& as.numeric(rownames(tie_L4))<=13480, select=c(Time, PC.sono, PC.EMG, 
	HT.sono, HT.EMG))
# check it
plot(c.tie_L4$HT.sono~c.tie_L4$Time, type="l")
par(new=T)
plot(c.tie_L4$PC.sono~c.tie_L4$Time, col="blue", type="l")

# Bird Yyy Landing 
head(yyy_L1)
plot(yyy_L1$HT.sono~yyy_L1$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points for subsetting in vivo HT strain data
yyy_i1<-identify(yyy_L1$Time, yyy_L1$HT.sono, labels= seq_along(yyy_L1$HT.sono))
# subset to crop noise and first and last wingbeats from strain
c.yyy_L1<-subset(yyy_L1, as.numeric(rownames(yyy_L1)) >=4110 
	& as.numeric(rownames(yyy_L1))<=14950, select=c(Time, PC.sono, PC.EMG, 
	HT.sono, HT.EMG))
# check it
plot(c.yyy_L1$HT.sono~c.yyy_L1$Time, type="l")
par(new=T)
plot(c.yyy_L1$PC.sono~c.yyy_L1$Time, col="blue", type="l")

# Bird Yyy Take-off (x4) all on same panel
head(yyy_T1)
plot(yyy_T1$HT.sono~yyy_T1$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points for subsetting in vivo HT strain data
yyy_i2<-identify(yyy_T1$Time, yyy_T1$HT.sono, labels= seq_along(yyy_T1$HT.sono))
# subset to crop noise and first and last wingbeats from strain
c.yyy_T1<-subset(yyy_T1, as.numeric(rownames(yyy_T1)) >=5620 
	& as.numeric(rownames(yyy_T1))<=15590, select=c(Time, PC.sono, PC.EMG, 
	HT.sono, HT.EMG))
# check it
plot(c.yyy_T1$HT.sono~c.yyy_T1$Time, type="l")
par(new=T)
plot(c.yyy_T1$PC.sono~c.yyy_T1$Time, col="blue", type="l")

head(yyy_T2)
plot(yyy_T2$HT.sono~yyy_T2$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points for subsetting in vivo HT strain data
yyy_i3<-identify(yyy_T2$Time, yyy_T2$HT.sono, labels= seq_along(yyy_T2$HT.sono))
# subset to crop noise and first and last wingbeats from strain
c.yyy_T2<-subset(yyy_T2, as.numeric(rownames(yyy_T2)) >=4510 
	& as.numeric(rownames(yyy_T2))<=15760, select=c(Time, PC.sono, PC.EMG, 
	HT.sono, HT.EMG))
# check it
plot(c.yyy_T2$HT.sono~c.yyy_T2$Time, type="l")
par(new=T)
plot(c.yyy_T2$PC.sono~c.yyy_T2$Time, col="blue", type="l")

head(yyy_T3)
plot(yyy_T3$HT.sono~yyy_T3$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points for subsetting in vivo HT strain data
yyy_i4<-identify(yyy_T3$Time, yyy_T3$HT.sono, labels= seq_along(yyy_T3$HT.sono))
# subset to crop noise and first and last wingbeats from strain
c.yyy_T3<-subset(yyy_T3, as.numeric(rownames(yyy_T3)) >=3290 
	& as.numeric(rownames(yyy_T3))<=16050, select=c(Time, PC.sono, PC.EMG, 
	HT.sono, HT.EMG))
# check it
plot(c.yyy_T3$HT.sono~c.yyy_T3$Time, type="l")
par(new=T)
plot(c.yyy_T3$PC.sono~c.yyy_T3$Time, col="blue", type="l")

head(yyy_T4)
plot(yyy_T4$HT.sono~yyy_T4$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points for subsetting in vivo HT strain data
yyy_i5<-identify(yyy_T4$Time, yyy_T4$HT.sono, labels= seq_along(yyy_T4$HT.sono))
# subset to crop noise and first and last wingbeats from strain
c.yyy_T4<-subset(yyy_T4, as.numeric(rownames(yyy_T4)) >=3800 
	& as.numeric(rownames(yyy_T4))<=14805, select=c(Time, PC.sono, PC.EMG, 
	HT.sono, HT.EMG))
# check it
plot(c.yyy_T4$HT.sono~c.yyy_T4$Time, type="l")
par(new=T)
plot(c.yyy_T4$PC.sono~c.yyy_T4$Time, col="blue", type="l")

# new subfolder to write new csvs that contain fully processed data and are ready for analysis is written
# into filepath names ('in vivo strain data for analysis')

#### Write csvs for cropped strains ####
write.csv(c.www_L1, "in vivo strain data for analysis/BirdWww1 in vivo_HT strain.csv", row.names=FALSE)
write.csv(c.www_T1, "in vivo strain data for analysis/BirdWww2 in vivo_HT strain.csv", row.names=FALSE)
write.csv(c.www_T2, "in vivo strain data for analysis/BirdWww3 in vivo_HT strain.csv", row.names=FALSE)
write.csv(c.www_T3, "in vivo strain data for analysis/BirdWww4 in vivo_HT strain.csv", row.names=FALSE)
write.csv(c.www_T4, "in vivo strain data for analysis/BirdWww5 in vivo_HT strain.csv", row.names=FALSE)

write.csv(c.xxx_L1, "in vivo strain data for analysis/BirdXxx1 in vivo_HT strain.csv", row.names=FALSE)
write.csv(c.xxx_T1, "in vivo strain data for analysis/BirdXxx2 in vivo_HT strain.csv", row.names=FALSE)
write.csv(c.xxx_T2, "in vivo strain data for analysis/BirdXxx3 in vivo_HT strain.csv", row.names=FALSE)
write.csv(c.xxx_T3, "in vivo strain data for analysis/BirdXxx4 in vivo_HT strain.csv", row.names=FALSE)
write.csv(c.xxx_T4, "in vivo strain data for analysis/BirdXxx5 in vivo_HT strain.csv", row.names=FALSE)

write.csv(c.uuu_T1, "in vivo strain data for analysis/BirdUuu1 in vivo_HT strain.csv", row.names=FALSE)
write.csv(c.uuu_T2, "in vivo strain data for analysis/BirdUuu2 in vivo_HT strain.csv", row.names=FALSE)
write.csv(c.uuu_L1, "in vivo strain data for analysis/BirdUuu3 in vivo_HT strain.csv", row.names=FALSE)
write.csv(c.uuu_L2, "in vivo strain data for analysis/BirdUuu4 in vivo_HT strain.csv", row.names=FALSE)
write.csv(c.uuu_L3, "in vivo strain data for analysis/BirdUuu5 in vivo_HT strain.csv", row.names=FALSE)

write.csv(c.tie_T1, "in vivo strain data for analysis/BirdTie1 in vivo_HT strain.csv", row.names=FALSE)
write.csv(c.tie_L1, "in vivo strain data for analysis/BirdTie2 in vivo_HT strain.csv", row.names=FALSE)
write.csv(c.tie_L2, "in vivo strain data for analysis/BirdTie3 in vivo_HT strain.csv", row.names=FALSE)
write.csv(c.tie_L3, "in vivo strain data for analysis/BirdTie4 in vivo_HT strain.csv", row.names=FALSE)
write.csv(c.tie_L4, "in vivo strain data for analysis/BirdTie5 in vivo_HT strain.csv", row.names=FALSE)

write.csv(c.yyy_L1, "in vivo strain data for analysis/BirdYyy1 in vivo_HT strain.csv", row.names=FALSE)
write.csv(c.yyy_T1, "in vivo strain data for analysis/BirdYyy2 in vivo_HT strain.csv", row.names=FALSE)
write.csv(c.yyy_T2, "in vivo strain data for analysis/BirdYyy3 in vivo_HT strain.csv", row.names=FALSE)
write.csv(c.yyy_T3, "in vivo strain data for analysis/BirdYyy4 in vivo_HT strain.csv", row.names=FALSE)
write.csv(c.yyy_T4, "in vivo strain data for analysis/BirdYyy5 in vivo_HT strain.csv", row.names=FALSE)

