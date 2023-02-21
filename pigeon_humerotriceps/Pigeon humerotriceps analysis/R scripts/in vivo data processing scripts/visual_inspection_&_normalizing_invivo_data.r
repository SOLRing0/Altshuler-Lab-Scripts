# Script to visually inspect and normalize distances of Robertson and Biewener strain data 
# for HT and Pec to resting length 
## NOTE: All ranges for each trial's resting length was chosen using visual assessment and 
## the identify function 

# get and set working directory 
getwd()
## select folder containing all in vivo strain data for processing 
setwd(choose.dir())

# read in csvs for Bird Www

www_L1<-read.csv("Bird Www/WwwMMr15L_EMGfiltdSONOFiltd.csv")
www_T1<-read.csv("Bird Www/WwwMMr07T_EMGfiltdSONOfixedFiltd.csv")
www_T2<-read.csv("Bird Www/WwwMMr08T_EMGfiltdSONOfixedFiltd.csv")
www_T3<-read.csv("Bird Www/WwwMMr09T_EMGfiltdSONOfixedFiltd.csv")
www_T4<-read.csv("Bird Www/WwwMMr10T_EMGfiltdSONOfixedFiltd.csv")

# read in csvs for Bird Xxx

xxx_L1<-read.csv("Bird Xxx/XxxMMr36L_EMGfiltdSONOfiltd.csv")
xxx_T1<-read.csv("Bird Xxx/XxxMMr30T_EMGfiltdSONOfixedFiltd.csv")
xxx_T2<-read.csv("Bird Xxx/XxxMMr35T_EMGfiltdSONOfixedFiltd.csv")
xxx_T3<-read.csv("Bird Xxx/XxxMMr37T_EMGfiltdSONOfixedFiltd.csv")
xxx_T4<-read.csv("Bird Xxx/XxxMMr39T_EMGfiltdSONOfiltd.csv")

# read in csvs for Bird Uuu

uuu_T1<-read.csv("Bird Uuu/TuuuMM06axo_EMGfiltdSONOfixedFiltd.csv")
uuu_T2<-read.csv("Bird Uuu/TuuuMM07axo_EMGfiltdSONOfixedFiltd.csv")
uuu_L1<-read.csv("Bird Uuu/LuuuMM06axo_EMGfiltdSONOfixedFiltd.csv")
uuu_L2<-read.csv("Bird Uuu/LuuuMM07axo_EMGfiltdSONOfixedFiltd.csv")
uuu_L3<-read.csv("Bird Uuu/LuuuMM08axo_EMGfiltdSONOfixedFiltd.csv")

# read in csvs for Bird Tie

tie_T1<-read.csv("Bird Tie/TtieMM07axo_EMGfiltdSONOFiltd.csv")
tie_L1<-read.csv("Bird Tie/LtieMM04axo_EMGfiltdSONOfixedFiltd.csv")
tie_L2<-read.csv("Bird Tie/LtieMM07axo_EMGfiltdSONOfixedFiltd.csv")
tie_L3<-read.csv("Bird Tie/LtieMM10axo_EMGfiltdSONOFiltd.csv")
tie_L4<-read.csv("Bird Tie/LtieMM12axo_EMGfiltdSONOFiltd.csv")

# read in csvs for Bird Yyy

yyy_L1<-read.csv("Bird Yyy/YyyMMr13L_EMGfiltdSONOfixedFiltd.csv")
yyy_T1<-read.csv("Bird Yyy/YyyMMr08T_EMGfiltdSONOfixedFiltd.csv")
yyy_T2<-read.csv("Bird Yyy/YyyMMr10T_EMGfiltdSONOfixedFiltd.csv")
yyy_T3<-read.csv("Bird Yyy/YyyMMr16T_EMGfiltdSONOfixedFiltd.csv")
yyy_T4<-read.csv("Bird Yyy/YyyMMr18T_EMGfiltdSONOfixedFiltd.csv")

# shift timing to start at zero for all strains to facilitate cropping

# Bird Www
www_L1$Time<-www_L1$Time-(min(www_L1$Time))
www_T1$Time<-www_T1$Time-(min(www_T1$Time))
www_T2$Time<-www_T2$Time-(min(www_T2$Time))
www_T3$Time<-www_T3$Time-(min(www_T3$Time))
www_T4$Time<-www_T4$Time-(min(www_T4$Time))

# Bird Xxx
xxx_L1$Time<-xxx_L1$Time-(min(xxx_L1$Time))
xxx_T1$Time<-xxx_T1$Time-(min(xxx_T1$Time))
xxx_T2$Time<-xxx_T2$Time-(min(xxx_T2$Time))
xxx_T3$Time<-xxx_T3$Time-(min(xxx_T3$Time))
xxx_T4$Time<-xxx_T4$Time-(min(xxx_T4$Time))

# Bird Uuu
uuu_T1$Time<-uuu_T1$Time-(min(uuu_T1$Time))
uuu_T2$Time<-uuu_T2$Time-(min(uuu_T2$Time))
uuu_L1$Time<-uuu_L1$Time-(min(uuu_L1$Time))
uuu_L2$Time<-uuu_L2$Time-(min(uuu_L2$Time))
uuu_L3$Time<-uuu_L3$Time-(min(uuu_L3$Time))


# Bird Tie
tie_T1$Time<-tie_T1$Time-(min(tie_T1$Time))
tie_L1$Time<-tie_L1$Time-(min(tie_L1$Time))
tie_L2$Time<-tie_L2$Time-(min(tie_L2$Time))
tie_L3$Time<-tie_L3$Time-(min(tie_L3$Time))
tie_L4$Time<-tie_L4$Time-(min(tie_L4$Time))

# Bird Yyy
yyy_L1$Time<-yyy_L1$Time-(min(yyy_L1$Time))
yyy_T1$Time<-yyy_T1$Time-(min(yyy_T1$Time))
yyy_T2$Time<-yyy_T2$Time-(min(yyy_T2$Time))
yyy_T3$Time<-yyy_T3$Time-(min(yyy_T3$Time))
yyy_T4$Time<-yyy_T4$Time-(min(yyy_T4$Time))

##### PLOTTING #####
### HT ###

# Bird Www Landing 
head(www_L1)
plot(www_L1$HT.sono~www_L1$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points to getting resting length for normalizing distances
www_i1<-identify(www_L1$Time, www_L1$HT.sono, labels= seq_along(www_L1$HT.sono))
# normalize strain values in each file using their 'resting' distances
www_L1$HT.sono<-www_L1$HT.sono-(mean(www_L1$HT.sono[c(70:2115,2250:4230)]))

# Bird Www Take-off (x4) all on same panel
head(www_T1)
plot(www_T1$HT.sono~www_T1$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points to getting resting length for normalizing distances
www_i2<-identify(www_T1$Time, www_T1$HT.sono, labels= seq_along(www_T1$HT.sono))
# normalize strain values in each file using their 'resting' distances
www_T1$HT.sono<-www_T1$HT.sono-(mean(www_T1$HT.sono[1:875]))

head(www_T2)
plot(www_T2$HT.sono~www_T2$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points to getting resting length for normalizing distances
www_i3<-identify(www_T2$Time, www_T2$HT.sono, labels= seq_along(www_T2$HT.sono))
# normalize strain values in each file using their 'resting' distances
www_T2$HT.sono<-www_T2$HT.sono-(mean(www_T2$HT.sono[1:2035]))

head(www_T3)
plot(www_T3$HT.sono~www_T3$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points to getting resting length for normalizing distances
www_i4<-identify(www_T3$Time, www_T3$HT.sono, labels= seq_along(www_T3$HT.sono))
# normalize strain values in each file using their 'resting' distances
www_T3$HT.sono<-www_T3$HT.sono-(mean(www_T3$HT.sono[c(1:290,1130:2320)]))

head(www_T4)
plot(www_T4$HT.sono~www_T4$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points to getting resting length for normalizing distances
www_i5<-identify(www_T4$Time, www_T4$HT.sono, labels= seq_along(www_T4$HT.sono))
# normalize strain values in each file using their 'resting' distances
www_T4$HT.sono<-www_T4$HT.sono-(mean(www_T4$HT.sono[c(1:770,960:1570)]))

# Bird Xxx Landing 
head(xxx_L1)
plot(xxx_L1$HT.sono~xxx_L1$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points to getting resting length for normalizing distances
xxx_i1<-identify(xxx_L1$Time, xxx_L1$HT.sono, labels= seq_along(xxx_L1$HT.sono))
# normalize strain values in each file using their 'resting' distances
xxx_L1$HT.sono<-xxx_L1$HT.sono-(mean(xxx_L1$HT.sono[1:3054]))

# Bird Xxx Take-off (x4) all on same panel
head(xxx_T1)
plot(xxx_T1$HT.sono~xxx_T1$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points to getting resting length for normalizing distances
xxx_i2<-identify(xxx_T1$Time, xxx_T1$HT.sono, labels= seq_along(xxx_T1$HT.sono))
# normalize strain values in each file using their 'resting' distances
xxx_T1$HT.sono<-xxx_T1$HT.sono-(mean(xxx_T1$HT.sono[1:4335]))

head(xxx_T2)
plot(xxx_T2$HT.sono~xxx_T2$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points to getting resting length for normalizing distances
xxx_i3<-identify(xxx_T2$Time, xxx_T2$HT.sono, labels= seq_along(xxx_T2$HT.sono))
# normalize strain values in each file using their 'resting' distances
xxx_T2$HT.sono<-xxx_T2$HT.sono-(mean(xxx_T2$HT.sono[1:3225]))

head(xxx_T3)
plot(xxx_T3$HT.sono~xxx_T3$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points to getting resting length for normalizing distances
xxx_i4<-identify(xxx_T3$Time, xxx_T3$HT.sono, labels= seq_along(xxx_T3$HT.sono))
# normalize strain values in each file using their 'resting' distances
xxx_T3$HT.sono<-xxx_T3$HT.sono-(mean(xxx_T3$HT.sono[1:6895]))

head(xxx_T4)
plot(xxx_T4$HT.sono~xxx_T4$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points to getting resting length for normalizing distances
xxx_i5<-identify(xxx_T4$Time, xxx_T4$HT.sono, labels= seq_along(xxx_T4$HT.sono))
# normalize strain values in each file using their 'resting' distances
xxx_T4$HT.sono<-xxx_T4$HT.sono-(mean(xxx_T4$HT.sono[1:3230]))

# Bird Uuu Take-off (x2)
head(uuu_T1)
plot(uuu_T1$HT.sono~uuu_T1$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points to getting resting length for normalizing distances
uuu_i1<-identify(uuu_T1$Time, uuu_T1$HT.sono, labels= seq_along(uuu_T1$HT.sono))
# normalize strain values in each file using their 'resting' distances
uuu_T1$HT.sono<-uuu_T1$HT.sono-(mean(uuu_T1$HT.sono[c(30:135,350:580,640:730,805:1170)]))

head(uuu_T2)
plot(uuu_T2$HT.sono~uuu_T2$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points to getting resting length for normalizing distances
uuu_i2<-identify(uuu_T2$Time, uuu_T2$HT.sono, labels= seq_along(uuu_T2$HT.sono))
# normalize strain values in each file using their 'resting' distances
uuu_T2$HT.sono<-uuu_T2$HT.sono-(mean(uuu_T2$HT.sono[c(40:265,330:475,584:900,990:1110,1360:1780)]))

# Bird Uuu Landing (x3) all on same panel
head(uuu_L1)
plot(uuu_L1$HT.sono~uuu_L1$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points to getting resting length for normalizing distances
uuu_i3<-identify(uuu_L1$Time, uuu_L1$HT.sono, labels= seq_along(uuu_L1$HT.sono))
# normalize strain values in each file using their 'resting' distances
uuu_L1$HT.sono<-uuu_L1$HT.sono-(mean(uuu_L1$HT.sono[c(55:200,530:690,
	1105:1355,1830:2080,2440:2660,2810:3015,3205:3450)]))

head(uuu_L2)
plot(uuu_L2$HT.sono~uuu_L2$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points to getting resting length for normalizing distances
uuu_i4<-identify(uuu_L2$Time, uuu_L2$HT.sono, labels= seq_along(uuu_L2$HT.sono))
# normalize strain values in each file using their 'resting' distances
uuu_L2$HT.sono<-uuu_L2$HT.sono-(mean(uuu_L2$HT.sono[c(30:110,245:370,
	450:550,640:705,885:980,1090:1165)]))

head(uuu_L3)
plot(uuu_L3$HT.sono~uuu_L3$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points to getting resting length for normalizing distances
uuu_i5<-identify(uuu_L3$Time, uuu_L3$HT.sono, labels= seq_along(uuu_L3$HT.sono))
# normalize strain values in each file using their 'resting' distances
uuu_L3$HT.sono<-uuu_L3$HT.sono-(mean(uuu_L3$HT.sono[c(20:75,135:240,
	355:415,555:660,725:760,835:865,920:1010,1080:1195)]))

# Bird Tie Take-off 
head(tie_T1)
plot(tie_T1$HT.sono~tie_T1$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
tie_i1<-identify(tie_T1$Time, tie_T1$HT.sono, labels= seq_along(tie_T1$HT.sono))
# normalize strain values in each file using their 'resting' distances
tie_T1$HT.sono<-tie_T1$HT.sono-(mean(tie_T1$HT.sono[25:2450]))

# Bird Tie Landing (x4) all on same panel
head(tie_L1)
plot(tie_L1$HT.sono~tie_L1$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
tie_i2<-identify(tie_L1$Time, tie_L1$HT.sono, labels= seq_along(tie_L1$HT.sono))
# normalize strain values in each file using their 'resting' distances
tie_L1$HT.sono<-tie_L1$HT.sono-(mean(tie_L1$HT.sono[55:5530]))

head(tie_L2)
plot(tie_L2$HT.sono~tie_L2$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
tie_i3<-identify(tie_L2$Time, tie_L2$HT.sono, labels= seq_along(tie_L2$HT.sono))
# normalize strain values in each file using their 'resting' distances
tie_L2$HT.sono<-tie_L2$HT.sono-(mean(tie_L2$HT.sono[55:6570]))

head(tie_L3)
plot(tie_L3$HT.sono~tie_L3$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
tie_i4<-identify(tie_L3$Time, tie_L3$HT.sono, labels= seq_along(tie_L3$HT.sono))
# normalize strain values in each file using their 'resting' distances
tie_L3$HT.sono<-tie_L3$HT.sono-(mean(tie_L3$HT.sono[110:8075]))

head(tie_L4)
plot(tie_L4$HT.sono~tie_L4$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
tie_i5<-identify(tie_L4$Time, tie_L4$HT.sono, labels= seq_along(tie_L4$HT.sono))
# normalize strain values in each file using their 'resting' distances
tie_L4$HT.sono<-tie_L4$HT.sono-(mean(tie_L4$HT.sono[85:4060]))

# Bird Yyy Landing 
head(yyy_L1)
plot(yyy_L1$HT.sono~yyy_L1$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points to getting resting length for normalizing distances
yyy_i1<-identify(yyy_L1$Time, yyy_L1$HT.sono, labels= seq_along(yyy_L1$HT.sono))
# normalize strain values in each file using their 'resting' distances
yyy_L1$HT.sono<-yyy_L1$HT.sono-(mean(yyy_L1$HT.sono[c(1:600,695:1100,
	1275:1755,1910:2345,2470:3220)]))

# Bird Yyy Take-off (x4) all on same panel
head(yyy_T1)
plot(yyy_T1$HT.sono~yyy_T1$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points to getting resting length for normalizing distances
yyy_i2<-identify(yyy_T1$Time, yyy_T1$HT.sono, labels= seq_along(yyy_T1$HT.sono))
# normalize strain values in each file using their 'resting' distances
yyy_T1$HT.sono<-yyy_T1$HT.sono-(mean(yyy_T1$HT.sono[c(1:1875,2015:2420,
	2590:2990,3115:4710)]))

head(yyy_T2)
plot(yyy_T2$HT.sono~yyy_T2$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points to getting resting length for normalizing distances
yyy_i3<-identify(yyy_T2$Time, yyy_T2$HT.sono, labels= seq_along(yyy_T2$HT.sono))
# normalize strain values in each file using their 'resting' distances
yyy_T2$HT.sono<-yyy_T2$HT.sono-(mean(yyy_T2$HT.sono[50:3165]))

head(yyy_T3)
plot(yyy_T3$HT.sono~yyy_T3$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points to getting resting length for normalizing distances
yyy_i4<-identify(yyy_T3$Time, yyy_T3$HT.sono, labels= seq_along(yyy_T3$HT.sono))
# normalize strain values in each file using their 'resting' distances
yyy_T3$HT.sono<-yyy_T3$HT.sono-(mean(yyy_T3$HT.sono[375:2845]))

head(yyy_T4)
plot(yyy_T4$HT.sono~yyy_T4$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points to getting resting length for normalizing distances
yyy_i5<-identify(yyy_T4$Time, yyy_T4$HT.sono, labels= seq_along(yyy_T4$HT.sono))
# normalize strain values in each file using their 'resting' distances
yyy_T4$HT.sono<-yyy_T4$HT.sono-(mean(yyy_T4$HT.sono[c(915:1050,1685:1735,
	1820:2275,2400:2560)]))

### PC ###

# Bird Www Landing 
head(www_L1)
plot(www_L1$PC.sono~www_L1$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points to getting resting length for normalizing distances
www_in1<-identify(www_L1$Time, www_L1$PC.sono, labels= seq_along(www_L1$PC.sono))
# normalize strain values in each file using their 'resting' distances
www_L1$PC.sono<-www_L1$PC.sono-(mean(www_L1$PC.sono[10:4200]))

# Bird Www Take-off (x4) all on same panel
head(www_T1)
plot(www_T1$PC.sono~www_T1$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points to getting resting length for normalizing distances
www_in2<-identify(www_T1$Time, www_T1$PC.sono, labels= seq_along(www_T1$PC.sono))
# normalize strain values in each file using their 'resting' distances
www_T1$PC.sono<-www_T1$PC.sono-(mean(www_T1$PC.sono[1:1360]))

head(www_T2)
plot(www_T2$PC.sono~www_T2$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points to getting resting length for normalizing distances
www_in3<-identify(www_T2$Time, www_T2$PC.sono, labels= seq_along(www_T2$PC.sono))
# normalize strain values in each file using their 'resting' distances
www_T2$PC.sono<-www_T2$PC.sono-(mean(www_T2$PC.sono[1:3110]))

head(www_T3)
plot(www_T3$PC.sono~www_T3$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points to getting resting length for normalizing distances
www_in4<-identify(www_T3$Time, www_T3$PC.sono, labels= seq_along(www_T3$PC.sono))
# normalize strain values in each file using their 'resting' distances
www_T3$PC.sono<-www_T3$PC.sono-(mean(www_T3$PC.sono[1:2445]))

head(www_T4)
plot(www_T4$PC.sono~www_T4$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points to getting resting length for normalizing distances
www_in5<-identify(www_T4$Time, www_T4$PC.sono, labels= seq_along(www_T4$PC.sono))
# normalize strain values in each file using their 'resting' distances
www_T4$PC.sono<-www_T4$PC.sono-(mean(www_T4$PC.sono[1:1720]))

# Bird Xxx Landing 
head(xxx_L1)
plot(xxx_L1$PC.sono~xxx_L1$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points to getting resting length for normalizing distances
xxx_in1<-identify(xxx_L1$Time, xxx_L1$PC.sono, labels= seq_along(xxx_L1$PC.sono))
# normalize strain values in each file using their 'resting' distances
xxx_L1$PC.sono<-xxx_L1$PC.sono-(mean(xxx_L1$PC.sono[1:2820]))

# Bird Xxx Take-off (x4) all on same panel
head(xxx_T1)
plot(xxx_T1$PC.sono~xxx_T1$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points to getting resting length for normalizing distances
xxx_in2<-identify(xxx_T1$Time, xxx_T1$PC.sono, labels= seq_along(xxx_T1$PC.sono))
# normalize strain values in each file using their 'resting' distances
xxx_T1$PC.sono<-xxx_T1$PC.sono-(mean(xxx_T1$PC.sono[1:4630]))

head(xxx_T2)
plot(xxx_T2$PC.sono~xxx_T2$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points to getting resting length for normalizing distances
xxx_in3<-identify(xxx_T2$Time, xxx_T2$PC.sono, labels= seq_along(xxx_T2$PC.sono))
# normalize strain values in each file using their 'resting' distances
xxx_T2$PC.sono<-xxx_T2$PC.sono-(mean(xxx_T2$PC.sono[1:3420]))

head(xxx_T3)
plot(xxx_T3$PC.sono~xxx_T3$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points to getting resting length for normalizing distances
xxx_in4<-identify(xxx_T3$Time, xxx_T3$PC.sono, labels= seq_along(xxx_T3$PC.sono))
# normalize strain values in each file using their 'resting' distances
xxx_T3$PC.sono<-xxx_T3$PC.sono-(mean(xxx_T3$PC.sono[1:7155]))

head(xxx_T4)
plot(xxx_T4$PC.sono~xxx_T4$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points to getting resting length for normalizing distances
xxx_in5<-identify(xxx_T4$Time, xxx_T4$PC.sono, labels= seq_along(xxx_T4$PC.sono))
# normalize strain values in each file using their 'resting' distances
xxx_T4$PC.sono<-xxx_T4$PC.sono-(mean(xxx_T4$PC.sono[1:3150]))

# Bird Uuu Take-off (x2)
head(uuu_T1)
plot(uuu_T1$PC.sono~uuu_T1$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points to getting resting length for normalizing distances
uuu_in1<-identify(uuu_T1$Time, uuu_T1$PC.sono, labels= seq_along(uuu_T1$PC.sono))
# normalize strain values in each file using their 'resting' distances
uuu_T1$PC.sono<-uuu_T1$PC.sono-(mean(uuu_T1$PC.sono[45:1480]))

head(uuu_T2)
plot(uuu_T2$PC.sono~uuu_T2$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points to getting resting length for normalizing distances
uuu_in2<-identify(uuu_T2$Time, uuu_T2$PC.sono, labels= seq_along(uuu_T2$PC.sono))
# normalize strain values in each file using their 'resting' distances
uuu_T2$PC.sono<-uuu_T2$PC.sono-(mean(uuu_T2$PC.sono[35:1770]))

# Bird Uuu Landing (x3) all on same panel
head(uuu_L1)
plot(uuu_L1$PC.sono~uuu_L1$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points to getting resting length for normalizing distances
uuu_in3<-identify(uuu_L1$Time, uuu_L1$PC.sono, labels= seq_along(uuu_L1$PC.sono))
# normalize strain values in each file using their 'resting' distances
uuu_L1$PC.sono<-uuu_L1$PC.sono-(mean(uuu_L1$PC.sono[80:4005]))

head(uuu_L2)
plot(uuu_L2$PC.sono~uuu_L2$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points to getting resting length for normalizing distances
uuu_in4<-identify(uuu_L2$Time, uuu_L2$PC.sono, labels= seq_along(uuu_L2$PC.sono))
# normalize strain values in each file using their 'resting' distances
uuu_L2$PC.sono<-uuu_L2$PC.sono-(mean(uuu_L2$PC.sono[50:3550]))

head(uuu_L3)
plot(uuu_L3$PC.sono~uuu_L3$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points to getting resting length for normalizing distances
uuu_in5<-identify(uuu_L3$Time, uuu_L3$PC.sono, labels= seq_along(uuu_L3$PC.sono))
# normalize strain values in each file using their 'resting' distances
uuu_L3$PC.sono<-uuu_L3$PC.sono-(mean(uuu_L3$PC.sono[35:900]))

# Bird Tie Take-off 
head(tie_T1)
plot(tie_T1$PC.sono~tie_T1$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
tie_in1<-identify(tie_T1$Time, tie_T1$PC.sono, labels= seq_along(tie_T1$PC.sono))
# normalize strain values in each file using their 'resting' distances
tie_T1$PC.sono<-tie_T1$PC.sono-(mean(tie_T1$PC.sono[45:2765]))

# Bird Tie Landing (x4) all on same panel
head(tie_L1)
plot(tie_L1$PC.sono~tie_L1$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
tie_in2<-identify(tie_L1$Time, tie_L1$PC.sono, labels= seq_along(tie_L1$PC.sono))
# normalize strain values in each file using their 'resting' distances
tie_L1$PC.sono<-tie_L1$PC.sono-(mean(tie_L1$PC.sono[c(60:3675,4875:5750)]))

head(tie_L2)
plot(tie_L2$PC.sono~tie_L2$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
tie_in3<-identify(tie_L2$Time, tie_L2$PC.sono, labels= seq_along(tie_L2$PC.sono))
# normalize strain values in each file using their 'resting' distances
tie_L2$PC.sono<-tie_L2$PC.sono-(mean(tie_L2$PC.sono[c(95:3770,4310:6460)]))

head(tie_L3)
plot(tie_L3$PC.sono~tie_L3$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
tie_in4<-identify(tie_L3$Time, tie_L3$PC.sono, labels= seq_along(tie_L3$PC.sono))
# normalize strain values in each file using their 'resting' distances
tie_L3$PC.sono<-tie_L3$PC.sono-(mean(tie_L3$PC.sono[c(390:2850,3360:8590)]))

head(tie_L4)
plot(tie_L4$PC.sono~tie_L4$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
tie_in5<-identify(tie_L4$Time, tie_L4$PC.sono, labels= seq_along(tie_L4$PC.sono))
# normalize strain values in each file using their 'resting' distances
tie_L4$PC.sono<-tie_L4$PC.sono-(mean(tie_L4$PC.sono[70:2170]))

# Bird Yyy Landing 
head(yyy_L1)
plot(yyy_L1$PC.sono~yyy_L1$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points to getting resting length for normalizing distances
yyy_in1<-identify(yyy_L1$Time, yyy_L1$PC.sono, labels= seq_along(yyy_L1$PC.sono))
# normalize strain values in each file using their 'resting' distances
yyy_L1$PC.sono<-yyy_L1$PC.sono-(mean(yyy_L1$PC.sono[1:3385]))

# Bird Yyy Take-off (x4) all on same panel
head(yyy_T1)
plot(yyy_T1$PC.sono~yyy_T1$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points to getting resting length for normalizing distances
yyy_in2<-identify(yyy_T1$Time, yyy_T1$PC.sono, labels= seq_along(yyy_T1$PC.sono))
# normalize strain values in each file using their 'resting' distances
yyy_T1$PC.sono<-yyy_T1$PC.sono-(mean(yyy_T1$PC.sono[1:4950]))

head(yyy_T2)
plot(yyy_T2$PC.sono~yyy_T2$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points to getting resting length for normalizing distances
yyy_in3<-identify(yyy_T2$Time, yyy_T2$PC.sono, labels= seq_along(yyy_T2$PC.sono))
# normalize strain values in each file using their 'resting' distances
yyy_T2$PC.sono<-yyy_T2$PC.sono-(mean(yyy_T2$PC.sono[35:3995]))

head(yyy_T3)
plot(yyy_T3$PC.sono~yyy_T3$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points to getting resting length for normalizing distances
yyy_in4<-identify(yyy_T3$Time, yyy_T3$PC.sono, labels= seq_along(yyy_T3$PC.sono))
# normalize strain values in each file using their 'resting' distances
yyy_T3$PC.sono<-yyy_T3$PC.sono-(mean(yyy_T3$PC.sono[160:2920]))

head(yyy_T4)
plot(yyy_T4$PC.sono~yyy_T4$Time, bty="n", xlab="Time", ylab="delta L (cm)", 
	type="l")
# find points to getting resting length for normalizing distances
yyy_in5<-identify(yyy_T4$Time, yyy_T4$PC.sono, labels= seq_along(yyy_T4$PC.sono))
# normalize strain values in each file using their 'resting' distances
yyy_T4$PC.sono<-yyy_T4$PC.sono-(mean(yyy_T4$PC.sono[c(150:270,465:565,
	715:900,1300:1455,1605:1780,1865:2375)]))

#### Write csvs for normalized distances ####
write.csv(www_L1, "Bird Www/WwwMMr15L_EMGfiltdSONOFiltd_N.csv", row.names=FALSE)
write.csv(www_T1, "Bird Www/WwwMMr07T_EMGfiltdSONOfixedFiltd_N.csv", row.names=FALSE)
write.csv(www_T2, "Bird Www/WwwMMr08T_EMGfiltdSONOfixedFiltd_N.csv", row.names=FALSE)
write.csv(www_T3, "Bird Www/WwwMMr09T_EMGfiltdSONOfixedFiltd_N.csv", row.names=FALSE)
write.csv(www_T4, "Bird Www/WwwMMr10T_EMGfiltdSONOfixedFiltd_N.csv", row.names=FALSE)

write.csv(xxx_L1, "Bird Xxx/XxxMMr36L_EMGfiltdSONOfiltd_N.csv", row.names=FALSE)
write.csv(xxx_T1, "Bird Xxx/XxxMMr30T_EMGfiltdSONOfixedFiltd_N.csv", row.names=FALSE)
write.csv(xxx_T2, "Bird Xxx/XxxMMr35T_EMGfiltdSONOfixedFiltd_N.csv", row.names=FALSE)
write.csv(xxx_T3, "Bird Xxx/XxxMMr37T_EMGfiltdSONOfixedFiltd_N.csv", row.names=FALSE)
write.csv(xxx_T4, "Bird Xxx/XxxMMr39T_EMGfiltdSONOfiltd_N.csv", row.names=FALSE)

write.csv(uuu_T1, "Bird Uuu/TuuuMM06axo_EMGfiltdSONOfixedFiltd_N.csv", row.names=FALSE)
write.csv(uuu_T2, "Bird Uuu/TuuuMM07axo_EMGfiltdSONOfixedFiltd_N.csv", row.names=FALSE)
write.csv(uuu_L1, "Bird Uuu/LuuuMM06axo_EMGfiltdSONOfixedFiltd_N.csv", row.names=FALSE)
write.csv(uuu_L2, "Bird Uuu/LuuuMM07axo_EMGfiltdSONOfixedFiltd_N.csv", row.names=FALSE)
write.csv(uuu_L3, "Bird Uuu/LuuuMM08axo_EMGfiltdSONOfixedFiltd_N.csv", row.names=FALSE)

write.csv(tie_T1, "Bird Tie/TtieMM07axo_EMGfiltdSONOFiltd_N.csv", row.names=FALSE)
write.csv(tie_L1, "Bird Tie/LtieMM04axo_EMGfiltdSONOfixedFiltd_N.csv", row.names=FALSE)
write.csv(tie_L2, "Bird Tie/LtieMM07axo_EMGfiltdSONOfixedFiltd_N.csv", row.names=FALSE)
write.csv(tie_L3, "Bird Tie/LtieMM10axo_EMGfiltdSONOFiltd_N.csv", row.names=FALSE)
write.csv(tie_L4, "Bird Tie/LtieMM12axo_EMGfiltdSONOFiltd_N.csv", row.names=FALSE)

write.csv(yyy_L1, "Bird Yyy/YyyMMr13L_EMGfiltdSONOfixedFiltd_N.csv", row.names=FALSE)
write.csv(yyy_T1, "Bird Yyy/YyyMMr08T_EMGfiltdSONOfixedFiltd_N.csv", row.names=FALSE)
write.csv(yyy_T2, "Bird Yyy/YyyMMr10T_EMGfiltdSONOfixedFiltd_N.csv", row.names=FALSE)
write.csv(yyy_T3, "Bird Yyy/YyyMMr16T_EMGfiltdSONOfixedFiltd_N.csv", row.names=FALSE)
write.csv(yyy_T4, "Bird Yyy/YyyMMr18T_EMGfiltdSONOfixedFiltd_N.csv", row.names=FALSE)

