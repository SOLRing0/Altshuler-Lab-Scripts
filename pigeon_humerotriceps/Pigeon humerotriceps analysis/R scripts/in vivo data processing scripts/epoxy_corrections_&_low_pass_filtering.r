# Script to apply epoxy distance corrections and low pass filtering to Robertson and Biewener 
# strain data for HT and Pec -- file naming convention indicates which ones need further processing
## load package for filtering signals (see paper for filtering parameters)
library(signal)

# get and set working directory 
getwd()
## select folder containing all in vivo data for processing 
setwd(choose.dir())

# read in csvs for Bird Www

www_L1<-read.csv("Bird Www/WwwMMr15L_EMGfiltd.csv")
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

uuu_T1<-read.csv("Bird Uuu/TuuuMM06axo_EMGfiltdSONOfixed.csv")
uuu_T2<-read.csv("Bird Uuu/TuuuMM07axo_EMGfiltdSONOfixed.csv")
uuu_L1<-read.csv("Bird Uuu/LuuuMM06axo_EMGfiltdSONOfixed.csv")
uuu_L2<-read.csv("Bird Uuu/LuuuMM07axo_EMGfiltdSONOfixed.csv")
uuu_L3<-read.csv("Bird Uuu/LuuuMM08axo_EMGfiltdSONOfixed.csv")

# read in csvs for Bird Tie

tie_T1<-read.csv("Bird Tie/TtieMM07axo_EMGfiltd.csv")
tie_L1<-read.csv("Bird Tie/LtieMM04axo_EMGfiltdSONOfixed.csv")
tie_L2<-read.csv("Bird Tie/LtieMM07axo_EMGfiltdSONOfixed.csv")
tie_L3<-read.csv("Bird Tie/LtieMM10axo_EMGfiltd.csv")
tie_L4<-read.csv("Bird Tie/LtieMM12axo_EMGfiltd.csv")

# read in csvs for Bird Yyy

yyy_L1<-read.csv("Bird Yyy/YyyMMr13L_EMGfiltdSONOfixedFiltd.csv")
yyy_T1<-read.csv("Bird Yyy/YyyMMr08T_EMGfiltdSONOfixedFiltd.csv")
yyy_T2<-read.csv("Bird Yyy/YyyMMr10T_EMGfiltdSONOfixed.csv")
yyy_T3<-read.csv("Bird Yyy/YyyMMr16T_EMGfiltdSONOfixed.csv")
yyy_T4<-read.csv("Bird Yyy/YyyMMr18T_EMGfiltdSONOfixedFiltd.csv")

# 250 Hz low pass filter strain data that is unfiltered
# and add correction for epoxy forming the crystals: 0.82mm for 2mm crystals 
# (PC), and 0.16mm for 1mm crystals (HT)
# butter(n,W,type) -- n=order, use 1 to prevent overshooting
# W is cutoff frequency expressed as a fraction of the Nyquist frequency,
# which is 1/2 of sampling frequency (in this case 5000/2)
# W=0.1 to get 250Hz low pass frequency
bf <- butter(1, 0.1, type="low") 

# Bird Www HT
www_L1$HT.sono <- (filter(bf, www_L1$HT.sono))+(0.16/10)

# Bird Www PC
www_L1$PC.sono <- (filter(bf, www_L1$PC.sono))+(0.82/10)

# write new csv for filtered SONO -- Bird Www
write.csv(www_L1, "Bird Www/WwwMMr15L_EMGfiltdSONOFiltd.csv", row.names=FALSE)

# Bird Uuu HT
uuu_T1$HT.sono <- (filter(bf, uuu_T1$HT.sono))+(0.16/10)
uuu_T2$HT.sono <- (filter(bf, uuu_T2$HT.sono))+(0.16/10)
uuu_L1$HT.sono <- (filter(bf, uuu_L1$HT.sono))+(0.16/10)
uuu_L2$HT.sono <- (filter(bf, uuu_L2$HT.sono))+(0.16/10)
uuu_L3$HT.sono <- (filter(bf, uuu_L3$HT.sono))+(0.16/10)

# Bird Uuu PC
uuu_T1$PC.sono <- (filter(bf, uuu_T1$PC.sono))+(0.82/10)
uuu_T2$PC.sono <- (filter(bf, uuu_T2$PC.sono))+(0.82/10)
uuu_L1$PC.sono <- (filter(bf, uuu_L1$PC.sono))+(0.82/10)
uuu_L2$PC.sono <- (filter(bf, uuu_L2$PC.sono))+(0.82/10)
uuu_L3$PC.sono <- (filter(bf, uuu_L3$PC.sono))+(0.82/10)

# write new csvs for filtered SONO -- Bird Uuu
write.csv(uuu_T1, "Bird Uuu/TuuuMM06axo_EMGfiltdSONOfixedFiltd.csv", row.names=FALSE)
write.csv(uuu_T2, "Bird Uuu/TuuuMM07axo_EMGfiltdSONOfixedFiltd.csv", row.names=FALSE)
write.csv(uuu_L1, "Bird Uuu/LuuuMM06axo_EMGfiltdSONOfixedFiltd.csv", row.names=FALSE)
write.csv(uuu_L2, "Bird Uuu/LuuuMM07axo_EMGfiltdSONOfixedFiltd.csv", row.names=FALSE)
write.csv(uuu_L3, "Bird Uuu/LuuuMM08axo_EMGfiltdSONOfixedFiltd.csv", row.names=FALSE)

# Bird Tie HT
tie_T1$HT.sono <- (filter(bf, tie_T1$HT.sono))+(0.16/10)
tie_L1$HT.sono <- (filter(bf, tie_L1$HT.sono))+(0.16/10)
tie_L2$HT.sono <- (filter(bf, tie_L2$HT.sono))+(0.16/10)
tie_L3$HT.sono <- (filter(bf, tie_L3$HT.sono))+(0.16/10)
tie_L4$HT.sono <- (filter(bf, tie_L4$HT.sono))+(0.16/10)

# Bird Tie PC
tie_T1$PC.sono <- (filter(bf, tie_T1$PC.sono))+(0.82/10)
tie_L1$PC.sono <- (filter(bf, tie_L1$PC.sono))+(0.82/10)
tie_L2$PC.sono <- (filter(bf, tie_L2$PC.sono))+(0.82/10)
tie_L3$PC.sono <- (filter(bf, tie_L3$PC.sono))+(0.82/10)
tie_L4$PC.sono <- (filter(bf, tie_L4$PC.sono))+(0.82/10)

# write new csvs for filtered SONO -- Bird Tie
write.csv(tie_T1, "Bird Tie/TtieMM07axo_EMGfiltdSONOFiltd.csv", row.names=FALSE)
write.csv(tie_L1, "Bird Tie/LtieMM04axo_EMGfiltdSONOfixedFiltd.csv", row.names=FALSE)
write.csv(tie_L2, "Bird Tie/LtieMM07axo_EMGfiltdSONOfixedFiltd.csv", row.names=FALSE)
write.csv(tie_L3, "Bird Tie/LtieMM10axo_EMGfiltdSONOFiltd.csv", row.names=FALSE)
write.csv(tie_L4, "Bird Tie/LtieMM12axo_EMGfiltdSONOFiltd.csv", row.names=FALSE)

# Bird Yyy HT
yyy_T2$HT.sono <- (filter(bf, yyy_T2$HT.sono))+(0.16/10)
yyy_T3$HT.sono <- (filter(bf, yyy_T3$HT.sono))+(0.16/10)

# Bird Yyy PC
yyy_T2$PC.sono <- (filter(bf, yyy_T2$PC.sono))+(0.82/10)
yyy_T3$PC.sono <- (filter(bf, yyy_T3$PC.sono))+(0.82/10)

# write new csvs for filtered SONO -- Bird Yyy
write.csv(yyy_T2, "Bird Yyy/YyyMMr10T_EMGfiltdSONOfixedFiltd.csv", row.names=FALSE)
write.csv(yyy_T3, "Bird Yyy/YyyMMr16T_EMGfiltdSONOfixedFiltd.csv", row.names=FALSE)

