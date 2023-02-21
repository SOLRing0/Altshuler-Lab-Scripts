########################################################################
########################################################################
# Last edited by B.Goller on 07 June 2016
########################################################################
########################################################################
#Clear the workspace
rm(list = ls())
########################################################################
########################################################################
if(! "circular" %in% .packages(all=TRUE))
  install.packages("circular")
if(! "ggplot2" %in% .packages(all=TRUE))
  install.packages("ggplot2")

library(circular)
library(ggplot2)
########################################################################
########################################################################
#DEFINE FUNCTIONS
########################################################################
########################################################################


#load each cell
#get the preferred direction
#keep the data as well
#plot them and fit a curve with the peak set by the prefdir
#compare the cell curves between cells and species

quick.mean <- function(x){return(mean(x,na.rm = TRUE))}
########################################################################
########################################################################
extract.num.from.string <- function(s) {
  s.split <- strsplit(s, "cell")
  s.id <- as.numeric(unlist(strsplit(s.split[[1]][1], "[^[:digit:]]")))
  s.id <- s.id[!is.na(s.id)][1:3]
  
  s.cell <- as.numeric(unlist(strsplit(s.split[[1]][2], "[^[:digit:]]")))
  s.cell <- s.cell[!is.na(s.cell)][1]
  return(c(s.id, s.cell))
}
########################################################################
########################################################################
response.to.density <- function(dir, firing) {
  #Use the firing to create repeats of the direction value
  output <- NULL
  
  num.dirs <- unique(dir)
  rep.min <- min(as.numeric(lapply(num.dirs, FUN = function(x) length(which(dir == x)))))
  
  for(i in 1:(length(num.dirs)*rep.min)) {
    output <- c(output,rep(dir[i],round(10*firing[i])))
  }
  return(output)
}
########################################################################
########################################################################
simple.tuning <- function(dir, response, threshold){
  result <- rep(0, length(dir))
  for(i in 1:length(dir)){
    if(response[i]>threshold){result[i] <- 1}
  }
  return(result)
}
########################################################################
########################################################################
dir.by.mean2 <- function(dir, firing){
  v <- rep(0,8)
  h <- rep(0,8)
  
  unique.dir <- unique(dir)
  rep.min <- min(as.numeric(lapply(unique.dir, FUN = function(x) length(which(dir == x)))))
  #print(rep.min)
  
  for(i in 1:length(unique.dir)){
    #Vertical component of firing vector
    v[i] <- sum(firing[which(dir == unique.dir[i])[1:rep.min]]*sin(-unique.dir[i]*pi/180))
    #Horizontal component of firing vector
    h[i] <- sum(firing[which(dir == unique.dir[i])[1:rep.min]]*cos(-unique.dir[i]*pi/180))
    #rep.test[i] <- length(firing[which(dir == unique.dir[i])])
  }
  #calculate preferred direction as arctangent of the sums of vertical and horizontal components
  pd <- atan(sum(v)/sum(h))*(180/pi)
  
  ifelse(sum(v)>0, ifelse(sum(h)>0, pd.adj<-360-pd, pd.adj <- 180+abs(pd)),
         ifelse(sum(h)>0,pd.adj <- abs(pd),pd.adj <- 180-abs(pd)))
  
  return(c(sum(v),sum(h),pd, pd.adj))
}
########################################################################
########################################################################
cell.mean.bydir <- function(dir, firing){
  unique.dir <- unique(dir)
  cell.resp.means <- matrix(0,length(unique.dir),2)
  
  for(i in 1:length(unique.dir)){
    cell.resp.means[i,] <- c(unique.dir[i],quick.mean(firing[which(dir == unique.dir[i])]))
  }
  
  return(cell.resp.means)
}
########################################################################
########################################################################
get.preferred.direction <- function(file.list){
  cell.ids <- matrix(NA,length(file.list),4)
  rayleigh.test <- matrix(NA,length(file.list),2)
  prefdir.sum <- matrix(NA,length(file.list),4)
  prefdir.mle <- matrix(NA,length(file.list),2)
  
  for(i in 1:length(file.list))
  {
    prefdir.mle.results <- NULL
    print(i)
    onecell.data <- read.csv(file.list[i], header = FALSE)
    onecell.data <- cbind(apply(onecell.data[,3:22],1,quick.mean),onecell.data)
    cell.ids[i,] <- extract.num.from.string(file.list[i])
    
    onecell.meanbase <- quick.mean(onecell.data[is.na(onecell.data[,2]),1])
    
    if(max(abs(onecell.data[!is.na(onecell.data[,2]),1] - onecell.meanbase)) > (onecell.meanbase * 0.2)) {
      onecell.relfiring <- onecell.data[!is.na(onecell.data[,2]),1]-onecell.meanbase + abs(min(onecell.data[!is.na(onecell.data[,2]),1]-onecell.meanbase))
      for.circ.density <- response.to.density(onecell.data[!is.na(onecell.data[,2]),2], onecell.relfiring)

      ray.test <- rayleigh.test(circular(for.circ.density,type = "angles",units = "degrees",rotation = "clock"))
      rayleigh.test[i,] <- c(round(ray.test$statistic,4), round(ray.test$p.value,4))
      
      if(rayleigh.test[i,2] <= 0.05){
        prefdir.mle.results <- mle.vonmises(circular(for.circ.density,type = "angles",units = "degrees",rotation = "clock"))
        
        prefdir.mle[i,] <- c(ifelse(as.numeric(unlist(prefdir.mle.results["mu"])) < 0, 
                                    360+as.numeric(unlist(prefdir.mle.results["mu"])), 
                                    as.numeric(unlist(prefdir.mle.results["mu"]))),
                             as.numeric(unlist(prefdir.mle.results["kappa"])))
        
        print(sprintf("Mu: %f and 1/Kappa: %f",as.numeric(prefdir.mle.results["mu"]),1/as.numeric(prefdir.mle.results["kappa"])))
        
        prefdir.sum[i,] <- dir.by.mean2(onecell.data[!is.na(onecell.data[,2]),2],(onecell.data[!is.na(onecell.data[,2]),1]-onecell.meanbase))
      }
    }
    else {
      prefdir.sum[i,] <- rep(NA,4)
      prefdir.mle[i] <- NA
      }
  }
  
  
  #col1 is the sum method, col2 is the mle method
  return(cbind(cell.ids,rayleigh.test,prefdir.sum[,4],prefdir.mle))
}
########################################################################
########################################################################
dir <- INSERT PATH TO "2015-CALAN-PSTH-direction/"
setwd(dir)
calan.file.list <- list.files(pattern = "CALAN.+\\.csv$", ignore.case = TRUE)
calan.prefdirs <- get.preferred.direction(calan.file.list)
colnames(calan.prefdirs) <- c("bird", "trk", "site", "cell", "rayleigh.stat", "rayleigh.pval", "pd.sum", "pd.mle", "pd.kappa")
plot(calan.prefdirs[,8],calan.prefdirs[,7], xlab = "hb mle method", ylab = "hb sum method")
########################################################################
########################################################################
#Set the working directory (must contain data file)
dir <- INSERT PATH TO "2015-ZF-PSTH-direction/"
setwd(dir)
zf.file.list <- list.files(pattern = "ZF.+\\.csv$", ignore.case = TRUE)
zf.prefdirs <- get.preferred.direction(zf.file.list)
colnames(zf.prefdirs) <- c("bird", "trk", "site", "cell", "rayleigh.stat", "rayleigh.pval", "pd.sum", "pd.mle", "pd.kappa")
plot(zf.prefdirs[,7],zf.prefdirs[,8], xlab = "zf mle method", ylab = "zf sum method")
########################################################################
########################################################################
#Pigeon files need extra work to process because they do not have a set number of directions tested
#Set the working directory (must contain data file)
dir <- INSERT PATH TO "PIGEON-PSTH-direction/"
setwd(dir)
pg.file.list <- list.files(pattern = "\\.csv$", ignore.case = TRUE)
pg.prefdirs <- get.preferred.direction(pg.file.list)
colnames(pg.prefdirs) <- c("bird", "trk", "site", "cell", "rayleigh.stat", "rayleigh.pval", "pd.sum", "pd.mle", "pd.kappa")
plot(pg.prefdirs[,7],pg.prefdirs[,8], xlab = "pg mle method", ylab = "pg sum method")
########################################################################
########################################################################
#generate a plot of the von Mises distributions
step <- seq(-180,180, by = 5)

cell.vonmises.curves <- matrix(NA,(nrow(calan.prefdirs)+nrow(zf.prefdirs)+nrow(pg.prefdirs)),(length(step)+1))

#plot(f1(circular(d1, units = 'degrees', rotation = 'clock'))+off1, type = 'l', ylim = c(0.3,1.1), col = 2)

for(i in 1:dim(calan.prefdirs)[1]){
  if(!is.na(calan.prefdirs[i,8])){
    fx <- function(x) dvonmises(x, circular(calan.prefdirs[i,8], units = 'degrees', rotation = 'clock'), calan.prefdirs[i,9])
    dx <- rep(calan.prefdirs[i,8], length(step))+step
    offx <- 1-max(fx(circular(dx, units = 'degrees', rotation = 'clock')))
    cell.vonmises.curves[i,] <- c(2,fx(circular(dx, units = 'degrees', rotation = 'clock')) + offx)
    #lines(fx(circular(dx, units = 'degrees', rotation = 'clock')) + offx, type = 'l', col = 2)
  }
}

for(j in 1:dim(zf.prefdirs)[1]){
  if(!is.na(zf.prefdirs[j,8])){
    fx <- function(x) dvonmises(x, circular(zf.prefdirs[j,8], units = 'degrees', rotation = 'clock'), zf.prefdirs[j,9])
    dx <- rep(zf.prefdirs[j,8], length(step))+step
    offx <- 1-max(fx(circular(dx, units = 'degrees', rotation = 'clock')))
    cell.vonmises.curves[(j+dim(calan.prefdirs)[1]),] <- c(1,fx(circular(dx, units = 'degrees', rotation = 'clock')) + offx)
    #lines(fx(circular(dx, units = 'degrees', rotation = 'clock')) + offx, type = 'l', col = 1)
  }
}

for(k in 1:dim(pg.prefdirs)[1]){
  if(!is.na(pg.prefdirs[k,8])){
    fx <- function(x) dvonmises(x, circular(pg.prefdirs[k,8], units = 'degrees', rotation = 'clock'), pg.prefdirs[k,9])
    dx <- rep(pg.prefdirs[k,8], length(step))+step
    offx <- 1-max(fx(circular(dx, units = 'degrees', rotation = 'clock')))
    cell.vonmises.curves[(k+dim(calan.prefdirs)[1]),] <- c(3,fx(circular(dx, units = 'degrees', rotation = 'clock')) + offx)
    #lines(fx(circular(dx, units = 'degrees', rotation = 'clock')) + offx, type = 'l', col = 1)
  }
}

vonmises.df <- data.frame(cell.vonmises.curves[!is.na(cell.vonmises.curves[,1]),])
vonmises.df$X1 <- factor(vonmises.df$X1,levels=c(1,2,3), labels=c("zf","hb","pg"))

# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# To use for fills, add
#scale_fill_manual(values=cbPalette)
# To use for line and point colors, add
#scale_colour_manual(values=cbPalette)

plot(as.numeric(vonmises.df[1,(2:dim(vonmises.df)[2])]), type = 'l', col = cbbPalette[as.numeric(vonmises.df[i,1])], ylim = c(0.3,1.1))
for(i in 2:dim(vonmises.df)[1]){
  lines(as.numeric(vonmises.df[i,(2:dim(vonmises.df)[2])]), type = 'l', col = cbbPalette[as.numeric(vonmises.df[i,1])])
}


dp1 <- qplot(X74, data=vonmises.df, geom="density", fill=X1, alpha=I(.5), 
             main="Preferred Direction +- 180deg", xlab="Intercept", ylab="Relative Proportion", xlim = c(0.3,1.1)) +
             scale_fill_manual(values = cbbPalette)
dp1 + coord_flip()

dp2 <- qplot(X56, data=vonmises.df, geom="density", fill=X1, alpha=I(.5), 
             main="Preferred Direction +- 90deg", xlab="Intercept", ylab="Relative Proportion", xlim = c(0.3,1.1)) +
  scale_fill_manual(values = cbbPalette)
dp2 + coord_flip()

dp3 <- qplot(X47, data=vonmises.df, geom="density", fill=X1, alpha=I(.5), 
             main="Preferred Direction +- 45deg", xlab="Intercept", ylab="Relative Proportion", xlim = c(0.3,1.1)) +
  scale_fill_manual(values = cbbPalette)
dp3 + coord_flip()

dp4 <- qplot(X41, data=vonmises.df, geom="density", fill=X1, alpha=I(.5), 
             main="Preferred Direction +- 15deg", xlab="Intercept", ylab="Relative Proportion", xlim = c(0.3,1.1)) +
  scale_fill_manual(values = cbbPalette)
dp4 + coord_flip()