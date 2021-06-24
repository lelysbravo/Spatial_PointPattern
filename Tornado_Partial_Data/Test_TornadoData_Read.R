library(spatstat)
library(raster)
library(rasterVis)
library(ncdf4)
library(reshape2)
library(dplyr)
library(plyr)
library(maptools)
library(maps)
library(sf)
library(tigris)
library(astsa)
library(TSA)
library(mapproj)
library(viridis)
## Testing load .RData files
#Load STP files Year 1987, Quarter 2
STPXXX=get(load("VNc1987Q2.RData",.GlobalEnv))
#STPXXX=load("VNc1987Q2.RData")
#Plot STP files
#plot(STPXXX,col=viridis(10),breaks=c(1,5,10,15,20,25,30,35,40,45,50))
plot(STPXXX,col=terrain.colors(10),breaks=c(1,5,10,15,20,25,30,35,40,45,50),xlab="Longitude",ylab="Latitude")
data("wrld_simpl")
map("state",add=TRUE)
#Load Tornado Reports
REPXXX=get(load("reports-1987.RData",envir = parent.frame(), verbose = FALSE))
# Plot second quarter reports
plot(REPXXX[[2]],add=TRUE,pch=c(3),cex=0.7)

## Testing load .RData files
#Load STP files Year 2011, Quarter 3
STPYYY=get(load("VNc2011Q2.RData"))
#Plot STP files
#plot(STPXXX,col=viridis(10),breaks=c(1,5,10,15,20,25,30,35,40,45,50))
plot(STPYYY,col=terrain.colors(10),breaks=c(1,5,10,15,20,25,30,35,40,45,50),xlab="Longitude",ylab="Latitude")
data("wrld_simpl")
map("state",add=TRUE)
#Load Tornado Reports Year 2011
REPYYY=get(load("reports-2011.RData"))
#Plot third quaerter reports
plot(REPYYY[[2]],add=TRUE,pch=c(3),cex=0.7)

##### Estimate rhohat
## 1987 Quarter 2
dev.off()
rho1987Q2=rhohat(unmark(REPXXX[[2]]),as.im(STPXXX),smoother="local",method="ratio") 
plot(rho1987Q2,xlab="",ylab="",main= " Year 1987 Quarter 2",ylim=c(0,2),
     cex.main=1.7,cex.axis=1.1,tcl=-0.2,cex.lab=1,legend=FALSE,font.lab=2,cex=1)
ylab.text = expression(rho*(Z(u)))
mtext(ylab.text,side=2, line =2.0,font.lab=2,cex=1)
xlab.text="sumax(STP)"
mtext(xlab.text,side=1, line =2.0,font.lab=2,cex=1)
## CDF test
cdft=cdf.test(unmark(REPXXX[[2]]),as.im(STPXXX), test="ad")
plot(cdf.test(unmark(REPXXX[[2]]),as.im(STPXXX), test="ad"),xlab="sumax(STP)",
     main=paste("Year 1987 Q2, P-Value AD test:", round(cdft$p.value,6),sep=" "),
     do.legend=FALSE,
     cex.main=1.2,cex.axis=1.5,ylab="Probability",tcl=-0.2,cex.lab=1.5,lwd = 2.2)

