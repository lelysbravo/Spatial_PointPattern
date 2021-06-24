#### Updated version of maptools needed
install.packages("maptools", repos = "http://R-Forge.R-project.org")
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
#
#
###########Read RDS files
STP1987<-readRDS("STPSUMAX1987.rds")
reports1987<-readRDS("reports-1987.rds")
STP2011<-readRDS("STPSUMAX2011.rds")
reports2011<-readRDS("reports-2011.rds")
#
#Plot STP files Year 1987 Quarter 2
#plot(STPXXX,col=viridis(10),breaks=c(1,5,10,15,20,25,30,35,40,45,50))
plot(STP1987[[2]],col=terrain.colors(10),breaks=c(1,5,10,15,20,25,30,35,40,45,50),xlab="Longitude",ylab="Latitude")
data("wrld_simpl")
map("state",add=TRUE)
# Plot second quarter reports
plot(reports1987[[2]],add=TRUE,pch=c(3),cex=0.7)

#Plot STP files Year 2011 Quarter 3
#plot(STPXXX,col=viridis(10),breaks=c(1,5,10,15,20,25,30,35,40,45,50))
plot(STP2011[[3]],col=terrain.colors(10),breaks=c(1,5,10,20,30,40,50,60,70,80),xlab="Longitude",ylab="Latitude")
data("wrld_simpl")
map("state",add=TRUE)
#Plot third quarter reports
plot(reports2011[[3]],add=TRUE,pch=c(3),cex=0.7)
##### Estimate rhohat
## 1987 Quarter 2
dev.off()
rho1987Q2=rhohat(unmark(reports1987[[2]]),as.im(STP1987[[2]]),smoother="local",method="ratio") 
plot(rho1987Q2,xlab="",ylab="",main= " Year 1987 Quarter 2",ylim=c(0,2),
     cex.main=1.7,cex.axis=1.1,tcl=-0.2,cex.lab=1,legend=FALSE,font.lab=2,cex=1)
ylab.text = expression(rho*(Z(u)))
mtext(ylab.text,side=2, line =2.0,font.lab=2,cex=1)
xlab.text="sumax(STP)"
mtext(xlab.text,side=1, line =2.0,font.lab=2,cex=1)
## CDF test
cdft=cdf.test(unmark(reports1987[[2]]),as.im(STP1987[[2]]), test="ad")
plot(cdf.test(unmark(reports1987[[2]]),as.im(STP1987[[2]]), test="ad"),xlab="sumax(STP)",
     main=paste("Year 1987 Q2, P-Value AD test:", round(cdft$p.value,6),sep=" "),
     do.legend=FALSE,
     cex.main=1.2,cex.axis=1.5,ylab="Probability",tcl=-0.2,cex.lab=1.5,lwd = 2.2)
#Initialize lists
berman1987=list()
roc1987=list()
auc1987=list()
#Berman test for Year 1987 all Quarters
for (j in 1:4){
  berman1987[[j]]= berman.test(unmark(reports1987[[j]]),as.im(STP1987[[j]]),"Z2")
  print(berman1987[[j]])
}
berman1987
#ROC curves and AUC for Year 1987 all Quarters
for (j in 1:4){
  roc1987[[j]]=roc(unmark(reports1987[[j]]),as.im(STP1987[[j]]),high=FALSE)
  auc1987[[j]]=auc(unmark(reports1987[[j]]),as.im(STP1987[[j]]),high=FALSE)
}
