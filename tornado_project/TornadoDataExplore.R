#Tornado Data Analysis
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
#library(RcppCNPy)
setwd("~/tornado_project")
datator=read.csv("1950-2017_tor.csv")
attach(datator)
#Reading netcdf files with the netcdf4 library
nc_stp_max=nc_open("stp_diurnal_max_79_17.nc")
print(nc_stp_max)
attributes(nc_stp_max)$names
attributes(nc_stp_max$var)$names
attributes(nc_stp_max$dim)$names
nc_stp_sum=nc_open("stp_diurnal_sum_79_17.nc")
print(nc_stp_sum)
attributes(nc_stp_sum)$names
attributes(nc_stp_sum$var)$names
attributes(nc_stp_sum$dim)$names
ncatt_get(nc_stp_max, attributes(nc_stp_max$var)$names)
#Read Lat lon data
latlon=nc_open("narr_latlon.nc")
latsras=raster("narr_latlon.nc",varname="lats")
lonsras=raster("narr_latlon.nc",varname="lons")
dim(latsras)
dim(lonsras)
#Get netcdf files values
stpmax = ncvar_get(nc_stp_max, attributes(nc_stp_max$var)$names)
stpsum = ncvar_get(nc_stp_sum, attributes(nc_stp_sum$var)$names)
dim(stpmax)
dim(stpsum)
# Close netcdf files
nc_close(nc_stp_max)
nc_close(nc_stp_sum)
#Read US borderline and States shapefiles
XUS=st_read("US_outline.shp")
XUSstates=st_read("states.shp")
#Setting mising values to zero
stpmax[stpmax>1000000]=0
stpsum[stpsum>1000000]=0
#Permuting elements of stpmax and stpsum
stpsum=aperm(stpsum,c(2,1,4,3))
stpmax=aperm(stpmax,c(2,1,4,3))
dim(stpsum)
dim(stpmax)
#Creating annual brick files with the brick function from library raster
stpmaxbr=list()
for (i in 1:dim(stpmax)[4]){
  stpmaxbr[[i]]=brick(ncols=349,nrows=277,nl=365)
  values(stpmaxbr[[i]])=stpmax[,,,i]
  projection(stpmaxbr[[i]])="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
}
rm(stpmax)
stpsumbr=list()
for (i in 1:dim(stpsum)[4]){
  stpsumbr[[i]]=brick(ncols=349,nrows=277,nl=365)
  values(stpsumbr[[i]])=stpsum[,,,i]
  projection(stpsumbr[[i]])="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
}
rm(stpsum)
extent(stpsumbr[[1]])
crs(stpsumbr[[1]])
ncell(stpsumbr[[1]])
nlayers(stpsumbr[[1]])
#
#Clean tornado data
#Bounding Box Continental US
lamin=24.396308
lamax=49.384358
lomin=-124.848974
lomax=-66.885444
datator_c=datator %>% filter(slat<=lamax & slat>=lamin & slon <=lomax & slon >=lomin)
plot(datator_c$slon,datator_c$slat,cex=0.5,col="red",xlab="Longitude",ylab="Latitude")
data("wrld_simpl")
plot(wrld_simpl,add=TRUE)
dimtor=dim(datator_c)
inyear=datator_c$yr[1]
fnyear=datator_c$yr[dimtor[1]]
nyear=fnyear-inyear+1
#Save Tornado Data by quarters and years and convert to a ppp object
torpppyr=list()
torqt=list()
datappp=list()
qt=list(c(12,1,2),c(3,4,5),c(6,7,8),c(9,10,11))
for (i in 1:nyear){
  for (j in 1:4){
    k=i+inyear-1
    torqt[[j]]=datator_c %>% filter(yr==k & mo %in% qt[[j]])
  }
  for (j in 1:4){
    print(j)
    datappp[[j]]=ppp(torqt[[j]]$slon,torqt[[j]]$slat,
                     c(lomin,lomax),c(lamin,lamax),
                     marks=as.factor(torqt[[j]]$mag))
  }
  print(i)
  torpppyr[[i]]=datappp
}
#Calculate the intensity matrix
intmat=matrix(0,ncol=4,nrow=nyear)
for(i in 1:nyear){
  for(j in 1:4){
  intmat[i,j]=intensity(unmark(torpppyr[[i]][[j]]))  
  }
}
ts.plot(ts(apply(intmat,1,mean),start=c(1950)))
colnames(intmat)=quarter.name
row.names(intmat)=as.character(inyear:fnyear)
intmat
write.csv(intmat,"~/tornado_project/intmat.csv")
dev.off()
pdf("reports.pdf")
quarter.name=c("Dec-Feb","Mar-May","Jun-Aug","Sep-Nov")
par(mfrow=c(2,2))
for (i in 1:nyear){
 for (j in 1:4){
   plot(as.ppp(torpppyr[[i]][[j]]),main=paste(quarter.name[j]," ",inyear+i-1,sep=""))
   plot(wrld_simpl,add=TRUE)
 }
}
dev.off()
#
#Converting raster files lat lon to points
VX=rasterToPoints(lonsras,spatial=FALSE)
#Ordering 2nd column
VXo=VX[order(VX[,2]),]
#Converting rasters to points
VY=rasterToPoints(latsras,spatial=FALSE)
#Ordering 2nd column
VYo=VY[order(VY[,2]),]
testpts=cbind(VXo[,3],VYo[,3])
names(testpts)=c("x","y")
rtest=raster(ncols=349,nrows=277)
#inidayqt=c(1,91,182,274)
#findayqt=c(90,181,273,365)
daysinqt=list(c(335:365,1:59),c(60:151),c(152:243),c(244:334))
####################################################################
#Analisys with sum stp
###################################################################
stpsumyr=list()
stpsumqt=list()
#Attaching lat lon data to the stp files
for (i in 1:length(stpsumbr)){
  for (j in 1:4){
  QTSUM= sum(stack(stpsumbr[[i]],layers = daysinqt[[j]])) 
  VV=rasterToPoints(QTSUM)
  stpsumqt[[j]]=rasterize(testpts,rtest,VV[,3])
  }
  print(i)
 stpsumyr[[i]] = stpsumqt
}
auc.ts=ts(rep(0,length(stpsumbr)*4),start=c(1979,1),freq=4)
k=1
for (i in 1:length(stpsumbr)){
  for (j in 1:4){
   auc.ts[k]= auc(unmark(torpppyr[[28+i]][[j]]),as.im(stpsumyr[[i]][[j]]),high=FALSE)
   print(auc.ts[k])
         k=k+1
  }
}
#Plot auc index and fit a trend
ts.plot(auc.ts,ylab="AUC index",main="Sum STP",xlab="Year",ylim=c(0,0.3),type="n")
grid(lty=1, col=gray(.9)); lines(auc.ts,type="b",pch=21:24)
trauc=lm(auc.ts~time(auc.ts))
abline(trauc,col="red")
boxplot(auc.ts~season(auc.ts),ylab="AUC Index",main="Sum STP")
modelauc=lm(auc.ts~time(auc.ts)+season(auc.ts)-1)
summary(modelauc)
####################################
#Analisys with sum of max stp
###################################
stpsumaxyr=list()
stpsumaxqt=list()
#Attaching lat lon data to the stp files
for (i in 1:length(stpmaxbr)){
  for (j in 1:4){
    QTSUM= sum(stack(stpmaxbr[[i]],layers = daysinqt[[j]])) 
    VV=rasterToPoints(QTSUM)
    stpsumaxqt[[j]]=rasterize(testpts,rtest,VV[,3])
  }
  print(i)
  stpsumaxyr[[i]] = stpsumaxqt
}
aucmax.ts=ts(rep(0,length(stpmaxbr)*4),start=c(1979,1),freq=4)
k=1
for (i in 1:length(stpmaxbr)){
  for (j in 1:4){
    aucmax.ts[k]= auc(unmark(torpppyr[[28+i]][[j]]),as.im(stpsumaxyr[[i]][[j]]),high=FALSE)
    print(aucmax.ts[k])
    k=k+1
  }
}
ts.plot(aucmax.ts,ylab="AUC index",main="sumax(STP)",xlab="Year",ylim=c(0,0.3),type="n")
grid(lty=1, col=gray(.9)); lines(aucmax.ts,type="b",pch=21:24)
traucmax=lm(aucmax.ts~time(aucmax.ts))
abline(traucmax,col="red")
boxplot(aucmax.ts~season(aucmax.ts),names=quarter.name,ylab="AUC Index",main="sumax(STP)")
modelaucmax=lm(aucmax.ts~time(aucmax.ts)+season(aucmax.ts)-1)
summary(modelaucmax)
############################################
#Plotting quarterly stp with reports
###########################################
dev.off()
pdf("reports-stpsum.pdf")
par(mfrow=c(2,2))
######################################
for (i in 1:length(stpsumbr)){
  for (j in 1:4){
    plot(crop(stpsumyr[[i]][[j]],XUS))
    plot(as.ppp(torpppyr[[28+i]][[j]]),add=TRUE)
    data("wrld_simpl")
    plot(wrld_simpl,add=TRUE)
    title(paste(quarter.name[j]," ",inyear+i-1,sep=""))
  }
}
pdf("reports-stpsumax.pdf")
par(mfrow=c(2,2))
######################################
for (i in 1:length(stpmaxbr)){
  for (j in 1:4){
    plot(crop(stpmaxbr[[i]][[j]],XUS))
    plot(as.ppp(torpppyr[[28+i]][[j]]),add=TRUE)
    data("wrld_simpl")
    plot(wrld_simpl,add=TRUE)
    title(paste(quarter.name[j]," ",inyear+i-1,sep=""))
  }
}
dev.off()
##############################################
#Quarterly medians of stpsumyr and stpsumaxyr
#############################################
qtsumstp=list()
qtsumstpM=list()
for(j in 1:4){
  qtsumstp[[j]]=stpsumyr[[1]][[j]]
  for(i in 2:length(stpsumbr)){
  qtsumstp[[j]]=addLayer(qtsumstp[[j]],stpsumyr[[i]][[j]]) 
    }
  qtsumstpM[[j]]=calc(qtsumstp[[j]],function(x){median(x)})
  }
qtsumaxstp=list()
qtsumaxstpM=list()
for(j in 1:4){
  qtsumaxstp[[j]]=stpsumaxyr[[1]][[j]]
  for(i in 2:length(stpmaxbr)){
    qtsumaxstp[[j]]=addLayer(qtsumaxstp[[j]],stpsumaxyr[[i]][[j]]) 
  }
qtsumaxstpM[[j]]=calc(qtsumaxstp[[j]],function(x){median(x)}) 
}
#########Saving raster files for Victor######################
Fig1Q1=qtsumaxstpM[[1]]
Fig1Q2=qtsumaxstpM[[2]]
Fig1Q3=qtsumaxstpM[[3]]
Fig1Q4=qtsumaxstpM[[4]]
save(Fig1Q1,file="Fig1Q1.RData")
save(Fig1Q2,file="Fig1Q2.RData")
save(Fig1Q3,file="Fig1Q3.RData")
save(Fig1Q4,file="Fig1Q4.RData")
##############################################################
##############Median Maps#####################################
par(mfrow=c(2,2),mar=c(2,2,2,1))
for(j in 1:4){
plot(crop(qtsumstpM[[j]],XUS))
data(wrld_simpl)
#plot(wrld_simpl,add=TRUE)
plot(states(),add=TRUE)
title(paste("Median sum(STP) ",quarter.name[j],sep=""))
}
par(mfrow=c(2,2),mar=c(2,2,2,1))
for(j in 1:4){
  plot(crop(qtsumaxstpM[[j]],XUS),col=viridis(10),breaks=c(1,5,10,15,20,25,30,35,40,45,50))
  data("wrld_simpl")
  #plot(map("state", resolution=0, add=TRUE))
  map("state",add=TRUE)
  #ColorBar("r", pal='jet')
  #plot(wrld_simpl,add=TRUE)
  title(paste("Median sumax(STP) ",quarter.name[j],sep=""))
}
####HISTOGRAMS#############
########################
hval=list(4)
par(mfrow=c(2,2),mar=c(2,2,2,1))
for(j in 1:4){
  hist(crop(qtsumstpM[[j]],XUS),main="")
  title(paste("Median sum(STP) ,",quarter.name[j],sep=""))
}
#par(mfrow=c(2,2),mar=c(2,2,2,1))
par(mfrow=c(2,2))
for(j in 1:4){
  hval[[j]]=hist(crop(qtsumaxstpM[[j]],XUS),main="",plot=FALSE)
  hval[[j]]$counts[hval[[j]]$counts==0]=1
}
#Using log scale for counts
for (j in 1:4){
  barplot(hval[[j]]$counts,names.arg=hval[[j]]$breaks[-1],log="y",
          space=0,xlab="Median sumax(STP)",ylab="Counts",cex.names=0.9,
          cex.lab=1.1,mgp=c(1,0,0),tcl=-0.1, pt.cex = 3,cex=1.3)
  #title(paste("Median sumax(STP) ",main="",quarter.name[j],sep=""))
legend("top",quarter.name[j],bty="n",pt.cex = 3,cex=1.3)
}
############################################
####Delete stpmaxbr and stpsumbr#############
rm(stpmaxbr)
rm(stpsumbr)
#############################################
#Analysis of Point pattern Dependence on stp 
#############################################
#Year 1987
###########################################
#Plot reports vs sumax stp quantiles
#Read US borderline and States shapefiles
XUS=st_read("US_outline.shp")
XUSstates=st_read("states.shp")
dev.off()
par(mfrow=c(2,2),mar=c(2,2,2,1))
#Initialize lists
VNc1987=list()
b1987=list()
Zcut1987=list()
Vt1987=list()
Qc1987=list()
Qt1987=list()
lam1987=list()
roc1987=list()
auc1987=list()
mod1987=list()
for (j in 1:4){
  VNc1987[[j]]=crop(stpsumaxyr[[9]][[j]],XUS)
   #VNc1987[[j]]=crop(stpsumaxyr[[9]][[j]],XUS) + 1 
   #b1987[[j]]=quantile(log(VNc1987[[j]]),probs=c(0.,.45,.65,.85,1.),type=2)
   b1987[[j]]=quantile(VNc1987[[j]],probs=c(.50,.60,.70,.80,.9,1.),type=2)
   Zcut1987[[j]]=cut(VNc1987[[j]],breaks=b1987[[j]])
   Vt1987[[j]]=tess(image=Zcut1987[[j]])
   plot(Vt1987[[j]],main=paste(quarter.name[j]),axes=TRUE,xlim=c(-125.3295, -66.53295),ylim=c(20, 55),
        ribside="right",ribsep=.01,ribscale=0.3,ribwid=0.015,ribargs=list(cex.axis=0.7,mgp=c(1,0,0)),
        col=topo.colors,cex.names=0.9,
        cex.lab=2,mgp=c(2,1,0),tcl=-0.25)
   plot(torpppyr[[37]][[j]],add=TRUE,pch=c(1),cex=0.4,xlab="Longitude",ylab="Latitude")
   data("wrld_simpl")
   #plot(map("state", resolution=0, add=TRUE))
   map("state",add=TRUE)
   #plot(XUSstates$geometry,add=TRUE)
   #plot(states(),add=TRUE)
   #plot(wrld_simpl[wrld_simpl$NAME=="United States",],add=TRUE)
}
#################Saving files for Victor#########################
Fig3Q1=Vt1987[[1]]$image
Fig3Q2=Vt1987[[2]]$image
Fig3Q3=Vt1987[[3]]$image
Fig3Q4=Vt1987[[4]]$image
save(Fig3Q1,file="Fig3Q1.RData")
save(Fig3Q2,file="Fig3Q2.RData")
save(Fig3Q3,file="Fig3Q3.RData")
save(Fig3Q4,file="Fig3Q4.RData")
list1987=torpppyr[[37]]
save(list1987,file="reports-1987.RData")
#
#
###################################################################
#Plot Intensity vs quantiles
#######################################################################
dev.off()
par(mfrow=c(2,2),mar=c(2,2,2,1))
for (j in 1:4){
   Qc1987[[j]] = quadratcount(unmark(torpppyr[[37]][[j]]), tess=Vt1987[[j]])
   Qt1987[[j]] = quadrat.test(unmark(torpppyr[[37]][[j]]), tess=Vt1987[[j]])
   lam1987[[j]] = intensity(Qc1987[[j]])
   #barplot(lam1987[[j]], main=paste("Intensity vs. Sumax(STP)","Quarter",j,sep=" "),ylim=c(0,1.0))
   barplot(lam1987[[j]], main=paste(quarter.name[j]),ylim=c(0,1.2),
           space=0.01,xlab="sumax(STP) percentile groups",ylab=expression(lambda(s)),cex.names=1,
           cex.lab=1,mgp=c(1,0,0),tcl=-0.2,font.lab=2 )
}
#Intensity estimation as a function of sumax stp
dev.off()
#par(mfrow=c(2,2),mar=c(2,2,2,1))
par(mfrow=c(2,2),mai=c(0.65,0.7,0.6,0.1))
rho1987=list()
for (j in 1:4){
  rho1987[[j]]=rhohat(unmark(torpppyr[[37]][[j]]),as.im(VNc1987[[j]]),smoother="local",method="ratio") 
  plot(rho1987[[j]],xlab="",ylab="",main=paste(quarter.name[j]),ylim=c(0,2),
       cex.main=1.7,cex.axis=1.1,tcl=-0.2,cex.lab=1,legend=FALSE,font.lab=2,cex=1)
       ylab.text = expression(rho*(Z(u)))
       mtext(ylab.text,side=2, line =2.0,font.lab=2,cex=1)
       xlab.text="sumax(STP)"
       mtext(xlab.text,side=1, line =2.0,font.lab=2,cex=1)
}
dev.off()
#par(mfrow=c(2,2),mar=c(3,2.5,2.5,2))
par(mfrow=c(2,2),mai=c(0.65,0.7,0.6,0.2))
# Test based on the exact values of a covariate (Anderson-Darling test)
for (j in 1:4){
cdft=cdf.test(unmark(torpppyr[[37]][[j]]),as.im(VNc1987[[j]]), test="ad")
plot(cdf.test(unmark(torpppyr[[37]][[j]]),as.im(VNc1987[[j]]), test="ad"),xlab="sumax(STP)",
     main=paste(quarter.name[j], "P-Value AD test:", round(cdft$p.value,6),sep=" "),
     do.legend=FALSE,
     cex.main=1.2,cex.axis=1.5,ylab="Probability",tcl=-0.2,cex.lab=1.5,lwd = 2.2)
legend("bottomright",legend=c("observed","expected"),col=c("black","red"),lty=c(1,2),lwd = 2.2)
}
#Berman test
berman1987=list()
for (j in 1:4){
 berman1987[[j]]= berman.test(unmark(torpppyr[[37]][[j]]),as.im(VNc1987[[j]]),"Z2")
 print(berman1987[[j]])
}
#ROC curves and AUC
for (j in 1:4){
roc1987[[j]]=roc(unmark(torpppyr[[37]][[j]]),as.im(VNc1987[[j]]),high=FALSE)
auc1987[[j]]=auc(unmark(torpppyr[[37]][[j]]),as.im(VNc1987[[j]]),high=FALSE)
}
#####Inhomogeneus Poisson model fit#############
##########################################
#Xim=as.im(stpsumaxyr[[9]][[1]])
#mod19870=ppm(unmark(torpppyr[[37]][[1]])~1)
#mod19871=ppm(unmark(torpppyr[[37]][[1]])~Xim)
#mod19872=ppm(unmark(torpppyr[[37]][[1]])~Xim+I(Xim^2))
#anova(mod19871,mod19872,test="LR")
#plot(effectfun(mod19871, "Xim", se.fit=TRUE))
#plot(effectfun(mod19872, "Xim", se.fit=TRUE))
dev.new()
par(mfrow=c(2,2),mar=c(2,2,2,1))
for (j in 1:4){
Z=as.im(stpsumaxyr[[9]][[j]])
mod1987[[j]]=ppm(unmark(torpppyr[[37]][[j]])~Z)
summary(mod1987[[j]])
plot(effectfun(mod1987[[j]], "Z", se.fit=TRUE),main=paste(quarter.name[j]),ylim=c(0,10))
}
#############################################
#Analysis of Point pattern Dependence on stp 
############################################
#Year 2011
############################################
#Plot reports vs sumax stp quantiles
dev.new()
par(mfrow=c(2,2),mar=c(2,2,2,1))
VNc2011=list()
b2011=list()
Zcut2011=list()
Vt2011=list()
Qc2011=list()
Qt2011=list()
lam2011=list()
roc2011=list()
auc2011=list()
mod2011=list()
for (j in 1:4){
  VNc2011[[j]]=crop(stpsumaxyr[[33]][[j]],XUS) 
  b2011[[j]]=quantile(VNc2011[[j]],probs=c(.50,.60,.70,.8,.9,1.),type=2)
  Zcut2011[[j]]=cut(VNc2011[[j]],breaks=b2011[[j]])
  Vt2011[[j]]=tess(image=Zcut2011[[j]])
  plot(Vt2011[[j]],main=paste(quarter.name[j]),axes=TRUE,xlim=c(-125.3295, -66.53295),ylim=c(20, 55),
       ribside="right",ribsep=.01,ribscale=0.3,ribwid=0.015,ribargs=list(cex.axis=0.7,mgp=c(1,0,0)),
       col=topo.colors,cex.names=0.9,
       cex.lab=2,mgp=c(2,1,0),tcl=-0.25)
  plot(torpppyr[[61]][[j]],add=TRUE,pch=c(1),cex=0.5)
  data("wrld_simpl")
  #plot(map("state", resolution=0, add=TRUE))
  map("state",add=TRUE)
  #plot(states(),add=TRUE)
  #plot(wrld_simpl[wrld_simpl$NAME=="United States",],add=TRUE)
}
#####################Save Files for Victor#########################
Fig4Q1=Vt2011[[1]]$image
Fig4Q2=Vt2011[[2]]$image
Fig4Q3=Vt2011[[3]]$image
Fig4Q4=Vt2011[[4]]$image
save(Fig4Q1,file="Fig4Q1.RData")
save(Fig4Q2,file="Fig4Q2.RData")
save(Fig4Q3,file="Fig4Q3.RData")
save(Fig4Q4,file="Fig4Q4.RData")
save(Vt2011,file="sumaxstp-tess-2011.RData")
list2011=torpppyr[[61]]
save(unlist(list2011),file="reports-2011.RData")
###################################################################
#Plot Intensity vs quantiles
dev.new()
par(mfrow=c(2,2),mar=c(2,2,2,1))
for (j in 1:4){
  Qc2011[[j]] = quadratcount(unmark(torpppyr[[61]][[j]]), tess=Vt2011[[j]])
  Qt2011[[j]] = quadrat.test(unmark(torpppyr[[61]][[j]]), tess=Vt2011[[j]])
  lam2011[[j]] = intensity(Qc2011[[j]])
  #barplot(lam2011[[j]], main=paste("Intensity vs. Sumax(STP)","Quarter",j,sep=" "),ylim=c(0,2.0))
  barplot(lam2011[[j]], main=paste(quarter.name[j]),ylim=c(0,1.8),
          space=0.01,xlab="sumax(STP) percentile groups",ylab=expression(lambda(s)),cex.names=1,
          cex.lab=1,mgp=c(1,0,0),tcl=-0.2,font.lab=2 )
}

#Intensity estimation as a function of sumax stp
dev.new()
#par(mfrow=c(2,2),mar=c(2,2,2,1))
par(mfrow=c(2,2),mai=c(0.65,0.7,0.6,0.1))
rho2011=list()
for (j in 1:4){
  rho2011[[j]]=rhohat(unmark(torpppyr[[61]][[j]]),as.im(VNc2011[[j]]),smoother="local",method="ratio") 
  plot(rho2011[[j]],xlab="",ylab="",main=paste(quarter.name[j]),ylim=c(0,2.5),
       cex.main=1.7,cex.axis=1.1,tcl=-0.2,cex.lab=1,legend=FALSE,font.lab=2,cex=1)
  ylab.text = expression(rho*(Z(u)))
  mtext(ylab.text,side=2, line =2.0,font.lab=2,cex=1)
  xlab.text="sumax(STP)"
  mtext(xlab.text,side=1, line =2.0,font.lab=2,cex=1)
}
dev.off()
#par(mfrow=c(2,2),mar=c(2,2,2,1))
par(mfrow=c(2,2),mai=c(0.65,0.7,0.6,0.2))
# Test based on the exact values of a covariate (Anderson-Darling test)
for (j in 1:4){
  cdft=cdf.test(unmark(torpppyr[[61]][[j]]),as.im(VNc2011[[j]]), test="ad")
  plot(cdf.test(unmark(torpppyr[[61]][[j]]),as.im(VNc2011[[j]]), test="ad"),xlab="sumax(STP)",
       main=paste(quarter.name[j], "P-Value AD test:", round(cdft$p.value,6),sep=" "),
       do.legend=FALSE,
       cex.main=1.2,cex.axis=1.5,ylab="Probability",tcl=-0.2,cex.lab=1.5,lwd = 2.2)
  legend("bottomright",legend=c("observed","expected"),col=c("black","red"),lty=c(1,2),lwd = 2.2)
}
#Berman test
berman2011=list()
for (j in 1:4){
  berman2011[[j]]= berman.test(unmark(torpppyr[[61]][[j]]),as.im(VNc2011[[j]]),"Z2")
  print(berman2011[[j]])
}
#ROC curves and AUC
for (j in 1:4){
  roc2011[[j]]=roc(unmark(torpppyr[[61]][[j]]),as.im(VNc1987[[j]]),high=FALSE)
  auc2011[[j]]=auc(unmark(torpppyr[[61]][[j]]),as.im(VNc1987[[j]]),high=FALSE)
}
#####Inhomogeneus Poisson model fit#############
################################################
dev.off()
par(mfrow=c(2,2),mar=c(2,2,2,1))
for (j in 1:4){
  Z=as.im(stpsumaxyr[[33]][[j]])
  mod2011[[j]]=ppm(unmark(torpppyr[[61]][[j]])~Z)
  summary(mod2011[[j]])
  plot(effectfun(mod2011[[j]], "Z", se.fit=TRUE),main=paste("Quarter",j,sep=" "),ylim=c(0,15))
}
