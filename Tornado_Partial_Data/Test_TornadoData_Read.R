## Testing load .RData files
#Load STP files Year 1987, Quarter 2
STPXXX=get(load("~/tornado_project/VNc1987Q2.RData"))
#Plot STP files
#plot(STPXXX,col=viridis(10),breaks=c(1,5,10,15,20,25,30,35,40,45,50))
plot(STPXXX,col=terrain.colors(10),breaks=c(1,5,10,15,20,25,30,35,40,45,50),xlab="Longitude",ylab="Latitude")
data("wrld_simpl")
map("state",add=TRUE)
#Load Tornado Reports
REPXXX=get(load("~/tornado_project/reports-1987.RData"))
plot(REPXXX[[2]],add=TRUE,pch=c(3),cex=0.7)

## Testing load .RData files
#Load STP files Year 2011, Quarter 3
STPYYY=get(load("~/tornado_project/VNc2011Q3.RData"))
#Plot STP files
#plot(STPXXX,col=viridis(10),breaks=c(1,5,10,15,20,25,30,35,40,45,50))
plot(STPXXX,col=terrain.colors(10),breaks=c(1,5,10,15,20,25,30,35,40,45,50),xlab="Longitude",ylab="Latitude")
data("wrld_simpl")
map("state",add=TRUE)
#Load Tornado Reports
REPYYY=get(load("~/tornado_project/reports-2011.RData"))
plot(REPYYY[[3]],add=TRUE,pch=c(3),cex=0.7)