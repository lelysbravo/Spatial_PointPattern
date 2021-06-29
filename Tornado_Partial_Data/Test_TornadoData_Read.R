#### Updated version of maptools needed ----
install.packages("maptools", repos = "http://R-Forge.R-project.org")

# required libraries ----
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
# Read RDS files ----
STP1987<-readRDS("Tornado_Partial_Data/STPSUMAX1987.rds")
reports1987<-readRDS("Tornado_Partial_Data/reports-1987.rds")
STP2011<-readRDS("Tornado_Partial_Data/STPSUMAX2011.rds")
reports2011<-readRDS("Tornado_Partial_Data/reports-2011.rds")
#
# Plot STP files Year 1987 Quarter 2 ----
#plot(STPXXX,col=viridis(10),breaks=c(1,5,10,15,20,25,30,35,40,45,50))
plot(STP1987[[2]],col=terrain.colors(10),breaks=c(1,5,10,15,20,25,30,35,40,45,50),xlab="Longitude",ylab="Latitude")
data("wrld_simpl")
map("state",add=TRUE)
# Plot second quarter reports
plot(reports1987[[2]],add=TRUE,pch=c(3),cex=0.7)

# Plot STP files Year 2011 Quarter 3 ----
#plot(STPXXX,col=viridis(10),breaks=c(1,5,10,15,20,25,30,35,40,45,50))
plot(STP2011[[3]],col=terrain.colors(10),breaks=c(1,5,10,20,30,40,50,60,70,80),xlab="Longitude",ylab="Latitude")
data("wrld_simpl")
map("state",add=TRUE)
#Plot third quarter reports
plot(reports2011[[3]],add=TRUE,pch=c(3),cex=0.7)
##### Estimate rhohat
## 1987 Quarter 2 ----
dev.off()
rho1987Q2=rhohat(unmark(reports1987[[2]]),as.im(STP1987[[2]]),smoother="local",method="ratio") 
plot(rho1987Q2,xlab="",ylab="",main= " Year 1987 Quarter 2",ylim=c(0,2),
     cex.main=1.7,cex.axis=1.1,tcl=-0.2,cex.lab=1,legend=FALSE,font.lab=2,cex=1)
ylab.text = expression(rho*(Z(u)))
mtext(ylab.text,side=2, line =2.0,font.lab=2,cex=1)
xlab.text="sumax(STP)"
mtext(xlab.text,side=1, line =2.0,font.lab=2,cex=1)
## CDF test ----
cdft=cdf.test(unmark(reports1987[[2]]),as.im(STP1987[[2]]), test="ad")
plot(cdf.test(unmark(reports1987[[2]]),as.im(STP1987[[2]]), test="ad"),xlab="sumax(STP)",
     main=paste("Year 1987 Q2, P-Value AD test:", round(cdft$p.value,6),sep=" "),
     do.legend=FALSE,
     cex.main=1.2,cex.axis=1.5,ylab="Probability",tcl=-0.2,cex.lab=1.5,lwd = 2.2)
# Initialize lists ----
berman1987=list()
roc1987=list()
auc1987=list()
# Berman test for Year 1987 all Quarters ----
for (j in 1:4){
  berman1987[[j]]= berman.test(unmark(reports1987[[j]]),as.im(STP1987[[j]]),"Z2")
  print(berman1987[[j]])
}
berman1987[[1]]$p.value
berman1987[[1]]$statistic
# ROC curves and AUC for Year 1987 all Quarters ----
for (j in 1:4){
  roc1987[[j]]=roc(unmark(reports1987[[j]]),as.im(STP1987[[j]]),high=FALSE)
  auc1987[[j]]=auc(unmark(reports1987[[j]]),as.im(STP1987[[j]]),high=FALSE)
}
auc1987
str(roc1987[[1]])
# Plot Intensity vs quantiles
# dev.new()
# par(mfrow=c(2,2),mar=c(2,2,2,1))
# for (j in 1:4){
#   Qc2011[[j]] = quadratcount(unmark(torpppyr[[61]][[j]]), tess=Vt2011[[j]])
#   Qt2011[[j]] = quadrat.test(unmark(torpppyr[[61]][[j]]), tess=Vt2011[[j]])
#   lam2011[[j]] = intensity(Qc2011[[j]])
#   #barplot(lam2011[[j]], main=paste("Intensity vs. Sumax(STP)","Quarter",j,sep=" "),ylim=c(0,2.0))
#   barplot(lam2011[[j]], main=paste(quarter.name[j]),ylim=c(0,1.8),
#           space=0.01,xlab="sumax(STP) percentile groups",ylab=expression(lambda(s)),cex.names=1,
#           cex.lab=1,mgp=c(1,0,0),tcl=-0.2,font.lab=2 )
# }

perform_tests = function (reports, stp, year) { 
  # if cdf test needs VNC, add parameter to function definition
  # change cdf test call in this function, pass in VNC to this function
  
  # number of quarters in a year 
  n_quarters = length(reports)
  
  # data frame to hold the results of each hypothesis test run by quarter 
  results = data.frame(
    #bearman = I(vector('list', n_quarters)), # store bearman test object
    #roc = I(vector('list', n_quarters)), 
    year = year,
    quarter = seq_len(n_quarters),
    berman_Z2 = rep(NA, n_quarters),
    berman_pvalue = rep(NA, n_quarters),
    cdf_stat = rep(NA, n_quarters),
    cdf_pvalue = rep(NA, n_quarters),
    auc = rep(NA, n_quarters)
  )
  
  # perform hypothesis test for each quarter 
  for (j in seq_len(n_quarters)) {
    
    # compute hypothesis test values 
    berman_results = berman.test(unmark(reports[[j]]), as.im(stp[[j]]), "Z2")
    cdf_results = cdf.test(unmark(reports[[j]]),as.im(stp[[j]]), test="ad")
    
    # save hypothesis results into data frame 
    results$berman_Z2[[j]] = berman_results$statistic
    results$cdf_stat[[j]] = cdf_results$statistic
    results$berman_pvalue[[j]] = berman_results$p.value
    results$cdf_pvalue[[j]] = cdf_results$p.value
    
    # if we need to obtain ROC use this line of code 
    #results$roc[[j]] =
      #roc(unmark(reports[[j]]), as.im(stp[[j]]), high = FALSE)
    
    # compute the area under the curve 
    results$auc[[j]] =
      auc(unmark(reports[[j]]), as.im(stp[[j]]), high = FALSE)
  }
  
  # return data frame containing hypothesis testing results 
  results
}

# apply hypothesis tests to 1987 and 2011 
# combine data frames into one to make a table 
df_results = rbind(
  perform_tests(reports1987, STP1987, 1987),
  perform_tests(reports2011, STP2011, 2011)
)
df_results

#install.packages("gt")
#install.packages("gtsummary")
library("gt")
library("gtsummary")



df_results %>% 
  group_by(year) %>%
  gt() %>% 
  tab_header(
    title = "Spatial Analysis of Tornadoes",
    subtitle = "Based on 1987 and 2011 Data"
  ) %>% 
  tab_spanner(
    label = "Berman Test", 
    columns = matches("berman")
  ) %>%
  tab_spanner(
    label = "Cdf Test",
    columns = matches("cdf")
  ) %>%
  cols_label(
    berman_Z2 = "Statistic",
    berman_pvalue = "p-value",
    cdf_stat = "Statistic",
    cdf_pvalue = "p-value", 
    quarter = "Quarter", 
    auc = "AUC"
  ) %>%
  fmt_number(
    columns = c(berman_Z2, cdf_stat, auc),
    decimals = 3
  ) %>% 
  fmt(
    columns = c(berman_pvalue, cdf_pvalue), 
    fns = gtsummary::style_pvalue # style function from gtsummary package 
  ) 


