# A simple function for full subsets multiple regression in ecology with R
# 
# R. Fisher
# S.K. Wilson
# S.M. Sin
# A.C. Lee
# T.J. Langlois

# Reproducible example for:
# Case Study 2: The role of large reef-associated predators in structuring adjacent soft-sediment communities

# A re-analysis of data presented in:
#   Langlois, T. J., M. J. Anderson, and R. C. Babcock. 2005. Reef-associated predators influence adjacent soft-sediment communities. Ecology 86: 1508–1519.

# Script information----
# This script is designed to work with long format data - where response variables are stacked one upon each other (see http://tidyr.tidyverse.org/)
# There are two random factors, Site and NTR location
# We have used a Tweedie error distribution to account for the high occurence of zero values in the dataset.
# We have implemented the ramdom effects and Tweedie error distribution using the mgcv() package

# librarys----
detach("package:plyr", unload=TRUE)#will error - don't worry
library(tidyr)
library(dplyr)
options(dplyr.width = Inf) #enables head() to display all coloums
library(mgcv)
library(MuMIn)
library(car)
library(doBy)
library(gplots)
library(RColorBrewer)
library(doParallel)
library(gamm4)
library(RCurl) #needed to download data from GitHub

rm(list=ls())
study<-"Clams"

# Source functions----
function_full_subsets_gam <- getURL("https://raw.githubusercontent.com/beckyfisher/FSSgam/master/function_full_subsets_gam_v1.11.R?token=AOSO6tZYAozKTAZ1Kt-aqlQIsiKuxONjks5ZZCtiwA%3D%3D", ssl.verifypeer = FALSE)
eval(parse(text = function_full_subsets_gam))

function_check_correlations <- getURL("https://raw.githubusercontent.com/beckyfisher/FSSgam/master/function_check_correlations_v1.00.R?token=AOSO6uxF2ON3UFyXj10uqm_N_94ZSEM3ks5ZZCyCwA%3D%3D", ssl.verifypeer = FALSE)
eval(parse(text = function_check_correlations))

# Bring in and format the data----
name<-"clams"

# Load the dataset
dat <-read.csv(text=getURL("https://raw.githubusercontent.com/beckyfisher/FSSgam/master/case_study2_dataset.csv?token=AOSO6uyYhat9-Era46nbjALQpTydsTskks5ZY3vhwA%3D%3D"))%>%
  rename(response=Abundance)%>%
  #   Transform variables
  mutate(sqrt.X4mm=sqrt(X4mm))%>%
  mutate(sqrt.X2mm=sqrt(X2mm))%>%
  mutate(sqrt.X1mm=sqrt(X1mm))%>%
  mutate(sqrt.X500um=sqrt(X500um))%>%
  na.omit()
head(dat,2)

# Set predictor variables---
pred.vars=c("depth","X4mm","X2mm","X1mm","X500um","X250um","X125um","X63um",
            "fetch","org","snapper","lobster") 

# predictor variables Removed at first pass---
# broad.Sponges and broad.Octocoral.Black and broad.Consolidated , "InPreds","BioTurb" are too rare

# Check for correalation of predictor variables- remove anything highly correlated (>0.95)---
round(cor(dat[,pred.vars]),2)
# nothing is highly correlated 

pdf(file=paste(name,"predictor_plots.pdf",sep = "_"),onefile=T)
for(p in 1:length(pred.vars)){
  par(mfrow=c(2,1))
  plot.dat=dat
  hist(plot.dat[,pred.vars[p]],main=pred.vars[p])
  plot(plot.dat[,pred.vars[p]])
}
dev.off()

# Review of individual predictors - we have to make sure they have an even distribution---
#If the data are squewed to low numbers try sqrt>log or if squewed to high numbers try ^2 of ^3
# Decided that X4mm, X2mm, X1mm and X500um needed a sqrt transformation
#Decided Depth, x63um, InPreds and BioTurb were not informative variables. 

# # Re-set the predictors for modeling----
pred.vars=c("sqrt.X4mm","sqrt.X2mm","sqrt.X1mm","sqrt.X500um",
            "fetch","org","snapper","lobster") 

# Check to make sure Response vector has not more than 80% zeros----
unique.vars=unique(as.character(dat$Taxa))
unique.vars.use=character()
for(i in 1:length(unique.vars)){
  temp.dat=dat[which(dat$Taxa==unique.vars[i]),]
  if(length(which(temp.dat$response==0))/nrow(temp.dat)<0.8){
    unique.vars.use=c(unique.vars.use,unique.vars[i])}
}
unique.vars.use     
write.csv(unique.vars.use,file=paste(name,"unique.vars.use.csv",sep = "_"))

# Run the full subset model selection----
resp.vars=unique.vars.use
use.dat=dat
factor.vars=c("Status")# Status as a Factor with two levels
out.all=list()
var.imp=list()

# Loop through the FSS function for each Taxa----
for(i in 1:length(resp.vars)){
  use.dat=dat[which(dat$Taxa==resp.vars[i]),]
  
  Model1=gam(response~s(lobster,k=3,bs='cr')+ s(Location,Site,bs="re"),
             family=tw(),  data=use.dat)
  out.list=full.subsets.gam(use.dat=use.dat,
                            test.fit=Model1,
                            pred.vars.cont=pred.vars,
                            pred.vars.fact=factor.vars,
                            linear.vars="Distance",
                            k=3,
                            null.terms="s(Location,Site,bs='re')",
                            max.models=600,
                            parallel=T)  
  names(out.list)
  
  out.list$failed.models # examine the list of failed models
  mod.table=out.list$mod.data.out  # look at the model selection table
  mod.table=mod.table[order(mod.table$AICc),]
  mod.table$cumsum.wi=cumsum(mod.table$wi.AICc)
  out.i=mod.table[which(mod.table$delta.AICc<=3),]
  out.all=c(out.all,list(out.i))
  # var.imp=c(var.imp,list(out.list$variable.importance$aic$variable.weights.raw)) #Either raw importance score
  var.imp=c(var.imp,list(out.list$variable.importance$aic$variable.weights.r2.scaled)) #Or importance score weighted by r2
  
  # plot the best models
  for(m in 1:nrow(out.i)){
   best.model.name=as.character(out.i$modname[m])

   png(file=paste(name,m,resp.vars[i],"mod_fits.png",sep="_"))
   if(best.model.name!="null"){
    par(mfrow=c(3,1),mar=c(9,4,3,1))
    best.model=out.list$success.models[[best.model.name]]
    plot(best.model,all.terms=T,pages=1,residuals=T,pch=16)
    mtext(side=2,text=resp.vars[i],outer=F)}  
   dev.off()
  }
}

# Model fits and importance---
names(out.all)=resp.vars
names(var.imp)=resp.vars
all.mod.fits=do.call("rbind",out.all)
all.var.imp=do.call("rbind",var.imp)
write.csv(all.mod.fits,file=paste(name,"all.mod.fits.csv",sep="_"))
write.csv(all.var.imp,file=paste(name,"all.var.imp.csv",sep="_"))

# Generic importance plots-
pdf(file=paste(name,"var_importance_heatmap.pdf",sep="_"),onefile=T)
heatmap.2(all.var.imp,notecex=0.4,  dendrogram ="none",
          col=colorRampPalette(c("white","yellow","red"))(10),
          trace="none",key.title = "",keysize=2,
          notecol="black",key=T,
          sepcolor = "black",margins=c(12,8), lhei=c(4,15),Rowv=FALSE,Colv=FALSE)
dev.off()
# END OF MODEL---