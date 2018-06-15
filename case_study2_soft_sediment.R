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

# Part 1-FSS modeling----
# This script is designed to work with long format data - where response variables are stacked one upon each other (see http://tidyr.tidyverse.org/)
# There are two random factors, Site and NTR location
# We have used a Tweedie error distribution to account for the high occurence of zero values in the dataset.
# We have implemented the ramdom effects and Tweedie error distribution using the mgcv() package

# Part 2 - custom plot of importance scores----
# using ggplot2()

# Part 3 - plots of the most parsimonious models----
# here we use plots of the raw response variables and fitted relationships - to allow for the plotting of interactions between continous predictor variables and factors with levels again using ggplot2()

# Part 1-FSS modeling----

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
  var.imp=c(var.imp,list(out.list$variable.importance$aic$variable.weights.raw)) #Or importance score weighted by r2
  
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


# Part 2 - custom plot of importance scores----


# Load the dataset
dat.taxa <-read.csv(text=getURL("https://raw.githubusercontent.com/beckyfisher/FSSgam/master/case_study2_all.var.imp.csv?token=AOSO6ma3MdgxTWJgICEtKgUVUGiZkRW0ks5ZbagowA%3D%3D"))%>%
  rename(resp.var=X)%>%
  gather(key=predictor,value=importance,2:ncol(.))
head(dat.taxa,5)


# Plotting defaults----
library(ggplot2)
# Theme-
Theme1 <-
  theme( # use theme_get() to see available options
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill="white"),
    legend.key = element_blank(), # switch off the rectangle around symbols in the legend
    legend.text = element_text(size=8),
    legend.title = element_text(size=8, face="bold"),
    legend.position = "top",
    legend.direction="horizontal",
    text=element_text(size=10),
    strip.text.y = element_text(size = 10,angle = 0),
    axis.title.x=element_text(vjust=0.3, size=10),
    axis.title.y=element_text(vjust=0.6, angle=90, size=10),
    axis.text.x=element_text(size=10,angle = 90, hjust=1,vjust=0.5),
    axis.text.y=element_text(size=10,face="italic"),
    axis.line.x=element_line(colour="black", size=0.5,linetype='solid'),
    axis.line.y=element_line(colour="black", size=0.5,linetype='solid'),
    strip.background = element_blank())


# colour ramps-
re <- colorRampPalette(c("mistyrose", "red2","darkred"))(200)

# Labels-
legend_title<-"Importance"

# Annotations-
dat.taxa.label<-dat.taxa%>%
  mutate(label=NA)%>%
  mutate(label=ifelse(predictor=="Distance"&resp.var=="BDS","X",ifelse(predictor=="Status"&resp.var=="BDS","X",ifelse(predictor=="sqrt.X500um"&resp.var=="BDS","X",label))))%>%
  mutate(label=ifelse(predictor=="lobster"&resp.var=="BMS","X",label))%>%
  mutate(label=ifelse(predictor=="sqrt.X4mm"&resp.var=="CPN","X",ifelse(predictor=="lobster"&resp.var=="CPN","X",label)))
head(dat.taxa.label,2)

# Plot gg.importance.scores ----
gg.importance.scores <- ggplot(dat.taxa.label, aes(x=predictor,y=resp.var,fill=importance))+
  geom_tile(show.legend=T) +
  scale_fill_gradientn(legend_title,colours=c("white", re), na.value = "grey98",
                       limits = c(0, max(dat.taxa.label$importance)))+
  scale_x_discrete(limits=c("Distance",
                            "Status",
                            "lobster",
                            "snapper",
                            "fetch",
                            "org",
                            "sqrt.X4mm",
                            "sqrt.X2mm",
                            "sqrt.X1mm",
                            "sqrt.X500um"),
               labels=c(
                 "Distance",
                 "Status",
                 "Lobster",
                 "Snapper",
                 "Fetch (km)",
                 "Organic content",
                 "Grain size: 4mm",
                 "            2mm",
                 "            1mm",
                 "            500um"
               ))+
scale_y_discrete(limits = c("CPN",
                            "BMS",
                            "BDS"),
                 labels=c("P. novizelandiae",
                          "M. striata",
                          "D. subrosea"))+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  Theme1+
  geom_text(aes(label=label))
gg.importance.scores

ggsave("Langlois.importance.scores.pdf",width = 15, height = 7,units = "cm")


# Part 3 - plots of the most parsimonious models----

### now  make a nice plot of the most interesting models-----
library(gridExtra)
library(grid)
# Theme-
Theme1 <-
  theme( # use theme_get() to see available options
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # legend.background = element_rect(fill="white"),
    legend.background = element_blank(),
    legend.key = element_blank(), # switch off the rectangle around symbols in the legend
    legend.text = element_text(size=15),
    legend.title = element_blank(),
    legend.position = c(0.2, 0.8),
    text=element_text(size=15),
    strip.text.y = element_text(size = 15,angle = 0),
    axis.title.x=element_text(vjust=0.3, size=15),
    axis.title.y=element_text(vjust=0.6, angle=90, size=15),
    axis.text.x=element_text(size=15),
    axis.text.y=element_text(size=15),
    axis.line.x=element_line(colour="black", size=0.5,linetype='solid'),
    axis.line.y=element_line(colour="black", size=0.5,linetype='solid'),
    strip.background = element_blank())

# Bring in and format the data----
name<-"clams"

dat <-read.csv(text=getURL("https://raw.githubusercontent.com/beckyfisher/FSSgam/master/case_study2_dataset.csv?token=AOSO6uyYhat9-Era46nbjALQpTydsTskks5ZY3vhwA%3D%3D"))%>%
  rename(response=Abundance)%>%
  #   Transform variables
  mutate(sqrt.X4mm=sqrt(X4mm))%>%
  mutate(sqrt.X2mm=sqrt(X2mm))%>%
  mutate(sqrt.X1mm=sqrt(X1mm))%>%
  mutate(sqrt.X500um=sqrt(X500um))%>%
  mutate(distance=as.numeric(as.character(Distance)))%>%
  na.omit()
head(dat,2)

# Manually make the most parsimonious GAM models for each taxa ----
# MODEL Bivalve.Dosina.subrosea 500um + distance x Status ----
dat.bds<-dat%>%filter(Taxa=="BDS")
gamm=gam(response~s(sqrt.X500um,k=3,bs='cr')+s(distance,k=1,bs='cr', by=Status)+ s(Location,Site,bs="re")+ Status, family=tw(),data=dat.bds)

# predict - status from MODEL Bivalve.Dosina.subrosea----
mod<-gamm
testdata <- expand.grid(distance=mean(mod$model$distance),
                        sqrt.X500um=mean(mod$model$sqrt.X500um),
                        Location=(mod$model$Location),
                        Site=(mod$model$Site),
                        Status = c("Fished","No-take"))
head(testdata)
fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)
#head(fits,2)
predicts.bds.status = testdata%>%data.frame(fits)%>%
  group_by(Status)%>% #only change here
  summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()
write.csv(predicts.bds.status,"predicts.csv") #there is some BUG in dplyr - that this fixes
predicts.bds.status<-read.csv("predicts.csv")
head(predicts.bds.status,5)

# predict - distance.x.status from MODEL Bivalve.Dosina.subrosea----
mod<-gamm
testdata <- expand.grid(distance=seq(min(dat$distance),max(dat$distance),length.out = 20),
                        sqrt.X500um=mean(mod$model$sqrt.X500um),
                        Location=(mod$model$Location),
                        Site=(mod$model$Site),
                        Status = c("Fished","No-take"))

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)
#head(fits,2)
predicts.bds.distance.x.status = testdata%>%data.frame(fits)%>%
  group_by(distance,Status)%>% #only change here
  # group_by(sqrt.X500um,Status)%>% #only change here
  summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()
write.csv(predicts.bds.distance.x.status,"predicts.csv") #there is some BUG in dplyr - that this fixes
predicts.bds.distance.x.status<-read.csv("predicts.csv")
head(predicts.bds.distance.x.status,5)

# predict 500um from MODEL Bivalve.Dosina.subrosea----
mod<-gamm
testdata <- expand.grid(sqrt.X500um=seq(min(dat$sqrt.X500um),max(dat$sqrt.X500um),length.out = 20),
                        distance=mean(mod$model$distance),
                        Location=(mod$model$Location),
                        Site=(mod$model$Site),
                        Status = c("Fished","No-take"))

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)
#head(fits,2)
predicts.bds.500um = testdata%>%data.frame(fits)%>%
  group_by(sqrt.X500um)%>% #only change here
  # group_by(sqrt.X500um,Status)%>% #only change here
  summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()
write.csv(predicts.bds.500um,"predicts.csv") #there is some BUG in dplyr - that this fixes
predicts.bds.500um<-read.csv("predicts.csv")
head(predicts.bds.500um,5)

# MODEL Bivalve.Myadora.striata  Lobster----
dat.bms<-dat%>%filter(Taxa=="BMS")
head(dat.bms,2)
gamm=gam(response~s(lobster,k=3,bs='cr')+ s(Location,Site,bs="re"), family=tw(),data=dat.bms)

# predict - lobster from model for Bivalve.Myadora.striata ----
mod<-gamm
testdata <- expand.grid(lobster=seq(min(dat$lobster),max(dat$lobster),length.out = 20),
                        Location=(mod$model$Location),
                        Site=(mod$model$Site),
                        Status = c("Fished","No-take"))

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)
#head(fits,2)
predicts.bms.lobster = testdata%>%data.frame(fits)%>%
  group_by(lobster)%>% #only change here
  # group_by(sqrt.X500um,Status)%>% #only change here
  summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()
write.csv(predicts.bms.lobster,"predicts.csv") #there is some BUG in dplyr - that this fixes
predicts.bms.lobster<-read.csv("predicts.csv")
head(predicts.bms.lobster,5)

# MODEL Decapod.P.novazelandiae 4mm + Lobster----
dat.cpn<-dat%>%filter(Taxa=="CPN")
head(dat.cpn,2)
gamm=gam(response~s(sqrt.X4mm,k=3,bs='cr')+s(lobster,k=3,bs='cr')+ s(Location,Site,bs="re"), family=tw(),data=dat.cpn)

# predict - sqrt.X4mm from model for Decapod.P.novazelandiae ----
mod<-gamm
testdata <- expand.grid(sqrt.X4mm=seq(min(dat$sqrt.X4mm),max(dat$sqrt.X4mm),length.out = 20),
                        lobster=mean(mod$model$lobster),
                        Location=(mod$model$Location),
                        Site=(mod$model$Site),
                        Status = c("Fished","No-take"))

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)
head(fits,2)
predicts.cpn.4mm = testdata%>%data.frame(fits)%>%
  group_by(sqrt.X4mm)%>% #only change here
  summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()
write.csv(predicts.cpn.4mm,"predicts.csv") #there is some BUG in dplyr - that this fixes
predicts.cpn.4mm<-read.csv("predicts.csv")
head(predicts.cpn.4mm,5)

# predict - lobster from model for Decapod.P.novazelandiae ----
mod<-gamm
testdata <- expand.grid(lobster=seq(min(dat$lobster),max(dat$lobster),length.out = 20),
                        sqrt.X4mm=mean(mod$model$sqrt.X4mm),
                        Location=(mod$model$Location),
                        Site=(mod$model$Site),
                        Status = c("Fished","No-take"))

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)
#head(fits,2)
predicts.cpn.lobster = testdata%>%data.frame(fits)%>%
  group_by(lobster)%>% #only change here
  # group_by(sqrt.X500um,Status)%>% #only change here
  summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()
write.csv(predicts.cpn.lobster,"predicts.csv") #there is some BUG in dplyr - that this fixes
predicts.cpn.lobster<-read.csv("predicts.csv")
head(predicts.cpn.lobster,5)

# PLOTS for Bivalve.Dosina.subrosea 500um + distance x Status ----
ggmod.bds.status<- ggplot(aes(x=Status,y=response,fill=Status,colour=Status), data=predicts.bds.status) +
  ylab(" ")+
  xlab('Status')+
  #   ggtitle(substitute(italic(name)))+
  scale_fill_manual(labels = c("Fished", "No-take"),values=c("red", "black"))+
  scale_colour_manual(labels = c("Fished", "No-take"),values=c("red", "black"))+
  scale_x_discrete(limits = rev(levels(predicts.bds.status$Status)))+
  geom_bar(stat = "identity")+
  geom_errorbar(aes(ymin = response-se.fit,ymax = response+se.fit),width = 0.5) +
  theme_classic()+
  Theme1+
  annotate("text", x = -Inf, y=Inf, label = "(a)",vjust = 1, hjust = -.1,size=5)+
  annotate("text", x = -Inf, y=Inf, label = "   Dosinia subrosea",vjust = 1, hjust = -.1,size=5,fontface="italic")
ggmod.bds.status

ggmod.bds.distance.x.status<- ggplot(aes(x=distance,y=response,colour=Status), data=dat.bds) +
  ylab(" ")+
  xlab('Distance (m)')+
  #   ggtitle(substitute(italic(name)))+
  scale_color_manual(labels = c("Fished", "No-take"),values=c("red", "black"))+
  geom_jitter(width = 0.25,height = 0,alpha=0.75, size=2,show.legend=FALSE)+
  # geom_point(alpha=0.75, size=2)+
  geom_line(data=predicts.bds.distance.x.status,show.legend=FALSE)+
  geom_line(data=predicts.bds.distance.x.status,aes(y=response - se.fit),linetype="dashed",show.legend=FALSE)+
  geom_line(data=predicts.bds.distance.x.status,aes(y=response + se.fit),linetype="dashed",show.legend=FALSE)+
  theme_classic()+
  Theme1+
  annotate("text", x = -Inf, y=Inf, label = "(b)",vjust = 1, hjust = -.1,size=5)
ggmod.bds.distance.x.status

ggmod.bds.500um<- ggplot() +
  ylab(" ")+
  xlab('Grain size: 500 um (sqrt)')+
#   ggtitle(substitute(italic(name)))+
  scale_color_manual(labels = c("Fished", "No-take"),values=c("red", "black"))+
  #   geom_jitter(width = 0.25,height = 0)+
  geom_point(data=dat.bds,aes(x=sqrt.X500um,y=response,colour=Status),  alpha=0.75, size=2,show.legend=FALSE)+
  geom_line(data=predicts.bds.500um,aes(x=sqrt.X500um,y=response),alpha=0.5)+
  geom_line(data=predicts.bds.500um,aes(x=sqrt.X500um,y=response - se.fit),linetype="dashed",alpha=0.5)+
  geom_line(data=predicts.bds.500um,aes(x=sqrt.X500um,y=response + se.fit),linetype="dashed",alpha=0.5)+
  theme_classic()+
  Theme1+
  annotate("text", x = -Inf, y=Inf, label = "(c)",vjust = 1, hjust = -.1,size=5)
ggmod.bds.500um

# PLOTS Bivalve M.striata lobster ----
ggmod.bms.lobster<- ggplot() +
  ylab("Abundance")+
  xlab(bquote('Density of legal lobster (no./25' *m^-2*')'))+
  scale_color_manual(labels = c("Fished", "SZ"),values=c("red", "black"))+
  geom_point(data=dat.bms,aes(x=lobster,y=response,colour=Status),  alpha=0.75, size=2,show.legend=FALSE)+
  geom_line(data=predicts.bms.lobster,aes(x=lobster,y=response),alpha=0.5)+
  geom_line(data=predicts.bms.lobster,aes(x=lobster,y=response - se.fit),linetype="dashed",alpha=0.5)+
  geom_line(data=predicts.bms.lobster,aes(x=lobster,y=response + se.fit),linetype="dashed",alpha=0.5)+
  theme_classic()+
  Theme1+
  annotate("text", x = -Inf, y=Inf, label = "(d)",vjust = 1, hjust = -.1,size=5)+
  annotate("text", x = -Inf, y=Inf, label = "   Myadora striata",vjust = 1, hjust = -.1,size=5,fontface="italic")+
  geom_blank(data=dat.bms,aes(x=lobster,y=response*1.05))#to nudge data off annotations
ggmod.bms.lobster

# PLOTS Decapod.P.novazelandiae 4mm + lobster ----
ggmod.cpn.lobster<- ggplot() +
  ylab(" ")+
  xlab(bquote('Density of legal lobster (no./25' *m^-2*')'))+
  scale_color_manual(labels = c("Fished", "SZ"),values=c("red", "black"))+
  geom_point(data=dat.cpn,aes(x=lobster,y=response,colour=Status),  alpha=0.75, size=2,show.legend=FALSE)+
  geom_line(data=predicts.cpn.lobster,aes(x=lobster,y=response),alpha=0.5)+
  geom_line(data=predicts.cpn.lobster,aes(x=lobster,y=response - se.fit),linetype="dashed",alpha=0.5)+
  geom_line(data=predicts.cpn.lobster,aes(x=lobster,y=response + se.fit),linetype="dashed",alpha=0.5)+
  theme_classic()+
  Theme1+
  annotate("text", x = -Inf, y=Inf, label = "(e)",vjust = 1, hjust = -.1,size=5)+
  annotate("text", x = -Inf, y=Inf, label = "  Pagurus novizelandiae",vjust = 1, hjust = -.1,size=5,fontface="italic")+
  geom_blank(data=dat.cpn,aes(x=lobster,y=response*1.05))#to nudge data off annotations
ggmod.cpn.lobster

ggmod.cpn.4mm<- ggplot() +
  ylab(" ")+
  xlab('Grain size: 4 mm (sqrt)')+
  scale_color_manual(labels = c("Fished", "No-take"),values=c("red", "black"))+
  geom_point(data=dat.cpn,aes(x=sqrt.X4mm,y=response,colour=Status),  alpha=0.75, size=2,show.legend=FALSE)+
  geom_line(data=predicts.cpn.4mm,aes(x=sqrt.X4mm,y=response),alpha=0.5)+
  geom_line(data=predicts.cpn.4mm,aes(x=sqrt.X4mm,y=response - se.fit),linetype="dashed",alpha=0.5)+
  geom_line(data=predicts.cpn.4mm,aes(x=sqrt.X4mm,y=response + se.fit),linetype="dashed",alpha=0.5)+
  theme_classic()+
  Theme1+
  annotate("text", x = -Inf, y=Inf, label = "(f)",vjust = 1, hjust = -.1,size=5)+
  annotate("text", x = -Inf, y=Inf, label = " ",vjust = 1, hjust = -.1,size=5,fontface="italic")
ggmod.cpn.4mm

# combined.plot using grid() and gridExtra()------
blank <- grid.rect(gp=gpar(col="white"))

# To see what they will look like use grid.arrange()
grid.arrange(ggmod.bds.status,ggmod.bds.distance.x.status,ggmod.bds.500um,
             ggmod.bms.lobster,blank,blank,
             ggmod.cpn.lobster,ggmod.cpn.4mm,blank,nrow=3,ncol=3)

# Use arrangeGrob ONLY - as we can pass this to ggsave! Note use of raw ggplot's
combine.plot<-arrangeGrob(ggmod.bds.status,ggmod.bds.distance.x.status,ggmod.bds.500um,
                          ggmod.bms.lobster,blank,blank,
                          ggmod.cpn.lobster,ggmod.cpn.4mm,blank,nrow=3,ncol=3)

ggsave(combine.plot,file="Langlois_gamm.plot.png", width = 30, height = 30,units = "cm")
























