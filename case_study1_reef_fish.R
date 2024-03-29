#    Copyright 2020 Australian Institute of Marine Science
#
#    Licensed under the Apache License, Version 2.0 (the "License");
#    you may not use this file except in compliance with the License.
#    You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#    Unless required by applicable law or agreed to in writing, software
#    distributed under the License is distributed on an "AS IS" BASIS,
#    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#    See the License for the specific language governing permissions and
#    limitations under the License.

# A simple function for full subsets multiple regression in ecology with R
#
# R. Fisher
# S.K. Wilson
# S.M. Sin
# A.C. Lee
# T.J. Langlois

# Reproducible example for:
# Case Study 1: The relative influence of management and habitat on fish abundance and biomass

# note this example was updated on the 11th Oct 2018 to demonstrate useage of the replacement functions
# generate.model.set and fit.model.set that have now superced full.subsets.gam in package FSSgam
# Between them these functions carry out the same analysis, take the same arguments and return the same
# outputs as full.subsets.gam with the only difference being that the model set generation and model
# fitting procedures are separated into two steps. This was done to make the function easier to use,
# because the model set can be interrogated, along with the correlation matrix of the predictors before model
# fitting is even attempted.

# Source functions----
devtools::install_github("beckyfisher/FSSgam_package")
library(FSSgam)
library(RCurl)
# load data
download.file("https://raw.githubusercontent.com/beckyfisher/FSSgam/master/case_study1_dataset.csv", destfile = "c1dat.csv")
dat <-read.csv("c1dat.csv")
dim(dat)

require(mgcv)
require(MuMIn)
require(doParallel)
require(plyr)

colnames(dat)

dat$SQRTSA=dat$SA
dat$sqrt.rug=sqrt(dat$rugosity)
dat$sqrtLC=sqrt(dat$LC)
dat$sqrtHC=sqrt(dat$HC)
dat$sqrtMacro=sqrt(dat$macro)

cat.preds="ZONE"
null.vars=c("depth","site","SQRTSA") # use as random effect and null model
cont.preds=c("complexity","sqrt.rug","sqrtLC","sqrtHC","sqrtMacro",
             "SCORE1","SCORE2") # use as continuous predictors.
cor(dat[,cont.preds])
# have a look at the distribution of the continuous predictors
pdf(file="pred_vars.pdf",onefile=T)
for(p in 1:length(cont.preds)){
par(mfrow=c(2,1))
 hist(dat[,cont.preds[p]],main=cont.preds[p])
 plot(jitter(dat[,cont.preds[p]]))
 }
dev.off()

# remove extreme outliers
dat$Piscivore.abundance[which(dat$Piscivore.abundance>150)]=NA
dat$Piscivore.biomass[which(dat$Piscivore.biomass>40000)]=NA
dat$Invertivore.biomass[which(dat$Invertivore.biomass>40000)]=NA
dat$site <- as.factor(dat$site)

resp.vars.fams=list("Herbivore.abundance"=tw(),
                    "Invertivore.abundance"=tw(),
                    "Piscivore.abundance"=tw(),
                    "Planktivore.abundance"=tw(),
                    "Herbivore.biomass"=tw(),
                    "Invertivore.biomass"=tw(),
                    "Piscivore.biomass"=tw(),
                    "Planktivore.biomass"=tw())
resp.vars=names(resp.vars.fams)

# take a look at the response variables
pdf(file="resp_vars.pdf",onefile=T)
for(r in 1:length(resp.vars)){
par(mfrow=c(2,1))
 hist(dat[,resp.vars[r]],main=resp.vars[r])
 plot(jitter(dat[,resp.vars[r]]))
 }
dev.off()

### now fit the models ---------------------------------------------------------
i=1
out.all=list()
var.imp=list()
fss.all=list()
top.all=list()
pdf(file="mod_fits_functional_biomass.pdf",onefile=T)
for(i in 1:length(resp.vars)){
 use.dat=na.omit(dat[,c(null.vars,cont.preds,cat.preds,resp.vars[i])])
 use.dat$response=use.dat[,resp.vars[i]]
 Model1=gam(response~s(complexity,k=4,bs='cr')+
                    +s(SQRTSA,bs='cr',k=4) +s(site,bs="re"),
                    family=tw(),
                    data=use.dat)

 model.set=generate.model.set(use.dat=use.dat,max.predictors=2,   # limit size here because null model already complex
                             test.fit=Model1,k=3,
                             pred.vars.cont=cont.preds,
                             pred.vars.fact=cat.preds,
                             null.terms="s(SQRTSA,bs='cr',k=3)+s(site,bs='re')+s(depth,bs='cr',k=3)")

 out.list=fit.model.set(model.set, parallel = TRUE)
 ?d#names(out.list)
 # examine the list of failed models
 #out.list$failed.models
 #out.list$success.models
 fss.all=c(fss.all,list(out.list))
 mod.table=out.list$mod.data.out
 mod.table=mod.table[order(mod.table$AICc),]
 out.i=mod.table
 out.all=c(out.all,list(out.i))
 var.imp=c(var.imp,list(out.list$variable.importance$aic$variable.weights.raw))
 all.less.2AICc=mod.table[which(mod.table$delta.AICc<2),]
 top.all=c(top.all,list(all.less.2AICc))

 # plot the all best models
 par(oma=c(1,1,4,1))
 for(r in 1:nrow(all.less.2AICc)){
 best.model.name=as.character(all.less.2AICc$modname[r])
 best.model=out.list$success.models[[best.model.name]]
 if(best.model.name!="null"){
  plot(best.model,all.terms=T,pages=1,residuals=T,pch=16)
  mtext(side=3,text=resp.vars[i],outer=T)}
 }
}
dev.off()

names(out.all)=resp.vars
names(var.imp)=resp.vars
names(top.all)=resp.vars
names(fss.all)=resp.vars

all.mod.fits=do.call("rbind",out.all)
all.var.imp=do.call("rbind",var.imp)
top.mod.fits=do.call("rbind",top.all)

require(car)
require(doBy)
require(gplots)
require(RColorBrewer)

pdf(file="var_importance_heatmap_functional_biomass.pdf",height=5,width=7,pointsize=10)
heatmap.2(all.var.imp,notecex=0.4,  dendrogram ="none",
                     col=colorRampPalette(c("white","yellow","orange","red"))(30),
                     trace="none",key.title = "",keysize=2,
                     notecol="black",key=T,
                     sepcolor = "black",margins=c(12,14), lhei=c(3,10),lwid=c(3,10),
                     Rowv=FALSE,Colv=FALSE)
dev.off()

write.csv(all.mod.fits[,-2],"all_model_fits_functional_biomass.csv")
write.csv(top.mod.fits[,-2],"top_model_fits_functional_biomass.csv")
write.csv(model.set$predictor.correlations,"predictor_correlations.csv")

#### pretty plots of best models -----------------------------------------------
zones=levels(dat$ZONE)
pdf("best_top_model_quick_plots.pdf",height=8,width=7,pointsize=12)
par(mfcol=c(4,2),mar=c(4,4,0.5,0.5),oma=c(2,0.5,0.5,0.5),bty="l")
for(r in 1:length(resp.vars)){
      tab.r=out.all[[resp.vars[r]]]
      top.mods.r=tab.r[1,]
      mod.r.m=as.character(top.mods.r[1,"modname"])
      mod.m=fss.all[[resp.vars[r]]]$success.models[[mod.r.m]]
      mod.vars=unique(unlist(strsplit(unlist(strsplit(mod.r.m,split="+",fixed=T)),
                       split=".by.")))
      # which continuous predictor is the variable included?
      plot.var=as.character(na.omit(mod.vars[match(cont.preds,mod.vars)]))
      # plot that variables, with symbol colours for zone
      plot(dat[,plot.var],dat[,resp.vars[r]],pch=16,
         ylab=resp.vars[r],xlab=plot.var,col=dat$ZONE)
      legend("topleft",legend=paste("(",LETTERS[r],")",sep=""),
             bty="n")
      range.v=range(dat[,plot.var])
      seq.v=seq(range.v[1],range.v[2],length=20)
      newdat.list=list(seq.v,# across the range of the included variable
                       mean(use.dat$depth), # for a median depth
                       mean(use.dat$SQRTSA),# for a median SQRTSA
                       "MANGROVE", # pick the first site, except don't predict on
                               # this by setting terms=c(plot.var,"ZONE")
                       zones)  # for each zone
      names(newdat.list)=c(plot.var,"depth","SQRTSA","site","ZONE")
      pred.vals=predict(mod.m,newdata=expand.grid(newdat.list),
                     type="response",se=T,exclude=c("site","SQRTSA","depth"))
      for(z in 1:length(zones)){
       zone.index=which(expand.grid(newdat.list)$ZONE==zones[z])
       lines(seq.v,pred.vals$fit[zone.index],col=z)
       lines(seq.v,pred.vals$fit[zone.index]+pred.vals$se[zone.index]*1.96,lty=3,col=z)
       lines(seq.v,pred.vals$fit[zone.index]-pred.vals$se[zone.index]*1.96,lty=3,col=z)}
}
legend("bottom",legend= zones,bty="n",ncol=2,col=c(1,2),pch=c(16,16),
   inset=-0.61,xpd=NA,cex=.8)
dev.off()

