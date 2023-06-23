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

# Source the package
devtools::install_github("beckyfisher/FSSgam_package")
library(FSSgam)
library(RCurl)

################################################################################
### Example showing use of uGamm to allow fitting with gamm4  ##################
# load data  coral data set
#dat <-read.csv(text=getURL("https://raw.githubusercontent.com/beckyfisher/FSSgam/master/extra_examples_coral_data.csv"))
dat <- read.csv("extra_examples_coral_data.csv")
colnames(dat)
head(dat)
str(dat)

cat.preds=c("Survey","bleach.pres","dredge.pres","dhw.fact")
null.vars=c("Site")
cont.preds=c("av.wave","Depth")

# get rid of NA's and unused columns
use.dat=na.omit(dat[,c(null.vars,cat.preds,cont.preds,"allcoral","totalpoints")])
use.dat$successes=use.dat$allcoral
use.dat$failures=use.dat$totalpoints-use.dat$allcoral
use.dat$trials=use.dat$totalpoints

#test.fit model for all coral, with total points as trials
require(MuMIn)
Model1=uGamm(cbind(successes,failures)~s(Depth,k=4,bs='cr'),
              family=binomial(), random=~(1|Site),
             data=use.dat,
             lme4=TRUE)

model.set=generate.model.set(use.dat=use.dat,
                          test.fit=Model1,
                          pred.vars.cont=cont.preds,
                          pred.vars.fact=cat.preds)
out.list=fit.model.set(model.set, parallel = TRUE, #r2.type = "dev", 
                       report.unique.r2 = TRUE)
# examine the output
names(out.list)
out.list$failed.models
length(out.list$success.models)
mod.table=out.list$mod.data.out
mod.table=mod.table[order(mod.table$AICc),]
head(mod.table)

# check the predictor correlation matrix
model.set$predictor.correlations

# now run the same thing using the non.linear correlation matrix
model.set=generate.model.set(use.dat=use.dat,
                          test.fit=Model1,
                          pred.vars.cont=cont.preds,
                          pred.vars.fact=cat.preds,
                          non.linear.correlations=TRUE)
model.set$predictor.correlations
out.list=fit.model.set(model.set)
mod.table=out.list$mod.data.out
mod.table=mod.table[order(mod.table$AICc),]
head(mod.table)

#--- now an example running across a range of response variables  ------------
resp.vars=c("Acropora.spp.","Turbinaria.spp.","Pocillopora.spp.","Porites.spp.")
# get rid of NA's and unused columns
use.dat=na.omit(dat[,c(null.vars,cat.preds,cont.preds,resp.vars,"totalpoints")])

out.all=list()
var.imp=list()
fss.all=list()
top.all=list()
i=1
pdf(file="mod_fits_all.pdf",onefile=T)
for(i in 1:length(resp.vars)){
 use.dat$response=use.dat[,resp.vars[i]]
 #test.fit model for the particular coral i, with total points as trials
 Model1=uGamm(cbind(use.dat$response,use.dat$totalpoints-use.dat$response)~s(Depth,k=4,bs='cr'),
              family=binomial(), random=~(1|Site),
             data=use.dat,
             lme4=TRUE)

 model.set=generate.model.set(use.dat=use.dat,
                          test.fit=Model1,
                          pred.vars.cont=cont.preds,
                          pred.vars.fact=cat.preds)
 out.list=fit.model.set(model.set)
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
  plot(best.model$gam,all.terms=T,pages=1,residuals=T,pch=16)
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

pdf(file="var_importance.pdf",height=5,width=7,pointsize=10)
heatmap.2(all.var.imp,notecex=0.4,  dendrogram ="none",
                     col=colorRampPalette(c("yellow","orange","red"))(30),
                     trace="none",key.title = "",keysize=2,
                     notecol="black",key=T,
                     sepcolor = "black",margins=c(12,14), lhei=c(3,10),lwid=c(3,10),
                     Rowv=FALSE,Colv=FALSE)
dev.off()

write.csv(all.mod.fits,"all_model_fits.csv")
write.csv(top.mod.fits,"top_model_fits.csv")

################################################################################









