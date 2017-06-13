### 11 May 2017 ################################################################
# Gastropod reproduction case study.
# The analysis uses the full.subsets.function to fit a complete model set
# to explore temporal patterns in gonadal somatix index (GSI) in two species of
# gastropods. GSI is modelled as a Gamma distribution.
# Both lunar.date and month are entered into the model as continous cyclic smooths.
# Sex and Species are included as interaction terms with each other, and as
# interaction terms between the two smoothers.

source("C:/Users/rfisher/Dropbox/FSS_paper/FSSgam/function_full_subsets_gam_v1.11.r")
setwd("C:/Users/rfisher/Dropbox/FSS_paper/case_study_TMS_ACL")
dat=read.csv(file="gastropod_dataset.csv")
dim(dat)

str(dat)
# specify factors
dat$year=as.factor(dat$year)

# response variable exploration, what distribution?
hist(dat$GSI)
range(dat$GSI)
 # slightly skewed, continuous variable, does not include zero = gamma

 # specify predictors
cyclic.vars=c("lunar.date","month")
factor.vars=c("Sex","Species")
cont.vars=c("lunar.date","month")

require(mgcv)
require(MuMIn)
use.dat=dat
start.fit=gam(GSI~s(lunar.date,k=5,bs='cc'),
              family="Gamma",
              data=use.dat)
out.list=full.subsets.gam(use.dat=use.dat,
                         test.fit=start.fit,
                         pred.vars.cont=cont.vars,
                         pred.vars.fact=factor.vars,
                         cyclic.vars=cyclic.vars,k=5,
                         parallel=T,factor.interactions=T,size=4)
names(out.list)

write.csv(out.list$predictor.correlations,"predictor_correlations.csv")
out.list$predictor.correlations

# examine the list of failed models
length(out.list$failed.models)
length(out.list$success.models)

# look at the model selection table
mod.table=out.list$mod.data.out
mod.table=mod.table[order(mod.table$AICc),]
head(mod.table)
write.csv(mod.table,"modfits.csv")

barplot(out.list$variable.importance$bic$variable.weights.per.mod,las=2,
        ylab="Relative variable importance")

write.csv(out.list$predictor.correlations,"predictor_correlations.csv")

# extract the best model
mod.table=mod.table[order(mod.table$AIC),]
head(mod.table)

best.model=out.list$success.models[[mod.table$modname[1]]]
plot(best.model,all.terms=T,pages=1)

gam.check(best.model)
summary(best.model)

# make a pretty plot
# get the data from the out.list object because this contains the interactions
model.dat=out.list$used.data
head(model.dat)
str(model.dat)

x.seq=seq(from=1,to=30,length=100)

best.mod.fact.vars=c("Sex","Species")
fact.grid=unique(model.dat[,best.mod.fact.vars])
sex.ltys=c(1,2)
sex.pchs=c(19,21)
sex.lvls=levels(model.dat$Sex)
spp.lvls=levels(model.dat$Species)
sex.cols=c("black","darkgrey")

# create some symbology and colour schemes for plotting
model.dat$pchs=16
for(a in 1:length(sex.lvls)){
    model.dat$cols[which(model.dat$Sex==sex.lvls[a])]=sex.cols[a]}

y.lim=c(0,ceiling(max(model.dat$GSI)))
lunar.lim=range(model.dat$lunar.date)
lunar.seq=seq(from=1,to=30,length=100)
month.lim=range(model.dat$month)
month.seq=seq(from=1,to=30,length=100)

month.labels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")

spp.labels=c("Monodonta labio", "Patelloida saccharina")
pdf(file="best_model.pdf",pointsize=12)
par(mfcol=c(2,2),mar=c(0,0.5,0.5,0.5),oma=c(5,4,0.5,2))
for(x in 1:length(spp.lvls)){
# lunar
spp=spp.lvls[x]
fact.grid.plot=fact.grid[grep(spp,fact.grid$Species),]
plot.dat=model.dat[grep(spp,model.dat$Species),]
plot(NA,ylim=y.lim,xlim=lunar.lim,ylab="",xlab="",main="",xpd=NA,yaxt="n",xaxt="n")
axis(side=2)
if(x==2){axis(side=1)}
for(r in 1:nrow(fact.grid.plot)){
 pred.vals=predict(best.model,newdata=data.frame(
                                      lunar.date=lunar.seq,
                                      month=mean(model.dat$month),
                                      Species=spp,
                                      Sex=fact.grid.plot[r,"Sex"]),
                                      type="response",se=T)
 lines(x.seq,pred.vals$fit,lwd=1.5,
       col=sex.cols[which(fact.grid.plot[r,"Sex"]==sex.lvls)])
 lines(x.seq,pred.vals$fit+1.96*pred.vals$se,lty=3,lwd=1.5,
       col=sex.cols[which(fact.grid.plot[r,"Sex"]==sex.lvls)])
 lines(x.seq,pred.vals$fit-1.96*pred.vals$se, lty=3,lwd=1.5,
       col=sex.cols[which(fact.grid.plot[r,"Sex"]==sex.lvls)])}
points(jitter(plot.dat$lunar.date),plot.dat$GSI,pch=16,col=plot.dat$cols)}
# month
for(x in 1:length(spp.lvls)){
spp=spp.lvls[x]
fact.grid.plot=fact.grid[grep(spp,fact.grid$Species),]
plot.dat=model.dat[grep(spp,model.dat$Species),]
plot(NA,ylim=y.lim,xlim=month.lim,ylab="",xlab="",main="",xpd=NA,yaxt="n",xaxt="n")
#if(x==1){axis(side=2)}
if(x==2){axis(side=1,at=c(2,4,6,8,10,12),labels=month.labels[c(2,4,6,8,10,12)])}
for(r in 1:nrow(fact.grid.plot)){
 pred.vals=predict(best.model,newdata=data.frame(
                                      lunar.date=mean(model.dat$lunar.date),
                                      month=month.seq,
                                      Species=spp,
                                      Sex=fact.grid.plot[r,"Sex"]),
                                      type="response",se=T)
 lines(x.seq,pred.vals$fit,lwd=1.5,
       col=sex.cols[which(fact.grid.plot[r,"Sex"]==sex.lvls)])
 lines(x.seq,pred.vals$fit+1.96*pred.vals$se,lty=3,lwd=1.5,
       col=sex.cols[which(fact.grid.plot[r,"Sex"]==sex.lvls)])
 lines(x.seq,pred.vals$fit-1.96*pred.vals$se, lty=3,lwd=1.5,
       col=sex.cols[which(fact.grid.plot[r,"Sex"]==sex.lvls)])}
points(jitter(plot.dat$month),plot.dat$GSI,pch=16,col=plot.dat$cols)}

mtext(side=1,text=c("Lunar date","Month of the year"),at=c(0.25,0.75),outer=T,line=2.5)
mtext(side=2,text="GSI",outer=T,line=2)
mtext(side=4,text=spp.labels,at=c(0.75,0.25),outer=T)
legend.dat=unique(model.dat[,c("Sex","cols")])
legend.dat=legend.dat[order(legend.dat$Sex),]

legend("bottom",legend=paste(legend.dat$Sex),bty="n",cex=1,
       inset=-0.275,ncol=2,lty=1,col=legend.dat$cols,pch=16,xpd=NA)
dev.off()


















