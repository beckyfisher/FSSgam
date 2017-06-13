
# set arg vals
use.dat
test.fit
pred.vars.cont
pred.vars.fact=NA
cyclic.vars=NA
linear.vars=NA
size=3
cov.cutoff=0.28
k=5,parallel=F
n.cores=4
weight.vec=rep(1,nrow(use.dat))
null.terms="",bs.arg="'cr'"
smooth.interactions=pred.vars.fact
factor.interactions=F
max.models=500
r2.type="r2.lm.est"
report.unique.r2=F

