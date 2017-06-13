
# set arg vals
use.dat
test.fit
pred.vars.cont
pred.vars.fact=NA
size=3
k=5
cyclic.vars=NA
linear.vars=NA
smooth.interactions=pred.vars.fact
factor.interactions=F
cov.cutoff=0.28
bs.arg="'cr'"
null.terms=""
weight.vec=rep(1,nrow(use.dat))
max.models=500
parallel=F
n.cores=4
r2.type="r2.lm.est"
report.unique.r2=F
