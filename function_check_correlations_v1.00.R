check.correlations=function(dat,parallel=F,n.cores=4){
  classes.dat=sapply(dat,class)
  fact.vars=names(which(classes.dat=="factor" | classes.dat=="character"))
  cont.vars=names(which(classes.dat=="integer" | classes.dat=="numeric"))
  if(length(cont.vars)>1){
   cor.mat=cor(dat[,cont.vars],use="pairwise.complete.obs")}else{
   cor.mat=matrix(1,ncol=1,nrow=1)
   colnames(cor.mat)=cont.vars
   rownames(cor.mat)=cont.vars}
  if(length(fact.vars)>0){
   if(length(cont.vars)>0){
    lm.grid=expand.grid(list(fact.var=fact.vars,cont.var=cont.vars))
    r.estimates=cbind(lm.grid,apply(lm.grid,MARGIN=1,FUN=function(x){
        sqrt(summary(lm(dat[,x[2]]~factor(dat[,x[1]])))$r.sq)}))

    fact.cont.upper.right=matrix(NA,ncol=length(fact.vars),nrow=length(cont.vars))
    colnames(fact.cont.upper.right)=fact.vars;rownames(fact.cont.upper.right)=cont.vars

    fact.cont.lower.left=matrix(NA,ncol=length(cont.vars),nrow=length(fact.vars))
    colnames(fact.cont.lower.left)=cont.vars;rownames(fact.cont.lower.left)=fact.vars

    fact.fact.lower.right=matrix(NA,ncol=length(fact.vars),nrow=length(fact.vars))
    colnames(fact.fact.lower.right)=fact.vars;rownames(fact.fact.lower.right)=fact.vars

    out.cor.mat=rbind(cbind(cor.mat,fact.cont.upper.right),
                      cbind(fact.cont.lower.left,fact.fact.lower.right))

    # assign the estimated r values to the upper right and lower left corners
    for(r in 1:nrow(r.estimates)){
       # upper right
       col.index=which(colnames(out.cor.mat)==r.estimates$fact.var[r])
       row.index=which(rownames(out.cor.mat)==r.estimates$cont.var[r])
       out.cor.mat[row.index,col.index]=r.estimates[r,3]
       # lower left
       col.index=which(colnames(out.cor.mat)==r.estimates$cont.var[r])
       row.index=which(rownames(out.cor.mat)==r.estimates$fact.var[r])
       out.cor.mat[row.index,col.index]=r.estimates[r,3]
     }
    }else{
        fact.fact.lower.right=matrix(NA,ncol=length(fact.vars),nrow=length(fact.vars))
    colnames(fact.fact.lower.right)=fact.vars;rownames(fact.fact.lower.right)=fact.vars
    out.cor.mat=fact.fact.lower.right}

  # estimate r values for fact-fact combinations
  lm.grid=expand.grid(list(fact.var1=fact.vars,fact.var2=fact.vars))
  require(nnet)
  if(parallel==T){
   require(doSNOW)
   cl=makePSOCKcluster(n.cores)
   registerDoSNOW(cl)
   out.cor.dat<-foreach(r = 1:nrow(lm.grid),.packages=c('nnet'),.errorhandling='pass')%dopar%{
    var.1=as.character(lm.grid[r,1])
    var.2=as.character(lm.grid[r,2])
    dat.r=na.omit(dat[,c(var.1,var.2)])
    fit <- try(summary(multinom(dat.r[,var.1] ~ dat.r[,var.2],trace=F))$deviance,silent=T)
    null.fit=try(summary(multinom(dat[,var.1] ~ 1,trace=F))$deviance,silent=T)
    if(class(fit)!="try-error"){
       if(round(fit,4)==round(null.fit,4)){r.est=0}else{
      r.est=sqrt(1-(fit/null.fit))}
      c(var.1,var.2,r.est)}}
   stopCluster(cl)
   #registerDoSEQ()
   }else{
    out.cor.dat=list()
    for(r in 1:nrow(lm.grid)){
          var.1=as.character(lm.grid[r,1])
          var.2=as.character(lm.grid[r,2])
          dat.r=na.omit(dat[,c(var.1,var.2)])
          fit <- try(summary(multinom(dat.r[,var.1] ~ dat.r[,var.2],trace=F))$deviance,silent=T)
          null.fit=try(summary(multinom(dat[,var.1] ~ 1,trace=F))$deviance,silent=T)
          out=NA
          if(class(fit)!="try-error"){
           if(round(fit,4)==round(null.fit,4)){r.est=0}else{
                   r.est=sqrt(1-(fit/null.fit))}
                   out=c(var.1,var.2,r.est)}
      out.cor.dat=c(out.cor.dat,list(out))}
      }

    for(r in 1:length(out.cor.dat)){
       out.cor.mat[which(colnames(out.cor.mat)==out.cor.dat[[r]][1]),
                   which(rownames(out.cor.mat)==out.cor.dat[[r]][2])]=
                   as.numeric(out.cor.dat[[r]][3])}}else{out.cor.mat=cor.mat}
  return(out.cor.mat)
}