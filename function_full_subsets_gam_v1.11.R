full.subsets.gam=function(use.dat,
                          test.fit,
                          pred.vars.cont=NA,
                          pred.vars.fact=NA,
                          cyclic.vars=NA,
                          linear.vars=NA,
                          factor.smooth.interactions=pred.vars.fact,
                          factor.factor.interactions=F,
                          smooth.smooth.interactions=F,
                          cov.cutoff=0.28,
                          cor.matrix=NA,
                          max.predictors=3,
                          k=5,
                          bs.arg="'cr'",
                          null.terms="",
                          max.models=500,
                          parallel=F,
                          n.cores=4,
                          r2.type="r2.lm.est",
                          report.unique.r2=F,
                          factor.interactions="previous.arg",
                          smooth.interactions="previous.arg",
                          size="previous.arg"){

  # manage previous version arguments
  if(factor.interactions!="previous.arg"){
     factor.factor.interactions=factor.interactions
     warning('Argument factor.interactions has been replaced with factor.factor.interactions.
              Please update your code as usage of factor.interactions will not be supported in
              future versions.')
     }
  if(is.na(smooth.interactions)==T){
     factor.smooth.interactions=smooth.interactions
     warning('Argument smooth.interactions has been replaced with factor.smooth.interactions.
              Please update your code as usage of smooth.interactions will not be supported in
              future versions.')}else{if(smooth.interactions!="previous.arg"){
     factor.smooth.interactions=smooth.interactions
     warning('Argument smooth.interactions has been replaced with factor.smooth.interactions.
              Please update your code as usage of smooth.interactions will not be supported in
              future versions.')
     }}
  if(size!="previous.arg"){
     size=max.predictors
     warning('Argument size has been replaced with max.predictors.
              Please update your code as usage of size will not be supported in
              future versions.')
     }

  # make an "intercept" term for the null model
  use.dat$intercept=1
  interaction.terms=NA
  linear.interaction.terms=NA
  all.predictors=unique(na.omit(c(pred.vars.cont,pred.vars.fact,linear.vars)))
  included.vars=all.predictors


  # check the null model will fit
  if(nchar(null.terms)>0){
    null.formula=as.formula(paste("~ intercept-1",null.terms,sep="+"))}else{
    null.formula=as.formula("~ intercept-1")}

  if(length(grep("dsm",class(test.fit)))>0){
    null.formula=as.formula(paste("~1",null.terms,sep="+"))
    null.fit=try(update(test.fit,formula=null.formula),silent=T)
  }else{
    null.fit=try(update(test.fit,formula=null.formula,data=use.dat),silent=T)}

  if(class(null.fit)[1]=="try-error"){
        stop(paste("Null model not successfully fitted, please check your inputs.
                   If there are no random effects try using 'gam' instead of 'uGamm' 
                   in your test.fit model call.",
                   " ",
                   "The following error message was provided:  ",
                   "  ",
                   null.fit, ""))}

  # check for missing predictor values
  if(max(is.na(use.dat[,all.predictors]))==1){
        stop("Predictor variables contain NA and AICc/BIC comparisons are invalid. 
        Remove rows with NA from the input data or interpolate missing predictors.")}

  # make the interaction terms vector
  interaction.terms=NA

  # if there are factors
  if(length(na.omit(pred.vars.fact))>0){
  # if there are two or more factors
  if(length(na.omit(pred.vars.fact))>1){
    # make all the interactions between factors
    if(class(factor.factor.interactions)=="logical"){
     if(factor.factor.interactions==T){
        if(length(pred.vars.fact)<2){
            stop("You have less than 2 factors. Please reset 'factor.factor.interactions' to 'False'")}
      factor.correlations=check.correlations(use.dat[,pred.vars.fact])
      fact.combns=list()
      fact.cmbns.max.predictors=max.predictors
      if(max.predictors>length(pred.vars.fact)){fact.cmbns.max.predictors=length(pred.vars.fact)}
      for(i in 2:fact.cmbns.max.predictors){
        if(i<=length(pred.vars.fact)){
        fact.combns=c(fact.combns,
         combn(pred.vars.fact,i,simplify=F)) }}
        # check which were correlated
        fact.combns=lapply(fact.combns,FUN=function(x){
                row.index=which(match(rownames(factor.correlations),x)>0)
                col.index=which(match(colnames(factor.correlations),x)>0)
                cor.mat.m=factor.correlations[row.index,col.index]
                out=x
                if(max(abs(cor.mat.m[upper.tri(cor.mat.m)]))>cov.cutoff){out=NA}
                return(out)})
        fact.combns[which(is.na(fact.combns))]=NULL
        tt=data.frame(lapply(fact.combns,FUN=function(x){
                   do.call("paste",as.list(use.dat[,x]))}))
        factor.interaction.terms=unlist(lapply(fact.combns,FUN=paste,collapse=".I."))
        colnames(tt)=factor.interaction.terms

      use.dat=cbind(use.dat,tt)
      pred.vars.fact=c(pred.vars.fact,factor.interaction.terms)
      }
    }
    # make only specified interactions between factors
    if(class(factor.factor.interactions)=="character"){
        if(length(factor.factor.interactions)<2){
            stop("You specified less than 2 factors as factor.factor.interactions.")}
        if(max(is.na(match(factor.factor.interactions,colnames(use.dat))))==1){
            stop("Not all specified factor.factor.interactions are supplied in use.dat")}
      factor.correlations=check.correlations(use.dat[,factor.factor.interactions])
      #if(min(factor.correlations,na.rm=T)>cov.cutoff){
      #   stop("All factors have a correlation higher than your cutoff value")}
      if(length(which(factor.correlations<cov.cutoff))>1){
        fact.combns=list()
        fact.cmbns.max.predictors=max.predictors
        if(max.predictors>length(factor.factor.interactions)){fact.cmbns.max.predictors=length(factor.factor.interactions)}
        for(i in 2:fact.cmbns.max.predictors){
          if(i<=length(factor.factor.interactions)){
          fact.combns=c(fact.combns,
           combn(factor.factor.interactions,i,simplify=F)) }}
          # check which were correlated
          fact.combns=lapply(fact.combns,FUN=function(x){
                  row.index=which(match(rownames(factor.correlations),x)>0)
                  col.index=which(match(colnames(factor.correlations),x)>0)
                  cor.mat.m=factor.correlations[row.index,col.index]
                  out=x
                  if(max(abs(cor.mat.m[upper.tri(cor.mat.m)]))>cov.cutoff){out=NA}
                  return(out)})
          fact.combns[which(is.na(fact.combns))]=NULL
          tt=data.frame(lapply(fact.combns,FUN=function(x){
                     do.call("paste",as.list(use.dat[,x]))}))
          factor.interaction.terms=unlist(lapply(fact.combns,FUN=paste,collapse=".I."))
          colnames(tt)=factor.interaction.terms

        use.dat=cbind(use.dat,tt)
        pred.vars.fact=c(pred.vars.fact,factor.interaction.terms)
      }
    }
   }
   # make sure the factors are factors
   for(f in 1:length(pred.vars.fact)){
       use.dat[,pred.vars.fact[f]]=factor(use.dat[,pred.vars.fact[f]])}

   # check which ones should be included as interactions with the smoothers
   factor.smooth.interactions=pred.vars.fact[which(unlist(lapply(strsplit(pred.vars.fact,
      split=".I."),function(x){
      max(is.na(match(x,factor.smooth.interactions)))}))==0)]

   # make the interaction terms between the factors and continuous predictors
   if(length(na.omit(factor.smooth.interactions))>0){
    all.interactions=expand.grid(setdiff(pred.vars.cont,linear.vars),factor.smooth.interactions)
    interaction.terms=paste(all.interactions$Var1,all.interactions$Var2,sep=".by.")

    # now interactions between linear continous predictors and factors
    if(length(na.omit(linear.vars))>0){
     linear.interactions=expand.grid(linear.vars,factor.smooth.interactions)
     linear.interaction.terms=paste(linear.interactions$Var1,linear.interactions$Var2,
                                sep=".t.")}
    }
   }

   # if we want smooth.smooth interactions
   smooth.smooth.interaction.terms=NA
    # for interactions amonst all continuous predictors
    if(class(smooth.smooth.interactions)=="logical"){
      if(smooth.smooth.interactions==T){
        if(length(pred.vars.cont)<2){
            stop("You have less than 2 continuous predictors you wish interactions for.
            Please reset 'smooth.smooth.interactions' to 'False'")}
      continuous.correlations=check.correlations(use.dat[,pred.vars.cont])
      cont.combns=list()
      cont.cmbns.max.predictors=max.predictors
      if(max.predictors>length(pred.vars.cont)){cont.cmbns.max.predictors=length(pred.vars.cont)}
      for(i in 2:cont.cmbns.max.predictors){
        if(i<=length(pred.vars.cont)){
        cont.combns=c(cont.combns,
         combn(pred.vars.cont,i,simplify=F)) }}
        # check which were correlated
        cont.combns=lapply(cont.combns,FUN=function(x){
                row.index=which(match(rownames(continuous.correlations),x)>0)
                col.index=which(match(colnames(continuous.correlations),x)>0)
                cor.mat.m=continuous.correlations[row.index,col.index]
                out=x
                if(max(abs(cor.mat.m[upper.tri(cor.mat.m)]))>cov.cutoff){out=NA}
                return(out)})
        cont.combns[which(is.na(cont.combns))]=NULL
        tt=data.frame(lapply(cont.combns,FUN=function(x){
                   do.call("paste",as.list(use.dat[,x]))}))
        smooth.smooth.interaction.terms=unlist(lapply(cont.combns,FUN=paste,collapse=".te."))
        colnames(tt)=smooth.smooth.interaction.terms
     }
    }
    # for only specific interactions amonst continuous predictors
    if(class(smooth.smooth.interactions)=="character"){
        if(length(smooth.smooth.interactions)<2){
            stop("You specified less than 2 variables as smooth.smooth.interactions.")}
        if(max(is.na(match(smooth.smooth.interactions,colnames(use.dat))))==1){
            stop("Not all specified smooth.smooth.interactions are supplied in use.dat")}
      continuous.correlations=check.correlations(use.dat[,smooth.smooth.interactions])
      cont.combns=list()
      cont.cmbns.max.predictors=max.predictors
      if(max.predictors>length(smooth.smooth.interactions)){cont.cmbns.max.predictors=length(smooth.smooth.interactions)}
      for(i in 2:cont.cmbns.max.predictors){
        if(i<=length(smooth.smooth.interactions)){
        cont.combns=c(cont.combns,
         combn(smooth.smooth.interactions,i,simplify=F)) }}
        # check which were correlated
        cont.combns=lapply(cont.combns,FUN=function(x){
                row.index=which(match(rownames(continuous.correlations),x)>0)
                col.index=which(match(colnames(continuous.correlations),x)>0)
                cor.mat.m=continuous.correlations[row.index,col.index]
                out=x
                if(max(abs(cor.mat.m[upper.tri(cor.mat.m)]))>cov.cutoff){out=NA}
                return(out)})
        cont.combns[which(is.na(cont.combns))]=NULL
        tt=data.frame(lapply(cont.combns,FUN=function(x){
                   do.call("paste",as.list(use.dat[,x]))}))
        smooth.smooth.interaction.terms=unlist(lapply(cont.combns,FUN=paste,collapse=".te."))
        colnames(tt)=smooth.smooth.interaction.terms
    }

  all.predictors=na.omit(unique(c(all.predictors,pred.vars.fact)))
  # calculate a correlation matrix between all predictors
  cc=check.correlations(use.dat[,all.predictors],parallel=parallel,n.cores=n.cores)
  if(length(cor.matrix)==1){
   cor.matrix=cc
   # replace NA's with zero.
   cor.matrix[which(cor.matrix=="NaN")]=0
   cor.matrix[which(is.na(cor.matrix)==T)]=0}else{
      # check if the user defined matrix has the same rownames and colnames
      check.predictors=c(match(all.predictors,colnames(cor.matrix)),
                         match(all.predictors,rownames(cor.matrix)))
      missing.predictors=unique(rep(all.predictors,2)[which(is.na(check.predictors))])
      if(length(missing.predictors)>0){
            stop(paste("Supplied cor.matrix is missing required predictors: ",
            paste(missing.predictors,collapse=", "),".",sep=""))}
  }

  # make all possible combinations
  if(length(na.omit(c(pred.vars.cont,
                      pred.vars.fact)))<max.predictors){
        stop("Model max.predictors is greater than the number of predictors.")}
  all.mods=list()
  for(i in 1:max.predictors){
    all.mods=c(all.mods,
     combn(na.omit(c(pred.vars.cont,pred.vars.fact,
                     interaction.terms,
                     linear.interaction.terms,
                     smooth.smooth.interaction.terms)),
                     i,simplify=F))
  }

  # remove redundant models
  use.mods=all.mods
  for(m in 1:length(all.mods)){
    mod.m=all.mods[[m]]
    mod.terms=unlist(strsplit(unlist(strsplit(unlist(strsplit(mod.m,
                               split=".by.",fixed=T)),
                               split=".t.",fixed=T)),
                               split=".te.",fixed=T))
    n.vars.m=unique(unlist(strsplit(unlist(strsplit(unlist(strsplit(unlist(strsplit(mod.m,
                                  split=".by.",fixed=T)),
                                  split=".I.",fixed=T)),
                                  split=".t.",fixed=T)),
                                  split=".te.",fixed=T)))
    cont.vars=na.omit(na.omit(c(pred.vars.cont,linear.vars))[match(mod.terms,
                                           na.omit(c(pred.vars.cont,linear.vars)))])
    fact.vars=unique(na.omit(pred.vars.fact[match(mod.terms,pred.vars.fact)]))

    # if there are factor vars
    if(length(fact.vars)>0){
      # check that any "by" factor vars are accompanied by a + term in its owns right
      if(max(is.na(match(fact.vars,mod.m)))==1){use.mods[[m]]=NA}}

    # remove the model if the predictors are correlated
    if(length(mod.terms)>1){
     row.index=which(match(rownames(cor.matrix),unique(mod.terms))>0)
     col.index=which(match(colnames(cor.matrix),unique(mod.terms))>0)
     cor.mat.m=cor.matrix[row.index,col.index]
     if(max(abs(cor.mat.m[upper.tri(cor.mat.m)]))>cov.cutoff){use.mods[[m]]=NA}
    }

    # remove the model if there are more than the number of terms specified in "max.predictors"
    if(length(n.vars.m)>max.predictors){use.mods[[m]]=NA}

    # remove the models if a continuous predictor occurs as a by, or a te, and as a single term
    if(length(cont.vars)>length(unique(cont.vars))){use.mods[[m]]=NA}

  }

  use.mods[which(is.na(use.mods))]=NULL

  # now make the models into gamm formula
  if(nchar(null.terms)==0){# if there is no bs='re' random effect random effect
                             # or other null term in the null model
    mod.formula=list(as.formula("~ intercept-1"))}
  if(nchar(null.terms)>0){# to add a bs='re' random effect
    mod.formula=list(null.formula)}

  for(m in 1:length(use.mods)){
     mod.m=use.mods[[m]]
     cont.smooths=mod.m[which(match(mod.m,setdiff(pred.vars.cont,linear.vars))>0)]
     by.smooths=mod.m[grep(".by.",mod.m,fixed=T)]
     te.smooths=mod.m[grep(".te.",mod.m,fixed=T)]
     factor.terms=mod.m[which(match(mod.m,pred.vars.fact)>0)]
     linear.terms=mod.m[which(match(mod.m,linear.vars)>0)]
     linear.interaction.terms=mod.m[grep(paste(linear.vars,".t.",sep=""),mod.m,fixed=T)]
     all.terms.vec=character()

     if(length(cont.smooths>0)){all.terms.vec=c(all.terms.vec,
                  paste("s(",cont.smooths,",k=",k,",bs=",bs.arg,")",sep=""))}
     if(length(by.smooths>0)){all.terms.vec=c(all.terms.vec,
         paste("s(",gsub(".by.",",by=",by.smooths,fixed=T),",k=",k,",bs=",bs.arg,")",sep=""))}
     if(length(te.smooths>0) & class(test.fit)[[1]]!="gamm4"){all.terms.vec=c(all.terms.vec,
         paste("te(",gsub(".te.",",",te.smooths,fixed=T),",k=",k,",bs=",bs.arg,")",sep=""))}
     if(length(te.smooths>0) & class(test.fit)[[1]]=="gamm4"){all.terms.vec=c(all.terms.vec,
         paste("t2(",gsub(".te.",",",te.smooths,fixed=T),",k=",k,",bs=",bs.arg,")",sep=""))}
      if(length(linear.interaction.terms>0)){all.terms.vec=c(all.terms.vec,
               gsub(".t.","*",linear.interaction.terms,fixed=T))}
     if(length(factor.terms>0)){all.terms.vec=c(all.terms.vec,factor.terms)}
     if(length(linear.terms>0)){all.terms.vec=c(all.terms.vec,linear.terms)}
     if(max(is.na(cyclic.vars))!=1){
       for(r in 1:length(cyclic.vars)){
           for(v in 1:length(all.terms.vec)){
             if(length(grep(cyclic.vars[r],all.terms.vec[v]))>0){
                  all.terms.vec[v]=gsub(paste("bs=",bs.arg,sep=""),"bs='cc'",all.terms.vec[v])
                  }}}}
     for(v in 1:length(all.terms.vec)){
         if(length(grep("te(",all.terms.vec[v],fixed=T))>0){
            bs.arg.v=c("","")
            smooth.vars.v=unlist(strsplit(gsub("te(","",all.terms.vec[v],fixed=T),split=","))[1:2]
            var.type.vec=unlist(lapply(smooth.vars.v,FUN=function(x){match(x,cyclic.vars)}))
            bs.arg.v[which(is.na(var.type.vec))]=bs.arg
            bs.arg.v[which(var.type.vec>0)]="'cc'"
            bs.arg.v=paste("bs=c(",paste0(bs.arg.v,collapse=","),")",sep="")
            all.terms.vec[v]=gsub("bs='cc'",bs.arg.v,all.terms.vec[v])
            all.terms.vec[v]=gsub(paste("bs=",bs.arg,sep=""),bs.arg.v,all.terms.vec[v])
         }
     }

     if(nchar(null.terms)==0){# if there is no bs='re' random effect
                                # or other null term in the null model
       formula.m=as.formula(paste("~",
               paste(all.terms.vec,collapse="+")))}
     if(nchar(null.terms)>0){#
       formula.m=as.formula(paste("~",
               paste(c(all.terms.vec,null.terms),collapse="+")))}
     mod.formula=c(mod.formula,list(formula.m))
  }

  names(mod.formula)=c("null",lapply(use.mods,FUN=paste,collapse="+"))

  # Is this too many models?
  n.mods=length(mod.formula)
  time.to.run=round(system.time(try(update(test.fit,formula=mod.formula[[n.mods]],data=use.dat),silent=T))[3]*n.mods/60)
  test.mod=try(update(test.fit,formula=mod.formula[[n.mods]],data=use.dat),silent=T)
  mod.gbs=round(object.size(test.mod)/1073741824*n.mods,1)
  if(n.mods>max.models){
        stop(paste("You have ",n.mods," models. If you want to fit all of these you need to
        increase 'max.models' from ",max.models,". Otherwise, if the model set
        is larger than you can realistically fit, try reducing the number of predictors,
        setting the covariance 'cov.cutoff' argument to less than ", cov.cutoff,
        "
        or setting 'factor.factor.interactions' to FALSE (if you have factors).",sep=""))
       }

  # now fit the models by updating the test fit (with or without parallel)
  if(parallel==T){
   require(doParallel)
   cl=makePSOCKcluster(n.cores)
   registerDoParallel(cl)
   out.dat<-foreach(l = 1:length(mod.formula),
                   .packages=c('mgcv','gamm4','MuMIn'),
                   .errorhandling='pass')%dopar%{
         if(length(grep("dsm",class(test.fit)))>0){
           out=update(test.fit,formula=mod.formula[[l]])}
         if(length(grep("dsm",class(test.fit)))==0){
        out=update(test.fit,formula=mod.formula[[l]],data=use.dat)}
   }
   stopCluster(cl)
   registerDoSEQ()
           }else{
      out.dat=list()
      for(l in 1:length(mod.formula)){
         if(length(grep("dsm",class(test.fit)))>0){
           out=try(update(test.fit,formula=mod.formula[[l]]),silent=T)}
         if(length(grep("dsm",class(test.fit)))==0){
        out=try(update(test.fit,formula=mod.formula[[l]],data=use.dat),silent=T)}
        out.dat=c(out.dat,list(out))}
  }
  names(out.dat)=names(mod.formula[1:n.mods])

  # find all the models that didn't fit and extract the error messages
#  model.success=lapply(lapply(out.dat,FUN=class),FUN=function(x){
#     x[1]!="try-error"})
  model.success=lapply(lapply(out.dat,FUN=class),FUN=function(x){
     length(grep("gam",x))>0})

  failed.models=mod.formula[which(model.success==F)]
  success.models=out.dat[which(model.success==T)]
  if(length(success.models)==0){
        stop("None of your models fitted successfully. Please check your input objects.")}

  # some functions for extracting model information
  require(MuMIn)
  wi<<-function(AIC.vals){# This function calculate the Aikaike weights:
   # wi=(exp(-1/2*AICc.vals.adj))/Sum.wi=1 to r (exp(-1/2*AICc.vals.adj))
   AICc.vals.adj=AIC.vals-min(na.omit(AIC.vals))
   wi.den=rep(NA,length(AICc.vals.adj))
   for(i in 1:length(AICc.vals.adj)){
    wi.den[i]=exp(-1/2*AICc.vals.adj[i])}
   wi.den.sum=sum(na.omit(wi.den))
   wi=wi.den/wi.den.sum
   return(wi)}

  # of the successful models, make a table indicating which variables are included
  var.inclusions=matrix(0,ncol=length(included.vars),length(success.models))
  colnames(var.inclusions)=c(included.vars)

  for(m in 1:length(success.models)){
        pred.vars.m=unique(
          unlist(strsplit(unlist(strsplit(unlist(strsplit(unlist(strsplit(unlist(strsplit(unlist(strsplit(names(success.models)[m],
          split="+",fixed=T)),
          split=".by.",fixed=T)),
          split=".I.",fixed=T)),
          split="*",fixed=T)),
          split=".t.",fixed=T)),
          split=".te.",fixed=T)))
        if(pred.vars.m[1]!="null"){var.inclusions[m,pred.vars.m]=1}}

  # now make a table of all the model summary data
  mod.data.out=data.frame("modname"=names(success.models))
  mod.data.out$formula=lapply(success.models,FUN=function(x){as.character(formula(x)[3])})
  mod.data.out$AICc=unlist(lapply(success.models,FUN=AICc))
  mod.data.out$BIC=unlist(lapply(success.models,FUN=BIC))
  mod.data.out$delta.AICc=round(mod.data.out$AICc-min(mod.data.out$AICc),3)
  mod.data.out$delta.BIC=round(mod.data.out$BIC-min(mod.data.out$BIC),3)
  mod.data.out$wi.AICc=round(wi(mod.data.out$AICc),3)
  mod.data.out$wi.BIC=round(wi(mod.data.out$BIC),3)
  mod.data.out$r2.vals=round(unlist(lapply(success.models,FUN=function(x){
        out=NA
        if(class(x)[1]=="gam" & r2.type=="dev"){out=summary(x)$dev.expl}
        if(class(x)[1]=="gam" & r2.type=="r2"){out=summary(x)$r.sq}
        if(class(x)[1]=="gam" & r2.type=="r2.lm.est"){
           out=summary(lm(x$y~predict(x)))$r.sq}
        if(class(x)[[1]]=="gamm4" & r2.type=="dev"){
           out=summary(x$gam)$dev.expl
           if(length(out)==0){out=NA}}
        if(class(x)[[1]]=="gamm4" & r2.type=="r2"){out=summary(x$gam)$r.sq}
        if(class(x)[[1]]=="gamm4" & r2.type=="r2.lm.est"){
           out=summary(lm(attributes(x$mer)$frame$y~
                        predict(x[[1]],re.form=NA,type="response")))$r.sq}
           if(is.null(out)){out=NA}
        return(out)})),3)

  # substract the null model r2 value from each model r2 value
  if(report.unique.r2==T){
  null.r2=mod.data.out$r2.vals[which(mod.data.out$modname=="null")]
  mod.data.out$r2.vals.unique=mod.data.out$r2.vals-null.r2}

  # now calculate the summed edf
  mod.data.out$edf=round(unlist(lapply(success.models,FUN=function(x){
        if(class(x)[1]=="gam"){
          edf.m=summary(x)$edf
          p.coeff.m=summary(x)$p.coeff}else{
           #edf.m=summary(x$gam)$edf
           #p.coeff.m=summary(x$gam)$p.coeff
           edf.m=x$gam$edf
           p.coeff.m=x$gam$p.coeff
           }
        edf.m[which(edf.m<1)]=1 # any edf<0 are reset to 1 to ensure proper
                                # parameter count when there is shrinkage (bs='cc')
        return(sum(c(edf.m,length(p.coeff.m))))})),2)
  # count the edf values less than 0.25 to check for serious shrinkage
  mod.data.out$edf.less.1=unlist(lapply(success.models,FUN=function(x){
        #if(class(x)[1]=="gam"){edf.m=summary(x)$edf}else{edf.m=summary(x$gam)$edf}
        if(class(x)[1]=="gam"){edf.m=summary(x)$edf}else{edf.m=x$gam$edf}
        return(length(which(edf.m<0.25)))}))
  # now add columns for the included predictors to the dataframe
  mod.data.out=cbind(mod.data.out,var.inclusions)

  # now calculate the variable importance
   # find the min number of models for each variable
  min.mods=min(colSums(mod.data.out[,included.vars]))
  # first for AICc
  var.weights=unlist(lapply(included.vars,FUN=function(x){
           sum(sort(mod.data.out$wi.AICc[which(mod.data.out[,x]==1)],decreasing=T)[1:min.mods])}))
  names(var.weights)=included.vars
  variable.weights.raw=var.weights
  #variable.weights.raw=colSums(mod.data.out[,included.vars]*mod.data.out$wi.AICc)
  aic.var.weights=list(variable.weights.raw=variable.weights.raw)
  # next for BIC
  var.weights=unlist(lapply(included.vars,FUN=function(x){
           sum(sort(mod.data.out$wi.BIC[which(mod.data.out[,x]==1)],decreasing=T)[1:min.mods])}))
  names(var.weights)=included.vars
  variable.weights.raw=var.weights
  #variable.weights.raw=colSums(mod.data.out[,included.vars]*mod.data.out$wi.BIC)
  bic.var.weights=list(variable.weights.raw=variable.weights.raw)
  # now return the list of outputs
  return(list(mod.data.out=mod.data.out,
              used.data=use.dat,
              predictor.correlations=cor.matrix,
              #mod.formula=mod.formula,
              failed.models=failed.models,
              success.models=success.models,
              variable.importance=
                 list(aic=aic.var.weights,bic=bic.var.weights)))
} #------------------ end function --------------------------------------------#