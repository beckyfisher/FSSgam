

cat.preds=c("survey","bleach_pres","dredge_pres","dhw_fact" )
null.vars=c("site")
cont.preds=c("Av_wave_height","Depth")

#model
use.dat=na.omit(dat[,c(null.vars,cat.preds,cont.preds,"coral")])
  Model1=uGamm(coral~s(Latitude,k=4,bs='cr'),
              family=gaussian(), random=~(1|site),
             data=use.dat,
             lme4=TRUE)

  out.list=full.subsets.gam(use.dat=use.dat,max.predictors=3,
                            test.fit=Model1,k=5,
                            pred.vars.cont=cont.preds,
                            pred.vars.fact=cat.preds,
                            factor.factor.interactions=TRUE)












