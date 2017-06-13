# FSSgam
# Consists of two functions: full.subsets.gam and check.correlations
# requires: doParallel, MuMIn, gamm4, mgcv, nnet

### arguments ################################################################
# use.dat - A dataframe, with columns matching those included in pred.vars.cont and 
#          pred.vars.fact, the response variable to be analysed and any other fields 
#          required for the analysis (such as random effects, see test.fit).
# test.fit - A gam model fitted via a call to gam (mgcv) or uGamm (MuMIn).
#          This can use any of the (preferably continuous) predictors in the call
#          and will be used as a model to update in the fitting of the model set.
#          This must contain the appropriate random effects and call to family
#          (if not gaussian) and (in the case of a uGamm cal, if gamm4 should be used, 
#          or gamm (see ? uGamm)
# pred.vars.cont - A character vector indicating the continuous predictors to use. 
#          By default all continuous predictors will be fitted using a smoother 
#          (but see argument linear.vars). These must match column names in use.dat exactly. 
#          More than one continuous smooth predictor must be supplied to run the 
#          full.subsets.gam function.
# pred.vars.fact - A character vector indicating the factor predictors to use.
#          These must match column names in use.dat exactly.
#          All pred.vars.fact variables will be coerced to a factor, with
#          unused levels simultaneously dropped via a call to factor. All
#          factor variables will be included as interactions with the continuous
#          smoothers via their specification as a by argument in the call to gam.
#          Factor interactions are not included. This may be possible by manually
#          creating interacting factors from two or more factor variables using paste.
#          This may cause issues with over specified models, but such models are not
#          expected to be successfully fit and will be removed as failed models.
#          inspection of the failed models list.
# cyclic.vars - NA if there are no cyclic predictors, or if there are cyclic predictors,
#          a character vector containing the names of any cyclic variables.
#          note that these must also be contained in the pred.vars.cont character
#          vector as a continuous predictor. Please also note there are
#          issues with bs='cc' and model selection as this uses by default
#          shrinkagee. With shrinkage, variables
#          are retained in models but with zero edf, which makes interpretation of
#          AICc and BIC confusing. To account for this always select only the
#          most parsinomious model (that with the fewest parameters), not
#          just that with the lowest AICc. Reported estimated degrees of freedom
#          in the model output table represent the sum of the edf of the smooth
#          terms plus the number of paramteric coefficients. When cyclic variables
#          are included and shrinkage is used, any estimated edf of the smooth terms
#          that are less than 1 are reset to 1 before summing to ensure the the
#          total number of predictors in the model is captured properly.
# linear.vars - NA if there are no continuous predictors you which to treat as
#          linear predictors (no smooth)
#          only use this where variables are clearly continuous in nature, but you are
#          confident a linear relationship is adequate. It may also be useful
#          for continuous predictors are not well distributed along the x-axis
#          (ie, sampling was conducted in clumped distances from a feature of
#          interest). Where this is necessary, transformations should be considered
#          where they can be used to theoretical "linearlize" response relationships
# size - An integer indicating the maximum model size to fit (ie the maximum
#          number of predictors to include in any one model. Default is 3.
# cov.cutoff - A numeric value between 0 and 1 indicating the covariance cutoff
#          value to use for excluding colinear models.
#          Defaults is 0.28 (see Graham MH (2003) CONFRONTING MULTICOLLINEARITY
#          IN ECOLOGICAL MULTIPLE REGRESSION. Ecology 84:2809-2815. It is highly
#          recommended to keep this value low, as correlation among predictors
#          can yield spurious results. Note that predictors with a correlation
#          greater than the specified value will still appear in the model set
#          but will never appear in the same model. Including highly correlated
#          predictors can make interpreting variable importance values difficult.
# k - An integer indicationg the dimension of the basis used to represent the
#          smooth term(see ?s).
#          The default value is 5. Higher values are not recommended unless
#          a complex trend between the response variable and the
#          continuous predictor variables is expected, and the data are sufficient
#          to support this. k can be reduced to as low
#          as 3 where there is trouble obtaining convergence, or sample size is low.
#          Note that this must be set to overide the default value, regardless
#          of what k is used in the test.fit
# parallel - A logical value indicating if parallel processing should be used.
#          The default is FALSE. Before changing to TRUE,
#          please make sure parallel processing works for you by runing the code:
#          cl=makePSOCKcluster(2)
#          registerDoParallel(cl)
# n.cores - The number of cores to use in the parallel processing. The default value
#          is 4.
# weight.vec - the weights vector to use for passing n trials in a binomial call
# null.terms -  A character vector indicating the form of any re smooths to be
#          included in gam [e.g. "s(site,bs='re')"] or any other fixed terms or
#          smooths that the user wants to include in the null model. Use of
#          bs="re" is an alternative way of fitting simple random structures
#          that avoids use of PQL and allows a the greater range of families
#          available in gam.mgcv to be used. see ?s and links therin.
#          Note: make sure you use "gam" instead of uGamm to make sure PQL is
#          not used. Other use of the null.terms argument are to pass
#          smooths for factors that the user want to ensure are contained in
#          all models, including the "null". For example, there may be known
#          effects of length on swimming speed, or depth on species abundances
#          that do not form part of the ecological question being addressed
#          but must be accounted for in model fits.
# bs.arg - Default is 'cr'. Can be used to choose a different smoother, see ?s for more
#          information on smoother provided in gam (mgcv)
# smooth.interactions - A character vector specify if/which factor variables
#          should be included as "by" argument interaction terms with the
#          continuous predictors. Default is the list of pred.vars.fact, meaning
#          that all will be included by default, otherwise only those included
#          in the character vector will be added. Note that specified factors must
#          also be included in pred.vars.fact. If specified as NA o factor
#          interactions will be include.
#          Note interactions are only included  where the cov.cutoff value is not
#          exceeded based on the calculated all predictor correlation matrix.
# factor.interactions - A logical value indicating if interactions between
#          factors should be included.
#          Default is 'F'. Note that this can massively increase the number of
#          models in the canidate set. Not recomemnded when
#          there are factors with many levels.
# max.models - The total number of models allowed to be fit. If the candidate set
#          is bigger than this value, an error message will be returned, asking
#          the user to reset this to a larger value if they are sure they want to
#          fit that many models. Default is 500,
# r2.type - The value to extract from the gam model fit to use as the R squared value.
#         Defaults to "r2.lm.est" which returns and estimated R squared
#         value based on a linear regression between the observed and
#         predicted values.
#         "r2" will return the adjusted R.sq as reported by gam, gamm or gamm4.
#         "dev" which is the deviance explained as reported by gam or
#         gamm. Note gamm4 does not currently return a deviance.
# report.unique.r2 - The estimate null model R2 is substracted from each
#         model R2 to give an idea of the unique variance explained. This can
#         be useful where null terms are included in the model set.

### known issues -------------------------------------------------------------
# This function assumes you know what you are doing. Non-gaussian mixed model
# gamm resorts to PQL meaning that AICc calls will not return the AIC of the
# actual model. Please thoroughly read the help files contained within the gamm4
# and mgcv packages, including information under details.
### -------------------------------------------------------------------------- 
