### FSSgam
Consists of two functions: full.subsets.gam and check.correlations
requires: doParallel, MuMIn, gamm4, mgcv, nnet

### Summary ################################################################
 
Full subsets information theoretic approaches are becoming an increasingly popular tool for exploring predictive power and variable importance where a wide range of candidate predictors are being considered.

This repository contains a simple function in the statistical programming language R that can be used to construct, fit and compare a complete model set of possible ecological or environmental predictors, given a response variable of interest. The function is based on Generalized Additive Models (GAM) and builds on the MuMIn package.

Advantages include the capacity to fit more predictors than there are replicates, automatic removal of collinear models, and model sets that include interactions between factors and smooth predictors. 

The function takes a range of arguments that allow control over the model set being constructed, including specifying cyclic and linear continuous predictors, specification of the smoothing algorithm used and the maximum complexity allowed for smooth terms. 

The use of the function is demonstrated via case studies that highlight how appropriate model sets can be easily constructed, and the broader utility of the approach for exploratory ecology.


### known issues
This function assumes you know what you are doing. Non-gaussian mixed model
gamm resorts to PQL meaning that AICc calls will not return the AIC of the
actual model. Please thoroughly read the help files contained within the gamm4
and mgcv packages, including information under details.



![alt text](https://user-images.githubusercontent.com/14978794/27793878-571b6bba-6032-11e7-9c3b-8e651238b116.png "Importance scores")
<small>Example of importance scores </small>
