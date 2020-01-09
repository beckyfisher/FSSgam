### FSSgam
Consists of two functions: full.subsets.gam and check.correlations

Requires: doSNOW, MuMIn, gamm4, mgcv, nnet
Last tested using R version 3.4.3

### Summary ################################################################
 
Full subsets information theoretic approaches are becoming an increasingly popular tool for exploring predictive power and variable importance where a wide range of candidate predictors are being considered.

This repository contains a simple function in the statistical programming language R that can be used to construct, fit and compare a complete model set of possible ecological or environmental predictors, given a response variable of interest. The function is based on Generalized Additive Models (GAM) and builds on the MuMIn package.

Advantages include the capacity to fit more predictors than there are replicates, automatic removal of models with correlated predictors, and model sets that include interactions between factors and smooth predictors, as all as smooth interactions with other smooths (via te). 

The function takes a range of arguments that allow control over the model set being constructed, including specifying cyclic and linear continuous predictors, specification of the smoothing algorithm used and the maximum complexity allowed for smooth terms. 

The use of the function is demonstrated via case studies that highlight how appropriate model sets can be easily constructed, and the broader utility of the approach for exploratory ecology.

### Argument and outputs
See Table S1 and Table S2 for argument details and an explanation of outputs returned by the function.

### known issues
This function assumes you know what you are doing. Non-gaussian mixed model
gamm resorts to PQL meaning that AICc calls will not return the AIC of the
actual model. Please thoroughly read the help files contained within the gamm4
and mgcv packages, including information under details.


### License

The code is released under the Apache License 2.0

    Copyright 2020 Australian Institute of Marine Science

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at 

       http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.