# Function Outputs

## Overview

This vignette is based on Appendix S2 of Fisher et al. (2018), which
documented the outputs returned by the original combined
[`full.subsets.gam()`](https://beckyfisher.github.io/FSSgam_package/reference/full.subsets.gam.html)
function. That function has since been split into
[`generate_model_set()`](https://beckyfisher.github.io/FSSgam_package/reference/generate_model_set.html),
which builds and returns the candidate model set, and
[`fit_model_set()`](https://beckyfisher.github.io/FSSgam_package/reference/fit_model_set.html),
which fits every model in that set and returns the fitted results.
[`full.subsets.gam()`](https://beckyfisher.github.io/FSSgam_package/reference/full.subsets.gam.html)
itself still exists and returns the same combined set of outputs
described here.

The tables below reproduce the original output descriptions, split
according to which of the two functions now returns each one — matching
how they are accessed in the case study vignettes
(e.g. `model.set$predictor.correlations` from
[`generate_model_set()`](https://beckyfisher.github.io/FSSgam_package/reference/generate_model_set.html)’s
return value, and `out.list$mod.data.out` from
[`fit_model_set()`](https://beckyfisher.github.io/FSSgam_package/reference/fit_model_set.html)’s).

## Outputs from `generate_model_set()`

| Output | Details |
|----|----|
| \$used.data | A data.frame which is identical to the data.frame initially supplied by the user, but with any hard coded interaction terms appended via cbind. |
| \$predictor.correlations | The matrix of estimated predictor correlations returned by the function check.correlations and used for model exclusion based on cov.cutoff (see the function-arguments vignette). |

## Outputs from `fit_model_set()`

| Output | Details |
|----|----|
| \$mod.data.out | A data.frame that contains the statistics associated with each model fit. This includes AICc and BIC, delta values (e.g. AICc-(min(AICc))), corresponding weight (omega_i) values (Burnham and Anderson 2002), an estimate of the model R2, and a column for each of the included predictor variables containing either 0 (variable not included in the model) or 1 (variable is present in the model). Use of BIC in information theoretic approaches has been heavily criticised because of the inherent assumption of BIC that there is a “true” model that is represented in the candidate set (Anderson and Burnham 2002). Rather than decide a-priori which model selection tool users should adopt, both AICc and BIC are supplied as part of the function outputs. To simplify output, only AICc and AICc based model weights, rather than AIC, are included as these are asymptotically equivalent at large sample sizes, and for small sample sizes AICc should be used in any case. Calculating R2 values is non-trivial for mixed models, especially non-gaussian cases. A range of methods for estimating R2 are supplied (see `r2.type` in the function-arguments vignette), as in our experience a single method rarely performs adequately across all scenarios. |
| \$failed.models | A list containing the try-error catch associated with models that failed to fit. Ideally the list of failed models should be empty, but when this is not the case interrogating failed.models provides a useful means of troubleshooting. Users can examine which models are not fitting and explore the reasons for this by fitting the failed models outside the call based on the listed formula. When a large number of models fail to fit properly it usually indicates poor specification of the initial test.fit or other arguments in the call to generate_model_set (such as the inclusion of factor interactions when there are few data within each level of the factor), or that inappropriate variables are being included in the model set. |
| \$success.models | A complete list of all successfully fitted models. This can be used for multimodel inference and creating model averaged predictions. |
| \$variable.importance | A list containing importance scores for each included predictor. To determine the relative importance of each predictor across the whole model set, the omega_i values are summed for all models containing each variable. The higher the combined weights for an explanatory parameter, the more important it is in the analysis (Burnham and Anderson 2002). An assumption of the use of summed model weights to infer variable importance is that the number of models in which the different predictors are present is uniform. As correlated predictors are removed from the model set, this is not always the case. To overcome this issue, the summed variable.importance scores are the summed weights for the best n models, where n is equal to the minimum number of models any one predictor is present in. |
| $`variable.importance`$aic | **\$aic** — Variable importance values based on summed AICc weights. |
| $`variable.importance`$bic | **\$bic** — Variable importance values based on summed BIC weights. |

## References

Anderson, D. R., and K. P. Burnham. 2002. “Avoiding Pitfalls When Using
Information-Theoretic Methods.” *Journal of Wildlife Management* 66:
912–18.

Burnham, K. P., and D. R. Anderson. 2002. *Model Selection and
Multimodel Inference: A Practical Information-Theoretic Approach*. 2nd
ed. Springer-Verlag.

Fisher, R., S. K. Wilson, T. M. Sin, A. C. Lee, and T. J. Langlois.
2018. “A Simple Function for Full-Subsets Multiple Regression in Ecology
with r.” *Ecology and Evolution* 8 (12): 6104–13.
<https://doi.org/10.1002/ece3.4134>.
