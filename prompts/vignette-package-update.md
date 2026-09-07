## Session: update vignettes for FSSgam 1.1.0.9000 (issues #19, #20)
Date: 2026-09-07
Model: Claude Opus 5 (1M context)

### Prompts and Responses

**User:** Can we tackle https://github.com/beckyfisher/FSSgam/issues/19 and
https://github.com/beckyfisher/FSSgam/issues/20. Ideally we would make as few changes as
possible on the outward website material, and update largely the vignettes.

**Claude:** Measured what the package update actually changes before editing anything.

#### Method

FSSgam 1.1.0.9000 (dev, sha b0e3b90) was installed into an isolated library so the existing
1.0.0 install -- the version the live site was built against -- stayed available. Each
vignette's code was purled and run twice, once against each version, on the same R 4.5.1 and
the same mgcv. Model tables were compared on term-set-normalised model names, because the
package now writes the terms within a name in byte order and a naive name match reports every
reordered model as a difference. Published values were extracted from the Appendix S3/S4/S5
.docx files in publication/ for comparison.

#### What the update changes

Three of the four analysis vignettes are unaffected. Over all candidates:

| vignette | candidates | max abs AICc difference | max abs wi.AICc difference |
| --- | --- | --- | --- |
| case-study-1 | 224 | 0.52 (0.02 on the published rows) | 0.003 |
| case-study-3 | 52 | 0 | 0 |
| extra-examples | 26 | 0 | 0 |
| case-study-2 | 139 per taxon | 126.7 | 0.987 |

Case study 2 changed substantially, and the change is a correction. Under 1.0.0 each candidate
was fitted with a Tweedie power parameter of about 1.01 rather than the value estimated for
that model; the test.fit itself had p = 1.26. A Tweedie with p near 1 approaches a Poisson, so
the assumed variance function was wrong, inflating both edf and the information criteria.

Isolated to a single model to confirm the cause, for Pagurus novizelandiae:

```
model  Distance+lobster+sqrt.X2mm, 95 rows, identical formula and data under both versions

direct mgcv::gam()   AICc 478.76  edf 19.53  Tweedie(p=1.242)
update(test.fit)     AICc 478.76  edf 19.53  Tweedie(p=1.242)
fit_model_set 1.0.0  AICc 605.50  edf 24.87  Tweedie(p=1.01)
fit_model_set new    AICc 478.76  edf 19.53  Tweedie(p=1.242)
```

The candidate set is identical either way (139 models, same formulae, same 95 rows, same
factor levels, same test.fit at edf 17.86). Only the fitted variance function differed, and
the corrected version reproduces a direct mgcv fit exactly.

The practical effect is on how many models sit close to the best one. For Pagurus
novizelandiae, 1.0.0 placed one model within 2 AICc units of the best; the corrected fit
places twelve, against eighteen in published Table A4.2, with edf of 17.7-21.1 against the
published 17.7-21.5. The corrected results are closer to the published analysis than what the
site previously showed. Dosinia subrosea and Myadora striata change less: two models within 2
AICc units becoming two and three.

#### Departures from the published tables that are not attributable to FSSgam

Case study 1 differs from Table A3.1 by a median of 0.66 and a maximum of 2.47 AICc across the
seventeen published rows, and 1.0.0 differs from 1.1.0.9000 by a median of 0.00 on those same
rows. The drift is therefore in mgcv and R since 2018, not in the full subsets procedure. The
top-ranked model is unchanged for all eight responses. Case study 3 reproduces Table A5.1's
best model exactly at published precision, with R2 and total edf reproducing for all six rows
and only the delta AICc values of rows 2-6 differing, enough to exchange rows 5 and 6.

Decision (RF): recompute live and record the gap in each case study, rather than freezing the
published numbers or printing both tables. Recorded in a 'Relationship to the published table'
section in each of the three case studies.

#### Changes made

- `function-arguments.Rmd`: documented `non.linear.correlations`, `null.cov.cutoff`,
  `save.model.fits`, `VI.mods`, `progress` and `logLik.fn`, which were missing; corrected the
  `max.models` default from 500 to 200; `check.correlations` to `check_correlations`.
- `function-outputs.Rmd`: documented `$n.mods`, `$mod.formula`, `$null.term.correlations`,
  `$included.vars` and `$test.fit`; named the `r2.vals.unique` and `edf.less.1` columns;
  recorded that `full_subsets_gam()` now exists alongside `full.subsets.gam()`.
- `case-study-1.Rmd`, `case-study-2.Rmd`, `case-study-3.Rmd`: model tables now carry the
  columns and labels of published Tables A3.1, A4.2 and A5.1 (BIC, delta BIC and BIC weights
  were previously computed but not shown); a 'Relationship to the published table' section in
  each; the byte-order model naming change recorded.
- `case-study-1.Rmd`: 'FFSgam' corrected to 'FSSgam' in three places; the garbled sentence on
  predictor transformation rewritten; a note added on why na.omit is applied per response
  rather than once.
- `case-study-2.Rmd`: the Tweedie correction appended to the existing dated note block.
- `faq.Rmd`: two dead links to case_study1_reef_fish.R and extra_examples.R repointed at the
  vignettes that replaced them; the rename entry extended to cover check_correlations,
  check_non_linear_correlations, build_inclusion_mat and extract_mod_dat, which were renamed
  without a deprecated alias and so fail with 'could not find function', and fit_mod_l, which
  is no longer exported.

No change was made to the analysis code, the model specifications, the data, or any figure.

#### Verification

All seven vignettes were rendered against FSSgam 1.1.0.9000, each in its own Rscript process,
with no errors, no deprecation warnings and no unresolved citations.

---

