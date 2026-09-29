## Session: cran-acceptance-install-note
Date: 2026-09-29
Model: Claude Opus 5.5 (claude-opus-5-5[1m])

### Prompts and Responses

**User:** fssgam has been accepted on cran. Can you update anything we need here, or on the docummentation package to accomodate that?

**Claude:** FSSgam 1.2.0 was published on CRAN on 2026-09-28. The README installation
section was rewritten to give `install.packages("FSSgam")` first, the GitHub development
version second, and to state that the site's vignettes are built against the package's
`dev` branch, which can be ahead of the CRAN release, with `master` holding the released
version. No vignette changed. `.github/workflows/pkgdown.yaml` still installs
`FSSgam_package`'s `dev` ref on every branch; whether the root site should instead build
against the CRAN release was left to the user, per Section 4 of this repository's
`CLAUDE.md`.

---

**User:** Help me decide on [whether the root site should build against the CRAN release]. / yes, go ahead I agree

**Claude:** Recommended, and on agreement implemented, building the root site against the
CRAN release and the `/dev/` site against `FSSgam_package`'s `dev` branch. Reason: a reader
who installs from CRAN should reproduce the numbers on the root site, which is not
guaranteed while it renders against unreleased code. The branch-conditional form removed on
2026-09-07 fell back to the package repository's default branch; the CRAN release is a fixed
target and does not have that defect. Cost: a vignette needing an unreleased feature stays
on `dev` until that feature is on CRAN.

`.github/workflows/pkgdown.yaml`:

``` r
if (Sys.getenv("DOCS_BRANCH") == "dev") {
  remotes::install_github("beckyfisher/FSSgam_package", ref = "dev")
} else {
  install.packages("FSSgam")
}
message("Building against FSSgam ", packageVersion("FSSgam"))
```

with `DOCS_BRANCH: ${{ github.ref_name }}`. The README installation text and `CLAUDE.md`
Sections 3 and 4 were revised to match. The workflow was parsed with `yaml::read_yaml()` and
the step's R code with `parse()`; the build itself is verified only by the CI runs on push.

---
