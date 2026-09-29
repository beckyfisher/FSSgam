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
