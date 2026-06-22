# Project: FSSgam (documentation companion repo)

This file provides **project-specific context only** — repo structure,
dependencies, the current task, and prompt logging. It does not impose
personal style preferences on collaborators.

**For Claude:** If a `CLAUDE.md` exists in the parent directory of this repo,
read it before proceeding — it may contain user-specific environment and
style preferences. This repo-level file takes precedence over any
parent-level file where they conflict, **except** where a parent-level
default is explicitly overridden below (Section 4 notes one such override).

---

## 1. Repo Type

**Documentation/pkgdown companion repo**, not a deployable R package in the
normal sense. `DESCRIPTION` declares `Package: FSSgam.docs` — that name is
deliberate and distinct from the real `FSSgam` package. This repo's R
package scaffolding (`DESCRIPTION`, `NAMESPACE`, `R/zzz.R`, empty `man/`)
exists solely so that `pkgdown::build_site()` has a namespace to load while
it renders the site; `R/zzz.R` exports nothing on purpose.

The real package lives in a **separate repo**:
https://github.com/beckyfisher/FSSgam_package — this repo's vignettes
`library(FSSgam)` and call into it, but its source code is out of scope
here. Do not edit, clone into, or suggest changes to that repo from this
session; if something here surfaces what looks like a bug in the package
itself, report it back to the user rather than patching it in this repo.

---

## 2. What This Repo Does

This is the permanent, citable documentation companion to:

> Fisher R, Wilson SK, Sin TM, Lee AC, Langlois TJ (2018). A simple function
> for full-subsets multiple regression in ecology with R. *Ecology and
> Evolution*, 8(12), 6104–6113. <https://doi.org/10.1002/ece3.4134>

It hosts the pkgdown site at <https://beckyfisher.github.io/FSSgam>: FAQs,
case-study vignettes reproducing the published analyses, and extra usage
examples. It is the reference researchers actually read and cite — treat
its case-study **scientific content** (numbers, figures, conclusions) as a
fixed historical record. Only the *syntax* used to call the package, and
clearly documentation-only text, should ever change here.

### Recently completed: v1.0.0 rename + appendix integration (2026-06-22)

`FSSgam_package` was modernised for CRAN submission: `generate.model.set()`
and `fit.model.set()` became deprecated aliases (each calls `.Deprecated()`
then forwards) for `generate_model_set()` and `fit_model_set()`.
`full.subsets.gam()` was **not** renamed. That rename has since been merged
into `FSSgam_package`'s `master`, so the default
`remotes::install_github("beckyfisher/FSSgam_package")` now picks it up —
no branch pin needed to get the new functions.

At the same time, all vignette code/prose was updated to the new names, the
richer Background/Methods/Results & Discussion write-ups from the original
published Appendix S1–S5 were folded into the case-study vignettes
(replacing the placeholder discussion that used to be there), and two new
reference vignettes were added (`function-arguments.Rmd`,
`function-outputs.Rmd`, from Appendix S1/S2). `references.bib` grew from 11
to 90 entries to support this. The site also gained a dual root/`dev`
pkgdown deployment (see Section 3) and a "Reference" navbar link out to
`FSSgam_package`'s own function reference site.

Full detail of how/why is in `prompts/vignette-integration-rename.md` and
the `dev` branch git history — check there before re-deriving any of this
from scratch. If a future package update reintroduces deprecation warnings
or renames something else, the same approach applies: re-grep the vignettes
for the old name, update call-sites and prose, then render and verify (see
Section 5).

---

## 3. Key Files and Structure

```
vignettes/
  case-study-1.Rmd          # Ningaloo reef fish (zoning/habitat) — calls generate_model_set/fit_model_set
  case-study-2.Rmd          # Langlois et al. soft-sediment reanalysis — same; has a dated "** note **"
                            # block tracking its own update history — append to it, don't replace it
  case-study-3.Rmd          # intertidal gastropod reproductive patterns — same
  extra-examples.Rmd        # binomial/uGamm example via gamm4 — same
  faq.Rmd                    # FAQ; includes a Q&A entry explaining the snake_case rename
  function-arguments.Rmd     # reference vignette built from the original Appendix S1, split by
                            # whether each argument belongs to generate_model_set() or fit_model_set()
  function-outputs.Rmd       # same, for Appendix S2 / function outputs
R/zzz.R                  # @noRd placeholder; package exports nothing, exists only for pkgdown
DESCRIPTION              # Package: FSSgam.docs (intentionally not "FSSgam") — do not rename.
                          # Version is intentionally a dev-style 4-part number (currently
                          # 1.11.0.9000) so pkgdown's `development: mode: auto` (see _pkgdown.yml)
                          # routes dev-branch builds to docs/dev/ instead of the root site — don't
                          # "fix" this back to a plain release-style version.
NAMESPACE                # auto-generated, currently has no exports — do not hand-edit
man/                      # empty — no exported functions to document
_pkgdown.yml              # site nav/config — do not restructure without being asked. Has
                          # `development: mode: auto`, and a custom "Reference" navbar component
                          # linking out to FSSgam_package's own pkgdown reference site.
.github/workflows/pkgdown.yaml  # builds + deploys the site. Triggers on push to main/master/dev.
                          # Installs FSSgam_package from GitHub; the install step picks
                          # FSSgam_package's `dev` ref only when *this* repo's branch is `dev` —
                          # currently redundant now that FSSgam_package's dev is merged to its
                          # master, but harmless; ask before removing it, since it exists to
                          # support the dual root/dev pkgdown-site setup, not just the rename.
                          # Deploys with `keep_files: true` so the root and /dev/ builds don't
                          # clobber each other on gh-pages.
references.bib            # bibtex entries cited from the case-study/function-reference vignettes.
                          # Keyed `{Surname}{Year}{ShortCamelCaseTitle}`, e.g. `Wilson2012Ningaloo`.
inst/CITATION              # paper citation metadata — do not change
README.md                 # repo summary, install instructions, citation, license. Now also has a
                          # one-sentence pointer near the top to FSSgam_package's function
                          # reference — the DOI/citation block itself is still off-limits.
superceded/                # manual pre-edit backups for untracked files (already in active use —
                          # see x_function_check_correlations_v1.00.R); follow this existing
                          # convention for any other untracked file you edit
ignore/                    # scratch workspace, gitignored
scratch/                   # scratch workspace, gitignored — holds working material like source
                          # appendix text or PR-description drafts; not part of the published site
docs/                      # pkgdown build output, gitignored — never hand-edit or commit into this
publication/                # supplementary files from the original paper — read-only, do not touch
```

---

## 4. Known Constraints or Decisions

- **Do not rename the `FSSgam.docs` package id**, and do not add exported
  functions to it. It exists purely so pkgdown has something to load.

- **Do not edit, clone into, or push to `beckyfisher/FSSgam_package`** from
  this repo's session. This repo only consumes that package.

- **The site installs the real package from GitHub at build time**
  (`.github/workflows/pkgdown.yaml`). The install step is branch-conditional:
  when *this* repo's branch is `dev` it installs `FSSgam_package`'s `dev`
  ref; otherwise it installs `FSSgam_package`'s default branch. Since the
  snake_case rename is now merged into `FSSgam_package`'s `master`, both
  paths currently resolve to equivalent functions — the conditional isn't
  load-bearing for the rename itself any more, but it's still what makes the
  dual root/`dev` pkgdown-site setup work (see `_pkgdown.yml`'s
  `development: mode: auto` and `DESCRIPTION`'s dev-style `Version`). Don't
  remove it without checking with the user first.

- **Do not alter historical results, figures, or scientific conclusions**
  in the case-study vignettes — a standing rule, not tied to any one task.
  Only function-call syntax and documentation-only text (prose, citations,
  cross-references) are ever in scope for changes here.

- **Vignettes here are R Markdown (`.Rmd`), not Quarto (`.qmd`) — this
  overrides the global CLAUDE.md preference for Quarto.** `DESCRIPTION`
  declares `VignetteBuilder: knitr` and `SystemRequirements: pandoc`; the
  whole pkgdown pipeline is built around knitr/rmarkdown. Do not convert
  these to `.qmd` — that's a much larger, unrequested structural change
  that would break the existing build.

- **Do not change `inst/CITATION`, the license text, or the DOI/citation
  block in `README.md`.**

- **Apache License 2.0** — retain existing licence; do not change it.

- **Line endings:** files in this repo can show as locally modified with no
  real content change (CRLF vs. the LF committed history) —
  `git diff --stat` showing close to equal insertions/deletions across
  every line of a file is the usual signature, but it isn't reliable on its
  own: a genuinely tiny one-line edit can *also* render as a whole-file
  diff in this environment without any CRLF involved. Confirm with a
  byte-level check (`cmp`/`xxd` on a sampled unrelated line, or
  `diff <(git show HEAD:file) file`) before assuming corruption either way.
  If real CRLF is found, strip it with
  `tr -d '\r' < file > tmp && mv tmp file` run as an **individual command
  per file** — `for`/`while` loops in Bash have not reliably written files
  in this environment (they execute and report success but silently no-op
  on the actual write); single standalone commands do work.

- **Render each vignette in its own `Rscript` process when verifying
  changes.** Calling `rmarkdown::render()` more than once inside the same R
  session has segfaulted in this environment.

- **Git awareness (same rule as the global CLAUDE.md):** tracked-and-clean
  → edit freely; tracked-but-uncommitted → warn which files are affected and
  wait for confirmation; untracked → back up to `superceded/` first
  (preserving the original filename) before editing — this repo already
  uses that exact convention (see `superceded/x_function_check_correlations_v1.00.R`).

- **Ask before installing any new R package.** Rendering the vignettes
  needs everything in `DESCRIPTION`'s `Imports`/`Suggests` (`doSNOW`,
  `MuMIn`, `gamm4`, `mgcv`, `nnet`, plus `car`, `doBy`, `dplyr`, `ggplot2`,
  `gplots`, `grid`, `gridExtra`, `purrr`, `RColorBrewer`, `tibble`, `tidyr`,
  `tidyverse`, `rmarkdown`, `knitr`, `R.rsp`) plus the `FSSgam` package
  itself (installed separately, see above) — check what's actually
  installed before assuming a missing-package error means something is
  broken.

- **Show the diff/summary of what changed before committing**, even when
  already asked to make the change — a confirmed preference, asked for
  explicitly more than once.

- **Never commit or push without being explicitly asked.**

---

## 5. Verifying changes when the package updates

When `FSSgam_package` changes in a way that affects this repo (a rename, a
new argument, a behavioural change), the workflow that worked for the
v1.0.0 rename:

1. Re-grep the vignettes for the old name/usage — don't trust a past
   session's call-site counts without re-checking.
2. Update call-sites and prose (bolded mentions, `?old.name` →
   `?new_name` cross-references), leaving argument names and local
   variable names alone unless the package itself changed them.
3. Check what's actually installed locally first
   (`packageVersion("FSSgam")`, `ls(getNamespace("FSSgam"))`) — a
   missing-function error usually just means the wrong branch/version is
   installed, not that something is broken in this repo.
4. Render each changed `.Rmd` in its **own** `Rscript` process (see
   Section 4) and confirm: no errors, no deprecation warnings, and — for
   any vignette with a `bibliography:` field — no unresolved `[@key]`
   citations in the rendered HTML.
5. If rendering ever produces **different numeric results** than before
   (not just a different warning), stop and flag it to the user rather
   than updating the published numbers — that indicates a behavioural
   regression in `FSSgam_package`, to be fixed there, not papered over
   here.
6. Summarise afterwards: which files changed and what moved (code vs.
   prose vs. new content), confirmation every changed vignette re-rendered
   cleanly, which package branch you tested against vs. what the committed
   CI workflow actually installs, and anything that looked like a genuine
   behavioural difference rather than a cosmetic rename.

---

## 6. Prompt Log

Session logs for this project are in `prompts/`. Create the folder if it
does not exist. Use a short kebab-case descriptor as the filename for each
session (e.g. `vignette-rename-update.md`). If a file with that name
already exists, append to it.

Log format for each session entry:

```
## Session: <descriptor>
Date: <YYYY-MM-DD>
Model: <model name and version>

### Prompts and Responses

**User:** <prompt text>

**Claude:** <summary of response — full code blocks where relevant, prose summarised>

---
```
