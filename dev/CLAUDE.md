# Project: FSSgam (documentation companion repo)

This file provides **project-specific context only** — repo structure,
dependencies, the current task, and prompt logging. It does not impose
personal style preferences on collaborators.

**For Claude:** If a `CLAUDE.md` exists in the parent directory of this
repo, read it before proceeding — it may contain user-specific
environment and style preferences. This repo-level file takes precedence
over any parent-level file where they conflict, **except** where a
parent-level default is explicitly overridden below (Section 5 notes one
such override).

------------------------------------------------------------------------

## 1. Repo Type

**Documentation/pkgdown companion repo**, not a deployable R package in
the normal sense. `DESCRIPTION` declares `Package: FSSgam.docs` — that
name is deliberate and distinct from the real `FSSgam` package. This
repo’s R package scaffolding (`DESCRIPTION`, `NAMESPACE`, `R/zzz.R`,
empty `man/`) exists solely so that
[`pkgdown::build_site()`](https://pkgdown.r-lib.org/reference/build_site.html)
has a namespace to load while it renders the site; `R/zzz.R` exports
nothing on purpose.

The real package lives in a **separate repo**:
<https://github.com/beckyfisher/FSSgam_package> — this repo’s vignettes
[`library(FSSgam)`](https://github.com/beckyfisher/FSSgam_package) and
call into it, but its source code is out of scope here. Do not edit,
clone into, or suggest changes to that repo from this session; if
something here surfaces what looks like a bug in the package itself,
report it back to the user rather than patching it in this repo.

------------------------------------------------------------------------

## 2. What This Repo Does

This is the permanent, citable documentation companion to:

> Fisher R, Wilson SK, Sin TM, Lee AC, Langlois TJ (2018). A simple
> function for full-subsets multiple regression in ecology with R.
> *Ecology and Evolution*, 8(12), 6104–6113.
> <https://doi.org/10.1002/ece3.4134>

It hosts the pkgdown site at <https://beckyfisher.github.io/FSSgam>:
FAQs, case-study vignettes reproducing the published analyses, and extra
usage examples. It is the reference researchers actually read and cite —
treat its case-study **scientific content** (numbers, figures,
conclusions) as a fixed historical record. Only the *syntax* used to
call the package should ever change here, and only when the package
itself changes.

### Current task: update vignettes for the FSSgam_package v1.0.0 rename

`FSSgam_package` was just modernised for CRAN submission. The headline
change for this repo: `generate.model.set()` and `fit.model.set()` are
now deprecated aliases (each calls
[`.Deprecated()`](https://rdrr.io/r/base/Deprecated.html) then forwards
to the new name) for `generate_model_set()` and `fit_model_set()`. The
old names still work, so nothing here is broken today — but every
vignette currently calls the deprecated names, which means re-rendering
any of them right now would emit
[`.Deprecated()`](https://rdrr.io/r/base/Deprecated.html) warnings into
the published output. `full.subsets.gam()` was **not** renamed; leave
every mention of it alone.

Files that need updating (counts are call-sites of the old names found
at the time this file was written — re-grep before relying on these):

| File                           | Old-name occurrences |
|--------------------------------|----------------------|
| `vignettes/case-study-1.Rmd`   | 2                    |
| `vignettes/case-study-2.Rmd`   | 4                    |
| `vignettes/case-study-3.Rmd`   | 4                    |
| `vignettes/extra-examples.Rmd` | 6                    |
| `vignettes/faq.Rmd`            | 2                    |

For each file:

1.  **Code chunks** — replace `generate.model.set(` →
    `generate_model_set(` and `fit.model.set(` → `fit_model_set(`.
    Arguments are unchanged (the package kept all public argument names,
    e.g. `use.dat`, `test.fit`, `pred.vars.cont`, exactly as they were —
    only the function names changed). Do not touch variable names like
    `model.set` or `out.list`; those are local to the vignettes, not
    part of the package API.
2.  **Prose** — update inline bolded mentions (`**generate.model.set**`,
    `**fit.model.set**`) and cross-reference links
    (`?generate.model.set` → `?generate_model_set`) to the new names.
    `case-study-2.Rmd` has a short narrative passage (near the top)
    explaining why `generate.model.set`/ `fit.model.set` superseded
    `full.subsets.gam` — update the function names there too, but
    preserve the narrative.
3.  **`faq.Rmd`** — beyond the mechanical rename, consider adding a
    short new Q&A entry explaining the rename itself (old names still
    work via deprecation warning; new code should call the snake_case
    names). This is genuinely useful content given the change, not just
    a find-and-replace.
4.  **Render and verify** — after editing, render each changed `.Rmd`
    (e.g. `rmarkdown::render("vignettes/case-study-1.Rmd")` or
    `devtools::build_rmd("vignettes/case-study-1.Rmd")`) and confirm: it
    renders without error, and the output contains **no**
    `'generate.model.set' is deprecated` /
    `'fit.model.set' is deprecated` warnings. A leftover warning means a
    call site was missed.
5.  **This requires the updated package to actually be installed.** See
    Section 4 — as of this writing the rename lives on the `dev` branch
    of `FSSgam_package`, not yet merged to `master`, so the default
    `remotes::install_github("beckyfisher/FSSgam_package")` will **not**
    pick it up. Install explicitly from the branch to test locally;
    don’t leave that branch pin in any committed workflow file without
    checking with the user first (see Section 4).

If rendering a case-study vignette ever produces **different numeric
results** than before (not just a different deprecation warning), stop
and flag it to the user rather than updating the published numbers —
that would indicate a behavioural regression in `FSSgam_package`, which
needs to be fixed in that repo, not papered over here.

------------------------------------------------------------------------

## 3. Key Files and Structure

    vignettes/
      case-study-1.Rmd      # Ningaloo reef fish (zoning/habitat) — uses generate.model.set/fit.model.set
      case-study-2.Rmd      # Langlois et al. soft-sediment reanalysis — same, plus narrative prose to update
      case-study-3.Rmd      # intertidal gastropod reproductive patterns — same
      extra-examples.Rmd    # binomial/uGamm example via gamm4 — same, 3 separate call sites
      faq.Rmd                # FAQ; prose-only mentions, good place to add a rename note
    R/zzz.R                  # @noRd placeholder; package exports nothing, exists only for pkgdown
    DESCRIPTION              # Package: FSSgam.docs (intentionally not "FSSgam") — do not rename
    NAMESPACE                # auto-generated, currently has no exports — do not hand-edit
    man/                      # empty — no exported functions to document
    _pkgdown.yml              # site nav/config — do not restructure without being asked
    .github/workflows/pkgdown.yaml  # installs FSSgam_package from GitHub, then builds + deploys the site
    inst/CITATION              # paper citation metadata — do not change
    README.md                 # repo summary, install instructions, citation, license — function names
                              # are not mentioned here, no change needed for this task
    superceded/                # manual pre-edit backups for untracked files (already in active use —
                              # see x_function_check_correlations_v1.00.R); follow this existing
                              # convention for any other untracked file you edit
    ignore/                    # scratch workspace, gitignored
    docs/                      # pkgdown build output, gitignored — never hand-edit or commit into this
    publication/                # supplementary files from the original paper — read-only, do not touch

------------------------------------------------------------------------

## 4. Known Constraints or Decisions

- **Do not rename the `FSSgam.docs` package id**, and do not add
  exported functions to it. It exists purely so pkgdown has something to
  load.

- **Do not edit, clone into, or push to `beckyfisher/FSSgam_package`**
  from this repo’s session. This repo only consumes that package.

- **The site installs the real package from GitHub at build time**
  (`.github/workflows/pkgdown.yaml` →
  `remotes::install_github("beckyfisher/FSSgam_package")`, no `ref=`
  pin, so it tracks that repo’s default branch). At the time of writing,
  the snake_case rename is on that repo’s `dev` branch, not yet merged
  to `master`. To render against the new functions locally, install with
  `remotes::install_github("beckyfisher/FSSgam_package", ref = "dev")`
  first. **Do not commit a `ref = "dev"` pin into the workflow file** —
  that’s a local testing step only; ask the user before changing what
  the CI workflow installs, since flipping it back at the wrong time
  would silently re-deploy docs built against the old function names.

- **Do not alter historical results, figures, or scientific
  conclusions** in the case-study vignettes. Only function-call syntax
  and clearly documentation-only text (FAQ prose, cross-references) are
  in scope for the current task.

- **Vignettes here are R Markdown (`.Rmd`), not Quarto (`.qmd`) — this
  overrides the global CLAUDE.md preference for Quarto.** `DESCRIPTION`
  declares `VignetteBuilder: knitr` and `SystemRequirements: pandoc`;
  the whole pkgdown pipeline is built around knitr/rmarkdown. Do not
  convert these to `.qmd` as part of this task or any other — that’s a
  much larger, unrequested structural change that would break the
  existing build.

- **Do not change `inst/CITATION`, the license text, or the DOI/citation
  block in `README.md`.**

- **Apache License 2.0** — retain existing licence; do not change it.

- **Line endings:** like `FSSgam_package`, files in this repo can show
  as locally modified with no real content change (CRLF vs. the LF
  committed history) — `git diff --stat` showing equal
  insertions/deletions across every line of a file is the signature. If
  you see this, check before editing: commit the line-ending
  normalization as its own separate commit first (same approach used in
  `FSSgam_package`), so it doesn’t get tangled with the content changes
  from this task.

- **Git awareness (same rule as the global CLAUDE.md):**
  tracked-and-clean → edit freely; tracked-but-uncommitted → warn which
  files are affected and wait for confirmation; untracked → back up to
  `superceded/` first (preserving the original filename) before editing
  — this repo already uses that exact convention (see
  `superceded/x_function_check_correlations_v1.00.R`).

- **Ask before installing any new R package.** Rendering the vignettes
  needs everything in `DESCRIPTION`’s `Imports`/`Suggests` (`doSNOW`,
  `MuMIn`, `gamm4`, `mgcv`, `nnet`, plus `car`, `doBy`, `dplyr`,
  `ggplot2`, `gplots`, `grid`, `gridExtra`, `purrr`, `RColorBrewer`,
  `tibble`, `tidyr`, `tidyverse`, `rmarkdown`, `knitr`, `R.rsp`) plus
  the `FSSgam` package itself (installed separately, see above) — check
  what’s actually installed before assuming a missing-package error
  means something is broken.

- **Never commit or push without being explicitly asked.**

------------------------------------------------------------------------

## 5. After finishing the vignette update

Summarise clearly: - Which files changed and what specifically moved
(code call-sites vs. prose vs. new FAQ content). - Confirmation that
every changed vignette re-rendered cleanly with no deprecation
warnings. - Whether you tested against the package installed from `dev`
or from `master` — and a reminder that the committed workflow still
installs from `master`, so the rendered docs won’t reflect these changes
on the live site until `FSSgam_package`’s `dev` branch is merged. -
Anything you found that looks like a genuine behavioural difference (not
just a renamed function) between the old and new package — flag, don’t
fix here.

------------------------------------------------------------------------

## 6. Prompt Log

Session logs for this project are in `prompts/`. Create the folder if it
does not exist. Use a short kebab-case descriptor as the filename for
each session (e.g. `vignette-rename-update.md`). If a file with that
name already exists, append to it.

Log format for each session entry:

    ## Session: <descriptor>
    Date: <YYYY-MM-DD>
    Model: <model name and version>

    ### Prompts and Responses

    **User:** <prompt text>

    **Claude:** <summary of response — full code blocks where relevant, prose summarised>

    ---
