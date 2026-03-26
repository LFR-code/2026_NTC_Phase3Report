# Todo: Results figures and tables

Handoff for Scarlett. All figures/tables are based on the
presentation in `refs/NTC_MPtests.pdf`. The report text is
largely written but contains no figures or tables yet.

## Model data (submodule)

Model outputs are in the `Models/` submodule
(`git@github.com:LFR-code/2026_NTC_Phase3Report-Models.git`,
private, LFS-tracked). Clone with:
```bash
git clone --recurse-submodules \
  git@github.com:LFR-code/2026_NTC_Phase3Report.git
```

### Data inventory

**OM history fits** (`Models/OMHistory/`):
- `fit_WCVI2023_OM_mle/` -- DDM OM
- `fit_1S_2024/` -- 1S-predM OM
- `fit_3S_2024/` -- 3S-predM OM

**Simulations** (`Models/Sims/`):
- `sim_DDM_OM-SISCA_DIM_minE0847_MP/` -- DDM minE MP sim
- `sim_1SpredMOM_minEMP_full/` -- 1S-predM minE MP sim
- `3S_allocTest/` -- 3S-predM with 3 allocation MPs:
  - `...-1S_DIM_MP1/` (SI allocation)
  - `...-1S_DIM_MP2/` (Catch allocation)
  - `...-1S_DIM_MP3/` (All PFMA 23 allocation)
- `3S_allocTest/statistics/fullStatTable.csv` --
  performance metrics

## Source code repos for figure functions

- **DDM sim plots**: adapt from
  `git@github.com:LFR-code/2026_NTC_BCSRIF-HerringMP-MS3.git`
- **PredM sim plots (1S and 3S)**: adapt from
  `git@github.com:LFR-code/2026_NTC_BCSRIF-HerringMP-ms3R-MICE.git`

## Figure code (new R scripts needed)

- [x] **Write `R/plotMethods.R`** -- methods section figures
  (yield curve, historical SSB comparison, HCR diagram).
  Uses OM history fit objects. Adapt DDM yield curve code
  from the MS3 repo.
- [x] **Write `R/plotResults.R`** -- results section figures
  (simulation projections, retrospective assessment fits,
  effective harvest rate scatters). Adapt DDM plots from
  MS3 repo, predM plots from ms3R-MICE repo. Use shared
  plotting functions across OMs where possible.
- [x] **Write `R/tableMaker.R`** -- allocation proportion
  tables and performance metrics tables. Uses
  `fullStatTable.csv` and sim objects.
- [x] **Source scripts in `index.Rmd`** -- update the TODO
  stubs (lines 101-109) to load the correct `.RData` files
  from `Models/` and source the new R scripts.

## Prerequisites (index.Rmd)

- [x] **Connect data objects** -- uncomment/update the
  `source()` and `load()` calls in `index.Rmd`
  (lines 101-109). The correct paths are now:
  - `Models/OMHistory/fit_WCVI2023_OM_mle/fit_WCVI2023_OM_mle.rds`
  - `Models/OMHistory/fit_1S_2024/fit_1S_2024.RData`
  - `Models/OMHistory/fit_3S_2024/fit_3S_2024.RData`
  - `Models/Sims/sim_DDM_OM-SISCA_DIM_minE0847_MP/sim_DDM_OM-SISCA_DIM_minE0847_MP.RData`
  - `Models/Sims/sim_1SpredMOM_minEMP_full/sim_1SpredMOM_minEMP_full.RData`
  - `Models/Sims/3S_allocTest/sim_parBattestAlloc-3S_predM_OM-1S_DIM_MP1/...RData`
  - `Models/Sims/3S_allocTest/sim_parBattestAlloc-3S_predM_OM-1S_DIM_MP2/...RData`
  - `Models/Sims/3S_allocTest/sim_parBattestAlloc-3S_predM_OM-1S_DIM_MP3/...RData`

## Methods figures (02-methods.Rmd)

- [x] **Yield curve figure** (slide 2) -- DDM yield curve
  showing Spawning Biomass vs Yield with labelled reference
  points ($0.4B_{MSY}$, LRP, $0.8B_{MSY}$, $B_{MSY}$,
  $B_0$). Place after the *SISCAH-DDM* OM subsection.

- [x] **Historical SSB comparison figure** (slide 3) -- Time
  series of aggregate SSB from all three OMs (DDM,
  1S-predM, 3S-predM) plus per-PFMA lines for 3S-predM,
  with $B_0$ and $B_{MSY|DDM}$ reference lines. Place at
  end of the Operating Models section or as a standalone
  comparison.

- [x] **HCR diagram** (slide 4) -- Two-panel plot: harvest
  rate ($U$) vs biomass (top) and TAC vs biomass (bottom),
  showing the closed/open regions and $B_{ref}$/$U_{ref}$
  labels. Place in the HCR subsection.

## Results figures -- DDM OM (03-results.Rmd)

- [x] **DDM simulation projections** (slide 6) -- Three-panel:
  (i) SSB with projection envelope + $B_0$, LRP, $B_{ref}$
  lines, (ii) commercial catch by fleet (reduction,
  seineRoe, gillnet), (iii) SOK/FSC catch. Place in
  "Reference case: DDM" subsection.

- [x] **DDM retrospective assessment fits** (slide 7) -- OM
  truth (red) vs SISCA final (blue) vs earlier estimates
  (grey) + Index/q points, with dashed line at projection
  start. Place after DDM simulation projections.

- [x] **DDM effective harvest rate scatter** (slide 8) --
  Two-panel: (i) $U$ vs $B$ scatter with HCR line overlay,
  (ii) TAC vs $B$ scatter with HCR line overlay. Place
  after DDM retrospective fits.

## Results figures -- 1S-predM OM (03-results.Rmd)

- [x] **1S-predM simulation projections** (slide 10) -- Same
  three-panel layout as DDM but for 1S-predM OM. Place in
  "Robustness test: 1S-Pred" subsection.

- [x] **1S-predM retrospective assessment fits** (slide 11) --
  Same layout as DDM version. Place after 1S-predM
  projections.

- [x] **1S-predM effective harvest rate scatter** (slide 12) --
  Same layout as DDM version. Place after 1S-predM
  retrospective fits.

## Results figures and tables -- 3S-predM OM (03-results.Rmd)

- [x] **Allocation proportions table(s)** (slides 14-15) --
  Table showing 2024 and long-term average allocation
  proportions (SI, Catch, All 23) by PFMA. Can be one
  combined table or two. Place at start of the 3S-predM
  spatial allocation subsection.

- [x] **3S-predM simulation projections** (slide 17) --
  Three-panel with per-PFMA SSB lines, commercial catch by
  PFMA, SOK/FSC by PFMA. Place in "Aggregate performance"
  subsection.

- [x] **3S-predM retrospective assessment fits** (slide 18) --
  Same layout as DDM/1S-predM versions. Place after
  3S-predM projections.

- [x] **3S-predM effective harvest rate scatter** (slide 19) --
  Same layout as DDM/1S-predM versions. Place after
  3S-predM retrospective fits.

- [x] **Performance metrics tables** (slides 20-21) -- Two
  tables showing probability of low biomass and probability
  of low catch by PFMA and allocation method (SI, Catch,
  All 23). Place in "Performance by PFMA" subsection.

## Cleanup

- [x] **Remove the TODO comment** at `03-results.Rmd` line 116
  once figures/tables are inserted.
