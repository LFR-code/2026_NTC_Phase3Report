# NTC-BCSRIF Phase 3 Working Paper Review

**Document**: ?aayaaqa (Herring) Spawn Dynamics -- Phase 3 Working Paper
(NTC-BCSRIF, Landmark Fisheries Research for Uu-a-thluk Fisheries)
**Reviewer**: Dr. Sarah Chen
**Date**: 2026-03-26
**Files reviewed**: `index.Rmd`, `00-execSummary.Rmd`, `01-intro.Rmd`,
`02-methods.Rmd`, `03-results.Rmd`, `04-discussion.Rmd`

---

## Overall Assessment

This is a well-structured working paper that clearly communicates the
Phase 3 MSE work. The narrative arc is coherent: the introduction
establishes the management problem, the methods describe a sensible
closed-loop framework, the results are presented systematically across
three OMs, and the discussion correctly identifies OM choice as the
dominant driver of performance. The writing is concise and technically
appropriate for a project working paper audience. The document is close
to ready for internal review but has several gaps and issues -- some
scientifically substantive, some presentational -- that should be
addressed before circulation.

---

## Critical Issues (Must Address)

### 1. The SISCAH-DDM $B_{MSY}$ value is inconsistently reported

In `02-methods.Rmd` (§Methods/SISCAH-DDM subsection), $B_{MSY}$ is
stated as "approximately 17.4 kt." In the HCR subsection (§Methods/
Step 1), $B_{ref}$ is described as "approximately 17.35 kt." In
`03-results.Rmd` (§Results/Reference case), $B_{MSY}$ appears as
"approximately 17.35 kt." In the figure caption for the HCR diagram
(`capHCRdiagram`), the value is 17.35 kt, and `plotHCRDiagram()` is
called with `Bref = 17.35`. In `plotResults.R`, the HCR scatter
caption also states $U_{ref} = 8.47\%$, but the methods text states
$U_{ref} = 8.5\%$.

These small numerical discrepancies -- 17.4 kt vs. 17.35 kt, and
8.5\% vs. 8.47\% -- create unnecessary confusion for readers
cross-checking figures against text. Settle on a single reported
precision for each parameter (I'd recommend two decimal places:
17.35 kt and 8.47\%) and apply it consistently throughout the text,
figure captions, and table headers.

### 2. Performance metrics are only computed for the 3S-Pred OM

The methods section (§Performance metrics) describes two metrics
(probability of low biomass, probability of low catch) and frames them
as the primary evaluation tool. However, in the results, these metrics
are only presented for the *SISCAH-3S-Pred* OM (the `makePerfTables()`
function reads from `fullStatTable.csv` in the `3S_allocTest`
directory). The *SISCAH-DDM* and *SISCAH-1S-Pred* results sections
contain no tables or quantitative metrics -- only qualitative
statements like "biomass is maintained above the LRP approximately
75\% of the time."

This is a meaningful inconsistency. Either: (a) compute the same
metrics for all three OMs so results are comparable, or (b) revise
the methods section to clearly state that these formal metrics are
only computed for the three-stock OM because that is the scenario
where Stat Area-level comparisons are relevant. Option (b) is
probably more defensible given the intent of the analysis, but the
methods text as written implies broader coverage.

### 3. Assessment model mismatch is introduced without prior setup

The MP uses a *SISCAH-DIM* (density-independent mortality) assessment
model, but this acronym appears first in `03-results.Rmd` (§Reference
case) without any earlier introduction. The methods section describes
the aggregate assessment step but never explicitly names the
assessment model used in the MP as *SISCAH-DIM*, nor explains why DIM
rather than DDM was chosen as the assessment model (despite the DDM OM
being the reference case). This is a meaningful design choice: using
DIM in the assessment while the reference OM is DDM means the models
are already structurally mismatched even in the reference case, which
bears directly on the observed over-estimation bias discussed in the
results and discussion.

Recommend adding one or two sentences to §Methods/Step 1 naming the
assessment model as *SISCAH-DIM* and briefly explaining the rationale
for this choice (e.g., DIM is the simpler model used for projection
because structural constraints on future mortality are avoided, or
because the intent was to test a simpler, more portable assessment).

### 4. The simulation projection figure for 3S-Pred shows per-area SSB
   but the text describes only aggregate performance

The caption for `fig-3Sprojections` (`cap3Sprojections`) describes a
figure with per-PFMA SSB lines, commercial catch by PFMA, and SOK/FSC
by PFMA (consistent with "slide 17" in the TODO file). However, the
`plotSimProj()` function in `plotResults.R` produces an *aggregate*
SSB panel (summing across stocks) with no per-area breakdown. The
figure caption therefore describes a figure the function does not
produce. Either the caption needs to be updated to match what
`plotSimProj()` actually generates, or a separate per-area plot
function is needed for the 3S-Pred OM to show Stat Area-level
dynamics. This is important because one of the central questions
under the 3S-Pred OM is whether individual Stat Areas are at risk
-- aggregate plots may obscure relevant Stat Area-level variation.

---

## Important Improvements (Strongly Recommended)

### Background / Introduction (`01-intro.Rmd`)

**I-1. The Phase 1 conceptual model is mentioned but never connected
to Phase 3.** The introduction references a "conceptual simulation
model... developed to represent WCVI herring population dynamics
across nine Nuu-chah-nulth (NCN) spawning areas" in Phase 1. The
*SISCAH-9S-Pred* model from Phase 2 is also mentioned. Neither
appears again in the document. Since Phase 3 focuses entirely on the
three-stock framework, the introduction creates a false impression
that the nine-area model might be relevant here. A single sentence
clarifying that the Phase 3 work builds exclusively on the three-stock
framework (per the Phase 2 recommendation) would tidy this up.

**I-2. The management context section should cite the specific DFO
reference point framework.** The statement "a re-opening of the
fishery would require a management procedure that defines how stock
status is assessed and how a TAC is determined" is accurate but
under-supported. DFO has a Precautionary Approach (PA) policy that
mandates specific reference points and HCR structures for commercial
fisheries. A citation (e.g., DFO 2009 PA policy, or the existing
WCVI herring IFMP citation `dfo2024_herringIFMP`) would strengthen
the framing that the proposed minimum escapement HCR is consistent
with the PA framework.

**I-3. The objectives section lists three OMs but does not state the
allocation methods as a co-equal objective.** The introduction lists
three OMs clearly but treats allocation testing as a subordinate
detail ("We test three allocation methods..."). The executive summary
gives allocation testing equal billing. Aligning these would
improve coherence.

### Methods (`02-methods.Rmd`)

**M-1. The SISCAH-DDM reference point description uses inconsistent
units.** The text states "unfished biomass is approximately 42 kt,
$B_{MSY}$ is approximately 17.4 kt... and the maximum sustainable
yield is approximately 1.6 kt." The "kt" abbreviation is used
throughout, but maximum sustainable yield for a stock of this size
is typically in the range of hundreds to low thousands of tonnes.
If the units are correct (1.6 kt = 1,600 t), this should be
confirmed by checking the yield curve figure; the executive summary
correctly states "approximately 1,600 t." No unit conversion
inconsistency exists, but the text would benefit from stating
"1.6 kt (1,600 t)" at first mention to reduce potential reader
confusion given that the SOK/FSC panel in the simulation figures
uses tonnes (t) rather than kt. Note that `plotResults.R` applies a
`sok_scale <- 1000` multiplier, reinforcing that units are mixed
across figures.

**M-2. The observation error model for simulated spawn indices is not
described.** The simulation design section states "the OM generates
data (spawn indices with observation error)" but gives no detail on
how observation error is implemented. Since this affects how well the
assessment tracks true biomass -- a result that is discussed at length
-- even a one-sentence description of the error model (log-normal,
estimated CV, etc.) is warranted.

**M-3. The catch allocation among commercial fleets is mentioned but
not quantified.** The methods state that "the allocated TAC is further
distributed among commercial fleets (seine roe, gillnet, and
spawn-on-kelp or bough) according to fixed fleet allocation
proportions." No values are provided for these proportions, nor is
there a reference to where they come from (historical data?
regulatory allocation?). This detail matters because the DDM results
section mentions "Commercial catches (split among roe fisheries)"
implying more than two roe fleets. Either provide the proportions in
a table or appendix, or add a citation to the source of these values.

**M-4. The number of stochastic replicates (100) should be justified
briefly.** In MSE applications, 100 replicates is at the lower end of
what is typically used. A brief statement that this was assessed as
adequate to characterise the distribution of outcomes (e.g., by
checking convergence of performance metric quantiles) would
pre-empt questions from readers.

**M-5. The "historical minimum" baseline for performance metrics is
ambiguous.** The metrics are defined relative to "historical minimum
observed biomass" and "historical minimum non-zero catch." It is not
clear whether "historical" refers to the OM conditioning period
(1951--2024) or some shorter reference period. The choice matters
because the minimum over a 73-year period that includes low-stock
years could be very low, making it easy to achieve. Clarify what
period defines the historical minimum.

**M-6. The HCR equation and caption use slightly different notation.**
The HCR equation uses $U_{ref}$ for the maximum harvest rate. The
figure caption (`capHCRdiagram`) sets "Parameters: $B_{ref} = 17.35$
kt, $U_{ref} = 8.5\%$" but -- as noted under Critical Issue 1 --
the code uses 8.47\%. Beyond the numerical inconsistency, stating
the parameter as a percentage in the caption but as a decimal in
the equation (0.0847) is fine but could be signposted: e.g.,
"$U_{ref} = 0.0847$ (8.47\%)."

### Results (`03-results.Rmd`)

**R-1. The 1S-Pred results section uses imprecise language.** The
statement "lots of simulated outcomes show closures when $B > B_{ref}$"
(§Robustness test: SISCAH-1S-Pred) is informal language that should
not appear in a technical working paper. Suggested revision:
"A substantial proportion of simulated outcomes result in fishery
closures even when true biomass exceeds $B_{ref}$, consistent with
the lower equilibrium biomass under the predation OM."

**R-2. The discussion of assessment model over-estimation belongs
partly in the methods or discussion, not results.** The results
section contains two paragraphs about why the *SISCAH-DIM* assessment
over-estimates biomass under predation OMs. This is an interpretation
of results, not a result itself. The interpretation is well-reasoned
and the discussion section (§Assessment model performance) covers
substantially the same material. Consider streamlining the results
section to report the observation (over-estimation occurs) and let the
discussion carry the mechanistic explanation.

**R-3. The allocation proportions reported in the text do not have a
clear source.** The text states "the 2024 allocation proportions based
on spawn indices are approximately 5\%, 38\%, and 57\% for Stat Areas
25, 24, and 23, respectively, while long-term average proportions are
approximately 11\%, 31\%, and 58\%." These values presumably come
from the `makeAllocTable()` function output. If so, a cross-reference
to the table (e.g., "Table \@ref(tab:tab-allocProportions)") would
tie the text to the data. Currently the table chunk label is
`tab-allocProportions` but there is no `\@ref()` in the text pointing
to it.

**R-4. The performance metrics table chunk may have a rendering
issue.** The code chunk `tab-perfMetrics` calls
`makePerfTables()` and then prints `perfTbls$biomass` and
`perfTbls$catch` as two separate kableExtra objects. When rendered
via bookdown/PDF, printing two kable objects sequentially from a
single chunk can produce placement issues or missing captions for
the second table. Consider splitting these into two separate chunks
with separate `results = "asis"` options and explicit `label=` values
in the kbl calls. Note that `tableMaker.R` already passes
`label = "perfBiomass"` and `label = "perfCatch"` to the kbl calls,
but the chunk-level label `tab-perfMetrics` only resolves one table
for cross-referencing. Having separate chunks (e.g.,
`tab-perfBiomass` and `tab-perfCatch`) will allow proper
`\@ref(tab:tab-perfBiomass)` cross-references in the text.

**R-5. The 1S-Pred section has no sub-headings.** The DDM section has
"Biological reference points" and "MP performance" sub-sections; the
1S-Pred section does not. For parallel structure and ease of
navigation, the 1S-Pred section should use the same sub-heading
structure (or the DDM sub-headings should be removed to match the
simpler structure used for the other two OMs).

### Discussion (`04-discussion.Rmd`)

**D-1. The conclusions section does not address the allocation
question explicitly.** The §Conclusions paragraph summarises the
primary finding (OM choice dominates, predation matters), but does
not provide a clear bottom line on the allocation question. Readers
coming specifically for guidance on how to allocate the TAC will not
find a direct answer. Suggest adding a sentence such as: "Among the
three allocation methods tested, differences in performance were
negligible under all operating models, suggesting that the choice of
allocation method is not a priority concern under current productivity
assumptions." This is implied by the discussion but not stated as a
conclusion.

**D-2. The limitation on flat projection uncertainty envelopes
deserves more specific guidance.** The limitations section notes that
"projection uncertainty envelopes appear flat in some simulations,
suggesting that the process error implementation may need
refinement." This is an important caveat that could affect confidence
in the reported performance statistics. The reader is left without
guidance on how flat is "flat" or which OM(s) exhibit this problem.
Is it a known numerical issue with the current simulation setup?
A brief indication of whether this affects the 3S-Pred OM, the
1S-Pred OM, or both -- and whether it tends to under- or
over-state uncertainty -- would help readers assess the severity
of the caveat.

**D-3. The predation modelling discussion would benefit from a
forward-looking sentence.** The discussion correctly highlights that
resolving uncertainty about predation mortality is key. However, it
does not say what kinds of data or analyses would help resolve this
uncertainty. A single sentence pointing toward empirical avenues
(e.g., continued predator diet sampling, bioenergetics model updates
as new abundance estimates become available) would make the
conclusion more actionable.

### Cross-Section Consistency

**X-1. The executive summary mentions age-dependent egg production;
the introduction and methods do not.** The executive summary
(`00-execSummary.Rmd`) mentions that Phase 2 models included
"age-dependent egg production." The introduction echoes this briefly
but the Phase 3 methods do not discuss whether age-dependent
fecundity is incorporated in the Phase 3 OMs. If it is, a brief
note in the operating model descriptions would complete the picture.
If it is not (i.e., it was deferred), that should be noted explicitly
since it was presented as a Phase 2 contribution.

**X-2. The introduction characterises SISCAH-3S-Pred as developed
in Phase 2, but Phase 3 also involves a 3S-Pred OM fit.** The index
loads `fit_3S_2024` as the conditioning fit, which implies the 3S-Pred
OM was re-fit or updated for 2024 data as part of Phase 3 work. The
introduction says the Phase 2 model is used "as an operating model"
in Phase 3, which is consistent if no structural changes were made.
But the data year in the filename (`_2024`) suggests the model was
re-conditioned. Clarify whether the Phase 3 OM is identical in
structure to the Phase 2 model but updated with 2024 data, or whether
structural changes were made. This matters for reproducibility and
for readers familiar with the Phase 2 report.

**X-3. The executive summary states the HCR "maintains biomass above
the limit reference point approximately 75\% of the time" (under DDM).
The results section echoes this but does not clarify whether this is
the LRP ($0.3 B_0 \approx 13.5$ kt) or $B_{ref}$ (17.35 kt).** These
are different thresholds, and the LRP is approximately 22\% below
$B_{ref}$. In the discussion, the statement "maintains biomass above
reference levels" does not disambiguate either. It appears from
context that the 75\% figure refers to the $B_{ref}$ threshold (not
the LRP), which is the more conservative metric. Consistently
specifying which threshold the 75\% applies to would remove ambiguity.

**X-4. The DDM results section mentions the MSY value is
"approximately 1.6 kt" -- but this conflicts slightly with the
executive summary's "approximately 1,600 t."** These are
numerically identical (1.6 kt = 1,600 t), but if DDM results
reference the MSY as a performance benchmark, it should be noted
that catches "near the MSY" of 1.6 kt are expected, not catches
near the 75th percentile of the projection. The results text says
"Commercial catches... are at levels similar to the historical
average" -- this could be made more precise by stating whether
projected catches approach or fall below the MSY benchmark.

---

## Minor Suggestions (Consider)

- **Bibliography completeness**: The `lfr2025_NTC` and `lfr2024_NTC`
  bib entries use "LFR" as the author name, which will render
  awkwardly in an author-year citation style. Consider using the
  full institutional author ("Landmark Fisheries Research") or
  the lead author names if authorship is known.

- **Terminology**: "spawn-on-kelp or bough" appears in the methods
  fleet description. This is an unusual formulation; the standard
  abbreviation "SOK" is used elsewhere. Consider "spawn-on-kelp
  (SOK) and bough" or simply "SOK/bough" for consistency with the
  figure labels, which use "SOK/FSC."

- **The `fig-yieldCurve` chunk uses a pre-rendered PDF figure**
  (`knitr::include_graphics("figure/DDMyieldCurve.pdf")`) rather
  than calling `plotYieldCurve(DDMfit)`. The `plotMethods.R` file
  defines `plotYieldCurve()` as a function. If `DDMfit` is loaded
  at build time, the figure should be generated live to ensure
  it updates automatically if model fits change. Using a static
  pre-rendered PDF is inconsistent with the dynamic approach used
  for all other figures and creates a reproducibility risk.

- **`plotHistSSB()` uses `pfmaCols <- c("tomato", "darkorange",
  "firebrick")`** for Stat Areas 25, 24, 23 respectively, and
  the legend lists them as "Stat Area 25", "Stat Area 24",
  "Stat Area 23." However, the aggregate 3S-Pred line is also
  plotted in "firebrick," the same colour as Stat Area 23. This
  could cause visual confusion. Consider using distinct colours
  for aggregate vs. per-area lines (e.g., a dark solid line for
  aggregate and lighter dashed variants for individual areas).

- **`plotResults.R` labels the SOK/FSC panel in t (tonnes) rather
  than kt.** The SSB and commercial catch panels use kt. While the
  unit difference is appropriate given the smaller SOK/FSC
  quantities, the axis label "SOK/FSC Catch (t)" should be
  mentioned in the figure caption so readers are not confused by
  the scale change between panels.

- **The `makePerfTables()` function formats probabilities to 3
  decimal places** (`digits = 3` in `kbl()`). Given that the
  reported range is 0 to 0.2\% (biomass) and 71\% to 100\%
  (catch), three decimal places for proportions (0.000 to 0.002
  and 0.71 to 1.000) is appropriate, but consider formatting as
  percentages for easier reading (0.0\% to 0.2\% and 71\% to
  100\%).

- **Line 94 in `03-results.Rmd`**: the phrase "lots of simulated
  outcomes" (flagged as Critical under issue R-1 above) appears
  here. Separately: "Assessment errors translate to some
  over-harvesting near the $B_{ref}$ level" -- the word "some"
  is vague. If quantified (e.g., "in approximately X\% of
  replicates"), this would be more defensible.

- **The `_bookdown.yml` does not specify a `language:` field.**
  For Canadian English spell-checking and hyphenation in LaTeX,
  adding `lang: en-CA` to `index.Rmd` YAML (or ensuring the
  LaTeX template includes `babel` with `[canadian]english`) would
  be consistent with the project convention.

- **Figure captions in `02-methods.Rmd` are defined inside their
  respective code chunks** (the `cap*` variable is assigned within
  the chunk, then the chunk's `fig.cap` argument references it).
  This pattern requires that the caption variable is assigned
  before any figure-generating call in the same chunk, which is
  satisfied here, but it is fragile. A cleaner pattern is to
  define all captions in a dedicated setup chunk in `index.Rmd`
  or at the top of each `.Rmd` file so they are available
  regardless of chunk evaluation order. This is a style issue,
  not a correctness issue.

---

## Positive Aspects

- The overall narrative is clear and well-scoped. The document
  does not overreach -- it describes exactly what was done and
  why, without inflating the results.

- The three-OM structure (reference case + two robustness tests
  of increasing complexity) is an appropriate and well-motivated
  MSE design. The rationale for each OM is explained clearly.

- The HCR formulation is described precisely with a proper
  mathematical equation, and the equation is consistent with the
  `plotHCRDiagram()` implementation (both use the same
  $\text{TAC} = B - B_{ref}$ ramp followed by a $U_{ref}$ cap
  formulation), which is good.

- The limitations section is frank and appropriately self-critical,
  particularly the acknowledgement of flat projection uncertainty
  envelopes and the absence of movement among Stat Areas.

- The R code infrastructure (separate `plotMethods.R`,
  `plotResults.R`, `tableMaker.R` scripts, the blob-patching
  `.patchBlob()` function in `index.Rmd`) is well-organised and
  will make it straightforward for collaborators to reproduce
  or update figures when model outputs change.

- The exec summary is concise and hits the key findings clearly.
  The bullet-point format with quantitative values (75\% closure
  rate, >50\% closure under 3S-Pred) gives the reader a useful
  quick read.

- The project correctly treats *SISCAH-3S-Pred* as a robustness
  OM rather than the reference case. The structural mismatch
  framing (single-stock MP, spatially structured OM) is clearly
  and consistently explained throughout.

---

*Reviewer: Dr. Sarah Chen*
*Review completed: 2026-03-26*
