# NTC-BCSRIF Phase 3 Working Paper: Pre-Circulation Review

**Document**: ?aayaaqa (Herring) Spawn Dynamics -- Phase 3 Working
Paper (NTC-BCSRIF, Landmark Fisheries Research for Uu-a-thluk
Fisheries)\
**Reviewer**: Dr. Sarah Chen\
**Date**: 2026-03-26\
**Files reviewed**: `index.Rmd`, `00-execSummary.Rmd`,
`01-intro.Rmd`, `02-methods.Rmd`, `03-results.Rmd`,
`04-discussion.Rmd`, `06-appendix.Rmd`\
**Review scope**: Final pre-circulation check -- broken prose,
terminology consistency, "we" framing, scientific accuracy,
figure/table cross-references, and remaining TODO comments.

---

## Overall Assessment

The document is in good shape for a pre-circulation draft. The
narrative is coherent, the three-OM structure is logically presented,
and the writing is technically appropriate for the intended audience.
Several editing artefacts remain in the source, and there are a small
number of substantive issues -- most importantly, three active TODO
comments that have not been resolved, an uncorrected word-duplication
error, and a systematic gap in figure cross-referencing -- that should
be addressed before the document goes to Uu-a-thluk reviewers.

---

## Issues Requiring Action Before Circulation

### 1. Three active TODO comments remain in source

All three are in `01-intro.Rmd` and will render as HTML comments in
the PDF (invisible to readers), but they represent unresolved
editorial instructions that need to be actioned, not just left in
place.

**TODO at line 3--4** (`01-intro.Rmd`):

```
<!-- TODO: Review framing of LFR role (Sean: LFR is a
contractor, not a collaborator) -->
```

The paragraph immediately following this TODO reads: "The ?aayaaqa
(Herring) Spawn Dynamics project is a multi-year project led by
Uu-a-thluk Fisheries, with **technical support** from Landmark
Fisheries Research." The framing "technical support" is already
reasonably consistent with a contractor role rather than a
collaborator role. However, the next substantive concern is that
line 92 reads:

> "*SISCAH-1S-Pred* and *SISCAH-3S-Pred* OMs **developed by
> Uu-a-thluk**"

These models were developed under contract by LFR (Landmark Fisheries
Research) and funded through the NTC-BCSRIF project -- describing
them as "developed by Uu-a-thluk" conflates the funding body /
project owner with the technical contractor. Suggested revision:

> "the *SISCAH-1S-Pred* and *SISCAH-3S-Pred* OMs developed under
> the NTC-BCSRIF project"

or, if attribution to LFR is appropriate here:

> "the *SISCAH-1S-Pred* and *SISCAH-3S-Pred* OMs developed by LFR
> for this project"

Once the attribution question is resolved, remove the TODO comment.

**TODO at line 36--37** (`01-intro.Rmd`):

```
<!-- TODO: Review Phase 3 goal statement (Sean: should be
front and centre, not buried at end of paragraph) -->
```

The Phase 3 goal statement currently appears at the end of a
paragraph that begins with Phase 1 and Phase 2 context:

> "One goal of Phase 3 is to evaluate potential re-opening procedures
> for the WCVI herring fishery."

This is indeed buried. The sentence reads "one goal" which implies
there are other goals not listed. Consider restructuring this
paragraph or opening the Introduction with the Phase 3 framing more
prominently, then providing the Phase 1/2 context as supporting
background. The "one goal" phrasing also needs attention -- if MSE
and the distribution methods testing are both goals, say so; if
testing the RP is the primary goal, say "the primary goal" or "the
goal."

**TODO at line 69--70** (`01-intro.Rmd`):

```
<!-- TODO: Review MSE definition (Sean: tailor to NTC
herring context; clean up incomplete edits) -->
```

The MSE definition paragraph (lines 71--85) has a clear editing
artefact: line 75 reads:

> "Objectives and associated **their** associated metrics are agreed
> upon..."

This is a word-duplication error -- "associated" appears twice and
"their" is a fragment from a prior revision. The sentence should
read either:

> "Objectives and their associated metrics are agreed upon..."

or

> "Objectives and associated performance metrics are agreed upon..."

Beyond this correction, the paragraph is fairly generic MSE boilerplate
and does not specifically connect to the WCVI herring context, the
Nuu-chah-nulth objectives, or the specific structure of this project's
MSE (three OMs, one RP, three distribution methods). Tailoring this
paragraph -- even briefly -- would address Sean's note. At minimum,
the final sentence of the paragraph could be modified to reference
"the WCVI herring RP evaluated here" rather than describing a generic
MSE framework.

---

### 2. Word-duplication error in `03-results.Rmd` line 96

The sentence reads:

> "which allows it to chase **spawn large spawn** indices towards the
> end of the historical period."

This is a clear editing artefact -- "spawn" appears twice with
"large" inserted between repetitions. The intended reading is most
likely one of:

> "which allows it to chase large spawn index values towards the end
> of the historical period."

or

> "which allows it to track large spawn indices towards the end of
> the historical period."

The word "chase" is also informal for a technical document (it
appears again in the same paragraph at line 91: "more flexibility to
chase high data points"). Both uses of "chase" should be replaced
with more precise language such as "fit to" or "track."

---

### 3. Figure cross-references are missing for several figures

Three figures are defined as code chunks in `03-results.Rmd` but are
never referenced in the prose via `\@ref()`:

- `fig-DDMretro` (retrospective SSB fits, DDM OM): defined at line 49,
  no `\@ref(fig:fig-DDMretro)` in text. The DDM results paragraph
  (lines 22--36) references the range
  `\@ref(fig:fig-DDMprojections)--\@ref(fig:fig-DDMhcrPhase)`, which
  skips the retrospective figure entirely.

- `fig-3SprojByArea` (per-Stat Area projections, 3S-Pred OM): defined
  at line 188, no `\@ref()` in text. This is the figure that shows
  individual Stat Area dynamics under the 3S-Pred OM, which is
  central to the spatial distribution results. Its absence from the
  prose text means readers may not notice it.

- `tab-perfBiomass` and `tab-perfCatch`: defined as separate chunks
  at lines 263 and 268 respectively, but neither chunk has a `\@ref()`
  call in the text. The prose in the "Performance by Stat Area"
  section describes the results in words but does not direct readers
  to the tables. Cross-references should be added, e.g., at the end
  of the paragraph beginning "The probability of low biomass..."
  (line 242), and again after "The probability of low catch..."
  (line 249).

The existing range-reference pattern
`\@ref(fig:fig-DDMprojections)--\@ref(fig:fig-DDMhcrPhase)` used in
the DDM and 3S-Pred sections will silently skip over the retrospective
figure. Given that retrospective assessment performance is discussed
substantively in the 1S-Pred section (and the issue of assessment
over-estimation is a key finding), each retrospective figure should
be individually referenced.

---

### 4. "reopening" (un-hyphenated) in `01-intro.Rmd` line 43

The document uses "re-opening" (hyphenated) consistently throughout,
which is correct for Canadian English. However, line 43 reads:

> "While policy requires that **reopening** procedures will perform
> well at the stock assessment region (SAR) level..."

This single un-hyphenated instance ("reopening") should be corrected
to "re-opening" for consistency.

---

### 5. Terminology slip: "management procedures" in `04-discussion.Rmd`
   line 184 and `02-methods.Rmd` line 114

The document has renamed "management procedure" to "re-opening
procedure (RP)" as per project convention. Two instances of the old
term remain:

**`04-discussion.Rmd` line 184**:
> "demonstrates the value of testing **management procedures** against
> robustness scenarios that include alternative ecological hypotheses."

Suggested revision:
> "demonstrates the value of testing **re-opening procedures** against
> robustness scenarios that include alternative ecological hypotheses."

**`02-methods.Rmd` line 114**:
> "The *SISCAH-DIM* model is the adopted assessment model for BC
> Herring SARs that have active fisheries with TACs set via
> **management procedures**."

This sentence is describing the DFO context (where "management
procedure" is the standard DFO terminology) rather than referring to
this project's RP, so retaining "management procedures" here is
arguably correct. However, it may cause confusion immediately after
the preceding sentence introduced "re-opening procedure (RP)." Consider
adding a parenthetical to clarify the distinction, e.g.:

> "The *SISCAH-DIM* model is the adopted assessment model for BC
> Herring SARs that have active fisheries with TACs set via DFO
> management procedures."

Similarly, `01-intro.Rmd` line 60 and `04-discussion.Rmd` line 172--174
also use "management procedures" in a DFO policy context, which is
appropriate, but check that the distinction is clear to readers
encountering these sentences after the RP terminology has been
introduced.

---

### 6. Incomplete sentence structure: "One goal of Phase 3 is..."
   (`01-intro.Rmd` line 38)

As noted under TODO issue 2 above, the phrase "One goal of Phase 3"
implies multiple goals without listing them. The Objectives section
(lines 101--129) covers both the RP evaluation and the distribution
method testing, but neither is explicitly framed as a separate goal
in the introduction. Either:

(a) revise to "The primary goal of Phase 3 is to evaluate potential
    re-opening procedures..."; or

(b) list both goals: "Phase 3 has two goals: (i) to evaluate a
    proposed re-opening procedure for the WCVI herring fishery under
    three operating models, and (ii) to test methods for distributing
    the resulting TAC among Stat Areas."

---

### 7. Misspelling: "intrepeted" in `03-results.Rmd` line 84

> "which is **intrepeted** as a large outlier by the *SISCAH-1S-Pred*
> OM."

Should be "interpreted."

---

### 8. "quite high" is informal (`03-results.Rmd` line 82)

> "Initial TACs are set **quite high** given the large recent spawn
> index from Stat Area 24"

"Quite high" is informal. Suggested revision:
> "Initial TACs are elevated given the large recent spawn index from
> Stat Area 24, which..."

---

## Substantive Issues

### S-1. The five aggregate performance metrics defined in Methods
   (`pLRP`, `MAC10`, `MAS10`, `Openings`, `p650`) do not appear in
   the Results

The Performance Metrics section of `02-methods.Rmd` (lines 219--239)
defines five quantitative metrics for aggregate performance evaluation
at the WCVI SAR scale. None of these metric names appear anywhere in
`03-results.Rmd`. The results text reports the pLRP value (96\%) and
notes that biomass is above the LRP, but uses no other metric names
and presents no summary table for aggregate performance under the DDM
or 1S-Pred OMs. The table `tab-omStats` (line 277) is referenced in
the text and presumably contains aggregate statistics, but without
seeing the rendered output it is not clear whether it covers the five
defined metrics.

If the five metrics are computed and presented in `tab-omStats`, the
Results text should name them when interpreting results (e.g., "pLRP
= 0.96, exceeding the 75\% conservation objective; MAC10 = X kt").
If the metrics are not all computed for all OMs, the methods section
should be revised to clarify which metrics apply to which OMs rather
than implying uniform coverage.

### S-2. The Stat Area-level performance metrics are discussed without
   a clear linkage to distribution methods under the DDM OM

The methods section defines the Stat Area-level metrics (probability
of low biomass, probability of low catch) and states they are only
relevant under the 3S-Pred OM. The results section correctly reports
these metrics under the 3S-Pred OM. However, the framing in the
introduction and executive summary suggests that the distribution
method comparison is a co-equal objective with the RP evaluation.
Readers who expect a full distribution method comparison under all
OMs may be surprised to find that the DDM and 1S-Pred results
sections contain no distribution method comparison at all. A single
sentence in the Results introduction clarifying that distribution
method comparisons are only meaningful under the 3S-Pred OM (because
the other OMs are single-stock) would pre-empt this confusion.

### S-3. The observation error model description is incomplete

The simulation design section (lines 195--202 in `02-methods.Rmd`)
states that observation errors are "parameterised from historical
residuals" but does not specify whether this means the standard
deviation of log-residuals, the full residual time series, or a CV.
This is not a fatal gap but is worth one sentence given that the
assessment over-estimation bias is a key finding -- readers will
want to know whether the observation errors are realistic relative
to the historical period.

---

## Minor Style and Consistency Issues

- **`00-execSummary.Rmd` line 67**: "Methods for distributing catch
  among Stat Areas were also tested." -- The word "Methods" here
  could be misread as referring to the Methods section of the
  document. Consider "Distribution methods for catch among Stat
  Areas were also tested" or "Three methods for distributing catch
  among Stat Areas were also evaluated."

- **`00-execSummary.Rmd` line 61**: "the standard Fisheries and
  Oceans, Canada (DFO) model" -- first use of the DFO acronym is
  in the executive summary, which is correct. Confirm the acronym
  is not re-defined in the Introduction. (It is re-introduced at
  `01-intro.Rmd` line 42 as "Fisheries and Oceans, Canada (DFO)"
  -- this second definition is redundant and the acronym should just
  be used directly.)

- **`02-methods.Rmd` lines 232--234**: indentation of the MAC10 and
  MAS10 sub-items is inconsistent. MAC10 text starts at column 1,
  MAS10 text starts at column 1, but the surrounding list items
  (pLRP, Openings, p650) have their continuation lines indented
  three spaces. Visual inconsistency only; will not affect rendering.

- **`02-methods.Rmd` line 187**: "spawn-on-kelp or bough" -- the
  "or" is ambiguous (are these the same gear with two names, or two
  separate gears?). The abbreviation "SOK/bough" used elsewhere in
  the figure captions, or simply "spawn-on-kelp/bough (SOKB)" at
  first mention, would be clearer. Note that the performance metric
  at line 232 uses "spawn-on-bough" (without "kelp") and the figure
  captions use "SOK/FSC." These should be reconciled: if FSC (food,
  social, and ceremonial) and SOKB are combined in the same panel,
  the caption should say so; if they are separate, "spawn-on-bough"
  and "SOK" should not be used interchangeably.

- **`01-intro.Rmd` lines 75--76**: the sentence continues after
  the word duplication artefact ("Objectives and associated their
  associated metrics") with another "associated" on the next line:
  "associated metrics are agreed upon." After fixing the duplication,
  confirm the surrounding prose reads naturally.

- **`04-discussion.Rmd` lines 147--155**: the paragraph on Hake
  predation reads:
  > "recent estimates of Pacific Hake biomass have reached very low
  > levels [@jtc2025hake] and it is unclear when biomass will recover
  > and to what level. Therefore, these simulations **likely
  > over-estimate** the impact of Hake predation on WCVI Herring."

  The logic here is inverted. If Hake biomass is very low, these
  simulations would likely *under-estimate* Hake predation if the
  simulations project Hake recovering to carrying capacity (which is
  what the preceding sentence says the surplus production models do).
  The argument appears to be: (i) surplus production models project
  Hake recovery to carrying capacity, (ii) real Hake biomass is
  currently very low and may not recover, therefore (iii) the
  simulations project more Hake predation than will actually occur
  -- which means the simulations over-estimate Hake impact, not
  under-estimate it. If that is the intended meaning, the connection
  between the low current Hake biomass observation and the conclusion
  needs to be spelled out more clearly. Suggest:

  > "Recent estimates of Pacific Hake biomass have reached very low
  > levels [@jtc2025hake] and recovery to carrying capacity is
  > uncertain. Because the surplus production model projects recovery,
  > these simulations may over-estimate the long-term impact of Hake
  > predation."

- **`06-appendix.Rmd` line 90**: The gear index description for $g$
  in the notation table contains a trailing comma after the
  description for gear (3) gillnet:
  > "(3) gillnet fisheries,"
  This is followed by "(4) Surface and (5) dive..." and continues
  as a single long multi-row description. The comma after
  "fisheries" followed by a new cell entry creates a minor
  grammatical awkwardness. Consider using a semicolon-separated
  list within the single description cell or restructuring as a
  sub-list.

- **`06-appendix.Rmd` line 427**: "**percent difference** in
  relative recruit per spawner" in `plotRPS2Cap` should be
  "\% difference" (using `\%` per project convention). Note that
  this chunk has `eval = FALSE` so it will not render, but if it
  is activated in future the caption will need the correction.

---

## Cross-Section Consistency

- **Executive summary vs. Introduction on age-dependent
  fecundity**: The exec summary does not mention age-dependent egg
  production; the Introduction mentions it in the Phase 1 and Phase
  2 context; the Discussion (§Fecundity-at-age limitation) notes
  that fecundity was applied to the 3S-Pred OM but not the others.
  The Methods section is silent on this. The limitation section in
  the Discussion is clear and honest about the inconsistency, but
  readers who only read the Methods will not know the OMs differ in
  their fecundity treatment. A single sentence in the operating
  model descriptions -- specifically the SISCAH-3S-Pred subsection
  -- noting the fecundity adjustment would make the methods complete.

- **Exec summary and Results on the "96\%" figure**: The exec
  summary states "the population stays above conservation thresholds
  96\% of the time" (line 25). The results section (line 30--31)
  says "biomass is maintained above the LRP ($0.3 B_0$)
  approximately 96\% of the time, exceeding the 75\% conservation
  objective." These are consistent and both specify the LRP as
  the threshold, which is correct. No issue here -- flagging only
  to confirm the exec summary should also specify the LRP as the
  threshold (currently it says "conservation thresholds" generally).

- **Methods performance metrics vs. Discussion conclusions**: The
  methods define `pLRP` with the criterion "pLRP $\geq$ 75\% over
  three herring generations." The conclusions section does not
  restate this criterion when reporting the 96\% figure. Adding
  "...exceeding the 75\% conservation threshold" to the conclusions
  paragraph (or confirming it is already present) would strengthen
  the link between the defined criterion and the reported result.

---

## Positive Aspects

- The numerical values for $B_{ref}$ (17.35 kt) and $U_{ref}$
  (8.47\%) are now fully consistent across all sections, figure
  captions, and code calls. This is a notable improvement from
  earlier drafts.

- The appendix model description is well-structured and the
  notation table is thorough. The two-change description (Humpback
  seasonal groups, coastwide predator distribution) clearly signals
  what distinguishes SISCAH-3S-Pred from SISCAH-1S-Pred for readers
  familiar with the Phase 2 report.

- The "management procedures" vs "re-opening procedure" terminology
  is well-maintained throughout except for the two instances noted
  above, both of which are defensible in context.

- The use of `\approx` to convey $B_{ref} \approx B_{MSY|DDM}$
  (rather than claiming exact equality) is scientifically precise
  and consistent with the context note from Sam -- this distinction
  is correctly handled throughout the document.

- The discussion of predation model limitations (Hake biomass
  uncertainty, Humpback stabilisation) is appropriately cautious
  and adds value for readers interested in the robustness of the
  3S-Pred projections.

- The HCR equation is correctly formulated and consistent between
  the mathematical statement and the prose description that follows
  it (the "lesser of surplus and $U_{ref} \cdot \hat{B}_t$"
  explanation matches the case-statement equation).

---

## Summary Priority List

Items to action before circulation, in priority order:

1. Fix the word-duplication error at `03-results.Rmd` line 96
   ("chase spawn large spawn indices").
2. Fix the misspelling at `03-results.Rmd` line 84 ("intrepeted").
3. Fix "reopening" to "re-opening" at `01-intro.Rmd` line 43.
4. Fix "Objectives and associated their associated metrics" at
   `01-intro.Rmd` line 75.
5. Add cross-references for `fig-DDMretro`, `fig-3SprojByArea`,
   `tab-perfBiomass`, and `tab-perfCatch` in `03-results.Rmd`.
6. Resolve and remove the three TODO comments in `01-intro.Rmd`
   (lines 3, 36, 69), including addressing "developed by Uu-a-thluk"
   attribution (line 92) and the "One goal" framing (line 38).
7. Replace "management procedures" with "re-opening procedures" in
   `04-discussion.Rmd` line 184 (and clarify the DFO-context uses
   as noted above).
8. Check the Hake over-estimate logic at `04-discussion.Rmd`
   lines 147--155 (possible logical inversion).
9. Verify that `tab-omStats` covers the five aggregate metrics
   defined in Methods, or revise the metrics description.

---

*Reviewer: Dr. Sarah Chen*\
*Review completed: 2026-03-26*
