# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code)
when working with code in this repository.

## Project overview

This is a bookdown working paper for Phase 3 of the
?aayaaqa (Herring) Spawn Dynamics project (NTC-BCSRIF).
It evaluates a WCVI herring re-opening management procedure
via closed-loop simulation (MSE) under three operating models
(*SISCAH-DDM*, *SISCAH-1S-Pred*, *SISCAH-3S-Pred*).
The output is a PDF report: `_book/NTC-BCSRIF-Phase3-WorkingPaper.pdf`.

## Build commands

```bash
make pdf      # render PDF via bookdown::render_book
make check    # lint: TODOs, line lengths >80, placeholders
make clean    # remove _book/, _bookdown_files/, LaTeX artefacts
```

Alternative render (PDF + Word): `Rscript renderDoc.R`

Requires: R, bookdown, knitr, kableExtra, tidyverse, pandoc.

## Document structure

- `index.Rmd` -- YAML front matter, knitr setup, package loading
- `00-execSummary.Rmd` -- Executive summary (unnumbered)
- `01-intro.Rmd` -- Introduction, management context, objectives
- `02-methods.Rmd` -- OMs, MP, HCR, allocation methods, sim design
- `03-results.Rmd` -- Results by OM, spatial allocation performance
- `04-discussion.Rmd` -- Implications, limitations, conclusions
- `05-references.Rmd` -- Auto-generated references (do not edit)
- `_bookdown.yml` -- Chapter ordering and output filename

## Key conventions

- **LaTeX template**: `templates/basic.latex` (primary);
  `templates/lfrDraft.latex` is also available
- **Citation style**: `csl/fisheries-research.csl`
- **Bibliography**: `bib/library.bib` (BibDesk-managed BibTeX)
- **Cross-references**: use bookdown syntax
  (`\@ref(sec:label)`, `\@ref(fig:label)`, `\@ref(tab:label)`)
- **Percentages**: use `\%` in Rmd files, never bare `%`
- **Math/symbols**: use LaTeX only ($B_{MSY}$, $\times$),
  never Unicode special characters
- **Line length**: 80 characters max in all Rmd/md files
- **Spelling**: Canadian English throughout
- **Terminology**: use "methods" not "methodology"
- **Figure/table captions**: keep caption character vectors on
  a single line
- **R code**: use named arguments in all function calls;
  use `drop = FALSE` when subsetting arrays with `[`
- **Plotting**: base plot only (no ggplot2)
- **Data scripts**: `index.Rmd` has TODO stubs for sourcing
  `ms3Rtools.r`, `ms3Rplots.r`, `ms3Rstats.r` and loading
  simulation `.Rdata` files (not yet connected)

## Domain terminology

- **OM**: operating model (simulated truth)
- **MP**: management procedure (assessment + HCR + allocation)
- **HCR**: harvest control rule (minimum escapement type)
- **PFMA**: Pacific Fishery Management Area (23, 24, 25)
- **TAC**: total allowable catch
- **SISCAH**: the herring model family used in this project
- **DDM/DIM**: density-dependent / density-independent mortality
- **MSE**: management strategy evaluation (closed-loop simulation)
- **$B_{ref}$**: minimum escapement reference = $B_{MSY|DDM}$
