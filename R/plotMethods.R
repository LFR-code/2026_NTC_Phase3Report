# plotMethods.R
#
# Methods-section figure wrappers for the NTC Phase 3 report.
# Calls utility functions in ms3Rplots.R where applicable.
# All figures use base plot (no ggplot2).
#
# Functions:
#   plotYieldCurve()  -- DDM equilibrium yield curve (slide 2)
#   plotHistSSB()     -- Historical SSB all three OMs (slide 3)
#   plotHCRDiagram()  -- Two-panel minimum escapement HCR (slide 4)


# plotYieldCurve()
# Wrapper around plotRefPtsExpC() from yieldCurveFun.R.
# Source yieldCurveFun.R before calling.
# Inputs:
#   empRC   -- saveEmpRefCurves object
#   histRpt -- DDMfit (history report object)
#   qProbs  -- quantile probabilities for ref pt bands
#   yrRange -- year indices for averaging (relative to end)
plotYieldCurve <- function( empRC   = saveEmpRefCurves,
                            histRpt = DDMfit,
                            Bref    = 17.35,
                            qProbs  = c(0.25, 0.5, 0.75),
                            yrRange = -c(150:50) )
{
  suppressWarnings(suppressMessages(
    plotRefPtsExpC(
      empRefCurves = empRC,
      histRpt      = histRpt,
      yrRange      = yrRange,
      qProbs       = qProbs,
      area         = "WCVI",
      Bref         = Bref)
  ))

} # END plotYieldCurve()


# plotHistSSB()
# Historical spawning biomass time series from all three OMs.
# Inputs:
#   DDMfit  -- DDMfit (readRDS from fit_WCVI2023_OM_mle.rds)
#   fit1S   -- fit1S  (loaded as 'reports' from fit_1S_2024.RData)
#   fit3S   -- fit3S  (loaded as 'reports' from fit_3S_2024.RData)
plotHistSSB <- function( DDMfit = DDMfit,
                         fit1S  = fit1S,
                         fit3S  = fit3S,
                         maxYear = 2023 )
{
  fYear <- DDMfit$fYear

  # DDM: single stock SB_pt (1 x nT)
  SB_DDM  <- DDMfit$repOpt$SB_pt[1, ]
  yrs_DDM <- seq(from = fYear, by = 1,
                 length.out = length(SB_DDM))
  keep_DDM <- yrs_DDM <= maxYear
  SB_DDM   <- SB_DDM[keep_DDM]
  yrs_DDM  <- yrs_DDM[keep_DDM]

  # 1S: single stock SB_pt (1 x nT)
  SB_1S  <- fit1S$repOpt$SB_pt[1, ]
  yrs_1S <- seq(from = fit1S$fYear, by = 1,
                length.out = length(SB_1S))
  keep_1S <- yrs_1S <= maxYear
  SB_1S   <- SB_1S[keep_1S]
  yrs_1S  <- yrs_1S[keep_1S]

  # 3S: three stocks SB_pt (nP x nT)
  SB_3S_mat <- fit3S$repOpt$SB_pt
  yrs_3S    <- seq(from = fit3S$fYear, by = 1,
                   length.out = ncol(SB_3S_mat))
  keep_3S   <- yrs_3S <= maxYear
  SB_3S_mat <- SB_3S_mat[, keep_3S, drop = FALSE]
  yrs_3S    <- yrs_3S[keep_3S]
  SB_3S_agg <- colSums(SB_3S_mat)

  # B0 reference values
  B0_DDM  <- DDMfit$repOpt$B0_p
  B0_1S   <- fit1S$repOpt$refPts$refCurves$SBeq_pf[1, 1]
  B0_3S   <- sum(fit3S$repOpt$refPts$refCurves$SBeq_pf[, 1])

  allYrs <- c(yrs_DDM, yrs_1S, yrs_3S)
  allSB  <- c(SB_DDM, SB_1S, SB_3S_agg, B0_DDM)

  par(mar = c(4, 5, 1, 1), oma = c(0, 0, 0, 0))
  plot(
    x    = range(allYrs),
    y    = c(0, max(allSB, na.rm = TRUE) * 1.1),
    type = "n",
    xlab = "Year",
    ylab = "Spawning Biomass (kt)",
    las  = 1
  )

  # B0 reference lines
  abline(h = B0_DDM, lty = 3, col = "grey50",    lwd = 1.5)
  abline(h = B0_1S,  lty = 3, col = "steelblue", lwd = 1.5)
  abline(h = B0_3S,  lty = 3, col = "firebrick",  lwd = 1.5)

  # Per-Stat Area lines for 3S-predM (dashed)
  pfmaCols <- c("tomato", "darkorange", "firebrick")
  for (p in seq_len(nrow(SB_3S_mat))) {
    lines(
      x   = yrs_3S,
      y   = SB_3S_mat[p, ],
      col = pfmaCols[p],
      lty = 2,
      lwd = 1
    )
  }

  # Aggregate lines
  lines(x = yrs_DDM, y = SB_DDM,
        col = "black", lwd = 2.5)
  lines(x = yrs_1S, y = SB_1S,
        col = "steelblue", lwd = 2.5)
  lines(x = yrs_3S, y = SB_3S_agg,
        col = "firebrick", lwd = 2.5)

  # B0 labels
  labX <- max(allYrs) - 1
  text(x = labX, y = B0_DDM * 1.03,
       labels = expression(B["0 | DDM"]),
       adj = c(1, 0), col = "grey40", cex = 0.85)
  text(x = labX, y = B0_1S * 1.03,
       labels = expression(B["0 | 1S"]),
       adj = c(1, 0), col = "steelblue", cex = 0.85)
  text(x = labX, y = B0_3S * 1.03,
       labels = expression(B["0 | 3S"]),
       adj = c(1, 0), col = "firebrick", cex = 0.85)

  legend(
    x      = "topright",
    bty    = "n",
    legend = c(
      "SISCAH-DDM",
      "SISCAH-1S-Pred",
      "SISCAH-3S-Pred (agg.)",
      paste0("  Stat Area ", c(25, 24, 23))
    ),
    col    = c("black", "steelblue", "firebrick",
               pfmaCols),
    lty    = c(1, 1, 1, 2, 2, 2),
    lwd    = c(2.5, 2.5, 2.5, 1, 1, 1),
    cex    = 0.8
  )
} # END plotHistSSB()


# plotHCRDiagram()
# Two-panel minimum escapement HCR diagram.
# Panel 1: harvest rate U vs assessed biomass.
# Panel 2: TAC vs assessed biomass.
# Parameters match the project HCR:
#   Bref = 17.35 kt, Uref = 0.0847
plotHCRDiagram <- function( Bref = 17.35, Uref = 0.0847 )
{
  Bmax <- 50
  Bseq <- seq(0, Bmax, length.out = 500)

  # Minimum escapement HCR: TAC = B - Bref (up to Uref cap)
  # Uref cap kicks in at Bcap = Bref / (1 - Uref)
  Bcap     <- Bref / (1 - Uref)
  Ufun     <- function(B) ifelse(B <= Bref, 0, pmin(1 - Bref / B, Uref))
  Ucurve   <- Ufun(Bseq)
  TACcurve <- Ucurve * Bseq
  Ymax_U   <- Uref * 1.2
  Ymax_T   <- max(TACcurve) * 1.2

  par(mfrow = c(2, 1), mar = c(2, 5, 1, 1), oma = c(3, 0, 0, 0))

  # ---- Panel 1: U vs B ----
  plot(
    x    = range(Bseq),
    y    = c(0, Ymax_U),
    type = "n",
    xlab = "",
    ylab = "Harvest Rate (U)",
    las  = 1,
    xaxt = "n"
  )
  axis(side = 1, labels = FALSE)
  # Red (closed), yellow (ramp), blue (full harvest) regions
  polygon(x = c(0, Bref, Bref, 0),
          y = c(0, 0, Ymax_U, Ymax_U),
          col = adjustcolor("red", alpha.f = 0.10), border = NA)
  polygon(x = c(Bref, Bcap, Bcap, Bref),
          y = c(0, 0, Ymax_U, Ymax_U),
          col = adjustcolor("goldenrod", alpha.f = 0.15), border = NA)
  polygon(x = c(Bcap, Bmax, Bmax, Bcap),
          y = c(0, 0, Ymax_U, Ymax_U),
          col = adjustcolor("steelblue", alpha.f = 0.12), border = NA)
  abline(v = Bref, lty = 2, col = "red",       lwd = 1.5)
  abline(v = Bcap, lty = 2, col = "steelblue", lwd = 1.2)
  abline(h = Uref, lty = 2, col = "grey50",    lwd = 1.2)
  lines(x = Bseq, y = Ucurve, lwd = 3, col = "black")
  panLab(x = 0.06, y = 0.88, txt = "Closed",   col = "red",       cex = 0.8)
  panLab(x = 0.38, y = 0.88, txt = expression(U == 1 - B[ref]/B),
         col = "goldenrod4", cex = 0.75)
  panLab(x = 0.72, y = 0.88, txt = expression(U == U[ref]),
         col = "steelblue",  cex = 0.8)
  text(x = Bref, y = Ymax_U * 0.97,
       labels = expression(B[ref]), adj = c(1.1, 1),
       col = "red", cex = 0.85)
  text(x = Bmax * 0.98, y = Uref * 1.04,
       labels = expression(U[ref]), adj = c(1, 0),
       col = "grey40", cex = 0.85)

  # ---- Panel 2: TAC vs B ----
  plot(
    x    = range(Bseq),
    y    = c(0, Ymax_T),
    type = "n",
    xlab = "",
    ylab = "TAC (kt)",
    las  = 1
  )
  polygon(x = c(0, Bref, Bref, 0),
          y = c(0, 0, Ymax_T, Ymax_T),
          col = adjustcolor("red", alpha.f = 0.10), border = NA)
  polygon(x = c(Bref, Bcap, Bcap, Bref),
          y = c(0, 0, Ymax_T, Ymax_T),
          col = adjustcolor("goldenrod", alpha.f = 0.15), border = NA)
  polygon(x = c(Bcap, Bmax, Bmax, Bcap),
          y = c(0, 0, Ymax_T, Ymax_T),
          col = adjustcolor("steelblue", alpha.f = 0.12), border = NA)
  abline(v = Bref, lty = 2, col = "red",       lwd = 1.5)
  abline(v = Bcap, lty = 2, col = "steelblue", lwd = 1.2)
  lines(x = Bseq, y = TACcurve, lwd = 3, col = "black")
  panLab(x = 0.38, y = 0.88,
         txt = expression(TAC == B - B[ref]),
         col = "goldenrod4", cex = 0.75)
  panLab(x = 0.72, y = 0.88,
         txt = expression(TAC == U[ref] %.% B),
         col = "steelblue",  cex = 0.8)
  text(x = Bref, y = Ymax_T * 0.97,
       labels = expression(B[ref]), adj = c(1.1, 1),
       col = "red", cex = 0.85)

  mtext(side = 1, outer = TRUE, text = "Assessed Biomass (kt)", line = 1.5)
} # END plotHCRDiagram()
