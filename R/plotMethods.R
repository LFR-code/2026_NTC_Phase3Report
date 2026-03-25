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
# Equilibrium yield curve from the DDM operating model fit,
# with labelled reference points.
# Style mirrors yieldCurveFun.R: smooth spline, rotated srt=90 labels at y=0.
# Inputs:
#   fit  -- DDMfit object (readRDS from fit_WCVI2023_OM_mle.rds)
plotYieldCurve <- function( fit = DDMfit )
{
  rp   <- fit$repOpt$refPts
  Beq  <- rp$refCurves$SBeq_pf[1, ]
  Yeq  <- rp$refCurves$Yeq_pf[1, ]
  Bmsy <- rp$FmsyRefPts$SBeqFmsy_p
  MSY  <- rp$FmsyRefPts$YeqFmsy_p
  B0   <- fit$repOpt$B0_p
  LRP  <- 0.3 * B0

  clrs <- RColorBrewer::brewer.pal(n = 8, name = "Dark2")

  # Append B0 point (F=0 gives B=B0, Y=0) so curve reaches unfished biomass
  Beq_full <- c(Beq, B0)
  Yeq_full <- c(Yeq, 0)

  fit_spl <- smooth.spline(x = Beq_full, y = Yeq_full, nknots = 10)

  adjustpos <- 1.7   # x-offset for rotated labels (mirrors yieldCurveFun.R)
  msy_x     <- 5     # x position for MSY label

  par(mar  = c(.25, .25, .25, .25),
      oma  = c(4, 4, 1, 1),
      xaxs = "r", yaxs = "r", las = 1)

  plot(
    x    = c(0, B0 * 1.05),
    y    = c(0, max(Yeq_full, na.rm = TRUE) * 1.2),
    axes = FALSE,
    type = "n"
  )
  axis(side = 2, las = 1, cex.axis = 1)
  axis(side = 1, cex.axis = 1)
  grid()
  box()
  lines(fit_spl, col = "black", lwd = 3)

  # B0
  points(x = B0, y = predict(fit_spl, B0)$y,
         pch = 16, col = clrs[1], cex = 1.5)
  text(x = B0 + adjustpos, y = 0,
       labels = expression(B[0]), pos = 4, srt = 90)

  # LRP (0.3 B0)
  segments(x0 = LRP, y0 = -1, x1 = LRP,
           y1 = predict(fit_spl, LRP)$y, lty = 2)
  points(x = LRP, y = predict(fit_spl, LRP)$y,
         pch = 16, col = clrs[2], cex = 1.5)
  text(x = LRP - adjustpos, y = 0,
       labels = expression("LRP, 0.3" * B[0]), pos = 4, srt = 90)

  # MSY horizontal line and label
  segments(x0 = -5, y0 = predict(fit_spl, Bmsy)$y,
           x1 = Bmsy, y1 = predict(fit_spl, Bmsy)$y, lty = 2)
  text(x = msy_x, y = predict(fit_spl, Bmsy)$y,
       labels = "MSY", pos = 3)

  # Bmsy
  segments(x0 = Bmsy, y0 = -1, x1 = Bmsy,
           y1 = predict(fit_spl, Bmsy)$y, lty = 2)
  points(x = Bmsy, y = predict(fit_spl, Bmsy)$y,
         pch = 16, col = clrs[3], cex = 1.5)
  text(x = Bmsy - adjustpos, y = 0,
       labels = expression(B[MSY]), pos = 4, srt = 90)

  # 0.4 Bmsy
  segments(x0 = 0.4 * Bmsy, y0 = -1, x1 = 0.4 * Bmsy,
           y1 = predict(fit_spl, 0.4 * Bmsy)$y, lty = 2)
  points(x = 0.4 * Bmsy, y = predict(fit_spl, 0.4 * Bmsy)$y,
         pch = 16, col = clrs[4], cex = 1.5)
  text(x = 0.4 * Bmsy - adjustpos, y = 0,
       labels = expression(0.4 * B[MSY]), pos = 4, srt = 90)

  # 0.8 Bmsy
  segments(x0 = 0.8 * Bmsy, y0 = -1, x1 = 0.8 * Bmsy,
           y1 = predict(fit_spl, 0.8 * Bmsy)$y, lty = 2)
  points(x = 0.8 * Bmsy, y = predict(fit_spl, 0.8 * Bmsy)$y,
         pch = 16, col = clrs[5], cex = 1.5)
  text(x = 0.8 * Bmsy - adjustpos, y = 0,
       labels = expression(0.8 * B[MSY]), pos = 4, srt = 90)

  mtext(side = 1, outer = TRUE, text = "Spawning Biomass (kt)",
        line = 2.5, font = 2, las = 0)
  mtext(side = 2, outer = TRUE, text = "Equilibrium Yield (kt)",
        line = 2.5, font = 2, las = 0)
} # END plotYieldCurve()


# plotHistSSB()
# Historical spawning biomass time series from all three OMs.
# Inputs:
#   DDMfit  -- DDMfit (readRDS from fit_WCVI2023_OM_mle.rds)
#   fit1S   -- fit1S  (loaded as 'reports' from fit_1S_2024.RData)
#   fit3S   -- fit3S  (loaded as 'reports' from fit_3S_2024.RData)
plotHistSSB <- function( DDMfit = DDMfit,
                         fit1S  = fit1S,
                         fit3S  = fit3S )
{
  fYear <- DDMfit$fYear
  lYear <- DDMfit$lYear

  # DDM: single stock SB_pt (1 x nT)
  SB_DDM  <- DDMfit$repOpt$SB_pt[1, ]
  yrs_DDM <- seq(from = fYear, by = 1, length.out = length(SB_DDM))

  # 1S: single stock SB_pt (1 x nT)
  SB_1S  <- fit1S$repOpt$SB_pt[1, ]
  yrs_1S <- seq(
    from   = fit1S$fYear,
    by     = 1,
    length.out = length(SB_1S)
  )

  # 3S: three stocks SB_pt (nP x nT); aggregate = rowSums
  SB_3S_mat <- fit3S$repOpt$SB_pt      # 3 x nT
  SB_3S_agg <- colSums(SB_3S_mat)
  yrs_3S    <- seq(
    from   = fit3S$fYear,
    by     = 1,
    length.out = ncol(SB_3S_mat)
  )

  # Reference lines from DDM
  B0_DDM  <- DDMfit$repOpt$B0_p
  rp_DDM  <- DDMfit$repOpt$refPts
  Fmsy    <- rp_DDM$FmsyRefPts$Fmsy_p
  Bmsy    <- rp_DDM$refCurves$SBeq_pf[
    1, which.min(abs(rp_DDM$refCurves$F - Fmsy))]

  allYrs <- c(yrs_DDM, yrs_1S, yrs_3S)
  allSB  <- c(SB_DDM, SB_1S, SB_3S_agg)

  par(mar = c(4, 5, 1, 1), oma = c(0, 0, 0, 0))
  plot(
    x    = range(allYrs),
    y    = c(0, max(allSB, na.rm = TRUE) * 1.1),
    type = "n",
    xlab = "Year",
    ylab = "Spawning Biomass (kt)",
    las  = 1
  )

  # Reference lines
  abline(h = B0_DDM, lty = 3, col = "grey50", lwd = 1.5)
  abline(h = Bmsy,   lty = 2, col = "blue",   lwd = 1.5)

  # Per-Stat Area lines for 3S-predM (dashed)
  pfmaCols <- c("tomato", "darkorange", "firebrick")
  pfmaLabs <- fit3S$stock
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
  lines(x = yrs_DDM, y = SB_DDM,  col = "black",      lwd = 2.5)
  lines(x = yrs_1S,  y = SB_1S,   col = "steelblue",  lwd = 2.5)
  lines(x = yrs_3S,  y = SB_3S_agg, col = "firebrick", lwd = 2.5)

  # Labels
  labX <- min(allYrs) + 1
  usr  <- par("usr")
  text(
    x      = max(allYrs) - 1,
    y      = B0_DDM * 1.03,
    labels = expression(B[0]),
    adj    = c(1, 0),
    col    = "grey40",
    cex    = 0.85
  )
  text(
    x      = max(allYrs) - 1,
    y      = Bmsy * 1.03,
    labels = expression(B[MSY * "|" * DDM]),
    adj    = c(1, 0),
    col    = "blue",
    cex    = 0.85
  )

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
