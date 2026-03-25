# plotResults.R
#
# Results-section figure wrappers for the NTC Phase 3 report.
# Calls utility functions in ms3Rplots.R where applicable.
# All figures use base plot (no ggplot2).
#
# Functions:
#   plotSimProj()   -- Simulation projection tulip plots (slides 6, 10, 17)
#   plotRetroFits() -- Retrospective assessment fits     (slides 7, 11, 18)
#   plotHCRPhase()  -- HCR phase scatter plots          (slides 8, 12, 19)


# plotSimProj()
# Three-panel simulation projection figure.
# Panel 1: SSB tulip with B0, LRP, Bref reference lines.
# Panel 2: Commercial catch tulip by fleet.
# Panel 3: SOK/FSC catch tulip.
# Inputs:
#   obj  -- simulation blob (e.g. simDDM, sim1S, sim3S_MP1)
plotSimProj <- function( obj = simDDM )
{
  # Dimensions and time
  fYear    <- obj$ctlList$opMod$fYear
  nT       <- obj$om$nT
  nS       <- obj$om$nS
  nP       <- obj$om$nP
  tMP      <- obj$om$tMP
  yrs      <- seq(from = fYear, by = 1, length.out = nT)
  goodReps <- obj$goodReps
  nReps    <- sum(goodReps)

  # Reference values
  rp1  <- obj$rp[[1]]
  B0   <- rp1$B0_sp[1]
  Bref <- obj$mp$hcr$Bref_ispt[1, 1, 1, tMP]
  LRP  <- 0.3 * B0

  # Fleet indices: exclude SOK/FSC fleets (type 2) from commercial panel
  sokF  <- which(obj$om$fleetType_f == 2)
  commF <- setdiff(obj$ctlList$opMod$commGears, sokF)
  fleetLabs <- obj$ctlList$opMod$fleets

  # Aggregate SSB across species (herring only: s=1) and areas
  SB_ispt <- obj$om$SB_ispt[goodReps, 1, , , drop = FALSE]  # i x 1 x P x T
  SB_ipt  <- apply(
    X      = SB_ispt,
    FUN    = sum,
    MARGIN = c(1, 3, 4),
    na.rm  = TRUE
  )
  SB_qpt <- apply(
    X      = SB_ipt,
    FUN    = quantile,
    MARGIN = c(2, 3),
    probs  = c(0.025, 0.5, 0.975),
    na.rm  = TRUE
  )

  # Commercial catch aggregate
  C_ispft  <- obj$om$C_ispft[goodReps, 1, , , , drop = FALSE] # i x 1 x P x F x T
  C_ipft   <- apply(
    X      = C_ispft[, , , commF, , drop = FALSE],
    FUN    = sum,
    MARGIN = c(1, 3, 4, 5),
    na.rm  = TRUE
  )                                                            # i x P x nCommF x T
  C_ipt_comm <- apply(
    X      = C_ipft,
    FUN    = sum,
    MARGIN = c(1, 2, 4),
    na.rm  = TRUE
  )                                                            # i x P x T
  C_qpt_comm <- apply(
    X      = apply(C_ipt_comm, c(1, 3), sum),                 # i x T
    FUN    = quantile,
    MARGIN = 2,
    probs  = c(0.025, 0.5, 0.975),
    na.rm  = TRUE
  )

  # SOK/FSC catch aggregate
  C_ipt_sok <- apply(
    X      = C_ispft[, , , sokF, , drop = FALSE],
    FUN    = sum,
    MARGIN = c(1, 2, 3, 5),
    na.rm  = TRUE
  )                                                            # i x 1 x P x T
  C_qpt_sok <- apply(
    X      = apply(C_ipt_sok, c(1, 3, 4), sum),               # i x P x T
    FUN    = quantile,
    MARGIN = c(2, 3),
    probs  = c(0.025, 0.5, 0.975),
    na.rm  = TRUE
  )

  traces <- sample(seq_len(nReps), size = min(3, nReps))

  op <- par(no.readonly = TRUE)
  on.exit(par(op))

  par(mfrow = c(3, 1),
      mar   = c(2, 6, 1, 1),
      oma   = c(3, 0, 0, 0))

  # ---- Panel 1: SSB aggregate ----
  plot(
    x    = range(yrs),
    y    = c(0, max(SB_qpt[3, 1, ], B0, na.rm = TRUE) * 1.1),
    type = "n",
    xlab = "",
    ylab = "Spawning Biomass (kt)",
    las  = 1,
    xaxt = "n"
  )
  axis(side = 1, labels = FALSE)
  polygon(
    x      = c(yrs, rev(yrs)),
    y      = c(SB_qpt[1, 1, ], rev(SB_qpt[3, 1, ])),
    col    = "grey75",
    border = NA
  )
  for (tIdx in traces)
    lines(x = yrs, y = SB_ipt[tIdx, 1, ], lwd = 0.5, col = "grey60")
  lines(x = yrs, y = SB_qpt[2, 1, ], lwd = 2.5)
  abline(v = yrs[tMP] - 0.5, lty = 2, col = "grey50", lwd = 1)
  abline(h = B0,   lty = 3, col = "grey40", lwd = 1.5)
  abline(h = Bref, lty = 2, col = "blue",   lwd = 2)
  abline(h = LRP,  lty = 2, col = "red",    lwd = 2)
  legend(
    x      = "topright",
    bty    = "n",
    legend = c(
      "Median SSB",
      expression(B[0]), expression(B[ref]), "LRP"
    ),
    col = c("black", "grey40", "blue", "red"),
    lty = c(1, 3, 2, 2),
    lwd = c(2.5, 1.5, 2, 2),
    cex = 0.75
  )

  # ---- Panel 2: commercial catch — stacked fleet bars + projection whiskers ----
  bw      <- 0.35   # half-width of each bar in year units
  projIdx <- tMP:nT # projection years only (for whiskers)
  plot(
    x    = range(yrs),
    y    = c(0, max(C_qpt_comm[3, ], na.rm = TRUE) * 1.1),
    type = "n",
    xlab = "",
    ylab = "Commercial Catch (kt)",
    las  = 1,
    xaxt = "n"
  )
  axis(side = 1, labels = FALSE)
  abline(v = yrs[tMP] - 0.5, lty = 2, col = "grey50", lwd = 1)
  # Stacked coloured fleet bars (median)
  if (length(commF) > 0) {
    fCols  <- RColorBrewer::brewer.pal(
      n = max(3, length(commF)), name = "Dark2")
    hasBar <- logical(length(commF))
    # Build per-fleet median matrix: nF x nT
    C_med_ft <- matrix(0, nrow = length(commF), ncol = nT)
    for (fi in seq_along(commF)) {
      f <- commF[fi]
      C_ft <- apply(C_ispft[, , , f, , drop = FALSE], c(1, 3, 5), sum)
      C_med_ft[fi, ] <- apply(apply(C_ft, c(1, 3), sum), 2, median)
    }
    # Draw stacked bars bottom-up
    bottom <- rep(0, nT)
    for (fi in seq_along(commF)) {
      top <- bottom + C_med_ft[fi, ]
      if (sum(top, na.rm = TRUE) > 0) {
        rect(xleft   = yrs - bw, xright  = yrs + bw,
             ybottom = bottom,    ytop    = top,
             col = adjustcolor(fCols[fi], alpha.f = 0.85),
             border = NA)
        hasBar[fi] <- TRUE
      }
      bottom <- top
    }
    # Whiskers on total stack — projection years only
    segments(x0  = yrs[projIdx], x1 = yrs[projIdx],
             y0  = C_qpt_comm[1, projIdx],
             y1  = C_qpt_comm[3, projIdx],
             lwd = 1, col = "grey20")
    if (any(hasBar))
      legend(
        x      = "topright",
        bty    = "n",
        legend = fleetLabs[commF][hasBar],
        fill   = fCols[seq_along(commF)][hasBar],
        border = NA,
        cex    = 0.7
      )
  }

  # ---- Panel 3: SOK/FSC catch — bars + projection whiskers (scaled kt -> t) ----
  sok_scale   <- 1000
  C_qpt_sok_t <- C_qpt_sok * sok_scale
  maxSOK <- max(C_qpt_sok_t[3, 1, ], na.rm = TRUE)
  if (maxSOK == 0 || is.na(maxSOK)) maxSOK <- 0.1
  plot(
    x    = range(yrs),
    y    = c(0, maxSOK * 1.1),
    type = "n",
    xlab = "",
    ylab = "SOK/FSC Catch (t)",
    las  = 1
  )
  abline(v = yrs[tMP] - 0.5, lty = 2, col = "grey50", lwd = 1)
  rect(xleft   = yrs - bw,           xright  = yrs + bw,
       ybottom = 0,                   ytop    = C_qpt_sok_t[2, 1, ],
       col = "steelblue3", border = NA)
  segments(x0  = yrs[projIdx], x1 = yrs[projIdx],
           y0  = C_qpt_sok_t[1, 1, projIdx],
           y1  = C_qpt_sok_t[3, 1, projIdx],
           lwd = 1, col = "steelblue4")

  mtext(side = 1, outer = TRUE, text = "Year", line = 1.5)
} # END plotSimProj()


# plotRetroFits()
# Retrospective spawning biomass fits: OM truth vs assessment fits.
# Self-contained: does not call plotRetroSB() because retroSB_itspt is all-NA
# in predM blobs. Uses retroSB_agg_itt when available (predM blobs), otherwise
# aggregates retroSB_itspt across species and stocks (DDM blob).
# Inputs:
#   obj   -- simulation blob
#   iRep  -- replicate to display (default 1)
plotRetroFits <- function( obj = simDDM, iRep = 1 )
{
  tMP   <- obj$om$tMP
  nT    <- obj$om$nT
  nP    <- obj$om$nP
  fYear <- obj$ctlList$opMod$fYear
  pT    <- nT - tMP + 1
  yrs   <- seq(from = fYear, by = 1, length.out = nT)

  # True OM SSB (herring only, s=1; aggregate across stocks)
  SB_pt   <- apply(
    X      = obj$om$SB_ispt[iRep, 1, , , drop = FALSE],
    FUN    = sum,
    MARGIN = c(3, 4),
    na.rm  = TRUE
  )                                             # nP x nT -> sum over p below
  SBtrue  <- colSums(SB_pt)                    # length nT

  # Retrospective assessment fits (pT x nT)
  # predM blobs store this in retroSB_agg_itt; DDM uses retroSB_itspt
  if (!is.null(obj$mp$assess$retroSB_agg_itt) &&
      length(dim(obj$mp$assess$retroSB_agg_itt)) == 3) {
    retroSB_tt <- obj$mp$assess$retroSB_agg_itt[iRep, , ]   # pT x nT
  } else {
    rSB_rep    <- obj$mp$assess$retroSB_itspt[iRep, , , , , drop = FALSE]
    retroSB_tt <- apply(rSB_rep, c(2, 5), sum, na.rm = FALSE)  # pT x nT
  }
  retroSB_tt[retroSB_tt < 0] <- NA

  # Spawn survey index (fleet 5 = dive survey, herring s=1, stock p=1)
  survIdx <- obj$mp$data$I_ispft[iRep, 1, 1, 5, ]   # length nT

  # Reference points
  Bref <- c(obj$ctlList$mp$hcr$Bref_p,
            obj$ctlList$mp$hcr$Bref_s)[1]
  B0   <- obj$rp[[1]]$B0_sp[1]
  LRP  <- 0.3 * B0

  ymax <- max(SBtrue, retroSB_tt, na.rm = TRUE) * 1.15

  par(mar = c(4, 5, 1, 1), oma = c(0, 0, 0, 0))
  plot(
    x    = range(yrs),
    y    = c(0, ymax),
    type = "n",
    xlab = "Year",
    ylab = "Spawning Biomass (kt)",
    las  = 1
  )

  # Reference lines
  abline(h = B0,   lty = 3, col = "grey50",  lwd = 1.2)
  abline(h = Bref, lty = 2, col = "blue",    lwd = 1.5)
  abline(h = LRP,  lty = 2, col = "red",     lwd = 1.2)
  abline(v = yrs[tMP] - 0.5, lty = 2, col = "black", lwd = 0.8)

  # Grey retrospective assessment lines
  for (tt in seq_len(pT))
    lines(x = yrs, y = retroSB_tt[tt, ], col = "grey60", lwd = 1)

  # Final assessment fit (most recent retro step)
  lines(x = yrs, y = retroSB_tt[pT, ], col = "steelblue", lwd = 2.5)

  # True OM SSB
  lines(x = yrs, y = SBtrue, col = "red", lwd = 2.5)

  # Scaled survey index points (non-NA only)
  ok <- !is.na(survIdx) & survIdx > 0
  if (any(ok)) {
    scale <- median(SBtrue[ok] / survIdx[ok], na.rm = TRUE)
    points(x = yrs[ok], y = survIdx[ok] * scale,
           pch = 16, col = "black", cex = 0.6)
  }

  legend(
    x      = "topleft",
    bty    = "n",
    legend = c("True SSB (OM)",
               "Final assessment",
               "Earlier assessments",
               "Survey index (scaled)"),
    col    = c("red", "steelblue", "grey60", "black"),
    lty    = c(1, 1, 1, NA),
    pch    = c(NA, NA, NA, 16),
    lwd    = c(2.5, 2.5, 1, NA),
    cex    = 0.8
  )
} # END plotRetroFits()


# plotHCRScatter()
# Two-panel HCR scatter plot: U vs true OM biomass (top) and TAC vs true OM
# biomass (bottom). Overlays the theoretical HCR curve so the reader can see
# where simulated outcomes fall relative to the rule.
# Axes are in raw units (kt), matching the reference slide style.
# Inputs:
#   obj  -- simulation blob (e.g. simDDM, sim1S, sim3S_MP1)
plotHCRScatter <- function( obj = simDDM )
{
  tMP      <- obj$om$tMP
  nT       <- obj$om$nT
  goodReps <- obj$goodReps
  pT       <- nT - tMP + 1

  # Bref: fixed HCR parameter (17.35 kt). DDM blobs store it as Bref_p;
  # predM blobs store it as Bref_s.
  Bref <- c(obj$ctlList$mp$hcr$Bref_p,
            obj$ctlList$mp$hcr$Bref_s)[1]
  Uref <- obj$ctlList$mp$hcr$inputF_s[1]

  # True OM spawning biomass (aggregate across stocks, herring species only)
  SB_ispt <- obj$om$SB_ispt[goodReps, 1, , , drop = FALSE]   # i x 1 x P x T
  B_it    <- matrix(
    apply(X = SB_ispt, FUN = sum, MARGIN = c(1, 4), na.rm = TRUE)[, tMP:nT],
    nrow = sum(goodReps), ncol = pT
  )

  # Aggregate TAC across all stocks (sum over p dimension)
  TAC_ispt_sub <- obj$mp$hcr$TAC_ispt[goodReps, 1, , , drop = FALSE]
  TAC_it <- matrix(
    apply(X = TAC_ispt_sub, FUN = sum, MARGIN = c(1, 4), na.rm = TRUE)[, tMP:nT],
    nrow = sum(goodReps), ncol = pT
  )
  U_it    <- TAC_it / (B_it + TAC_it)

  # Flatten and filter bad values
  Bvec <- as.vector(B_it)
  Uvec <- as.vector(U_it)
  Tvec <- as.vector(TAC_it)
  ok   <- is.finite(Bvec) & Bvec > 0 & is.finite(Uvec) & is.finite(Tvec)
  Bvec <- Bvec[ok];  Uvec <- Uvec[ok];  Tvec <- Tvec[ok]

  # HCR curve
  Bmax  <- max(Bvec, Bref * 3, na.rm = TRUE) * 1.05
  Bseq  <- seq(0, Bmax, length.out = 500)
  Ufun  <- function(B) ifelse(B <= Bref, 0, pmin(1 - Bref / B, Uref))
  Ucurve  <- Ufun(Bseq)
  TACcurve <- Ucurve * Bseq

  # Boundary between ramp and flat-U zones: Bcap = Bref / (1 - Uref)
  Bcap  <- Bref / (1 - Uref)
  ptCol <- adjustcolor("grey40", alpha.f = 0.25)

  # Y limits: cover scatter points plus a small margin
  Umax   <- max(Uvec, Uref, na.rm = TRUE) * 1.15
  TACmax <- max(TACcurve, Tvec, na.rm = TRUE) * 1.1
  if (!is.finite(TACmax) || TACmax <= 0) TACmax <- 1

  op <- par(no.readonly = TRUE)
  on.exit(par(op))

  par(mfrow = c(2, 1),
      mar   = c(2, 5, 1, 1),
      oma   = c(3, 0, 0, 0))

  # ---- Panel 1: U vs B ----
  plot(
    x    = c(0, Bmax),
    y    = c(0, Umax),
    type = "n",
    xlab = "",
    ylab = expression(paste("Harvest Rate  ", italic(U))),
    las  = 1,
    xaxt = "n"
  )
  axis(side = 1, labels = FALSE)
  # Coloured regions: red (closed), yellow (ramp), blue (full harvest)
  polygon(x = c(0, Bref, Bref, 0),
          y = c(0, 0, Umax, Umax),
          col = adjustcolor("red", alpha.f = 0.10), border = NA)
  polygon(x = c(Bref, Bcap, Bcap, Bref),
          y = c(0, 0, Umax, Umax),
          col = adjustcolor("goldenrod", alpha.f = 0.15), border = NA)
  polygon(x = c(Bcap, Bmax, Bmax, Bcap),
          y = c(0, 0, Umax, Umax),
          col = adjustcolor("steelblue", alpha.f = 0.12), border = NA)
  abline(v = Bref, lty = 2, col = "red",       lwd = 1.5)
  abline(v = Bcap, lty = 2, col = "steelblue", lwd = 1.2)
  abline(h = Uref, lty = 2, col = "grey50",    lwd = 1.2)
  points(x = Bvec, y = Uvec, pch = 16, cex = 0.45, col = ptCol)
  lines(x = Bseq, y = Ucurve, lwd = 2.5, col = "black")
  text(x = Bref, y = Umax * 0.97,
       labels = expression(B[ref]), adj = c(1.1, 1),
       col = "red", cex = 0.85)
  text(x = Bmax * 0.97, y = Uref * 1.04,
       labels = expression(U[ref]), adj = c(1, 0),
       col = "grey40", cex = 0.85)

  # ---- Panel 2: TAC vs B ----
  plot(
    x    = c(0, Bmax),
    y    = c(0, TACmax),
    type = "n",
    xlab = "",
    ylab = "TAC (kt)",
    las  = 1
  )
  polygon(x = c(0, Bref, Bref, 0),
          y = c(0, 0, TACmax, TACmax),
          col = adjustcolor("red", alpha.f = 0.10), border = NA)
  polygon(x = c(Bref, Bcap, Bcap, Bref),
          y = c(0, 0, TACmax, TACmax),
          col = adjustcolor("goldenrod", alpha.f = 0.15), border = NA)
  polygon(x = c(Bcap, Bmax, Bmax, Bcap),
          y = c(0, 0, TACmax, TACmax),
          col = adjustcolor("steelblue", alpha.f = 0.12), border = NA)
  abline(v = Bref, lty = 2, col = "red",       lwd = 1.5)
  abline(v = Bcap, lty = 2, col = "steelblue", lwd = 1.2)
  points(x = Bvec, y = Tvec, pch = 16, cex = 0.45, col = ptCol)
  lines(x = Bseq, y = TACcurve, lwd = 2.5, col = "black")
  text(x = Bref, y = TACmax * 0.97,
       labels = expression(B[ref]), adj = c(1.1, 1),
       col = "red", cex = 0.85)

  mtext(side = 1, outer = TRUE,
        text = "True Spawning Biomass (kt)", line = 1.5)
} # END plotHCRScatter()
