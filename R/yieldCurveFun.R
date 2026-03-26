# yield curve with SB on X axis and yield on Y
plotRefPtsExpC <- function( empRefCurves = SOGempRefCurves,
                            histRpt = histRpts[[1]],
                            yrRange = yrRange,
                            sIdx = 1, pIdx = 1, qProbs,
                            area = 'SoG',
                            Bref = NULL )
{
  
  nReps <- dim(saveEmpRefCurves$C_isptk)[1]

  # Compute USR from MCMC posteriors if available;
  # otherwise use dummy USR (USR line is not plotted)
  if (!is.null(histRpt$posts$SB_ipt)) {
    Dstar   <- histRpt$ctlList$hypo$Dstar - 1950
    Dscalar <- histRpt$ctlList$hypo$Dscalar
    SB_it   <- histRpt$posts$SB_ipt[, 1, ]
    nMCMC   <- dim(SB_it)[1]
    Busr_i  <- apply(
      X = Dscalar * SB_it[, Dstar],
      FUN = mean, MARGIN = 1)
    set.seed(1234)
    samples <- sample(
      x = seq_len(nMCMC),
      size = nReps,
      replace = FALSE)
    Busr_i <- Busr_i[samples[seq_len(nReps)]]
    USR_q  <- quantile(Busr_i, probs = qProbs)
  } else {
    # MLE fit: no posteriors, use zero USR
    USR_q <- rep(0, length(qProbs))
  }

  eqListQuant <- solveSimEqbriaQuant(
    empRefCurves = empRefCurves,
    yrRange      = yrRange,
    USR_q        = USR_q,
    qProbs       = qProbs)
  
  # ref curves
  C_qk <- eqListQuant$C_qk
  Y_qk <- eqListQuant$Y_qk
  B_qk <- eqListQuant$B_qk
  F_qk <- eqListQuant$F_qk
  U_qk <- eqListQuant$U_qk
  
  maxU <- max(U_qk[2,])
  maxC <- max(C_qk[2,])
  maxB <- max(B_qk[2,])
  
  # Ref pts
  Umsy_q  <- eqListQuant$Umsy_q
  Bmsy_q  <- eqListQuant$Bmsy_q
  MSY_q   <- eqListQuant$MSY_q
  
  
  Ucrash_q <- eqListQuant$Ucrash_q
  Ccrash_q <- eqListQuant$Ccrash_q
  Bcrash_q <- eqListQuant$Bcrash_q
  
  ##
  empUmsy_q <- round(Umsy_q,2)
  empBmsy_q <- round(Bmsy_q,2)
  empB0_q   <- round(B_qk[,1],2)
  empMSY_q  <- round(MSY_q,2)
  
  Uusr_q    <- eqListQuant$Uusr_q
  YeqUSR_q  <- eqListQuant$YeqUSR_q
  
  maxUidx <- length(U_qk[2,])
  
  maxU <- 1.3*max(Ucrash_q)
  
  qIdx <- which(qProbs %in% range(qProbs))
  medIdx <- which(qProbs == 0.5)
  
  B_itk <- empRefCurves$B_isptk[,1,1,,]
  
  if (area == 'SoG'){
    adjustpos <- 3
    msy_x = 10
    label <- c(expression(B["0 | OM 1"]),expression("LRP, 0.3"*B["0 | OM 1"]),expression(B['MSY | OM 1']),
               expression("0.4"*B['MSY | OM 1']),expression("0.8"*B['MSY | OM 1']),
               expression("Provisional USR"["OM 1"]))
    msyLab <- "MSY | OM 1"
  }
  else if (area == 'PRD'){
    adjustpos <- 1.7
    msy_x = 5
    label <- c(expression(B["0 | DDM OM"]),expression("LRP, 0.3"*B["0 | DDM OM"]),expression(B['MSY | DDM OM']),
               expression("0.4"*B['MSY | DDM OM']),expression("0.8"*B['MSY | DDM OM']),
               expression("Provisional USR"["DDM OM"]))
    msyLab <- "MSY | DDM OM"
  }
  else if (area == 'WCVI'){
    adjustpos <- 1.7
    msy_x = 5
    label <- c(expression(B["0 | DDM OM"]),expression("LRP, 0.3"*B["0 | DDM OM"]),expression(B['MSY | DDM OM']),
               expression("0.4"*B['MSY | DDM OM']),expression("0.8"*B['MSY | DDM OM']),
               expression("Provisional USR"["DDM OM"]))
    msyLab <- "MSY | DDM OM"
  }
  
  YLIM <- c(0,max(Y_qk[medIdx,])*1.2)
  XLIM <- c(0,max(B_qk[medIdx,])*1.1)
  
  
  
  plotLayout <- matrix( c(1),
                        ncol = 1, byrow = T)
  layout(plotLayout)
  
  par(mar = c(.25,.25,.25,.25),
      oma = c(4,4,1,1), xaxs = 'r', yaxs = 'r', las = 1)
  
  clrs <- RColorBrewer::brewer.pal(8,"Dark2")
  
  # plot yield
  plot( x = c(-0.3,maxB), y = c(0,0.25),ylim = YLIM,xlim = XLIM,
        axes = FALSE, type = "n" )
  axis(side = 2, las =1, cex.axis = 1)
  axis(side = 1, cex.axis = 1)
  grid()
  box()
  # lines( x = B_qk[medIdx,1:maxUidx], y = Y_qk[medIdx,1:maxUidx], col = "black", lwd = 3, lty = 1)
  
  # Trim noisy high-harvest-rate tail where B < 1 kt
  Braw <- B_qk[medIdx, 1:maxUidx]
  Yraw <- Y_qk[medIdx, 1:maxUidx]
  keep <- Braw > 1
  Braw <- Braw[keep]
  Yraw <- Yraw[keep]

  # Fit spline directly on B vs Y with higher spar
  fit <- smooth.spline(x = Braw, y = Yraw, spar = 0.7)
  # Extend curve to origin
  lines(x = c(0, fit$x), y = c(0, fit$y),
        col = "black", lwd = 3)
  
  # B0
  points(x = empB0_q[2], y = predict(fit, empB0_q[2])$y, pch = 16, col = clrs[1], cex = 1.5)
  text(x = empB0_q[2]+adjustpos, y = 0, label[1], pos = 4, srt = 90)
  
  # 0.3B0
  segments(x0 = predict(fit, 0.3*empB0_q[2])$x, y0 = -1,
           x1 = predict(fit, 0.3*empB0_q[2])$x, y1 = predict(fit, 0.3*empB0_q[2])$y, lty = 2)
  points(predict(fit, 0.3*empB0_q[2]), pch = 16, col = clrs[2], cex = 1.5)
  text(x = 0.3*empB0_q[2]-adjustpos, y = 0, label[2], pos = 4, srt = 90)
  
  # MSY and label
  segments(x0 = -5, y0 = predict(fit,empBmsy_q[2])$y, x1 = empBmsy_q[2], y1 = predict(fit, empBmsy_q[2])$y, lty = 2)
  text(x = msy_x, y = predict(fit,empBmsy_q[2])$y, msyLab, pos = 3)
  
  
  # Bmsy
  segments(x0 = empBmsy_q[2], y0 = -1, x1 = empBmsy_q[2], y1 = predict(fit,empBmsy_q[2])$y, lty = 2)
  points(x = empBmsy_q[2], y = predict(fit, empBmsy_q[2])$y, pch = 16, col = clrs[3], cex = 1.5)
  text(x = empBmsy_q[2]-adjustpos, y = 0, label[3], pos = 4, srt = 90)
  
  # 0.4Bmsy
  segments(x0 = 0.4*empBmsy_q[2], y0 = -1,
           x1 = 0.4*empBmsy_q[2], y1 = predict(fit, 0.4*empBmsy_q[2])$y, lty = 2)
  points(predict(fit, 0.4*empBmsy_q[2]), pch = 16, col = clrs[4], cex = 1.5)
  text(x = 0.4*empBmsy_q[2]-adjustpos, y = 0, label[4], pos = 4, srt = 90)
  
  # 0.8Bmsy
  segments(x0 = 0.8*empBmsy_q[2], y0 = -1,
           x1 = 0.8*empBmsy_q[2], y1 = predict(fit, 0.8*empBmsy_q[2])$y, lty = 2)
  points(predict(fit, 0.8*empBmsy_q[2]), pch = 16, col = clrs[5], cex = 1.5)
  text(x = 0.8*empBmsy_q[2]-adjustpos, y = 0, label[5], pos = 4, srt = 90)
  
  # USR
  # segments(x0 = USR_q[2], y0 = -1,
  #          x1 = USR_q[2], y1 = predict(fit, USR_q[2])$y, lty = 2)
  # points(x = USR_q[2], y = predict(fit, USR_q[2])$y, col = clrs[6], pch = 16, cex = 1.5)
  # text(x = USR_q[2]-adjustpos, y = 0, label[6], pos = 4, srt = 90, cex = 1)
  
  # if(history == T){
  #   # historic data: 2004 - 2023
  #   endYr <- 2023
  #   startYr <- endYr - histYrs +1
  #   B_it <- array(0,dim = c(200,histYrs))
  #   Y_it <- array(0,dim = c(200,histYrs))
  #   load(paste0("Outputs/",folder,"/sim_parBat",folder,1,"/sim_parBat",folder,1,".RData"))
  #   startIdx <- 73 - histYrs +1
  #   B_it <- blob$om$SB_ispt[,1,1,startIdx:73]
  #   Y_it <- blob$om$C_ispt[,1,1,startIdx:73]
  #
  #   B_t <- apply(B_it, MARGIN = 2, FUN = median)
  #   catch <- read.csv('history/catchData.csv')
  #   group <- group_by(catch, Year) %>% summarise(C = sum(Value)) %>%
  #     filter(Year >= startYr)
  #   divClrs <- paletteer_c('grDevices::Grays', n = length(B_t))
  #
  #   lines(B_t, group$C, col = 'gray', lwd = 2)
  #   points(B_t, group$C, col = 'black', bg = divClrs, pch = 21, )
  #   text(B_t[1],group$C[1], startYr, pos = 1)
  #   text(B_t[histYrs],group$C[histYrs], endYr, pos = 1)
  # }
  
  
  # Bref from external SAR (if supplied)
  if (!is.null(Bref)) {
    segments(x0 = Bref, y0 = -1,
             x1 = Bref,
             y1 = predict(fit, Bref)$y, lty = 2)
    points(x = Bref, y = predict(fit, Bref)$y,
           pch = 16, col = clrs[6], cex = 1.5)
    text(x = Bref - adjustpos, y = 0,
         expression(B["ref | 2023 DDM CSAS"]),
         pos = 4, srt = 90)
  }

  mtext( side = 1, outer = TRUE, text = "Spawning Biomass (kt)",
         line = 2.5, font = 2, las = 0)
  mtext(side = 2, text = "Yield (kt)", line = 2.5, font = 2, las = 0)


}# plotRefPtsExpC()

solveSimEqbriaQuant <- function(  empRefCurves,
                                  sIdx = 1, pIdx = 1,
                                  USR_q = 38, qProbs,
                                  reps = NULL,
                                  yrRange = -(150:50))
{
  
  
  
  # Pull empirical ref curves
  if(is.null(reps))
    reps <- 1:dim(empRefCurves$C_isptk)[1]
  
  pT <- dim(empRefCurves$C_isptk)[4]
  yrRange <- pT+yrRange
  yrRange <- round(yrRange)
  
  C_qtk <- apply(X = empRefCurves$C_isptk[reps,sIdx,pIdx,yrRange,,drop = FALSE], FUN = quantile, MARGIN = 4:5, probs = qProbs)
  Y_qtk <- apply(X = empRefCurves$Y_isptk[reps,sIdx,pIdx,yrRange,,drop = FALSE], FUN = quantile, MARGIN = 4:5, probs = qProbs)
  B_qtk <- apply(X = empRefCurves$B_isptk[reps,sIdx,pIdx,yrRange,,drop = FALSE], FUN = quantile, MARGIN = 4:5, probs = qProbs)
  F_qtk <- apply(X = empRefCurves$F_isptk[reps,sIdx,pIdx,yrRange,,drop = FALSE], FUN = quantile, MARGIN = 4:5, probs = qProbs)
  U_qtk <- apply(X = empRefCurves$U_isptk[reps,sIdx,pIdx,yrRange,,drop = FALSE], FUN = quantile, MARGIN = 4:5, probs = qProbs)
  M_qtk <- apply(X = empRefCurves$M_isptk[reps,sIdx,pIdx,yrRange,,drop = FALSE], FUN = quantile, MARGIN = 4:5, probs = qProbs)
  R_qtk <- apply(X = empRefCurves$R_isptk[reps,sIdx,pIdx,yrRange,,drop = FALSE], FUN = quantile, MARGIN = 4:5, probs = qProbs)
  
  C_qk <- apply(X = C_qtk, FUN = median, MARGIN = c(1,3))
  Y_qk <- apply(X = Y_qtk, FUN = median, MARGIN = c(1,3))
  B_qk <- apply(X = B_qtk, FUN = median, MARGIN = c(1,3))
  F_qk <- apply(X = F_qtk, FUN = median, MARGIN = c(1,3))
  U_qk <- apply(X = U_qtk, FUN = median, MARGIN = c(1,3))
  M_qk <- apply(X = M_qtk, FUN = median, MARGIN = c(1,3))
  R_qk <- apply(X = R_qtk, FUN = median, MARGIN = c(1,3))
  
  # Order
  Forder <- order(F_qk[which(qProbs == 0.5),])
  C_qk <- C_qk[,Forder]
  Y_qk <- Y_qk[,Forder]
  B_qk <- B_qk[,Forder]
  F_qk <- F_qk[,Forder]
  U_qk <- U_qk[,Forder]
  M_qk <- M_qk[,Forder]
  R_qk <- R_qk[,Forder]
  
  B0_q      <- c(0)
  Bmsy_q    <- c(0)
  Umsy_q    <- c(0)
  Bmsy_q    <- c(0)
  MSY_q     <- c(0)
  Ucrash_q  <- c(0)
  Ccrash_q  <- c(0)
  Bcrash_q  <- c(0)
  Mcrash_q  <- c(0)
  Rcrash_q  <- c(0)
  Uusr_q    <- c(0)
  YeqUSR_q  <- c(0)
  
  
  # Make spline
  for( qq in 1:length(qProbs))
  {
    
    
    CUspline <- splinefun(x = U_qk[qq,], y = C_qk[qq,])
    YUspline <- splinefun(x = U_qk[qq,], y = Y_qk[qq,])
    BUspline <- splinefun(x = U_qk[qq,], y = B_qk[qq,])
    MUspline <- splinefun(x = U_qk[qq,], y = M_qk[qq,])
    RUspline <- splinefun(x = U_qk[qq,], y = R_qk[qq,])
    
    BUSRspline <- splinefun(x = U_qk[qq,], y = B_qk[qq,] - USR_q[qq])
    
    
    maxUidx <- length(B_qk[qq,])
    
    # max(which(B_qk[qq,] > 1e-4))
    # cat("maxU = ", U_k[maxUidx],"\n")
    # if(U_k[maxUidx] < 1e-4)
    #   browser()
    # Note, the curve isn't smooth, so we need
    # to account for local optima when finding
    # these crash points - optimise and uniroot are
    # not super clever.
    minUcrash <- max(0,max(Ucrash_q))
    Useq <- seq(from = 0, to = U_qk[qq,maxUidx], length.out = 1000)
    Yseq <- YUspline(Useq,deriv = 1)
    
    Ucrash <- try(optimize( interval = c(minUcrash,U_qk[qq,maxUidx]), f=YUspline, deriv=1),
                   silent = TRUE)
    
    
    if( inherits(Ucrash, "try-error") ){
      Ucrash <- list()
      Ucrash$minimum <- 0
    }
    
    
    k <- 0
    
    
    while(any(Yseq < YUspline(Ucrash$minimum,deriv = 1)))
    {
      minUidx <- which.min(Yseq)
      if(minUidx == length(Yseq))
        break
      
      newInt <- c(Useq[minUidx - 3], Useq[minUidx + 3])
      
      Ucrash <- try(optimize( interval = newInt, f=YUspline, deriv=1),
                     silent = TRUE)
      k <- k+1
      
      if(k == 4)
        break
    }
    
    Ucrash <- Ucrash$minimum
    
    
    if(Ucrash > 0  )
    {
      minUmsy <- max(0.01,max(Umsy_q))
      
      if(sign(YUspline(minUmsy,deriv = 1)) != sign(YUspline(Ucrash,deriv = 1)))
      {
        Umsy  <- try(uniroot(   interval = c(minUmsy,Ucrash),
                                f = YUspline,
                                deriv = 1 )$root,
                     silent = TRUE)
        
      }
      else Umsy <- Ucrash
      
      # Need to add a check that the estimate of MSY from this is great
      
      Yseq <- YUspline(Useq)
      k <- 0
      
      while(any(Yseq > YUspline(Umsy)))
      {
        
        maxUidx <- which.max(Yseq)
        newInt <- c(Useq[maxUidx-5],Useq[maxUidx+5])
        
        Umsy  <- try(uniroot( interval = newInt,
                              f = YUspline,
                              deriv = 1 )$root,
                     silent = TRUE)
        k <- k+1
        
        if(k == 4)
          break
      }
      
    }
    else Umsy <- 0
    
    Uusr <- 0
    if(USR_q[qq] < B_qk[qq,1])
      Uusr <- try(uniroot(  interval = c(0.0,Umsy),
                            f = BUSRspline,
                            deriv = 0 )$root,
                   silent = TRUE)
    
    if( inherits(Umsy, "try-error") )
      Umsy <- 0
    
    if( inherits(Uusr, "try-error") )
    {
      Uusr <- 0
    }
    
    
    
    B0_q[qq]      <- BUspline(0)
    Bmsy_q[qq]    <- BUspline(Umsy)
    Umsy_q[qq]    <- Umsy
    MSY_q[qq]     <- YUspline(Umsy)
    Ucrash_q[qq]  <- Ucrash
    Ccrash_q[qq]  <- YUspline(Ucrash)
    Bcrash_q[qq]  <- BUspline(Ucrash)
    Mcrash_q[qq]  <- MUspline(Ucrash)
    Rcrash_q[qq]  <- RUspline(Ucrash)
    Uusr_q[qq]    <- Uusr
    YeqUSR_q[qq]  <- YUspline(Uusr)
    
  }
  
  
  outList <- list(  C_qk = C_qk,
                    Y_qk = Y_qk,
                    U_qk = U_qk,
                    F_qk = F_qk,
                    B_qk = B_qk,
                    M_qk = M_qk,
                    R_qk = R_qk,
                    B0_q      = B0_q,
                    Bmsy_q    = Bmsy_q,
                    Umsy_q    = Umsy_q,
                    MSY_q     = MSY_q,
                    Ucrash_q  = Ucrash_q,
                    Ccrash_q  = Ccrash_q,
                    Bcrash_q  = Bcrash_q,
                    Mcrash_q  = Mcrash_q,
                    Rcrash_q  = Rcrash_q,
                    Uusr_q    = Uusr_q,
                    YeqUSR_q  = YeqUSR_q )
  
  
  outList
}