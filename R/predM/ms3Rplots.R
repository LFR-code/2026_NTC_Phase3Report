# <><><><><><><><><><><><><><><><><><><><><><><><><><><>
# ms3Rplots.R
#
# Plots for ms3R.
#
# Author: SDN Johnson
# Date: October 8, 2019
#
# Last Update: July 9, 2020
#
# <><><><><><><><><><><><><><><><><><><><><><><><><><><>

library(ClimProjDiags)

compareTACs <- function(  obj1 = blob,
                          obj2 = blob,
                          iRep = 1,
                          sIdx = 1,
                          pIdx = 1,
                          xAxLab = "Blob1",
                          yAxLab = "Blob2" )
{
  # Get time range of projections
  tMP <- obj1$om$tMP
  nT  <- obj1$om$nT

  years <- seq(from = 1951, by = 1, length.out = nT)

  # Get TAC series
  C1_t   <- obj1$om$C_ispt[iRep,sIdx,pIdx,tMP:nT]
  C2_t   <- obj2$om$C_ispt[iRep,sIdx,pIdx,tMP:nT]

  TAC1_t <- obj1$mp$hcr$TAC_ispt[iRep,sIdx,pIdx,tMP:nT]
  TAC2_t <- obj2$mp$hcr$TAC_ispt[iRep,sIdx,pIdx,tMP:nT]

  pT <- nT - tMP + 1

  colFunc <- colorRampPalette(c("grey20","grey75"))
  colVec  <- colFunc(pT)

  plot( x = c(0,max(TAC1_t)),
        y = c(0,max(TAC2_t)),
        type = "n",
        xlab = xAxLab,
        ylab = yAxLab)

    points(x = TAC1_t, y = TAC2_t, col = colVec, pch = 16)
    points(x = C1_t, y = C2_t, col = colVec, pch = 21, lwd = 2)
    abline( a = 0, b = 1, lwd = 2 )
    abline( a = 0, b = 2, lwd = 1, lty = 2 )
    legend( x =  "topleft", bty = "n",
            legend = c(years[tMP],years[nT],"y = x","y = 2x","TAC","Catch"),
            col = c(colVec[c(1,pT)],"black","black","black","black"),
            pch = c(16,16,NA,NA,16,21),
            lwd = c(NA,NA,2,1,NA,NA),
            lty = c(NA,NA,1,2,NA,NA),
            pt.lwd = c(NA,NA,NA,NA,NA,2) )

}



# plotGridTulipBtCt()
# Plots a grid of simulation envelopes
# for biomass and legal catch (TAC)
plotGridTulipBtCtUt <- function(  blobList  = wtdSims,
                                  labels    = NULL,
                                  dep       = FALSE,
                                  yLimB     = c(0,35),
                                  yLimC     = c(0,2),
                                  yLimHR    = c(0,.1),
                                  traces    = 3,
                                  traceSeed = 1234,
                                  proj      = TRUE,
                                  maxProjT  = NULL,
                                  goodReps  = NULL )
{
  # Load blobs
  nSims <- length(blobList)
  simLabs <- c()
  nTs     <- c()
  for(simIdx in 1:nSims)
  {
    simLabs[simIdx] <- blobList[[simIdx]]$ctlList$ctl$mpName
    names(blobList)[simIdx] <- simLabs[simIdx]
    nTs[simIdx] <- blob$om$nT
  }

  if(!is.null(labels))
    plotLabs <- labels
  else plotLabs <- simLabs


  # Now set up plotting area
  par(  mfcol = c(3,nSims),
        mar   = c(.1,.1,.1,.1),
        oma   = c(5,6,3,4) )

  fYear <- 1951
  nT <- blobList[[1]]$om$nT
  tMP <- blobList[[1]]$om$tMP

  yrs <- seq(from = fYear, length.out = max(nTs))
  if(is.null(maxProjT))
    maxProjT <- max(nTs) - tMP

  if(proj)
    xLim <- range(yrs[(tMP-3):(tMP + maxProjT)])
  else xLim <- range(yrs)

  if(is.null(goodReps))
    goodReps <- blobList[[1]]$goodReps

  nReps <- sum(goodReps)
  nS    <- dim(blobList[[1]]$om$SB_ispt)[2]
  nP    <- dim(blobList[[1]]$om$SB_ispt)[3]
  nT    <- dim(blobList[[1]]$om$SB_ispt)[4]


  if(traces > 0)
  {
    set.seed(traceSeed)
    traceIdx <- sample(1:nReps, size = traces )
  }


  for(simIdx in 1:nSims)
  {
    blob <- blobList[[simIdx]]

    yrs <- seq(from = fYear, length.out = nTs[simIdx])
    # cat(simLabs[simIdx],"\n")

    # Get model states
    SB_ispt   <- blobList[[simIdx]]$om$SB_ispt[which(goodReps),,,,drop = FALSE]

    nF <- blob$om$nF
    C_ispt      <- blobList[[simIdx]]$om$C_ispt[which(goodReps),,,,drop = FALSE]
    uC_ispt     <- C_ispt


    if( nF > 5 )
    {


      uC_ispft    <- blobList[[simIdx]]$om$C_ispft
      P_ispft     <- blobList[[simIdx]]$om$P_ispft

      uC_ispft[,,,6,] <- P_ispft[,,,6,]

      C_ispt  <- apply(X = uC_ispft[,,,1:6,,drop = FALSE], FUN = sum, MARGIN = c(1,2,3,5), na.rm = T)

      uC_ispt <- apply(X = uC_ispft[,,,1:6,,drop = FALSE], FUN = sum, MARGIN = c(1,2,3,5), na.rm = T)
      uC_ispt <- uC_ispt[which(goodReps),,,,drop = FALSE]
    }

    # TAC_ispt  <- blobList[[simIdx]]$mp$hcr$TAC_ispt[1:nReps,,,,drop = FALSE]

    U_ispt    <- uC_ispt
    U_ispt    <- uC_ispt/(SB_ispt + uC_ispt)

    # Load ref pts
    B0_isp     <- array(NA,dim = c(nReps,nS,nP))
    Bmsy_isp   <- array(NA,dim = c(nReps,nS,nP))
    USR_isp    <- array(NA,dim = c(nReps,nS,nP))
    Bcrash_isp <- array(NA,dim = c(nReps,nS,nP))
    Ucrash_isp <- array(NA,dim = c(nReps,nS,nP))
    MSY_isp    <- array(NA,dim = c(nReps,nS,nP))
    # Umsy_isp  <- array(0,dim = c(nReps,nS,nP))
    # MSY_isp   <- array(0,dim = c(nReps,nS,nP))

    # Pull out ref pts, use ms3Rstats functions



    goodRepIdx <- which(goodReps)
    for( i in 1:nReps)
    {
      # j <- goodRepIdx[i]
      # Pull out Dstar and Dscalar
      tdxUSR        <- blob$ctlList$opMod$histCtl$hypo$Dstar - 1950
      scalarUSR     <- blob$ctlList$opMod$histCtl$hypo$Dscalar
      #
      # message("i = ", i, ", j = ", j)
      USR_isp[i,,]  <- apply(X = scalarUSR * SB_ispt[i,,,tdxUSR,drop = FALSE], FUN = mean, MARGIN = 2:3 )
      #
    }

    # Busr_q <- quantile(USR_isp, probs = c(0.25, 0.5, 0.75), na.rm = T)
    # eqListQuant <- solveSimEqbriaQuant( empRefCurves = empRefCurves,
    #                                     maxXspline = maxXspline,
    #                                     yrRange = yrRange,
    #                                     USR_q = Busr_q, reps = 1:nReps,
    #                                     qProbs = c(0.25, 0.5, 0.75) )


    # B0_isp[,1,1]     <- eqListQuant$B0_q[2]
    # Bcrash_isp[,1,1] <- eqListQuant$Bcrash_q[2]
    # Ucrash_isp[,1,1] <- eqListQuant$Ucrash_q[2]
    # MSY_isp[,1,1]    <- eqListQuant$MSY_q[2]
    # Bmsy_isp[,1,1]   <- eqListQuant$Bmsy_q[2]

    LRP_isp <- 0.3*B0_isp
    # USR_isp <- 0.8*Bmsy_isp

    # message("USR = ", Busr_q[2], ", B0 = ", B0_isp[1,1,1])

    sbylab <- "Spawning Biomass (kt)"

    if(dep)
    {
      sbylab <- expression(SB[t]/SB[0])
      for( t in 1:nTs[simIdx] )
      {
        SB_ispt[,,,t] <- SB_ispt[,,,t] / B0_isp
      }
      LRP_isp <- LRP_isp/B0_isp
      Bcrash_isp <- Bcrash_isp/B0_isp
      Bmsy_isp <-  Bmsy_isp/B0_isp
      USR_isp <- USR_isp/B0_isp
      B0_isp  <- B0_isp/B0_isp
      # USR_isp <- USR_isp/Bmsy_isp
      # TRP_isp <- TRP_isp/Bmsy_isp
    }

    B0 <- median(B0_isp)
    LRP <- median(LRP_isp)
    Bcrash <- median(Bcrash_isp)
    Ucrash <- median(Ucrash_isp)
    Bmsy <- median(Bmsy_isp)
    USR <- median(USR_isp)
    MSY <- median(MSY_isp)
    # TRP <- mean(TRP_isp)

    SB_qspt   <- apply(X = SB_ispt, FUN = quantile, MARGIN = 2:4, na.rm = T,
                        probs = c(0.025, 0.5, .975 ) )

    C_qspt   <- apply(X = C_ispt, FUN = quantile, MARGIN = 2:4, na.rm = T,
                        probs = c(0.025, 0.5, .975 ) )

    U_qspt   <- apply(X = U_ispt, FUN = quantile, MARGIN = 2:4, na.rm = T,
                        probs = c(0.025, 0.5, .975 ) )


    plot( x = xLim, y = yLimB,
          type = "n", axes = FALSE )
      mfg <- par("mfg")
      mtext(side = 3, text = plotLabs[simIdx], font = 2, cex = 1.5, line = 1)

      if(mfg[2] == 1)
      {
        axis(side = 2, las = 1, cex.axis = 1.5 )
        mtext( side = 2, text = sbylab, line = 4 , cex = 1.5)
      }
      # if(mfg[2] == mfg[4]){
      #   axis(side = 4, las = 1)
      #   legend("topright", bty = "n", lty = c(2,2,5,2,5), col=c("grey30","darkgreen","purple","orange","red"),
      #          c("B0","USR","Bmsy","LRP","Bcrash"), lwd = c(2,2,2,2,2))
      # }
      grid()
      box()

      # polygon( x = c(yrs,rev(yrs)),
      #          y = c(SB_qspt[1,1,1,],rev(SB_qspt[3,1,1,])),
      #          col = "grey60", border = NA )
      lines( x = yrs, y = SB_qspt[2,1,1,], lwd = 3 )

      abline(v = yrs[tMP] - 0.5, lty = 2, lwd =.9 )
      abline(v = 2052, lty = 2, lwd =1 )

      if(mfg[2] == mfg[4]){
        axis(side = 4, las = 1, cex.axis = 1.5 )
        # legend("topright", bty = "n", lty = c(2,2,5,2,5), col=c("grey30","darkgreen","purple","orange","red"),
        #        legend=c(bquote(SB[0]),"USR",bquote(B[MSY]),"LRP",bquote(B[crash])), lwd = c(1,2,2,2,2), cex=1.2)
      }
      abline(h = B0, lty = 2, col = "grey30" )
      # abline(h = TRP, lty = 2, col = "darkgreen" )
      abline(h = USR, lty = 2, col = "darkgreen", lwd=2)
      abline(h = Bcrash, lty = 5, col = "red", lwd=2)
      abline(h = Bmsy, lty = 5, col = "purple", lwd=2)
      abline(h = LRP, lty = 2, col = "orange", lwd=2 )
      if(traces > 0)
        for( t in traceIdx )
          lines( x = yrs, y = SB_ispt[t,1,1,], lwd = 1, col = "black")
      # abline(h = 2.6, lty = 2, col = "red" )

    plot( x = xLim, y = yLimC,
          type = "n", axes = FALSE )
      mfg <- par("mfg")

      if( mfg[2] == 1 )
      {
        axis(side = 2, las = 1, cex.axis = 1.5 )
        mtext( side = 2, text = "Catch (kt)", line = 4 , cex = 1.5)
      }
      if( mfg[2] == mfg[4] ){
        # legend("topright", bty = "n",lty = 5, col="purple", "MSY", lwd=2, cex=1.2)
        axis(side = 4, las = 1, cex.axis = 1.5 )
      }

      grid()
      box()

      # polygon( x = c(yrs,rev(yrs)),
      #          y = c(C_qspt[1,1,1,],rev(C_qspt[3,1,1,])),
      #          col = "grey60", border = NA )
      lines( x = yrs, y = C_qspt[2,1,1,], lwd = 3 )
      # abline( h = mean(MSY_isp[,1,1]), lty = 2, col = "darkgreen")
      abline(v = yrs[tMP] - 0.5, lty = 2, lwd = .9)
      abline(h = MSY, col = "purple", lty = 5, lwd=2)

      if(traces > 0)
        for( t in traceIdx )
          lines( x = yrs, y = C_ispt[t,1,1,], lwd = 1, col = "black")

    plot( x = xLim, y = yLimHR,
          type = "n", axes = FALSE )
      mfg <- par("mfg")

      if( mfg[2] == 1 )
      {
        axis(side = 2, las = 1, cex.axis = 1.5 )
        mtext( side = 2, text = "Harvest Rate", line = 4 , cex = 1.5)
      }

      if( mfg[2] == mfg[4] ){
        axis(side = 4, las = 1, cex.axis = 1.5 )
        # legend("topright", lty = 5,bty = "n", col="red", legend=bquote(U[crash]), lwd=2, cex=1.2)
      }
      axis(side = 1 , cex.axis = 1.5 )


      grid()
      box()

      # polygon( x = c(yrs,rev(yrs)),
      #          y = c(U_qspt[1,1,1,],rev(U_qspt[3,1,1,])),
      #          col = "grey60", border = NA )
      lines( x = yrs, y = U_qspt[2,1,1,], lwd = 3 )
      if(simLabs[simIdx] != "NoFish"){
        abline( h = blobList[[simIdx]]$ctlList$mp$hcr$Uref_p[1], lty = 2, col = "grey30")
        abline(h = Ucrash, lty = 5, col = "red", lwd=2)
      }
      abline(v = yrs[tMP] - 0.5, lty = 2, lwd =.9 )

      if(traces > 0)
        for( t in traceIdx )
          lines( x = yrs, y = U_ispt[t,1,1,], lwd = 1, col = "black")

  }

  mtext(side = 1, text=  "Year", outer = TRUE, line = 2, cex = 1.5)
}



# genSimInfo()
genSimInfo <- function(simFolder = here('mserproject'))
{

  # Read in sims
  sims <- list.files(file.path(simFolder))
  sims <- sims[grepl("sim",sims)]

  readInfoFile <- function( sim )
  {
    infoPath <- file.path(simFolder,sim,'infoFile.txt' )
    info <- lisread(infoPath)
    info.df <- as.data.frame(info)
    info.df$simLabel <- sim

    info.df
  }

  # Read in info files, sort by  scenarios
  info.df <- lapply( X = sims, FUN = readInfoFile )
  info.df <- do.call( "rbind", info.df )

}  #END genSimInfo()


plotMtBt <- function( obj = blob,
                      iRep = 1 )
{

  # Check if depM
  opMod     <- obj$ctlList$opMod
  repObj    <- obj$ctlList$opMod$histRpt
  depM      <- repObj$densityDepM
  juveMage  <- repObj$juveMage

  # Dimensions
  nS <- obj$om$nS
  nP <- obj$om$nP
  nT <- obj$om$nT

  M_spt <- array(0,dim = c(nS,nP,nT))
  B_spt <- array(0,dim = c(nS,nP,nT))

  # Get Bt and Mt series
  M_spt[1:nS,1:nP,] <- obj$om$M_iaxspt[iRep,juveMage+1,1,,,]
  B_spt[1:nS,,]     <- obj$om$B_ispt[iRep,,,]

  # Pull m1 and Mb values
  Mb_sp             <- obj$rp[[iRep]]$Mb_sp
  m1_sp             <- obj$rp[[iRep]]$m1_sp

  # Pull totB0
  totB0_sp <- obj$rp[[iRep]]$totB0_sp

  Bseq    <- seq(from = 1, to = 1.2*max(B_spt,na.rm = T), length.out = 1000)
  Mseq_sp <- array(0,dim = c(nS,nP,1000))

  M0_sp   <- array(0,dim = c(nS,nP))

  for(s in 1:nS )
    for(p in 1:nP )
    {
      Mseq_sp[s,p,] <- Mb_sp[s,p] + exp(-m1_sp[s,p] * Bseq/totB0_sp[s,p])
      M0_sp[s,p]    <- Mb_sp[s,p] + exp(-m1_sp[s,p])
    }

  par(mfrow = c(nS,nP),
      oma   = c(3,3,1,1),
      mar   = c(1,1,1,1) )
  for( s in 1:nS )
    for( p in 1:nP )
    {
      plot( x = c(0,1.2*max(B_spt,na.rm = T)),
            y = c(0,1.2*max(M_spt,na.rm = T)),
            axes = FALSE, type = "n", xaxs = "i", yaxs = "i")
      axis( side = 1)
      axis( side = 2, las = 1)
      grid()
      box()
      lines(x  = Bseq, y = Mseq_sp[s,p,], col = "salmon", lwd = 3)
      points( x = B_spt[s,p,], y = M_spt[s,p,], pch = 16 )
      segments( x0 = totB0_sp[s,p], y0 = 0, y1 = M0_sp[s,p], lty = 2)
      segments( x0 = 0, x1 = totB0_sp[s,p], y0 = M0_sp[s,p], lty = 2)
    }

  mtext( side = 1, outer = T, text = "Age 2+ Biomass (kt)", line = 2)
  mtext( side = 2, outer = T, text = "Natural Mortalty (/yr)", line = 2)

}



# multiSimBtIt()
# Plot BtIt for same rep across different mp/scenario combinations
multiSimBtIt <- function( iRep = 1,
                           simFolder,
                           OMs = c("rw5_loPondM, rw15_hiPondM"),
                           MPs = c("NoFish"),
                           yLimB     = c(0,22))


{


  info.df <- genSimInfo(simFolder=simFolder) %>%
              mutate_if( is.factor, as.character) %>%
              filter( scenario %in% OMs,
                      mp %in% MPs )

  sims    <- expand.grid(OMs, MPs)
  names(sims) <- c('OM','MP')


  # Set up plot environment
  par(mfcol = c(3,nrow(info.df) ), mar = c(1.5,1.5,0.5,1),
      oma = c( 3,4,1,0))

  # Loop over MPs
  for( i in 1:nrow(info.df) )
  {
    subDF <-  info.df %>%
              filter( mp == sims$MP[i] &
                      scenario == sims$OM[i])

    simID <- subDF[1,"simLabel"]

    simFile <- paste(simID,".RData",sep = "")
    simPath <- file.path(simFolder,simID,simFile)

    # Load blob
    load(simPath)

    # Plot BtIt
    plotBtIt_p(obj = blob, iRep = 1, f=5, addCatch=TRUE,
               parArg=FALSE, YlabOn=FALSE, legdOn=FALSE)

    mtext(side = 3, text = paste(sims$OM[i],sims$MP[i]) , line = 33, cex = 1)

  }

  mtext( side =2, outer = TRUE, text = "Spawning biomass kt)", line = 1)



} # END multiSimBtCtRt()

plotMultiHCR <- function(LCPs, UCPs, hcrNames, xlab)
{
  # number of plots
  nI <- length(LCPs)

  par(mfrow=c(1,nI), mgp=c(1.4,0.6,0), mar=c(2,3,1,1), oma = c(2,1,0,0))
  for(i in 1:nI)
    plotHockeyStickHCR(LCP = LCPs[i], UCP = UCPs[i],
                        hcrName=hcrNames[i])

  mtext( side = 1, text = xlab, outer = TRUE)

}

# plotHockeyStickHCR
plotHockeyStickHCR <- function( LCP = .5, UCP = .6, refHR = .1,
                                refHRaxis = NULL, hcrName='',
                                refHRlab='', mpLab=NULL,
                                refB = 1, yAXT='n',
                                yLim = NULL,
                                xLab='', yLab='Harvest Rate')
                                # xLab = expression(paste("Stock Status (",B[t]/B[MSY],")",sep = "")),
                                # yLab = "Target Fishing Mortality (F)" )
{
  if(is.null(yLim))
    yLim <- c(0,1.5*refHR)


  plot( x = c(0,refB*1.2), y = yLim, type = "n", xlab = "", ylab = "", las = 1,
        yaxt=yAXT, main=hcrName )
    segments( x0 = 0, x1 = LCP * refB, y0 = 0, lwd = 1.5 )
    segments( x0 = LCP*refB, x1 = UCP*refB, y0 = 0, y1 = refHR, lwd = 1.5 )
    segments( x0 = UCP*refB, x1 = refB*1.2, y0 = refHR, y1 = refHR, lwd = 1.5 )
    abline( v = c(UCP*refB, LCP * refB), lwd = .8, lty = 2)
    abline( h = refHR, lty = 3, lwd = .8)
    mtext( side = 1, text = xLab, line = 2 )
    mtext( side = 2, text = yLab, line = 2 )

    if(is.null(refHRaxis))
      refHRaxis = refHR

    if(!is.null(mpLab))
      mtext( side = 4, text = mpLab, line = 0.5, cex=0.7 )


    axis(side=2, at=c(0,refHR), tick=TRUE, labels=c('0',refHRaxis), las=1)
    text(x=0.2,y=refHR+0.004, refHRlab)


} # END plotHockeyStickHCR()

# plotTACallocationError()
# A multi-panel for understanding the
# level of error in TAC allocations under
# aggregate MPs. Not sure of the best way
# to look at this...
plotTACallocationError <- function( obj = blob,
                                    iRep = 1 )
{
  goodReps <- obj$goodReps

  # Get biomass arrays
  SB_ispt        <- obj$om$SB_ispt[goodReps,,,]
  VB_ispt        <- obj$om$vB_ispft[goodReps,,,2,]
  totB_ispt      <- obj$om$B_ispt[goodReps,,,]
  retroSB_itspt  <- obj$mp$assess$retroSB_itspt[goodReps,,,,]

  ctlList <- obj$ctlList

  retroSB_itspt[retroSB_itspt < 0] <- NA
  nReps     <- sum(goodReps)

  # Model dims
  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP
  nT      <- obj$om$nT



}


makeAMREs <- function(  obj = blob,
                        sIdx = 1, pIdx = 1,
                        pt = 1 )
{
  # Pull model dims
  nS    <- obj$om$nS
  nP    <- obj$om$nP
  nX    <- obj$om$nX
  nF    <- obj$om$nF
  tMP   <- obj$om$tMP
  nT    <- obj$om$nT

  # Get goodreps
  goodReps <- which(obj$goodReps)
  nReps    <- length(goodReps)

  repObj  <- obj$ctlList$opMod$histRpt

  # Biological pars
  B0_i      <- rep(0,nReps)
  R0_i      <- rep(0,nReps)
  totB0_i   <- rep(0,nReps)
  h_i       <- rep(0,nReps)
  M_i       <- array(0,dim = c(nReps))
  M0_i      <- array(0,dim = c(nReps))
  m1_i      <- rep(0,nReps)
  M_it      <- array(0,dim = c(nReps,tMP - 1))

  # Nuisance pars
  q_if      <- array(0,dim = c(nReps,nF))
  tauObs_if <- array(0,dim = c(nReps,nF))

  # Selectivity parameters
  SelA_ift    <- array(0,dim = c(nReps,nF,tMP - 1))
  SelB_ift    <- array(0,dim = c(nReps,nF,tMP - 1))

  estSelA_ift    <- array(0,dim = c(nReps,nF,tMP - 1))
  estSelB_ift    <- array(0,dim = c(nReps,nF,tMP - 1))



  # Pull true values
  for( idx in 1:length(goodReps) )
  {
    i <- goodReps[idx]

    B0_i[idx]       <- obj$rp[[i]]$B0_sp[sIdx,pIdx,1]
    R0_i[idx]       <- obj$rp[[i]]$R0_sp[sIdx,pIdx]
    totB0_i[idx]    <- obj$rp[[i]]$totB0_sp[sIdx,pIdx]
    h_i[idx]        <- repObj$rSteepness_p[pIdx]
    M_i[idx]        <- obj$rp[[i]]$Mb_sp[sIdx,pIdx]
    m1_i[idx]       <- obj$rp[[i]]$m1_sp[sIdx,pIdx]
    M0_i[idx]       <- M_i[idx] + exp(-m1_i[idx])
    q_if[idx,]      <- repObj$qComb_pg[pIdx,]
    tauObs_if[idx,] <- repObj$tauComb_pg[pIdx,]

    SelA_ift[idx,,1:(tMP-1)]   <- repObj$SelAlpha_pgt[1,,]
    SelB_ift[idx,,1:(tMP-1)]   <- SelA_ift[idx,,1:(tMP-1)] + repObj$SelBeta_pgt[1,,]

  }

  # And now adjust the estimates as well
  estSelA_ift[1:nReps,,]   <- obj$mp$assess$retroSelA_itspft[goodReps,pt,sIdx,pIdx,,1:(tMP-1)]
  estSelB_ift[1:nReps,,]   <- obj$mp$assess$retroSelA_itspft[goodReps,pt,sIdx,pIdx,,1:(tMP-1)] + obj$mp$assess$retroSelB_itspft[goodReps,pt,sIdx,pIdx,,1:(tMP-1)]


  if(!is.null(obj$mp$assess$retroM0_itsp))
    errM0_i     <- -1*(M0_i - obj$mp$assess$retroM0_itsp[goodReps,pt,sIdx,pIdx])/M0_i

  errB0_i     <- -1*(B0_i - obj$mp$assess$retroB0_itsp[goodReps,pt,sIdx,pIdx])/B0_i
  errtotB0_i  <- -1*(totB0_i - obj$mp$assess$retrototB0_itsp[goodReps,pt,sIdx,pIdx])/totB0_i
  errR0_i     <- -1*(R0_i - obj$mp$assess$retroR0_itsp[goodReps,pt,sIdx,pIdx])/R0_i
  errh_i      <- -1*(h_i - obj$mp$assess$retroh_itsp[goodReps,pt,sIdx,pIdx])/h_i
  errM_i      <- -1*(M_i - obj$mp$assess$retroM_itsp[goodReps,pt,sIdx,pIdx])/M_i
  errm1_i     <- -1*(m1_i - obj$mp$assess$retrom1_itsp[goodReps,pt,sIdx,pIdx])/m1_i
  errq_if     <- -1*(q_if - obj$mp$assess$retroq_itspf[goodReps,pt,sIdx,pIdx,])/q_if
  errtau_if   <- -1*(tauObs_if - obj$mp$assess$retrotauObs_itspf[goodReps,pt,sIdx,pIdx,])/tauObs_if

  if(is.null(obj$mp$assess$retroM0_itsp))
  {
    retroM0_itsp <- obj$mp$assess$retroM_itsp + exp(-1 * obj$mp$assess$retrom1_itsp)
    errM0_i       <- -1*(M0_i - retroM0_itsp[goodReps,pt,sIdx,pIdx])/M0_i
  }

  # Get errors for historical period only - will update to
  # future once we have tv sel in the projections. Right now,
  # we're good as all tvSel happens in the history and the
  # historical estimates for the newRV will be captured here.
  errSelA_ift   <- -1*(SelA_ift - estSelA_ift)/SelA_ift
  errSelB_ift   <- -1*(SelB_ift - estSelB_ift)/SelB_ift

  errB0_q   <- quantile(errB0_i, probs = c(0.025,0.5,0.975),na.rm = T)
  errtotB0_q<- quantile(errtotB0_i, probs = c(0.025,0.5,0.975),na.rm = T)
  errR0_q   <- quantile(errR0_i, probs = c(0.025,0.5,0.975),na.rm = T)
  errh_q    <- quantile(errh_i, probs = c(0.025,0.5,0.975), na.rm = T)
  errM_q    <- quantile(errM_i, probs =c(0.025,0.5,0.975), na.rm = T)
  errM0_q   <- quantile(errM0_i, probs =c(0.025,0.5,0.975), na.rm = T)
  errm1_q   <- quantile(errm1_i, probs =c(0.025,0.5,0.975), na.rm = T)
  errq_qf   <- apply(X = errq_if, FUN = quantile, probs =c(0.025,0.5,0.975), MARGIN = 2, na.rm = T)
  errtau_qf <- apply(X = errtau_if, FUN = quantile, probs =c(0.025,0.5,0.975), MARGIN = 2, na.rm = T)

  errSelA_qft <- apply(X = errSelA_ift, FUN = quantile, probs =c(0.025,0.5,0.975), MARGIN = 2:3, na.rm = T)
  errSelB_qft <- apply(X = errSelB_ift, FUN = quantile, probs =c(0.025,0.5,0.975), MARGIN = 2:3, na.rm = T)

  errList <- list(  errB0_i       = errB0_i,
                    errtotB0_i    = errtotB0_i,
                    errR0_i       = errR0_i,
                    errM0_i       = errM0_i,
                    errh_i        = errh_i,
                    errM_i        = errM_i,
                    errm1_i       = errm1_i,
                    errq_if       = errq_if,
                    errtau_if     = errtau_if,
                    errSelA_ift   = errSelA_ift,
                    errSelB_ift   = errSelB_ift )

  quantList <- list(  errB0_q       = errB0_q,
                      errtotB0_q    = errtotB0_q,
                      errR0_q       = errR0_q,
                      errM0_q       = errM0_q,
                      errh_q        = errh_q,
                      errM_q        = errM_q,
                      errm1_q       = errm1_q,
                      errq_qf       = errq_qf,
                      errtau_qf     = errtau_qf,
                      errSelA_qft   = errSelA_qft,
                      errSelB_qft   = errSelB_qft )

  outList <- list(  errList = errList,
                    quantList = quantList )

  outList
} # END makeAMREs


# plotAMREs()
plotAMREs <- function(  obj = blob,
                        sIdx = 1, pIdx = 1,
                        pt = 1, plotFleets = c(4:5),
                        fleetNames = c("Surface","Dive") )
{
  # Pull model dims
  nReps <- obj$nSims
  nS    <- obj$om$nS
  nP    <- obj$om$nP
  nX    <- obj$om$nX
  nF    <- obj$om$nF
  tMP   <- obj$om$tMP
  nT    <- obj$om$nT

  repObj  <- obj$ctlList$opMod$histRpt


  quantErrors <- makeAMREs( obj = obj, sIdx = sIdx,
                            pIdx = pIdx, pt = pt)$quantList

  errB0_q       <- quantErrors$errB0_q
  errtotB0_q    <- quantErrors$errtotB0_q
  errR0_q       <- quantErrors$errR0_q
  errM0_q       <- quantErrors$errM0_q
  errh_q        <- quantErrors$errh_q
  errM_q        <- quantErrors$errM_q
  errm1_q       <- quantErrors$errm1_q
  errq_qf       <- quantErrors$errq_qf
  errtau_qf     <- quantErrors$errtau_qf
  errSelA_qft   <- quantErrors$errSelA_qft
  errSelB_qft   <- quantErrors$errSelB_qft

  labels <- c("B0","totB0","R0","h","Mb","m1","M0")

  nErr <- length(labels) + 2*length(plotFleets)



  for(f in 1:length(plotFleets))
    labels <- c(labels,paste0("q_",fleetNames[f]),paste0("tau_",fleetNames[f]))

  plot( x = c(1,nErr), xlab = "", ylab = "",
        y = c(-1,1), axes = FALSE, type = "n" )
    axis( side = 1, at = 1:nErr, labels = labels )
    axis(side = 2, las =1)
    grid()
    box()

    # B0 error
    segments(x0 = 1, y0 = errB0_q[1], y1 = errB0_q[3], lwd = 3, col = "grey60" )
    points(x = 1, y = errB0_q[2], pch = 16, col = "grey60" )

    # totB0 error
    segments(x0 = 2, y0 = errtotB0_q[1], y1 = errtotB0_q[3], lwd = 3, col = "grey60" )
    points(x = 2, y = errtotB0_q[2], pch = 16, col = "grey60" )


    # R0 error
    segments(x0 = 3, y0 = errR0_q[1], y1 = errR0_q[3], lwd = 3, col = "grey60" )
    points(x = 3, y = errR0_q[2], pch = 16, col = "grey60" )


    # h error
    segments(x0 = 4, y0 = errh_q[1], y1 = errh_q[3], lwd = 3, col = "grey60" )
    points(x = 4, y = errh_q[2], pch = 16, col = "grey60" )

    # Mb error
    segments(x0 = 5, y0 = errM_q[1], y1 = errM_q[3], lwd = 3, col = "grey60" )
    points(x = 5, y = errM_q[2], pch = 16, col = "grey60" )

    # m1 error
    segments(x0 = 6, y0 = errm1_q[1], y1 = errm1_q[3], lwd = 3, col = "grey60" )
    points(x = 6, y = errm1_q[2], pch = 16, col = "grey60" )

    # M0 error
    segments(x0 = 7, y0 = errM0_q[1], y1 = errM0_q[3], lwd = 3, col = "grey60" )
    points(x = 7, y = errM0_q[2], pch = 16, col = "grey60" )

    nPlotFleets <- length(plotFleets)
    # Catchability error
    for(j in 1:length(plotFleets))
    {
      fIdx <- plotFleets[j]
      segments(x0 = 8 + (j-1) * 2, y0 = errq_qf[1,fIdx], y1 = errq_qf[3,fIdx], lwd = 3, col = "grey60" )
      points(x = 8 + (j-1) * 2, y = errq_qf[2,fIdx], pch = 16, col = "grey60" )

      segments(x0 = 9 + (j-1) * 2, y0 = errtau_qf[1,fIdx], y1 = errtau_qf[3,fIdx], lwd = 3, col = "grey60" )
      points(x = 9 + (j-1) * 2, y = errtau_qf[2,fIdx], pch = 16, col = "grey60" )
    }


    mtext( side = 2, text = "Relative Error in AM estimates", line = 3)
}


# plotAMREsel()
# Relative error distributions for selectivity
# parameters at a given projection time-step
plotAMREsel <- function(  obj = blob,
                          sIdx = 1, pIdx = 1,
                          pt = 1 )
{
  # Pull model dims
  nReps <- obj$nSims
  nS    <- obj$om$nS
  nP    <- obj$om$nP
  nX    <- obj$om$nX
  nF    <- obj$om$nF
  tMP   <- obj$om$tMP
  nT    <- obj$om$nT

  repObj  <- obj$ctlList$opMod$histRpt

  # Get fleet names
  fleetNames_f <- dimnames(repObj$C_pgt)[[2]]

  # Make error distributions
  quantErrors <- makeAMREs( obj = obj, sIdx = sIdx,
                            pIdx = pIdx, pt = pt)$quantList

  errList <- list()
  errList$errSelA_qft   <- quantErrors$errSelA_qft
  errList$errSelB_qft   <- quantErrors$errSelB_qft

  parNames <- c("L50A","L95A")

  fleetCols <- RColorBrewer::brewer.pal(nF,"Paired")

  if(pt > 4)
    nFleets <- nF
  else nFleets <- nF-1

  par(mfrow = c(4,1), oma = c(3,5,3,3), mar = c(.1,.1,.1,.1))
  for( j in 1:4)
  {
    # maxY <- max(abs(errList[[j]]),na.rm =T)
    plot( x = c(0,nFleets+1), xlab = "", ylab = "",
          y = c(-1,1), axes = FALSE, type = "n" )
      mfg <- par("mfg")
      if( mfg[1] == mfg[3])
        axis( side = 1, at = 1:(nFleets), labels = fleetNames_f[1:(nFleets)] )

      abline(h = 0, lty = 2, lwd = 1)

      axis(side = 2, las =1)
      legend(x = "topright", legend = parNames[j], bty= "n")
      grid()
      box()

      # Now, loop over fleets and plot the quantiles
      for( f in 1:(nFleets))
      {
        if( (j %in% c(3,4)) & (!f %in% c(1,7,8,9)) )
          next

        segments( x0 = f, y0 = errList[[j]][1,f,1], y1 = errList[[j]][3,f,1], col = fleetCols[f],
                  lwd = 3 )
        points( x = f, y = errList[[j]][2,f,1], col = fleetCols[f],
                pch = 16, cex = 2 )
      }
  }


  mtext( side = 2, text = "Relative Error in AM estimates", line = 3, outer = TRUE)
} # END plotAMREsel()


# plotAMREs()
plotMultiAMREs <- function( blobList = list(blob1 = blob),
                            sIdx = 1, pIdx = 1,
                            pt = 1, plotFleets = c(7,8),
                            fleetNames = c("RV","SFW","newRV") )
{
  obj <- blobList[[1]]

  # Pull model dims
  nReps <- obj$nSims
  nS    <- obj$om$nS
  nP    <- obj$om$nP
  nX    <- obj$om$nX
  nF    <- obj$om$nF

  repObj <- obj$ctlList$opMod$histRpt

  nPlotFleets <- length(plotFleets)
  nErr        <- 2+nX+2*nPlotFleets
  nBlobs      <- length(blobList)


  xJitter   <- seq( from = -.3, to = .3, length.out = nBlobs)
  blobCols  <- RColorBrewer::brewer.pal(nBlobs,"Dark2")



  labels <- c("B0","h","M_m","M_f")

  for(f in 1:length(plotFleets))
    labels <- c(labels,paste0("q_",fleetNames[f]),paste0("tau_",fleetNames[f]))


  plot( x = c(0.5,nErr+.5), xlab = "", ylab = "",
        y = c(-1,1), axes = FALSE, type = "n",
        xaxs = "i")
    axis( side = 1, at = 1:nErr, labels = labels )
    axis(side = 2, las =1)
    abline(v = 0.5 + 1:(nErr-1))
    abline(h = 0, lty = 2, col = "grey60")
    box()

    for( k in 1:nBlobs)
    {
      errList <- makeAMREs( blobList[[k]],
                            sIdx = sIdx, pIdx = pIdx,
                            pt = pt )$quantList
      # B0 error
      segments(x0 = xJitter[k] + 1, y0 = errList$errB0_q[1], y1 = errList$errB0_q[3], lwd = 3, col = blobCols[k] )
      points(x = xJitter[k] + 1, y = errList$errB0_q[2], pch = 16, col = blobCols[k] )

      # h error
      segments(x0 = xJitter[k] + 2, y0 = errList$errh_q[1], y1 = errList$errh_q[3], lwd = 3, col = blobCols[k] )
      points(x = xJitter[k] + 2, y = errList$errh_q[2], pch = 16, col = blobCols[k] )

      # M_m error
      segments(x0 = xJitter[k] + 3, y0 = errList$errM_qx[1,1], y1 = errList$errM_qx[3,1], lwd = 3, col = blobCols[k] )
      points(x = xJitter[k] + 3, y = errList$errM_qx[2,1], pch = 16, col = blobCols[k] )

      # M_f error
      segments(x0 = xJitter[k] + 4, y0 = errList$errM_qx[1,2], y1 = errList$errM_qx[3,2], lwd = 3, col = blobCols[k] )
      points(x = xJitter[k] + 4, y = errList$errM_qx[2,2], pch = 16, col = blobCols[k] )

      # Catchability error
      for(j in 1:nPlotFleets)
      {
        fIdx <- plotFleets[j]
        segments(x0 = xJitter[k] + 5 + (j-1)*2, y0 = errList$errq_qf[1,fIdx], y1 = errList$errq_qf[3,fIdx], lwd = 3, col = blobCols[k] )
        points(x = xJitter[k] + 5 + (j-1)*2, y = errList$errq_qf[2,fIdx], pch = 16, col = blobCols[k] )

        segments(x0 = xJitter[k] + 6 + (j-1)*2, y0 = errList$errtau_qf[1,fIdx], y1 = errList$errtau_qf[3,fIdx], lwd = 3, col = blobCols[k] )
        points(x = xJitter[k] + 6 + (j-1)*2, y = errList$errtau_qf[2,fIdx], pch = 16, col = blobCols[k] )
      }

    }
    legend( x = "topleft", bg ="white",
            legend = c( "Sim Lengths",
                        "Sim Index",
                        "Sim All Data"),
            lwd = 3, col = blobCols,
            pch = 16, pt.cex = 2 )

    mtext( side = 2, text = "Relative Error in AM estimates", line = 3)
}


# plotMultiTulipAssError
# A wrapper for plotTulipAssError that creates a grid
# over OMs and MPs
plotMultiTulipAssError <- function( groupFolder = "Outputs/tunedHRs" )
{

  # First, we'll go over the sims in the groupFolder
  # and arrange by scenario/MP. Use the same
  # routine we used to make the weighted blobs

  # Now read the infoFile in each sim folder
  dirList   <- list.dirs(groupFolder)
  dirList   <- dirList[grepl(x = dirList, pattern = "sim_")]

  infoList  <- file.path(dirList,"infoFile.txt" )
  infoList  <- lapply(X = infoList, FUN = lisread)
  infoList  <- lapply(X = infoList, FUN = as.data.frame)
  info.df   <- do.call(rbind, infoList)

  # Now pull MPs and scenarios
  scenarios <- unique(info.df$scenario)
  mps       <- unique(info.df$mp)
  mps       <- mps[order(mps)]

  nScen <- length(scenarios)
  nMP   <- length(mps)

  par(  mfrow = c(nScen, nMP),
        oma   = c(4,5,3,4),
        mar   = c(.1,.1,.1,.1))

  for(sIdx in 1:nScen)
    for(mIdx in 1:nMP )
    {

      blobRow <- info.df |>
                  filter( scenario == scenarios[sIdx],
                          mp == mps[mIdx] )

      folderName <- blobRow$simLabel


      load(file.path(groupFolder,folderName,paste0(folderName,".Rdata")))

      plotTulipAssError(  obj = blob,
                          sIdx = 1,
                          pIdx = 1, noPar = TRUE )

      if(sIdx == 1)
        mtext(side = 3, text = mps[mIdx], font = 2)

      mfg <- par("mfg")
      if(mfg[1] == mfg[3])
        axis(side = 1)
      if(mfg[2] == 1 )
        axis( side = 2, las = 1)
      if(mfg[2] == mfg[4])
        axis(side = 4, las = 1)

    }
  mtext(side = 1, text = "Year", outer = TRUE, line = 3)
  mtext( side = 2, text = "Relative error in SISCAH SSB forecast", outer = TRUE, line = 3)



} # plotMultiTulipAssError



# Envelopes of simulated assessment errors
plotTulipAssError <- function(  obj = blob, sIdx = 1, pIdx = 1, noPar = FALSE )
{

  goodReps <- obj$goodReps

  # Model dims
  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP
  nT      <- obj$om$nT
  pT      <- nT - tMP + 1

  # Get biomass arrays
  SB_it        <- obj$om$endSB_ispt[goodReps,sIdx,pIdx,1:nT]
  retroSB_itt  <- obj$mp$assess$retroSB_itspt[goodReps,1:pT,sIdx,pIdx,1:nT]

  ctlList <- obj$ctlList

  retroSB_itt[retroSB_itt < 0] <- NA
  nReps     <- sum(goodReps)


  SB_it[SB_it == 0]     <- NA

  assErr_it <- array(NA, dim = c(nReps,pT) )
  for(t in 1:pT)
    assErr_it[,t] <- (retroSB_itt[,t,tMP + t - 1] - SB_it[,tMP + t - 1])/SB_it[,tMP + t - 1]

  assErr_qt <- apply( X = assErr_it, na.rm = T,
                      FUN = quantile,
                      probs = c(0.025, 0.5, 0.975),
                      MARGIN = c(2) )

  years <- seq(from = 1951, by = 1, length.out = nT )

  traces <- sample( 1:nReps, 3)

  # Now plot
  plot( x = range(years[(tMP):nT]),
        y = range(assErr_qt,na.rm = T ),
        type = "n", xlab = "", ylab = "", axes = FALSE )

        mfg <- par("mfg")

        if(!noPar)
        {
          axis( side = 1 )
          axis( side = 2, las = 1 )
        }

        box()
        grid()

        # Plot baseline
        polygon(  x = c(years[tMP:nT],rev(years[tMP:nT])),
                  y = c(assErr_qt[1,],rev(assErr_qt[3,])),
                  border = NA, col = "grey60")
        lines( x = years[tMP:nT], y = assErr_qt[2,],
                col = "black", lwd = 2)
        for( traceIdx in traces )
          lines(  x = years[tMP:nT], y = assErr_it[traceIdx,],
                  col = "black", lwd = .8 )

        abline( h = 0, lty = 2, lwd = 1 )
  if( !noPar )
  {
    mtext(  side = 1, text = "Year", outer = TRUE, line = 2 )
    mtext(  side = 2, text = "Relative assessment error",
            line = 2, outer = TRUE )
  }


} # END plotTulipAssError

# Envelopes of simulated assessment errors
plotTulipAssErrorMort <- function(  obj = blob,
                                    pt = 1, sIdx = 1, pIdx = 1 )
{

  goodReps <- obj$goodReps

  # Model dims
  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP
  nT      <- obj$om$nT


  # Get biomass arrays
  M_it        <- obj$om$M_iaxspt[goodReps,2,1,sIdx,pIdx,1:(tMP-1)]
  retroM_it   <- obj$mp$assess$retroM_itaxspt[goodReps,pt,2,1,sIdx,pIdx,1:(tMP-1)]

  ctlList <- obj$ctlList

  retroM_it[retroM_it < 0] <- NA
  nReps     <- sum(goodReps)


  M_it[M_it == 0]     <- NA

  assErr_it <- array(NA, dim = c(nReps,tMP-1) )
  assErr_it[,1:(tMP-1)] <- (retroM_it[,1:(tMP-1)] - M_it[,1:(tMP-1)])/M_it[,1:(tMP-1)]

  assErr_qt <- apply( X = assErr_it, na.rm = T,
                      FUN = quantile,
                      probs = c(0.025, 0.5, 0.975),
                      MARGIN = c(2) )

  years <- seq(from = 1951, by = 1, length.out = tMP-1 )

  tIdx <- 1:(tMP-1)

  traces <- sample( 1:nReps, 3)

  # Now plot
  plot( x = range(years[tIdx]),
        y = range(assErr_qt[,tIdx],na.rm = T ),
        type = "n", xlab = "", ylab = "", axes = FALSE )
        mfg <- par("mfg")
        axis( side = 1 )
        axis( side = 2, las = 1 )
        box()
        grid()

        # Plot baseline
        polygon(  x = c(years,rev(years)),
                  y = c(assErr_qt[1,],rev(assErr_qt[3,])),
                  border = NA, col = "grey60")
        lines( x = years, y = assErr_qt[2,],
                col = "black", lwd = 2)
        for( traceIdx in traces )
          lines(  x = years, y = assErr_it[traceIdx,],
                  col = "black", lwd = .8 )

        abline( h = 0, lty = 2, lwd = 1 )
  mtext(  side = 1, text = "Year", outer = TRUE, line = 2 )
  mtext(  side = 2, text = "Relative assessment error",
          line = 2, outer = TRUE )


} # END plotTulipAssError

# plotTulipBtCtBaseSim()
# Overlays the biomass can catch tulips
# from the baseline omniscient sim (black) and
# stochastic simulation (red) so we can see
# how the different results are made
plotTulipBtCtBaseSim <- function( sim = 1,
                                  groupFolder = "diffCV_fixedF_longGrid",
                                  var         = "SB_ispt",
                                  save        = FALSE,
                                  fYear       = 1956,
                                  proj        = TRUE )
{
  # Load loss reports
  objLoss <- .loadLoss(sim = sim, groupFolder)

  # Load blob, save reference points
  .loadSim(sim = sim, groupFolder)

  rp <- blob$rp[[1]]

  scenName  <- blob$ctlList$ctl$scenarioName
  mpName    <- blob$ctlList$ctl$mpName

  stamp <- paste(scenName,":",mpName, sep = "")

  EmsyMSRefPts <- rp$EmsyMSRefPts
  FmsyRefPts   <- rp$FmsyRefPts


  simFolder <- objLoss$simFolder

  speciesNames    <- objLoss$speciesNames
  stockNames      <- objLoss$stockNames

  baseState_ispt  <- objLoss$baseStates[[var]]
  simState_ispt   <- objLoss$simStates[[var]]

  B0_isp          <- objLoss$baseStates[["SB_ispt"]][,,,1]


  # Pull model dimensions
  tMP <- objLoss$tMP
  nT  <- objLoss$nT
  nF  <- objLoss$nF
  nS  <- objLoss$nS
  nP  <- objLoss$nP
  pT  <- objLoss$pT

  if( proj )
    tIdx <- (tMP-1):nT
  else tIdx <- 1:nT

  years <- seq( from = fYear, by = 1, length.out = nT )

  nReps <- dim(simState_ispt)[1]

  traces <- sample(1:nReps, size = 3)

  # BeqFmsySS_sp <- array(NA, dim =c(nS+1, nP+1))
  # BeqEmsyMS_sp <- array(NA, dim =c(nS+1, nP+1))

  # BeqFmsySS_sp[1:nS,1:nP] <- FmsyRefPts$BeqFmsy_sp
  # BeqFmsySSS_sp[1:nS,p+1] <- apply( X = FmsyRefPts$BeqFmsy_sp,
  #                                   FUN =  )

  if( save )
  {
    graphics.off()
    filename <- paste("baseSimTulipOverlay_",var,".pdf")
    outFile <- file.path(simFolder,filename)

    pdf( outFile, width = 11, height = 8 )
  }

  # put on depletion scale
  if( var == "SB_ispt")
  {
    maxY_sp <- matrix(1, nrow = nS + 1, ncol = nP +1 )
    yLab = expression(paste("Spawning biomass relative to unfished (", B[t]/B[0], ")",sep = ""))
    for( p in 1:(nP+1) )
    {

      for( s in 1:(nS+1) )
        for( i in 1:nReps)
        {
          baseState_ispt[i,s,p,]  <- baseState_ispt[i,s,p,] / B0_isp[i,s,p]
          simState_ispt[i,s,p,]   <- simState_ispt[i,s,p,] / B0_isp[i,s,p]
        }

      maxY_sp[,p] <- max( baseState_ispt[,,p,tIdx], simState_ispt[,,p,tIdx])
    }


  } else {
    yLab    <- "Catch (kt)"
    maxY_sp <- matrix(NA, nrow = nS+1, ncol = nP+1 )

    for( s in 1:(nS+1))
      for( p in 1:(nP+1))
        maxY_sp[s,p] <- max( baseState_ispt[,s,p,tIdx], simState_ispt[,s,p,tIdx], na.rm = T )
  }

  baseState_qspt <- apply( X = baseState_ispt, FUN = quantile,
                            probs = c(0.025, 0.5, 0.975 ),
                            MARGIN = c(2,3,4) )

  simState_qspt <- apply( X = simState_ispt, FUN = quantile,
                            probs = c(0.025, 0.5, 0.975 ),
                            MARGIN = c(2,3,4) )

  par( mfcol = c((nP+1), nS+1),
        mar = c(0.1,0.1,0.1,0.1),
        oma = c(4,4.5,3,3) )

  if( var == "C_ispt" )
    par( mar = c(.1, 1.5, .1, 1.5) )


  basePolyCol <- "grey50"
  baseLineCol <- "black"
  simPolyCol <- scales::alpha("red",.3)
  simLineCol <- "red"

  for( s in 1:(nS+1) )
    for( p in 1:(nP+1) )
    {
      plot( x = range(years[tIdx]),
            y = c(0, maxY_sp[s,p] ),
            type = "n", xlab = "", ylab = "", axes = FALSE )
        mfg <- par("mfg")
        if( mfg[1] == mfg[3] )
          axis( side = 1 )
        if( mfg[2] == 1 )
          axis( side = 2, las = 1 )
        if( var == "C_ispt" & mfg[2] != 1 )
          axis( side = 2, las = 1)
        if( mfg[2] == mfg[4] )
          mtext( side = 4, text = stockNames[p], line = 0.5)
        if( mfg[1] == 1 )
          mtext( side = 3, text = speciesNames[s])
        box()
        grid()

        # Plot baseline
        polygon(  x = c(years,rev(years)),
                  y = c(baseState_qspt[1,s,p,],rev(baseState_qspt[3,s,p,])),
                  border = NA, col = basePolyCol)
        lines( x = years, y = baseState_qspt[2,s,p,],
                col = baseLineCol, lwd = 2)
        for( traceIdx in traces )
          lines(  x = years, y = baseState_ispt[traceIdx,s,p,],
                  col = baseLineCol, lwd = .8 )

        # Plot stochastic sim
        polygon(  x = c(years,rev(years)),
                  y = c(simState_qspt[1,s,p,],rev(simState_qspt[3,s,p,])),
                  border = NA, col = simPolyCol)
        lines( x = years, y = simState_qspt[2,s,p,],
                col = simLineCol, lwd = 2)
        for( traceIdx in traces )
          lines(  x = years, y = simState_ispt[traceIdx,s,p,],
                  col = simLineCol, lwd = .8 )

        # Plot some lines
        abline( v = years[tMP], lty = 2, lwd = 1 )
        if( s <= nS & p <= nP & var == "SB_ispt")
        {
          abline( h = FmsyRefPts$BeqFmsy_sp[s,p]/B0_isp[1,s,p],
                  lty = 3, lwd = 1 )
          abline( h = EmsyMSRefPts$BeqEmsy_sp[s,p]/B0_isp[1,s,p],
                  lty = 4, lwd = 1 )
        }


    }
  mtext( side = 1, text = "Year", outer = TRUE, line = 2 )
  mtext( side = 2, text = yLab, outer = TRUE, line = 2 )

  mtext( side = 1, text = stamp, outer = TRUE, line = 2.5,
          adj = .85, cex = .6, col = "grey75")

  if(save)
    dev.off()

}



# plotBatchConvergenceRate()
# Abandoned performance metric for simulations,
# basically measuring the usefulness of a model
# under lower data quality conditions - can be
# used to set a threshold for including
# results in the paper...
plotBatchConvergenceRate <- function( groupFolder = "diffCV_newObsCV_short",
                                      prefix = "diffCV",
                                      AMlabs = c( singleStock = "singleStock",
                                                  hierMultiStock = "hierMultiStock",
                                                  dataPooled = "dataPooled",
                                                  coastwide = "coastwide",
                                                  totalAgg = "totalAgg" ) )
{
  # First, read info files from the relevant
  # sims
  simFolderList <- list.dirs(here::here("Outputs",groupFolder))
  simFolderList <- simFolderList[grepl(prefix, simFolderList)]


  info.df <-  readBatchInfo( batchDir = here::here("Outputs",groupFolder) ) %>%
                filter( grepl( prefix, simLabel ))

  # Break up MP names into factor levels
  # MP labels are AM_Fsrce_eqbm
  splitMP <- function(  mpLab,
                        breaks = "_",
                        factorNames = c("AM","Fsrce","eqbm") )
  {
    splitMP <- stringr::str_split(mpLab, breaks)[[1]]

    outList <- vector(  mode = "list",
                        length = length(splitMP) )
    names(outList) <- factorNames
    for( k in 1:length(splitMP))
      outList[[k]] <- splitMP[k]

    outList
  }

  MPfactors <- lapply( X = info.df$mp, FUN = splitMP )
  MPfactors <- do.call(rbind, MPfactors)

  # Read in the performance tables
  mpTable <- cbind( info.df, MPfactors ) %>%
              mutate( perfPath = here::here("Outputs",groupFolder,simLabel,"simPerfStats.csv") )


  perfTableList <- lapply(  X = mpTable$perfPath,
                            FUN = read.csv,
                            header = TRUE,
                            stringsAsFactors = FALSE )

  names( perfTableList ) <- mpTable$simLabel

  # Summarise the perf tables into convergence rates
  summariseTable <- function( table )
  {
    table <- table %>%
              group_by(simLabel) %>%
              summarise(  obsCVmult = mean(projObsErrMult),
                          minConv   = min(pGoodReps),
                          maxConv   = max(pGoodReps),
                          meanConv  = mean(pGoodReps) )

    table
  }

  # Summarise and make table
  perfTableSummList <- lapply(  X = perfTableList,
                                FUN = summariseTable )

  convRateTable <- do.call(rbind, perfTableSummList )


  mpTable <- mpTable %>% left_join( convRateTable,
                                    by = "simLabel" )

  CVmults <- unique(mpTable$obsCVmult)
  CVmults <- CVmults[order(CVmults)]

  AMs       <- names(AMlabs)
  AMcols    <- RColorBrewer::brewer.pal(length(AMs), "Set2")
  AMjitter  <- seq( from = -.3, to = .3, length.out = length(AMs) )

  AMwidth <- AMjitter[2] - AMjitter[1]

  # Now set up plotting area
  plot( x = c(0,length(CVmults) + 1), y = c(0,1),
        xlab = "", ylab = "", axes=FALSE,
        type = "n" )
    axis( side = 1, at = 1:length(CVmults),
          labels = CVmults )
    axis( side = 2, las = 1  )
    box()
    for( aIdx in 1:length(AMs) )
    {
      amLab <- AMs[aIdx]
      SStable <- mpTable %>% filter( AM == amLab,
                                      grepl("SS", eqbm) ) %>%
                              dplyr::select( obsCVmult,
                                              meanConv )

      SStable <- SStable[order(SStable$obsCVmult),]

      MStable <- mpTable %>% filter( AM == amLab,
                                      grepl("MS", eqbm) ) %>%
                              dplyr::select( obsCVmult,
                                              meanConv )
      MStable <- MStable[order(MStable$obsCVmult),]

      rect( xleft = 1:3 + AMjitter[aIdx] - AMwidth/2,
            xright = 1:3 + AMjitter[aIdx],
            ybottom = 0,
            ytop = SStable$meanConv,
            border = "black",
            col = AMcols[aIdx] )

      rect( xleft = 1:3 + AMjitter[aIdx],
            xright = 1:3 + AMjitter[aIdx] + AMwidth/2,
            ybottom = 0,
            ytop = MStable$meanConv,
            border = "black",
            col = AMcols[aIdx] )



    }



} # END plotBatchConvergenceRate()


# plotBatchLossDists_CV()
# Multi-panel plot of relative and absolute
# loss under a batch of MPs for all stock/species
# combinations, including data-pooled and
# coast-wide aggregations.
plotBatchLossDists_CV <- function(  groupFolder = "diffCV_fixedF_longGrid",
                                    prefix = "sim_parBatfixedF",
                                    lossType = "abs",
                                    var = "C_ispt",
                                    period = 62:72,
                                    lossList = NULL,
                                    dim1 = 3,   # species (D,E,R, DP)
                                    dim2 = 1,   # stocks (H,Q,W, CW)
                                    qProbs = c(.05,.5,.95),
                                    refPts = "MSrefPts" )
{
  # First, read info files from the relevant
  # sims

  simFolderList <- list.dirs(here::here("Outputs",groupFolder),recursive = FALSE)
  simFolderList <- simFolderList[grepl(prefix, simFolderList)]

  infoFiles     <- lapply(  X = file.path(simFolderList,"infoFile.txt"),
                            FUN = lisread )


  if( lossType == "rel" )
  {
    lossArrayName <- "totRelLoss_isp"
    yLab <- "Total Relative Loss"
  }

  if( lossType == "abs" )
  {
    lossArrayName <- "totAbsLoss_isp"
    yLab <- "Total Absolute Loss (kt)"
  }

  # Break up MP names into factor levels
  # MP labels are AM_Fsrce_eqbm
  splitMP <- function(  mpLab,
                        breaks = "_",
                        factorNames = c("AM","Fsrce","eqbm") )
  {
    splitMP <- stringr::str_split(mpLab, breaks)[[1]]

    outList <- vector(  mode = "list",
                        length = length(splitMP) )
    names(outList) <- factorNames
    for( k in 1:length(splitMP))
      outList[[k]] <- splitMP[k]

    outList
  }



  for( lIdx in 1:length(infoFiles) )
  {
    infoFiles[[lIdx]] <- c(infoFiles[[lIdx]], splitMP(infoFiles[[lIdx]]$mp))
    infoFiles[[lIdx]] <- as.data.frame(infoFiles[[lIdx]])
  }

  info.df <- do.call(rbind, infoFiles) %>%
              mutate_if(is.factor, as.character)

  if( !is.null(refPts) )
    info.df <- info.df %>% filter( eqbm == refPts )




  splitCVmult <- function( scenName )
  {
    splitUnderscore <- stringr::str_split( scenName, "_" )[[1]][2]

    splitObsErr     <- stringr::str_split( splitUnderscore, "obsErr")[[1]][1]

    CV <- as.numeric(splitObsErr)

    CV
  }

  info.df<- info.df %>%
              mutate( CVmult = sapply( X = scenario, FUN = splitCVmult) )


  # Load loss files, place in a list
  simFolderNames  <- info.df$simLabel
  if( is.null(lossList) )
  {
    lossList        <- lapply(  X = simFolderNames,
                                FUN = .loadLoss,
                                folder = groupFolder )
    names(lossList) <- simFolderNames
  }


  # Calculate total loss for variable/period
  totLossList  <- lapply( X = lossList,
                          FUN = calcTotalLossPeriod,
                          var = var, period = period )
  names(totLossList) <- names(lossList)


  nS  <- lossList[[1]]$nS
  nP  <- lossList[[1]]$nP
  nT  <- lossList[[1]]$nT
  tMP <- lossList[[1]]$tMP
  pT  <- lossList[[1]]$pT

  # Get number of MPs
  nSims <- length(lossList)


  AMs       <- unique(info.df$AM)
  Fsources  <- unique(info.df$Fsrce)
  eqbm      <- unique(info.df$eqbm)
  MPs       <- unique(info.df$mp)
  CVs       <- unique(info.df$CVmult)
  CVs       <- CVs[order(CVs)]

  nAM       <- length(AMs)
  nSrce     <- length(Fsources)
  nEqbm     <- length(eqbm)
  nMP       <- length(MPs)



  # Some axis labels and translations from
  # simple axis label to info file AM name
  AMcodes <- c( SS = 1,
                MS = 2,
                DP = 3,
                CW = 4,
                TA = 5 )

  AMdictionary <- c(  singleStock     = "SS",
                      hierMultiStock  = "MS",
                      dataPooled      = "DP",
                      coastwide       = "CW",
                      totalAgg        = "TA" )

  AMcols          <- RColorBrewer::brewer.pal(nAM, "Dark2")
  names(AMcols)   <- names(AMdictionary)
  eqbmPCH  <- 15 + 1:nEqbm
  names(eqbmPCH) <- eqbm
  FsrceLTY   <- 1:nSrce
  names(FsrceLTY)  <- Fsources


  MPgrid <- list( AMs = names(AMdictionary),
                  Fsrce = Fsources,
                  eqbm = eqbm )

  MPgrid <- expand.grid(MPgrid)


  xJitter <- seq(from = -.4, to = .4, length.out = nMP )


  nReps         <- dim(lossList[[1]]$lossRel[[var]])[1]
  speciesNames  <- lossList[[1]]$speciesNames
  stockNames    <- lossList[[1]]$stockNames

  # Make an array to hold loss function values
  totLossArray_misp <- array(NA,  dim = c(nSims, nReps, nS+1, nP+1 ),
                                  dimnames = list(  simID = info.df$simLabel,
                                                    rep = 1:nReps,
                                                    species = speciesNames,
                                                    stock = stockNames ) )

  for( k in 1:nSims )
  {
    simID <- info.df$simLabel[k]
    totLossArray_misp[k,,,] <- totLossList[[simID]][[lossArrayName]]
  }

  totLossQuantiles_qmsp <- apply( X = totLossArray_misp,
                                  FUN = quantile,
                                  probs = qProbs,
                                  MARGIN = c(1,3,4) )

  # Start plotting windows
  nSpec   <- length(dim1)
  nStock  <- length(dim2)

  par(  mfrow = c(nSpec,nStock),
        mar = c(0.2,1.5,.2,1),
        oma = c(3,3,3,3) )

  for( sIdx in 1:nSpec )
    for( pIdx in 1:nStock )
    {
      s <- dim1[sIdx]
      p <- dim2[pIdx]

      lossRange <- max(abs(range(totLossQuantiles_qmsp[,,s,p])))

      plot( y = c(0,lossRange), x = c(0.5,3.5),
            type = "n", axes = FALSE  )

        mfg <- par("mfg")
        if( mfg[1] == mfg[3])
          axis( side = 1, at = 1:length(CVs), labels = CVs )

        axis( side = 2, las = 1 )

        if( mfg[1] == 1 )
          mtext( side = 3, text = speciesNames[s], font = 2, line = 1 )

        if( mfg[2] == mfg[4] )
          mtext( side = 4, text = stockNames[p], font = 2, line = 1 )

        grid()
        box()

        jitIdx <- 0
        for( amIdx in 1:nAM )
          for( eqIdx in 1:nEqbm )
          {
            jitIdx <- jitIdx + 1
            AMid    <- names(AMdictionary)[amIdx]
            eqbmID  <- eqbm[eqIdx]

            subInfo <- info.df %>%
                        filter( AM    == AMid,
                                eqbm  == eqbmID ) %>%
                        arrange( CVmult )

            simLabels <- subInfo$simLabel

            FsrceID <- subInfo$Fsrce[1]
            AMcode  <- AMdictionary[amIdx]

            MPlab <- paste( AMid, FsrceID, eqbmID, sep = "_")

            points( x = 1:3 + xJitter[jitIdx],
                    y = totLossQuantiles_qmsp[2,simLabels,s,p],
                    pch = eqbmPCH[eqbmID],
                    col = AMcols[AMid], cex = 1.5 )
            # lines(  x = 1:3 + xJitter[jitIdx],
            #         y = totLossQuantiles_qmsp[2,simLabels,s,p],
            #         lty = 1, lwd = .8, col = AMcols[AMid] )

            segments( y0 = totLossQuantiles_qmsp[1,simLabels,s,p],
                      y1 = totLossQuantiles_qmsp[3,simLabels,s,p],
                      x0 = 1:3 + xJitter[jitIdx],
                      lty = FsrceLTY[FsrceID],
                      col = AMcols[AMid], lwd = 2  )

            abline( v = 1:2 + .5,
                    col = "black", lty = 2, lwd = 3 )


          }

      }


  legend( x       = "bottomright",
          bty     = "n",
          legend  = c(AMdictionary,eqbm),
          pch     = c(rep(NA,5),eqbmPCH),
          lty     = c(rep(1,5),NA,NA),
          col     = c(AMcols,"black","black"),
          lwd     = 2 )


  mtext( side = 2, outer = TRUE, text = yLab, line = 2 )
  mtext( side = 1, outer = TRUE, text = "Observation Error SD multiplier ", line = 2 )


}


# plotBatchLossDists()
# Multi-panel plot of relative and absolute
# loss under a batch of MPs for all stock/species
# combinations, including data-pooled and
# coast-wide aggregations.
plotBatchLossDists <- function( groupFolder = "diffCV_fixedF_longGrid",
                                prefix = "sim_parBatfixedF",
                                lossType = "rel",
                                var = "C_ispt",
                                period = 62:72,
                                lossList = NULL )
{
  # First, read info files from the relevant
  # sims

  simFolderList <- list.dirs(here::here("Outputs",groupFolder))
  simFolderList <- simFolderList[grepl(prefix, simFolderList)]

  infoFiles     <- lapply(  X = file.path(simFolderList,"infoFile.txt"),
                            FUN = lisread )

  if( lossType == "rel" )
  {
    lossArrayName <- "totRelLoss_isp"
    yLab <- "Total Relative Loss"
  }

  if( lossType == "abs" )
  {
    lossArrayName <- "totAbsLoss_isp"
    yLab <- "Total Absolute Loss (kt)"
  }

  # Break up MP names into factor levels
  # MP labels are AM_Fsrce_eqbm
  splitMP <- function(  mpLab,
                        breaks = "_",
                        factorNames = c("AM","Fsrce","eqbm") )
  {
    splitMP <- stringr::str_split(mpLab, breaks)[[1]]

    outList <- vector(  mode = "list",
                        length = length(splitMP) )
    names(outList) <- factorNames
    for( k in 1:length(splitMP))
      outList[[k]] <- splitMP[k]

    outList
  }



  for( lIdx in 1:length(infoFiles) )
  {
    infoFiles[[lIdx]] <- c(infoFiles[[lIdx]], splitMP(infoFiles[[lIdx]]$mp))
    infoFiles[[lIdx]] <- as.data.frame(infoFiles[[lIdx]])
  }

  info.df <- do.call(rbind, infoFiles) %>%
              mutate_if(is.factor, as.character)


  # Load loss files, place in a list
  simFolderNames  <- info.df$simLabel
  if( is.null(lossList) )
  {
    lossList        <- lapply(  X = simFolderNames,
                                FUN = .loadLoss,
                                folder = groupFolder )
    names(lossList) <- simFolderNames
  }


  # Calculate total loss for variable/period
  totLossList  <- lapply( X = lossList,
                          FUN = calcTotalLossPeriod,
                          var = var, period = period )
  names(totLossList) <- names(lossList)


  nS  <- lossList[[1]]$nS
  nP  <- lossList[[1]]$nP
  nT  <- lossList[[1]]$nT
  tMP <- lossList[[1]]$tMP
  pT  <- lossList[[1]]$pT

  # Get number of MPs
  nSims <- length(lossList)


  AMs       <- unique(info.df$AM)
  Fsources  <- unique(info.df$Fsrce)
  eqbm      <- unique(info.df$eqbm)

  nAM       <- length(AMs)
  nSrce     <- length(Fsources)
  nEqbm     <- length(eqbm)

  # Some axis labels and translations from
  # simple axis label to info file AM name
  AMcodes <- c( SS = 1,
                MS = 2,
                DP = 3,
                CW = 4,
                TA = 5 )

  AMdictionary <- c(  singleStock     = "SS",
                      hierMultiStock  = "MS",
                      dataPooled      = "DP",
                      coastwide       = "CW",
                      totalAgg        = "TA" )

  AMcols          <- RColorBrewer::brewer.pal(nAM, "Dark2")
  names(AMcols)   <- names(AMdictionary)
  FsrcePCH  <- 15 + 1:nSrce
  names(FsrcePCH) <- Fsources
  eqbmLTY   <- 1:nEqbm
  names(eqbmLTY)  <- eqbm

  nReps         <- dim(lossList[[1]]$lossRel[[var]])[1]
  speciesNames  <- lossList[[1]]$speciesNames
  stockNames    <- lossList[[1]]$stockNames

  # Make an array to hold loss function values
  totLossArray_misp <- array(NA,  dim = c(nSims, nReps, nS+1, nP+1 ),
                                  dimnames = list(  simID = info.df$simLabel,
                                                    rep = 1:nReps,
                                                    species = speciesNames,
                                                    stock = stockNames ) )

  for( k in 1:nSims )
  {
    simID <- info.df$simLabel[k]
    totLossArray_misp[k,,,] <- totLossList[[simID]][[lossArrayName]]
  }

  totLossQuantiles_qmsp <- apply( X = totLossArray_misp,
                                  FUN = quantile,
                                  probs = c(0.25, 0.5, 0.75),
                                  MARGIN = c(1,3,4) )



  # Loop and plot - start by plotting ALL
  # combos
  yJitter <- seq(from = -.3, to = .3, length = nSrce * nEqbm )
  names(yJitter) <- apply(expand.grid(Fsources, eqbm), 1, paste, collapse="_")

  # Plot total loss for given period
  # on grid of Species/stocks
  par(  mfcol = c(nS+1,nP+1),
        mar = c(.2,1.5,.2,1),
        oma = c(3,3.5,3,3) )

  for( s in 1:(nS+1) )
    for( p in 1:(nP+1) )
    {
      lossRange <- max(abs(range(totLossQuantiles_qmsp[,,s,p])))

      plot( y = c(0,lossRange), x = c(0,6),
            type = "n", axes = FALSE  )

        mfg <- par("mfg")
        if( mfg[1] == mfg[3])
          axis( side = 1, at = AMcodes, labels = names(AMcodes) )

        axis( side = 2, las = 1 )

        if( mfg[1] == 1 )
          mtext( side = 3, text = speciesNames[s], font = 2, line = 1 )

        if( mfg[2] == mfg[4] )
          mtext( side = 4, text = stockNames[p], font = 2, line = 1 )

        grid()
        box()

        for( k in 1:nSims )
        {
          simID   <- info.df$simLabel[k]
          AMid    <- info.df$AM[k]
          FsrceID <- info.df$Fsrce[k]
          eqbmID  <- info.df$eqbm[k]
          AMcode  <- AMdictionary[AMid]

          jitY    <- yJitter[paste(FsrceID,eqbmID,sep = "_")]

          points( y   = totLossQuantiles_qmsp[2,simID,s,p],
                  x   = AMcodes[AMcode] + jitY,
                  pch = FsrcePCH[FsrceID],
                  col = AMcols[AMid], cex = 1.5 )
          segments( y0 = totLossQuantiles_qmsp[1,simID,s,p],
                    y1 = totLossQuantiles_qmsp[3,simID,s,p],
                    x0 = AMcodes[AMcode] + jitY,
                    lty = eqbmLTY[eqbmID],
                    col = AMcols[AMid], lwd = 2  )
        }
    } # END p loop
    # END s loop

  legend( x       = "bottomright",
          bty     = "n",
          legend  = c(Fsources,eqbm),
          pch     = c(FsrcePCH,NA,NA),
          lty     = c(NA,NA,eqbmLTY),
          lwd     = 2 )


  mtext( side = 2, outer = TRUE, text = yLab, line = 2 )
  mtext( side = 1, outer = TRUE, text = "AM Configuration", line = 2 )


} # END plotBatchLossDists()

# plotTotLossDists()
# Plots relative and absolute loss
# for a given simulation. Requires
# loss for a given baseline to have
# been calculated first, and saved
# into the sim folder
plotTotLossDists <- function( sim = 1,
                              groupFolder = "shortGrid",
                              lossType    = "rel",
                              var         = "SB_ispt",
                              save        = FALSE )
{
  # Load loss reports
  objLoss <- .loadLoss(sim = sim, groupFolder)
  simFolder <- objLoss$simFolder

  # Species/stock names and model dims
  speciesNames    <- objLoss$speciesNames
  stockNames      <- objLoss$stockNames
  tMP             <- objLoss$tMP
  nT              <- objLoss$nT
  nF              <- objLoss$nF
  nS              <- objLoss$nS
  nP              <- objLoss$nP
  pT              <- objLoss$pT
  simFolder       <- objLoss$simFolder

  loss_shortTerm  <- calcTotalLossPeriod(objLoss,var, period = tMP:(tMP+10))
  loss_medTerm    <- calcTotalLossPeriod(objLoss,var, period = tMP:(tMP+20))
  loss_longTerm   <- calcTotalLossPeriod(objLoss,var, period = tMP:nT)

  if( lossType == "rel" )
  {
    lossListName  <- "totRelLoss_isp"
    yLab          <- "Total relative loss (unitless)"
  }

  if( lossType == "abs" )
  {
    lossListName  <- "totAbsLoss_isp"
    yLab          <- "Total loss (kt)"
  }

  if( save )
    pdf( file = file.path(simFolder,paste(lossType,var,"LossCleveland.pdf",sep = "_")),
          width = 11, height = 8.5 )

  # Plot window - just for loss right now
  par(  mfcol = c(nS + 1, nP + 1),
        mar = c(.5,.5,.5,.5),
        oma = c(3,3,3,3) )

  for( s in 1:(nS+1) )
    for( p in 1:(nP + 1) )
    {
      lossShort_i <- loss_shortTerm[[lossListName]][,s,p]
      lossMed_i   <- loss_medTerm[[lossListName]][,s,p]
      lossLong_i  <- loss_longTerm[[lossListName]][,s,p]

      lossShort_q <- quantile(lossShort_i, probs = c(0.025, 0.5, 0.975) )
      lossMed_q   <- quantile(lossMed_i, probs = c(0.025, 0.5, 0.975) )
      lossLong_q  <- quantile(lossLong_i, probs = c(0.025, 0.5, 0.975) )

      maxLoss     <- max( abs(lossShort_q[is.finite(lossShort_q)]),
                          abs(lossMed_q[is.finite(lossMed_q)]),
                          abs(lossLong_q[is.finite(lossLong_q)]), na.rm = T )

      plot( x = c(0,4), y = c(0,maxLoss), type = "n",
            axes = FALSE )
        mfg <- par("mfg")
        if( mfg[1] == mfg[3] )
          axis( side = 1, at = 1:3, labels = c("S", "M", "L") )
        if( mfg[2] == 1 )
          axis( side = 2, las = 1)
        if( mfg[1] == 1 )
        {
          mtext( side = 3, text = speciesNames[s] )
        }
        if( mfg[2] == mfg[4] )
        {
          mtext( side = 4, text = stockNames[p] )
        }
        grid()
        box()

        segments( x0 = 1, y0 = lossShort_q[1], y1 = lossShort_q[3]  )
        points( x = 1, y = lossShort_q[2]  )
        segments( x0 = 2, y0 = lossMed_q[1], y1 = lossMed_q[3]  )
        points( x = 2, y = lossMed_q[2]  )
        segments( x0 = 3, y0 = lossLong_q[1], y1 = lossLong_q[3]  )
        points( x = 3, y = lossLong_q[2]  )

    }

  mtext( side = 1, text = "Time Period", outer = T, line = 2)
  mtext( side = 2, text = yLab, outer = T, line = 2)

  if( save )
    dev.off()
} # END plotLossDists()

# plotLossTulips()
# Plots relative and absolute loss
# for a given simulation. Requires
# loss for a given baseline to have
# been calculated first, and saved
# into the sim folder
plotLossTulip <- function(  sim = 1,
                            groupFolder = "shortGrid",
                            lossType    = "rel",
                            var         = "SB_ispt",
                            save        = FALSE )
{
  # Load loss reports
  objLoss <- .loadLoss(sim = sim, groupFolder)
  simFolder <- objLoss$simFolder

  if(lossType == "rel" )
  {
    loss_ispt <- objLoss$lossRel[[var]]
    yLab      <- "Relative Loss"
  }

  if(lossType == "abs" )
  {
    loss_ispt <- objLoss$lossRaw[[var]]
    yLab      <- "Absolute loss (kt)"
  }

  # Pull model dimensions
  tMP <- objLoss$tMP
  nT  <- objLoss$nT
  nF  <- objLoss$nF
  nS  <- objLoss$nS
  nP  <- objLoss$nP
  pT  <- objLoss$pT


  nReps <- dim(loss_ispt)[1]

  traces <- sample(1:nReps, size = 3)

  # Time steps
  fYear   <- objLoss$fYear
  tdxPlot <- tMP:nT
  years   <- seq(from = fYear, by = 1, length.out = nT)

  if( save )
  {
    graphics.off()
    pdf( file = file.path(simFolder,paste(lossType,var,"LossEnvelopes.pdf",sep = "_")),
          width = 11, height = 8.5 )
  }

  # Plot window - just for loss right now
  par(  mfcol = c(nS + 1, nP + 1),
        mar = c(.5,.5,.5,.5),
        oma = c(3,3,3,3) )

  speciesNames  <- c(objLoss$speciesNames,"data-pooled")
  stockNames    <- c(objLoss$stockNames,"coast-wide")

  loss_qspt <- apply( X = loss_ispt, FUN = quantile,
                      probs = c(0.025, 0.5, 0.975),
                      MARGIN = c(2,3,4) )

  plotYrs <- years[tdxPlot]


  # Loop and plot
  for( s in 1:(nS+1) )
    for( p in 1:(nP+1) )
    {
      maxLoss <- 1


      plot( x = range(years[tdxPlot]), y = c(-maxLoss,maxLoss),
            type = "n", axes = FALSE  )
      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1 )
      if( mfg[2] == 1 )
        axis( side = 2, las = 1)
      if( mfg[1] == 1 )
      {
        mtext( side = 3, text = speciesNames[s] )
      }
      if( mfg[2] == mfg[4] )
      {
        mtext( side = 4, text = stockNames[p] )
      }
      grid()
      box()

      abline( v = years[tMP], lty = 2, lwd = 3 )

      polygon(  x = c(years,rev(years)),
                y = c(loss_qspt[1,s,p,],rev(loss_qspt[3,s,p,])),
                border = NA, col = scales::alpha("grey10",.4) )
      lines(  x = years, y = loss_qspt[2,s,p,],
              col = "black", lwd = 3 )
      for( traceIdx in traces )
        lines( x = years, y = loss_ispt[traceIdx,s,p,],
                col = "black", lwd = 1 )
      grid()

      abline( h = 0, lty = 2, lwd = 1 )


    }
    mtext( side = 1, text = "Year", line = 2, outer = TRUE )
    mtext( side = 2, text = yLab, line = 2, outer = TRUE )

  if( save )
    dev.off()

} # END plotLoss()


# compareSimCondIdxData()
# Compares simulated index data to the data
# from the conditioning assessment.
compareSimCondIdxData <- function(  obj = blob, sIdx = 1,
                                    pIdx = 1, iRep = 1 )
{
  # get model dims and repObj
  nT <- obj$om$nT
  tMP <- obj$om$tMP
  nF <- obj$om$nF

  repObj <- obj$ctlList$opMod$histRpt

  # First, pull out the real data
  condI_ft              <- repObj$I_pgt[pIdx,,]
  condI_ft[5,1:(tMP-1)] <- repObj$combI_pt[pIdx,]
  simI_ft               <- obj$mp$data$I_ispft[iRep,sIdx,pIdx,,1:(tMP-1)]


  # Check which flets have data
  dataFleets <- c()
  for( f in 1:nF )
    if(any(condI_ft[f,] > 0))
      dataFleets <- c(dataFleets,f)

  condI_ft[condI_ft < 0] <- NA
  simI_ft[simI_ft < 0] <- NA

  fYear <- obj$ctlList$opMod$fYear
  yrs   <- seq( from = fYear, length.out = tMP + 3, by = 1)

  # Multi-panel plot overlaying
  #   1. Data
  #   2. Residuals - annually with distro at RHS
  #   3. Ratios of designs - if blended

  # Get standardised resids (obs errors)
  condResid_ft <- array(NA, dim = c(nF,tMP-1))
  condResid_ft[dataFleets,] <- repObj$zComb_pt
  simResid_ft               <- obj$om$errors$delta_ispft[iRep,sIdx,pIdx,,1:(tMP-1)]

  # get ratios
  condRatioI_ft <- repObj$rI_pgt[pIdx,,]
  simRatioI_ft  <- obj$om$rI_ispft[iRep,sIdx,pIdx,,]

  par(  mfrow = c(3,length(dataFleets)),
        oma   = c(3,3,2,2),
        mar   = c(.5,2,.5,2) )
  for( f in dataFleets )
  {
    # First plot data
    plot(x = range(yrs), y = c(0,max(condI_ft[f,], simI_ft[f,],na.rm = T)),
          type = "n", axes = FALSE )
      axis( side = 2, las = 1)
      mtext(side = 2, text= "Spawn Index (kt)", line = 2)
      grid()
      box()
      points( x = yrs[1:(tMP-1)], y = condI_ft[f,],  pch = 16, col = "grey40" )
      points( x = yrs[1:(tMP-1)], y = simI_ft[f,],   pch = 21, col = "salmon", bg = NA )

    # Now plot residuals
    plot(x = range(yrs), y = range(condResid_ft[f,], simResid_ft[f,],na.rm = T),
          type = "n", axes = FALSE )
      abline(h = 0, lty = 2)
      mtext(side = 2, text= "Standardised residual", line = 2)
      axis( side = 2, las = 1)
      grid()
      box()
      points( x = yrs[1:(tMP-1)], y = condResid_ft[f,], pch = 16, col = "grey40" )
      segments( x0 = yrs[1], x1 = yrs[tMP-1],
                y0 = mean(condResid_ft[f,],na.rm = T), lty = 2, col = "grey40" )
      segments( x0 = yrs[1], x1 = yrs[tMP-1], lty = 2, col = "salmon",
                y0 = mean(simResid_ft[f,],na.rm = T) )
      points( x = yrs[1:(tMP-1)], y = simResid_ft[f,], pch = 21, col = "salmon", bg = NA )

    # need to add violin plots here

    # Now ratios - although they are the same
    plot(x = range(yrs), y = c(0,1),
          type = "n", axes = FALSE )
      axis( side = 2, las = 1)
      axis( side = 1 )
      mtext(side = 2, text= "Design ratio", line = 2)
      grid()
      box()
      rect( xleft = yrs[1:(tMP-1)]-0.3,
            xright =yrs[1:(tMP-1)],
            ybottom = 0,
            ytop = condRatioI_ft[f,], col = "grey40" )
      rect( xleft = yrs[1:(tMP-1)],
            xright =yrs[1:(tMP-1)]+0.3,
            ybottom = 0,
            ytop = simRatioI_ft[f,], col = "salmon" )

  }
  mtext(side = 1, text= "Year", outer = TRUE, line = 2)
} # END compareSimCondIdxData()

# compareSimCondAgeData()
# Compares simulated age data to the data used
# to fit the conditioning assessment.
compareSimCondAgeData <- function(  obj = blob, sIdx = 1,
                                    pIdx = 1, iRep = 1,
                                    fIdx = 1 )
{
  # get model dims and repObj
  tMP <- obj$om$tMP
  nF  <- obj$om$nF
  nA  <- obj$om$nA

  repObj <- obj$ctlList$opMod$histRpt
  nT  <- repObj$nT
  yrs <- seq(from = 1951, length.out = nT)


  # Now, pull conditioned ages
  condA_at <- repObj$A_apgt[,pIdx,fIdx,]
  condA_at[condA_at < 0] <- NA
  sel_a <- repObj$sel_apgt[,pIdx,fIdx,1]

  for( a in 1:nA )
    if(sel_a[a] == 0)
      condA_at[a,] <- 0

  # Simulated ages
  simA_at  <- obj$mp$data$A_iaxspft[iRep,,1,sIdx,pIdx,fIdx,1:nT]

  # Now, we need to go through and sweep out the
  # sample sizes and turn into proportions, and
  # then plot
  condSampSize_t  <- apply(X = condA_at, FUN = sum, MARGIN = 2)
  simSampSize_t   <- apply(X = simA_at, FUN = sum, MARGIN = 2)

  for( t in 1:nT )
  {
    condA_at[,t]  <- condA_at[,t] / condSampSize_t[t]
    simA_at[,t]   <- simA_at[,t] / simSampSize_t[t]
  }

  # Check which years have data
  simDatIdx   <- which(!is.na(simSampSize_t))
  condDatIdx  <- which(!is.na(condSampSize_t))
  dataTimeIdx <- union(simDatIdx,condDatIdx)

  nPlotYrs  <- length(dataTimeIdx) + 2
  nCols     <- ceiling(sqrt(nPlotYrs))
  nRows     <- ceiling(sqrt(nPlotYrs))

  par(  mfcol = c(nRows,nCols),
        oma = c(4,4,3,3),
        mar = c(.1,.1,.1,.1)
      )

  for( tIdx in dataTimeIdx )
  {
    yr <- yrs[tIdx]

    plot(x = c(1,nA), y = c(0,1), type = "n",
          axes = FALSE )
      mfg <- par("mfg")
      if(mfg[1] == mfg[3] | tIdx == max(dataTimeIdx))
        axis(side = 1)
      if(mfg[1] == 1)
        axis(side = 3)
      if( mfg[2] == 1 )
        axis(side = 2, las = 1)
      if(mfg[2] == mfg[4])
        axis(side = 4, las = 1)
      grid()
      box()
      rect( xleft = 1:nA - .3,
            xright = 1:nA,
            ybottom = rep(0,nA),
            ytop = condA_at[,tIdx],
            col = "grey40", border = NA )
      rect( xleft = 1:nA,
            xright = 1:nA + .3,
            ybottom = rep(0,nA),
            ytop = simA_at[,tIdx],
            col = "salmon", border = NA )
      text( x = nA/2, y = 0.9,
            label = yr, font = 2 )
  }
  plot(x = c(1,nA), y = c(0,1), type = "n",
        axes = FALSE)
  legend( x = 2, y = 0.8, bty = "n",
          legend = c( "Real Data",
                      "Simulated"),
          pch = 15, cex = 1.5, col = c("grey40","salmon"))
  mtext( side = 1, text = "Age", outer= TRUE, line = 2)
  mtext( side = 2, text = "Proportion caught-at-age", outer= TRUE, line = 2)
}

plotSSvsMSrefPts <- function( obj = blob )
{
  # Pull reference points and curves
  rp            <- obj$rp[[1]]
  refCurves     <- rp$refCurves
  EmsyRefPts    <- rp$EmsyRefPts
  EmsyMSRefPts  <- rp$EmsyMSRefPts
  FmsyRefPts    <- rp$FmsyRefPts

  nP  <- obj$om$nP
  nS  <- obj$om$nS
  tMP <- obj$om$tMP
  nT  <- obj$om$nT

  # Pull qF
  qF_sp       <- blob$om$qF_ispft[1,,,2,nT]

  speciesNames <- blob$om$speciesNames
  stockNames   <- blob$om$stockNames

  # Want to compare Umsy and Bmsy between
  # multi-stock and single-stock models
  BmsySS_sp <- FmsyRefPts$BeqFmsy_sp
  BmsyMS_sp <- EmsyMSRefPts$BeqEmsy_sp

  UmsySS_sp <- FmsyRefPts$Umsy_sp
  UmsyMS_sp <- EmsyMSRefPts$Umsy_sp

  B0_sp     <- obj$ctlList$opMod$histRpt$B0_sp

  speciesCols <- RColorBrewer::brewer.pal(n = nS, "Set2")

  cols_sp   <- matrix(  speciesCols, nrow = nS, ncol = nP,
                        byrow = FALSE )
  pch_sp    <- matrix(  21:23, nrow = nS, ncol = nP,
                        byrow = TRUE )


  par(  mfrow = c(1,2),
        mar = c(1.5,4,0,0),
        oma = c(2,1,1,1) )

  # First plot Umsy
  plot( x = c(0,max(UmsySS_sp, UmsyMS_sp) ),
        y = c(0,max(UmsySS_sp, UmsyMS_sp) ),
        type = "n",
        xlab = "",
        ylab = "", las = 1 )
    mtext(side =1, text = "Single-species Umsy", line = 2 )
    mtext(side =2, text = "Multi-species Umsy", line = 3 )
    points( x = UmsySS_sp, y = UmsyMS_sp,
            cex = 3 * B0_sp / max(B0_sp),
            bg = cols_sp, pch = pch_sp )
    abline( a = 0, b = 1, lty = 2 )

  # First plot Umsy
  plot( x = c(0,max(BmsySS_sp, BmsyMS_sp) ),
        y = c(0,max(BmsySS_sp, BmsyMS_sp) ),
        type = "n", las = 1,
        xlab = "",
        ylab = "" )
    mtext(side =1, text = "Single-species Bmsy", line = 2 )
    mtext(side =2, text = "Multi-species Bmsy", line = 2 )
    points( x = BmsySS_sp, y = BmsyMS_sp,
            cex = 3 * B0_sp / max(B0_sp),
            bg = cols_sp, pch = pch_sp )
    abline( a = 0, b = 1, lty = 2 )


} # END plotSSvsMSrefPts


# plotEffYieldCurves()
# Function for plotting effort based
# yield curves for each species and
# the complex in a stock area.
plotEffYieldCurves <- function( obj = blob )
{
  # First, pull reference points and curves
  rp            <- obj$rp[[1]]
  refCurves     <- rp$refCurves
  EmsyRefPts    <- rp$EmsyRefPts
  EmsyMSRefPts  <- rp$EmsyMSRefPts
  FmsyRefPts    <- rp$FmsyRefPts

  nT  <- obj$om$nT
  nP  <- obj$om$nP
  nS  <- obj$om$nS
  tMP <- obj$om$tMP

  # Pull qF
  qF_sp       <- blob$om$qF_ispft[1,,,2,nT]

  speciesNames <- blob$om$speciesNames

  # Now compute MSY for SS and MS curves
  # MSY is still the same for SS, just need the effort
  # that corresponds to it, which is Fmsy/qF

  # get SS ref points
  EmsySS_sp <- rp$FmsyRefPts$Fmsy_sp / qF_sp
  MSYSS_sp  <- rp$FmsyRefPts$YeqFmsy_sp

  # Now get MS ref points
  Yeq_spe   <- rp$refCurves$EffCurves$Yeq_spe
  Yeq_spe[Yeq_spe < 0] <- NA
  Yeq_pe    <- apply(X = Yeq_spe, FUN = sum, MARGIN = c(2,3), na.rm = T )
  Yeq_pe[Yeq_pe == 0] <- NA
  Yeq_pe[,1] <- 0
  Eseq      <- rp$refCurves$EffCurves$E

  EmsyMS_p  <- EmsyMSRefPts$EmsyMS_p
  MSYMS_sp  <- EmsyMSRefPts$YeqEmsy_sp
  MSYMS_p   <- EmsyMSRefPts$YeqEmsy_p
  BmsyMS_sp <- EmsyMSRefPts$BeqEmsy_sp

  Yeq_e <- apply( X = Yeq_pe, FUN = sum, na.rm = T,
                  MARGIN = 2 )
  maxEval <- max(which(Yeq_e > 0))

  specCols <- RColorBrewer::brewer.pal(n = nS, "Dark2")

  par( mfrow = c(3,1), mar = c(1,1,1,1), oma = c(3,3,1,1) )
  for( p in 1:nP )
  {

    plot( x = c(0,Eseq[maxEval]), y = c(0, max(Yeq_pe[p,],na.rm = T) ),
          type = "n", xlab = "", ylab = "", axes = F )
      axis(side = 2, las = 1)
      box()
      grid()

      for( s in 1:nS )
      {
        lines( x = Eseq, y = Yeq_spe[s,p,],
               col = specCols[s], lty = 1, lwd = 2 )
        # lines( x = Eseq, y = Beq_spe[s,p,],
        #        col = specCols[s], lty = 2, lwd = 2 )

      }
      mfg <- par("mfg")
      if( mfg[1] == mfg[3])
        axis(side = 1)

      lines(  x = Eseq, y = Yeq_pe[p,],
              col = "black", lty = 1, lwd = 3 )

      segments( x0 = EmsyMS_p[p], col = "grey40",
                y0 = 0, y1 = MSYMS_p[p], lty = 2  )

      segments( x0 = 0, x1 = EmsyMS_p[p], col = "grey40",
                y0 = c(MSYMS_p[p],MSYMS_sp[,p]), lty = 2  )

      if(  p == 1 )
        legend( x = "topright", bty = "n",
                col = c(specCols,"black"),
                lwd = c(2,2,2,3),
                legend = c(speciesNames,"Complex") )

  }

  mtext( outer = TRUE, side = 1, text = "Total Effort Index", line = 2 )
  mtext( outer = TRUE, side = 2, text = "Equilibrium Yield (kt)", line = 2 )


  outList <- list(  EmsyMS_p = EmsyMS_p,
                    EmsySS_sp = EmsySS_sp,
                    MSYMS_p = MSYMS_p,
                    MSYSS_sp = MSYSS_sp,
                    MSYMS_sp = MSYMS_sp )
} # END plotEffYieldCurves()

# plotFYieldCurves()
# Function for plotting effort based
# yield curves for each species and
# the complex in a stock area.
plotFYieldCurves <- function( obj = blob,
                              maxF = 1,
                              xLim =c(0,1.2),
                              yLim = c(0,50) )
{
  # First, pull reference points and curves
  rp            <- obj$rp[[1]]
  refCurves     <- rp$refCurves
  FmsyRefPts    <- rp$FmsyRefPts

  nT  <- obj$om$nT
  nP  <- obj$om$nP
  nS  <- obj$om$nS
  tMP <- obj$om$tMP


  speciesNames <- blob$om$speciesNames

  # Now compute MSY for SS and MS curves
  # MSY is still the same for SS, just need the effort
  # that corresponds to it, which is Fmsy/qF

  # get SS ref points
  Fmsy_sp   <- rp$FmsyRefPts$Fmsy_sp
  MSY_sp    <- rp$FmsyRefPts$YeqFmsy_sp
  Bmsy_sp   <- rp$FmsyRefPts$BeqFmsy_sp


  # Now get MS ref points
  Yeq_spf   <- rp$refCurves$Yeq_spf
  Yeq_spf[Yeq_spf < 0] <- NA
  Beq_spf   <- rp$refCurves$Beq_spf
  Beq_spf[Yeq_spf < 0] <- NA
  Beq_spf[Beq_spf < 0] <- NA
  Fseq      <- rp$refCurves$F

  maxF <- Fseq[max(which(!is.na(Yeq_spf[1,1,])))]

  specCols <- RColorBrewer::brewer.pal(n = nS, "Dark2")

  par( mfrow = c(nP,nS), mar = c(1,1,1,1), oma = c(3,3,1,1) )
  for( s in 1:nS )
    for( p in 1:nP )
    {
      if(is.null(xLim))
        xLim <- c(0,maxF)
      if(is.null(yLim))
        yLim <- c(0, max(Yeq_spf[s,p,],Beq_spf[s,p,],na.rm = T) )

      plot( x = xLim, y = yLim,
            type = "n", xlab = "", ylab = "", axes = F,
            xaxs = "i" )
        axis(side = 2, las = 1)
        box()
        grid()


        mfg <- par("mfg")
        if( mfg[1] == mfg[3])
          axis(side = 1)

        lines(  x = Fseq, y = Yeq_spf[s,p,],
                col = "black", lty = 1, lwd = 3 )
        lines(  x = Fseq, y = Beq_spf[s,p,],
                col = "steelblue", lty = 1, lwd = 3 )


        segments( x0 = Fmsy_sp[s,p], col = "red",
                  y0 = 0, y1 = Bmsy_sp[s,p], lty = 2  )

        segments( x0 = 0, x1 = Fmsy_sp[p], col = "red",
                  y0 = c(MSY_sp[s,p],Bmsy_sp[s,p]), lty = 2  )


        legend( x = "topright",
                legend = c( "Yield",
                            paste0("MSY = ",round(MSY_sp[s,p],2)),
                            "Spawning Biomass",
                            paste0("Bmsy = ",round(Bmsy_sp[s,p],2)),
                            paste0("Fmsy = ",round(Fmsy_sp[s,p],2)) ),
                bty = "n", lwd = c(2,NA,2,NA,NA),
                col = c("black",NA,"steelblue",NA,NA) )

    }

  mtext( outer = TRUE, side = 1, text = "Fishing mortality (/yr)", line = 2 )
  mtext( outer = TRUE, side = 2, text = "Equilibrium Yield/Biomass (kt)", line = 2 )


  outList <- list(  Fmsy_sp = Fmsy_sp,
                    MSY_sp  = MSY_sp  )
} # END plotEffYieldCurves()


# Plot predator catchability vs prey biomass
plotPredqvB <- function(  obj = blob, showProj = FALSE, iRep = 1,
                          nTrace = 3, traceSeed = 8675209 )
{
  # Get pred gears
  predGears <- obj$ctlList$opMod$predGears
  predType  <- obj$ctlList$opMod$predType
  FRType    <- obj$ctlList$opMod$predFRtype

  # just the functional response ones
  predGears <- predGears[predType != 0 & FRType == 1]

  # Now pull pred effort
  tMP <- obj$om$tMP
  nP  <- obj$om$nP
  nF  <- obj$om$nF
  nT  <- obj$om$nT

  # age 2+ biomass
  B_t   <- obj$om$B_ispt[iRep,1,1,1:(tMP-1)]
  vB_ft <- obj$om$vB_ispft[iRep,1,1,,1:(tMP-1)]
  q_ft  <- obj$om$qF_ispft[iRep,1,1,,1:(tMP-1)]


  Brange <- c(0,max(B_t))

  gearNames <- names(obj$om$fleetType_f)

  qF_ift      <- obj$om$qF_ispft[,1,1,,tMP:nT]
  vBproj_ift  <- obj$om$vB_ispft[,1,1,,tMP:nT]

  eta_f     <- obj$om$eta_if[1,]
  h_f       <- obj$om$h_if[1,]
  lambda_f  <- obj$om$lambda_if[1,]



  nPred <- length(predGears)
  par(mfcol = c(nPred,1), oma = c(3,4.5,1,1),
        mar = c(.1,1,.1,1) )
  for( f in predGears)
  {
    qrange <- c(0,max(q_ft[f,],na.rm = T))

    qBtable <- data.frame(q = q_ft[f,], B = vB_ft[f,])

    # Smoother
    smoothline <- loess.smooth(x = vB_ft[f,], y = q_ft[f,])

    # Fit Disc
    Bseq <- seq( from = 0.001, to = max(vB_ft[f,]),
                  length.out = 200 )
    qseq <- eta_f[f] * Bseq^(lambda_f[f]-1) / (1 + h_f[f] * eta_f[f] * Bseq^lambda_f[f])

    qrange <- range(qseq,qrange)

    if(showProj)
    {
      qrange <- range(qrange,qF_ift[,f,])
    }



    plot(x = Brange, y = qrange, las = 1,
          type = "n", axes = FALSE )
      mfg <- par("mfg")
      if(mfg[1] == mfg[3])
        axis(side = 1)
      if(mfg[2] == 1)
        axis(side = 2, las = 1)
      grid()
      box()

      points( x = vB_ft[f,], y = q_ft[f,],
              col = "grey75", pch = 16, cex = 1.2 )

      lines(x = Bseq, y = qseq, lwd = 2, col = "darkgreen")
      lines(smoothline, lwd = 2, lty = 2, col = "black" )

      if(showProj)
        points( x = vBproj_ift[,f,],
                y = qF_ift[,f,] )


      legend( x = "topright",  bty ="n",
              legend = c( gearNames[f],
                          paste0("eta = ", signif(eta_f[f],2)),
                          paste0("h = ", signif(h_f[f],2)),
                          paste0("lambda = ", signif(lambda_f[f],2)),
                          "Loess Smoother",
                          "Disc equation" ),
              lty = c(NA,NA,NA,NA,2,1),
              lwd = c(NA,NA,NA,NA,2,2),
              col = c(NA,NA,NA,NA,"black","darkgreen") )

  }

  mtext(side = 1, outer = TRUE, text = "Vulnerable Herring Biomass (kt)", line = 2 )
  mtext( side = 2, outer = TRUE, text = "Predation Mortality per unit of predator", line = 3 )


}

solveSimEqbria <- function(empRefCurves,
                            maxXspline,
                            sIdx = 1, pIdx = 1)
{



  # Pull empirical ref curves
  C_k <- empRefCurves$C_spk[sIdx,pIdx,]
  B_k <- empRefCurves$B_spk[sIdx,pIdx,]
  F_k <- empRefCurves$F_spk[sIdx,pIdx,]
  U_k <- empRefCurves$U_spk[sIdx,pIdx,]
  M_k <- empRefCurves$M_spk[sIdx,pIdx,]
  R_k <- empRefCurves$R_spk[sIdx,pIdx,]

  # Order
  Forder <- order(F_k)
  C_k <- C_k[Forder]
  B_k <- B_k[Forder]
  F_k <- F_k[Forder]
  U_k <- U_k[Forder]
  M_k <- M_k[Forder]
  R_k <- R_k[Forder]

  # Make spline
  CUspline <- splinefun(x = U_k, y = C_k)
  BUspline <- splinefun(x = U_k, y = B_k)

  maxU <- max(U_k)
  maxC <- max(C_k)
  maxB <- max(B_k)
  maxM <- max(M_k)
  maxR <- max(R_k)

  Umsy  <- try(uniroot(  interval = c(0.01,maxXspline),
                        f = CUspline,
                        deriv = 1 )$root)
  if( class(Umsy) == "try-error")
    Umsy <- 0

  MSY  <- CUspline(Umsy)
  Bmsy <- BUspline(Umsy)



  outList <- list(  C_k = C_k,
                    U_k = U_k,
                    F_k = F_k,
                    B_k = B_k,
                    M_k = M_k,
                    R_k = R_k,
                    Umsy  = Umsy,
                    Bmsy  = Bmsy,
                    MSY   = MSY )

  outList
}

plotSimEqbria <- function(  folder = "CC_refPts",
                            saveEmpRefCurves = NULL,
                            maxXspline = 0.1,
                            sIdx = 1, pIdx = 1 )
{
  # Load ref curve list
  if(is.null(saveEmpRefCurves))
    load(here::here("Outputs",folder,"empRefCurves.RData"))

  simEqList <- solveSimEqbria(  empRefCurves = saveEmpRefCurves,
                                maxXspline = maxXspline)


  # ref curves
  C_k <- simEqList$C_k
  B_k <- simEqList$B_k
  F_k <- simEqList$F_k
  U_k <- simEqList$U_k
  M_k <- simEqList$M_k
  R_k <- simEqList$R_k

  maxU <- max(U_k)
  maxC <- max(C_k)
  maxB <- max(B_k)
  maxR <- max(R_k)
  maxM <- max(M_k)

  # Ref pts
  Umsy  <- simEqList$Umsy
  Bmsy  <- simEqList$Bmsy
  MSY   <- simEqList$MSY

  empUmsy <- round(Umsy,2)
  empBmsy <- round(Bmsy,2)
  empMSY  <- round(MSY,2)

  maxUidx <- max(which(B_k > 1e-2))
  Ucrash  <- (U_k[maxUidx])
  Bcrash  <- (B_k[maxUidx])
  Ccrash  <- (C_k[maxUidx])
  Mcrash  <- (M_k[maxUidx])
  Rcrash  <- (R_k[maxUidx])


  maxU <- 1.2*Ucrash

  par(mfrow = c(4,1),
      mar = c(.1,2,.1,2),
      oma = c(4,4,1,1) )


  # First, plot yield
  plot( x = c(0,maxU), y = c(0,1.2*maxC),
        axes = FALSE, type = "n" )
    axis(side = 2, las =1)
    grid()
    box()
    lines( x = U_k[1:maxUidx], y = C_k[1:maxUidx], col = "black", lwd = 3)
    points( x = empUmsy, y = empMSY, pch = 16, col = "grey40", cex = 2)
    points( x = Ucrash, y = Ccrash, pch = 21, bg = "red", cex = 2)
    # points( x = U_k[1:maxUidx], y = C_k[1:maxUidx], col = "black")
    mtext(side = 2, text = "Yield (kt)", line = 3)
    legend( "topright", bty = "n", cex = 2,
            legend = c( paste("Umsy   = ", empUmsy, sep = "" ),
                        paste(" MSY   = ", empMSY, sep = "") ,
                        paste("Ucrash = ", round(Ucrash,2),sep = "") ) )

  # Biomass
  plot( x = c(0,maxU), y = c(0,maxB),
        axes = FALSE, type = "n" )
    axis(side = 2, las =1)
    grid()
    box()
    lines( x = U_k[1:maxUidx], y = B_k[1:maxUidx], col = "black", lwd = 3)
    points( x = empUmsy, y = empBmsy, pch = 16, col = "grey40", cex = 2)
    points( x = Ucrash, y = Bcrash, pch = 21, bg = "red", cex = 2)
    legend( "topright", bty = "n", cex = 2,
            legend = c( paste("Bmsy = ", empBmsy, sep = "") ) )
    mtext(side = 2, text = "Spawning\nBiomass (kt)", line = 2.5)

  # Mortality
  plot( x = c(0,maxU), y = c(0,maxM),
        axes = FALSE, type = "n" )
    axis(side = 2, las =1)
    grid()
    box()
    lines( x = U_k[1:maxUidx], y = M_k[1:maxUidx], col = "black", lwd = 3)
    points( x = Ucrash, y = Mcrash, pch = 21, bg = "red", cex = 2)
    mtext(side = 2, text = "Natural\nMortality (/yr)", line = 2.5)

  # Recruitment
  plot( x = c(0,maxU), y = c(0,maxR),
        axes = FALSE, type = "n" )
    axis(side = 2, las =1)
    axis(side = 1)
    grid()
    box()
    lines( x = U_k[1:maxUidx], y = R_k[1:maxUidx], col = "black", lwd = 3)
    points( x = Ucrash, y = Rcrash, pch = 21, bg = "red", cex = 2)
    mtext(side = 2, text = "Recruitment (1e6)", line = 3)

    mtext(side = 1, text = "Harvest rate C / (C + B) (/yr)", line = 2.5)

}


# plotFATFRlogit()
# Plots SISCA outputs of predation mortality rates
# for selected predators and fits logistic Foraging Arena
# Theory style functional responses to those outputs
plotFATFRlogit <- function( obj = blob, iRep = 1 )
{
  # Get pred gears
  predGears <- obj$ctlList$opMod$predGears
  predType  <- obj$ctlList$opMod$predType
  FRType    <- obj$ctlList$opMod$predFRtype

  # just the functional response ones
  predGears <- predGears[predType != 0 & FRType == 2]
  gearNames <- names(obj$om$fleetType_f)

  # Now pull model dims
  tMP <- obj$om$tMP
  nP  <- obj$om$nP
  nF  <- obj$om$nF
  nT  <- obj$om$nT


  # Predator biomass
  P_ft   <- obj$om$E_ipft[iRep,1,,1:(tMP-1)]
  F_ft   <- obj$om$F_ispft[iRep,1,1,,1:(tMP-1)]
  q_ft   <- obj$om$qF_ispft[iRep,1,1,,1:(tMP-1)]

  qF_ispft     <- obj$om$qF_ispft[,1,1,,tMP:nT, drop=F]
  Pproj_ipft   <- obj$om$E_ipft[1:1,1,,tMP:nT, drop=F]

  Prange <- c(0,max(P_ft,Pproj_ipft, na.rm=T))

  eta_f     <- obj$om$eta_if[iRep,]
  h_f       <- obj$om$h_if[iRep,]
  lambda_f  <- obj$om$lambda_if[iRep,]
  qmin_f    <- obj$om$qmin_if[iRep,]

  # Now, we loop over predators and plot
  # the lines through the points.
  par(mfrow = c(length(predGears),2),
      mar = c(2,1,2,1),
      oma = c(3,4,1,4) )
  for( f in predGears )
  {
    maxP <- max(P_ft[f,],Pproj_ipft[,1,f,])
    Pseq <- seq(from = 0, to = maxP, length.out =1000)

    plot( x = c(0,maxP), y = c(0,1.5*max(F_ft[f,])),
          type = "n", axes = FALSE )
      mfg <- par("mfg")
      axis(side = 1)
      axis(side = 2, las = 1)
      grid()
      box()

      optMseq <- FATFRlogit(  P = Pseq,
                              v = eta_f[f],
                              k = h_f[f],
                              lambda = lambda_f[f],
                              qmin = qmin_f[f])

      lines(x = Pseq, y = optMseq, col = "royalblue", lwd = 4)
      points( x = P_ft[f,], y = F_ft[f,], pch = 16, col = "salmon" )

      legend( x = "bottomright", bty = "n",
              legend = c( gearNames[f],
                          paste0("qmin = ",signif(qmin_f[f],3)),
                          paste0("v = ",signif(eta_f[f],3)),
                          paste0("k = ",signif(h_f[f],3)),
                          paste0("lambda = ",signif(lambda_f[f],3))))

      plot( x = c(0,maxP), y = c(0,max(q_ft[f,])),
          type = "n", axes = FALSE )
      mfg <- par("mfg")
      axis(side = 1)
      axis(side = 4, las = 1)
      grid()
      box()
      points( x = P_ft[f,], y = q_ft[f,], pch = 16, col = "salmon" )

      optqseq <- optMseq/Pseq
      lines(x = Pseq, y = optqseq, col = "royalblue", lwd = 4)



  }
  mtext( side = 2, text = "Predation mortality rate (/yr)", outer = TRUE, line = 2.5)
  mtext( side = 4, text = "Per-capita predation mortality rate (yr N)^{-1}", outer = TRUE, line = 2.5)

  mtext(side = 1, outer = TRUE, line = 2, text = "Predator Abundance")

} # END plotFATFRlogit


saveBatchGridData <- function(sims = NULL, nHRs=NULL,
                              folder = "impPred_tvM_80Sr20Gill_perfU_new",
                              period = 50, sIdx=1, pIdx=1,
                              catchType='C+P')
{

  simFolders <- list.files(here('Outputs',folder),'sim_parBat')

  if(is.null(sims))
    sims <- 1:length(simFolders)

  nSims <- length(sims)
  blobList <- vector(mode = "list", length = nSims)

  .loadSim(sims[1], folder = folder)

  nS <- blob$om$nS
  nP <- blob$om$nP
  nF <- blob$om$nF
  nT <- blob$om$nT
  nG <- blob$om$nF
  nA <- blob$om$nA

  tMP <- blob$om$tMP

  goodReps <- blob$goodReps
  nReps    <- length(goodReps)

  # load in operating model to get weight and maturity at-age data
  load('history/fit_impPred_tvM/fit_impPred_tvM.RData')
  W_apt <- reports$data$W_apt
  meanWt_ap <- apply(X=W_apt, FUN=mean, MAR=c(1,2))
  W_iat <- array(NA, dim=c(nReps, nA,nT))
  W_iat[1,1:nA,1:(tMP-1)] <- W_apt[,pIdx,1:(tMP-1)]
  W_iat[1,1:nA,tMP:nT] <- meanWt_ap[,pIdx]

  mat_a <- reports$repOpt$mat_a
  mat_at <- array(NA, dim=c(nA,nT))
  mat_at[,1:nT] <-mat_a
  mat_iat <- array(NA, dim=c(nReps, nA,nT))
  mat_iat[1,,] <- mat_at

  if(nReps >1)
    for(i in 2:nReps)
    {
      W_iat[i,1:nA,1:nT] <- W_iat[1,1:nA,]
      mat_iat[i,1:nA,1:nT] <- mat_iat[1,1:nA,]
    }


  # pred gears
  gears         <- 1:nG
  predG         <- blob$ctlList$opMod$predGear
  fishG         <- gears[!(gears %in% predG)]
  sokFleets     <- which(blob$om$fleetType_f %in% c(2,3))
  nonSokFleets  <- which(blob$om$fleetType_f %in% c(0,1))
  sokPredG      <- predG[which(predG %in% sokFleets)]
  nonSokPredG   <- predG[which(!predG %in% sokFleets)]
  sokFishG      <- fishG[which(fishG %in% sokFleets)]
  nonSokFishG   <- fishG[which(!fishG %in% sokFleets)]

  # SOK conversion factor to converts dead ponded fish (kt) to egg yield (numbers)
  gamma_f <- reports$repOpt$gamma_g # Egg yield (numbers) = Catch (kt) / gamma


  # Arrays to hold empirical eqbm yields
  mnC_spk  <- array( NA, dim = c(nS,nP,nSims) ) # mean Fisheries/survey catch at eqlbm
  mdC_spk    <- array( NA, dim = c(nS,nP,nSims) ) # median Fisheries/survey catch at eqlbm

  mnP_spk  <- array( NA, dim = c(nS,nP,nSims) ) # mean ponded fish at eqlbm
  mdP_spk    <- array( NA, dim = c(nS,nP,nSims) ) # median ponded fish at eqlbm

  mnCp_spk <- array( NA, dim = c(nS,nP,nSims) ) # mean Predator consumption at eqlbm
  mdCp_spk   <- array( NA, dim = c(nS,nP,nSims) ) # mean Predator consumption at eqlbm

  mdSB_spk <- array( NA, dim = c(nS,nP,nSims) ) # median SB
  mnSB_spk <- array( NA, dim = c(nS,nP,nSims) ) # mean SB

  mdB_spk <- array( NA, dim = c(nS,nP,nSims) ) # median biomass
  mnB_spk <- array( NA, dim = c(nS,nP,nSims) ) # mean biomass

  mdStartSB_k <- rep(NA, nSims) # median Spawning biomass at start of year
  mnStartSB_k <- rep(NA, nSims)  #mean Spawning biomass at start of year
  mdStartB_k <- rep(NA, nSims) # median biomass at start of year
  mnStartB_k <- rep(NA, nSims)  #mean biomass at start of year

  mdVB_spfk <- array( NA, dim = c(nS,nP,nF,nSims) ) # median vulnerable biomass
  mnVB_spfk <- array( NA, dim = c(nS,nP,nF,nSims) ) # mean vulnerable biomass

  mdU1_spk <- array( NA, dim = c(nS,nP,nSims) ) # median HR using start of year SB
  mnU1_spk <- array( NA, dim = c(nS,nP,nSims) ) # mean HR using start of year SB
  mdU2_spk <- array( NA, dim = c(nS,nP,nSims) ) # median HR using start of year SB
  mnU2_spk <- array( NA, dim = c(nS,nP,nSims) ) # mean HR using end of year SB

  F_spk <- array( NA, dim = c(nS,nP,nSims) )
  E_spk <- array( NA, dim = c(nS,nP,nSims) )
  M_spk <- array( NA, dim = c(nS,nP,nSims) )
  R_spk <- array( NA, dim = c(nS,nP,nSims) )
  E_pk  <- array( NA, dim = c(nP,nSims) )

  # Biomass U, and catch time series
  mdSB_sptk <- array( NA, dim = c(nS,nP,nT,nSims) ) # median biomass
  mnSB_sptk <- array( NA, dim = c(nS,nP,nT,nSims) ) # mean biomass
  mdVB_spftk <- array( NA, dim = c(nS,nP,nF,nT,nSims) ) # median vulnerable biomass
  mnVB_spftk <- array( NA, dim = c(nS,nP,nF,nT,nSims) ) # mean vulnerable biomass

  mnC_sptk  <- array( NA, dim = c(nS,nP,nT,nSims) ) # mean Fisheries/survey catch at eqlbm
  mdC_sptk    <- array( NA, dim = c(nS,nP,nT,nSims) ) # median Fisheries/survey catch at eqlbm
  mnCp_sptk  <- array( NA, dim = c(nS,nP,nT,nSims) ) # mean predator catch at eqlbm
  mdCp_sptk    <- array( NA, dim = c(nS,nP,nT,nSims) ) # median predator catch at eqlbm
  mdU2_sptk <- array( NA, dim = c(nS,nP,nT,nSims) ) # median HR using end of year SB
  mnU2_sptk <- array( NA, dim = c(nS,nP,nT,nSims) ) # mean HR using end of year SB

  mdAlloc_spfk <- array( NA, dim = c(nS,nP,length(fishG),nSims) ) # mean catch/ponding allocation among fishing fleets

  mpNames_k <- vector('character',nSims)
  inputF_k  <- vector('numeric',nSims)

  qF_sp <- blob$om$qF_ispft[1,,,2,nT]

  for( x in sims )
  {
    .loadSim(x, folder = folder)
    # Calculate spawning biomass at start of the year from Numbers at-age

    startB_iat <- array( 0, dim = c(length(goodReps),nA,nT))
    startSB_iat <- array( 0, dim = c(length(goodReps),nA,nT))

    startB_iat[1:nReps,,]  <- blob$om$N_iaxspt[1:nReps,,1,sIdx,pIdx,] * W_iat[1:nReps,,]
    startSB_iat[1:nReps,,]  <- blob$om$N_iaxspt[,,1,sIdx,pIdx,] * W_iat[1:nReps,,] * mat_iat[1:nReps,,]
    startB_it  <- apply(X=startB_iat, FUN=sum, MAR=c(1,3) )
    startSB_it  <- apply(X=startSB_iat, FUN=sum, MAR=c(1,3) )

    mdStartB_k[x] <- median(startSB_it[goodReps,nT-0:period])
    mnStartB_k[x] <- mean(startSB_it[goodReps,nT-0:period])

    mdStartSB_k[x] <- median(startSB_it[goodReps,nT-0:period])
    mnStartSB_k[x] <- mean(startSB_it[goodReps,nT-0:period])

    C_ispft <- array( 0, dim = c(length(goodReps),nS,nP,nF, nT))
    C_ispft[,,,nonSokFleets,] <- blob$om$C_ispft[,,,nonSokFleets,]
    C_ispft[,,,sokFleets,]    <- blob$om$P_ispft[,,,sokFleets,]

    C_ispt <- apply(X = C_ispft[goodReps,,,fishG,,drop=FALSE], FUN = sum, MARGIN = c(1,2,3,5), na.rm = T)

    mdC_sptk[,,,x]  <- apply(X = C_ispt[goodReps,,,,drop=FALSE], FUN = median, MARGIN = c(2,3,4), na.rm = T )
    mnC_sptk[,,,x]  <- apply(X = C_ispt[goodReps,,,,drop=FALSE], FUN = mean, MARGIN = c(2,3,4), na.rm = T )

    mdC_spk[,,x]  <- apply(X = C_ispt[goodReps,,,nT-0:period,drop = FALSE], FUN = median, MARGIN = c(2,3), na.rm = T )
    mnC_spk[,,x]  <- apply(X = C_ispt[goodReps,,,nT-0:period,drop = FALSE], FUN = mean, MARGIN = c(2,3), na.rm = T )

    alloc_ispft <- array( 0, dim = c(length(goodReps),nS,nP,length(fishG), nT))
    for(fIdx in fishG)
        alloc_ispft[goodReps,1:nS,1:nP,fIdx,1:nT] <- C_ispft[goodReps,1:nS,1:nP,fIdx,1:nT]/C_ispt[goodReps,1:nS,1:nP,1:nT]

    mdAlloc_spft <- apply(X = alloc_ispft[goodReps,,,,,drop = FALSE], FUN = median, MARGIN = c(2,3,4,5), na.rm = T)

    mdAlloc_spfk[,,,x] <- apply(X = mdAlloc_spft[,,,nT-0:period,drop = FALSE], FUN = median, MARGIN = c(1,2,3), na.rm = T)

    if(!is.null(predG) & length(predG) > 0)
    {

      Cp_ispt  <- apply(X = C_ispft[goodReps,,,nonSokPredG,,drop=FALSE], FUN = sum, MARGIN = c(1,2,3,5), na.rm = T)

      mdCp_spk[,,x]  <- apply(X = Cp_ispt[goodReps,,,nT-0:period,drop = FALSE], FUN = median, MARGIN = c(2,3), na.rm = T )

      mnCp_spk[,,x]  <- apply(X = Cp_ispt[goodReps,,,nT-0:period,drop = FALSE], FUN = mean, MARGIN = c(2,3), na.rm = T )

      mdCp_spk[,,x]  <- apply(X = Cp_ispt[goodReps,,,nT-0:period,drop = FALSE], FUN = median, MARGIN = c(2,3), na.rm = T )
      mnCp_spk[,,x]  <- apply(X = Cp_ispt[goodReps,,,nT-0:period,drop = FALSE], FUN = mean, MARGIN = c(2,3), na.rm = T )
    }

    # SB is spawning biomass at spawn timing before ponded fish are added back in
    # endSB is spawning biomass at spawn timing after ponded fish are added back in

    mdSB_spk[,,x]  <- apply(X = blob$om$SB_ispt[goodReps,,,nT-0:period,drop = FALSE], FUN = median, MARGIN = c(2,3), na.rm = T )
    mnSB_spk[,,x]  <- apply(X = blob$om$SB_ispt[goodReps,,,nT-0:period,drop = FALSE], FUN = mean, MARGIN = c(2,3), na.rm = T )

    mdSB_sptk[,,,x]  <- apply(X = blob$om$SB_ispt[goodReps,,,,drop = FALSE], FUN = median, MARGIN = c(2,3,4), na.rm = T )
    mnSB_sptk[,,,x]  <- apply(X = blob$om$SB_ispt[goodReps,,,,drop = FALSE], FUN = mean, MARGIN = c(2,3,4), na.rm = T )

    mdB_spk[,,x]  <- apply(X = blob$om$B_ispt[goodReps,,,nT-0:period,drop = FALSE], FUN = median, MARGIN = c(2,3), na.rm = T )
    mnB_spk[,,x]  <- apply(X = blob$om$B_ispt[goodReps,,,nT-0:period,drop = FALSE], FUN = mean, MARGIN = c(2,3), na.rm = T )

    vB_ispft <- blob$om$vB_ispf - blob$om$C_ispft

    mdVB_spfk[,,,x]  <- apply(X = vB_ispft[goodReps,,,,nT-0:period,drop = FALSE], FUN = median, MARGIN = c(2,3,4), na.rm = T )
    mnVB_spfk[,,,x]  <- apply(X = vB_ispft[goodReps,,,,nT-0:period,drop = FALSE], FUN = mean, MARGIN = c(2,3,4), na.rm = T )

    mdU1_spk[,,x]  <- mdC_spk[,,x]/(mdStartSB_k[x]) # HR using start of year SB
    mnU1_spk[,,x]  <- mnC_spk[,,x]/(mnStartSB_k[x])

    mdU2_spk[,,x]  <- mdC_spk[,,x]/(mdC_spk[,,x] + mdSB_spk[,,x]) # HR using end of year SB + Catch for denominator
    mnU2_spk[,,x]  <- mnC_spk[,,x]/(mnC_spk[,,x] + mnSB_spk[,,x])

    mdU2_sptk[,,,x]  <- mdC_sptk[,,,x]/(mdC_sptk[,,,x] + mdSB_sptk[,,,x]) # HR using end of year SB + Catch for denominator
    mnU2_sptk[,,,x]  <- mnC_sptk[,,,x]/(mnC_sptk[,,,x] + mnSB_sptk[,,,x])

    # mdU2_spk[,,x]  <- mdC_spk[,,x]/(mdSB_spk[,,x]) # HR using end of year SB
    # mnU2_spk[,,x]  <- mnC_spk[,,x]/(mnSB_spk[,,x])

    F_spk[,,x]  <- apply(X = blob$om$F_ispft[1,,,fishG,nT,drop = FALSE], FUN = sum, MARGIN = c(2,3))
    M_spk[,,x]  <- apply(X = blob$om$M_iaxspt[goodReps,2,1,,,nT-0:period,drop = FALSE], FUN = median, MARGIN = c(2,3), na.rm = T )
    R_spk[,,x]  <- apply(X = blob$om$R_ispt[goodReps,,,nT-0:period,drop = FALSE], FUN = median, MARGIN = c(2,3), na.rm = T )
    E_pk[,x]    <- apply(X = blob$om$E_ipft[goodReps,,2,nT-0:period,drop = FALSE], FUN = median, MARGIN = c(2), na.rm = T )
    E_spk[,,x]  <- F_spk[,,x] / qF_sp

    mpNames_k[x]   <- blob$ctlList$ctl$mpName
    inputF_k[x]    <- blob$ctlList$mp$hcr$inputF

    # Clean up
    gc()
  }



    saveEmpRefCurves <- list( mdC_spk  = mdC_spk,
                              mdCp_spk = mdCp_spk,
                              mnC_spk  = mnC_spk,
                              mnCp_spk = mnCp_spk,
                              mdSB_spk  = mdSB_spk,
                              mnSB_spk  = mnSB_spk,
                              mdB_spk  = mdB_spk,
                              mnB_spk  = mnB_spk,
                              mdStartSB_k = mdStartSB_k,
                              mnStartSB_k = mnStartSB_k,
                              mdVB_spfk  = mdVB_spfk,
                              mnVB_spfk  = mnVB_spfk,
                              mdU1_spk  = mdU1_spk,
                              mnU1_spk  = mnU1_spk,
                              mdU2_spk  = mdU2_spk,
                              mnU2_spk  = mnU2_spk,
                              F_spk  = F_spk,
                              E_spk  = E_spk,
                              M_spk  = M_spk,
                              R_spk  = R_spk,
                              E_pk   = E_pk,
                              mdAlloc_spfk = mdAlloc_spfk,
                              mpNames_k = mpNames_k,
                              inputF_k = inputF_k,
                              mdSB_sptk = mdSB_sptk,
                              mdVB_spftk = mdVB_spftk,
                              mdU2_sptk = mdU2_sptk)

    save( saveEmpRefCurves, file = here::here("Outputs",folder,"empRefCurves.RData"))

    # sim grid with time series columns and U, SB (end of yr), VB (Seine Roe), VB (gillnet) for rows, ordered by input F for s=1, p=1

    actualOrder <- order(inputF_k)

    Ugrid_kt  <- mdU2_sptk[sIdx, pIdx, 73:nT, actualOrder] %>% t() %>% data.frame()
    row.names(Ugrid_kt) <- inputF_k[actualOrder]
    SBgrid_kt <- mdSB_sptk[sIdx, pIdx, 73:nT, actualOrder] %>% t() %>% data.frame()
    row.names(SBgrid_kt) <- inputF_k[actualOrder]
    VB2grid_kt <- mdVB_spftk[sIdx, pIdx,2, 73:nT, actualOrder] %>% t() %>% data.frame()
    row.names(VB2grid_kt) <- inputF_k[actualOrder]
    VB3grid_kt <- mdVB_spftk[sIdx, pIdx,3, 73:nT, actualOrder] %>% t() %>% data.frame()
    row.names(VB3grid_kt) <- inputF_k[actualOrder]

    write.csv(Ugrid_kt, here::here("Outputs",folder,"Ugrid_kt.csv"))
    write.csv(SBgrid_kt, here::here("Outputs",folder,"SBgrid_kt.csv"))
    write.csv(VB2grid_kt, here::here("Outputs",folder,"VB2grid_kt.csv"))
    write.csv(VB3grid_kt, here::here("Outputs",folder,"VB3grid_kt.csv"))

} # END saveBatchGridData()

# plotEmpYieldCurves
# Function to plot median simuated yield
# resulting from a grid of constant fishing
# mortality rates - used for checking the
# reference point calculations
plotEmpYieldCurves <- function( sims = 1:200, nHRs=NULL,
                                folder = "expPred_conM_Fgrid_noProcErr",
                                indepVar = "U",
                                maxXspline = 0.5,
                                xLim = c(0,2),
                                plotVB=FALSE,
                                plotYPRcurves = FALSE,
                                redoEmpRefCurves = FALSE,
                                period = 50, sIdx=1, pIdx=1,
                                curveType='mean',
                                hrType = 'endOfYrSB')
{
  nSims <- length(sims)
  blobList <- vector(mode = "list", length = nSims)

  .loadSim(sims[1], folder = folder)

  nS <- blob$om$nS
  nP <- blob$om$nP
  nF <- blob$om$nF
  nT <- blob$om$nT
  nG <- blob$om$nF
  nA <- blob$om$nA

  tMP <- blob$om$tMP

  goodReps <- blob$goodReps
  nReps    <- length(goodReps)

  # pred gears
  gears <- 1:nG
  predG <- blob$ctlList$opMod$predGear
  fishG <- gears[!(gears %in% predG)]


  qF_sp <- blob$om$qF_ispft[1,,,2,nT]

  if(!file.exists(here::here("Outputs",folder,"empRefCurves.RData")) |
      redoEmpRefCurves )
  {
    saveBatchGridData(folder = batchFolder)
  } else {
    # Load ref curve list
    load(here::here("Outputs",folder,"empRefCurves.RData"))

    mdC_spk  <- saveEmpRefCurves$mdC_spk
    mdCp_spk <- saveEmpRefCurves$mdCp_spk
    mdSB_spk  <- saveEmpRefCurves$mdSB_spk
    mdB_spk  <- saveEmpRefCurves$mdB_spk
    mdStartSB_k  <- saveEmpRefCurves$mdStartSB_k
    mdVB_spfk  <- saveEmpRefCurves$mdVB_spfk
    mdU1_spk  <- saveEmpRefCurves$mdU1_spk
    mdU2_spk  <- saveEmpRefCurves$mdU2_spk

    mnC_spk  <- saveEmpRefCurves$mnC_spk
    mnCp_spk <- saveEmpRefCurves$mnCp_spk
    mnSB_spk  <- saveEmpRefCurves$mnSB_spk
    mnB_spk  <- saveEmpRefCurves$mnB_spk
    mnStartSB_k  <- saveEmpRefCurves$mnStartSB_k
    mnVB_spfk  <- saveEmpRefCurves$mnVB_spfk
    mnU1_spk  <- saveEmpRefCurves$mnU1_spk
    mnU2_spk  <- saveEmpRefCurves$mnU2_spk

    F_spk  <- saveEmpRefCurves$F_spk
    E_spk  <- saveEmpRefCurves$E_spk
    E_pk   <- saveEmpRefCurves$E_pk

    inputF_k  <- saveEmpRefCurves$inputF_k
    mpNames_k <- saveEmpRefCurves$mpNames_k

  }

  # Pull F based ref curves from RP object
  refCurves       <- blob$rp[[1]]$refCurves
  Fvec            <- refCurves$F
  Uvec            <- refCurves$Ueq_spf
  BeqRefCurve_spf <- refCurves$Beq_spf
  YeqRefCurve_spf <- refCurves$Yeq_spf

  # Pull effort based ref points
  EmsyRefPts      <- blob$rp[[1]]$EmsyRefPts
  Evec            <- refCurves$EffCurves$E

  YeqRefCurve_spe <- refCurves$EffCurves$Yeq_spe
  BeqRefCurve_spe <- refCurves$EffCurves$Beq_spe

  Xmsy_sp <- array(0, dim = c(nS,nP))
  SBmsy_sp <- array(0, dim = c(nS,nP))
  MSY_sp <- array(0, dim = c(nS,nP))

  if(is.null(pIdx))
    pIdx <- 1:nP

  if(is.null(sIdx))
    sIdx <- 1:nS

  par(  mfcol = c(length(pIdx),length(sIdx)),
        mar = c(1,1.5,1,1.5),
        oma = c(5,5,3,3) )

  for( s in sIdx )
    for( p in pIdx )
    {
      if(curveType=='median')
      {
        C   <- mdC_spk[sIdx,pIdx,] # fisheries/survey catch
        Cp  <- mdCp_spk[sIdx,pIdx,] # predator catch
        B   <- mdB_spk[sIdx,pIdx,]
        SB1 <- mdStartSB_k
        SB2 <- mdSB_spk[sIdx,pIdx,]
        VB2 <- mdVB_spfk[sIdx,pIdx,2,]
        VB3 <- mdVB_spfk[sIdx,pIdx,3,]
        VB6 <- mdVB_spfk[sIdx,pIdx,6,]
        F   <- F_spk[sIdx,pIdx,]
        E   <- E_spk[sIdx,pIdx,]

        if(hrType == 'startOfYrSB')
          U   <- mdU1_spk[sIdx,pIdx,]

        if(hrType == 'endOfYrSB')
          U   <- mdU2_spk[sIdx,pIdx,]

        if(hrType == 'startOfYrB')
          U   <- C/B

        if(hrType == 'VB_SR')
          U   <- C/(VB2+C)

      }

      if(curveType=='mean')
      {
        C   <- mnC_spk[sIdx,pIdx,] # fisheries/survey catch
        Cp  <- mnCp_spk[sIdx,pIdx,] # predator catch
        B   <- mnB_spk[sIdx,pIdx,]
        SB1 <- mnStartSB_k
        SB2 <- mnSB_spk[sIdx,pIdx,]
        VB2 <- mnVB_spfk[sIdx,pIdx,2,]
        VB3 <- mnVB_spfk[sIdx,pIdx,3,]
        VB6 <- mdVB_spfk[sIdx,pIdx,6,]
        F   <- F_spk[sIdx,pIdx,]
        E   <- E_spk[sIdx,pIdx,]

        if(hrType == 'startOfYrSB')
          U   <- mnU1_spk[sIdx,pIdx,]

        if(hrType == 'endOfYrSB')
          U   <- mnU2_spk[sIdx,pIdx,]

        if(hrType == 'startOfYrB')
          U   <- C/B

        if(hrType == 'VB_SR')
          U   <- C/(VB2+C)
      }

      if( indepVar == "F" )
        maxX <- max(F)

      if( indepVar == "E" )
        maxX <- max(E_pk[p,],E_spk[sIdx,pIdx,])

      if( indepVar == "U" )
        maxX <- max(U)

      if(xLim[2] < maxX )
        maxX <- xLim[2]

        if( indepVar == "F" )
          idxOrder <- order(F)

        if( indepVar == "U" )
          idxOrder <- order(U)

        idxOrder <- order(inputF_k)
        # inputF_k <- inputF_k[idxOrder]

        if(!is.null(nHRs))
          idxOrder <- idxOrder[1:nHRs]

        if(any(U<0))
          idxOrder <- which(U[idxOrder]>0)

        C   <- C[idxOrder]
        Cp  <- Cp[idxOrder]
        B   <- B[idxOrder]
        SB1 <- SB1[idxOrder]
        SB2 <- SB2[idxOrder]
        VB2 <- VB2[idxOrder]
        VB3 <- VB3[idxOrder]
        U   <- U[idxOrder]
        F   <- F[idxOrder]
        E   <- E[idxOrder]

        if( indepVar == "F" )
        {
          X     <- F
          Xvec  <- Fvec
          Yeq   <- YeqRefCurve_spf[sIdx,pIdx,]
          Beq   <- BeqRefCurve_spf[sIdx,pIdx,]

          xLab  <- "Fishing Mortality"
        }

        if( indepVar == "U" )
        {
          X     <- U
          Xvec  <- Uvec
          Yeq   <- YeqRefCurve_spf[sIdx,pIdx,]
          Beq   <- BeqRefCurve_spf[sIdx,pIdx,]


          if(hrType == 'endOfYrSB')
            xLab  <- "Harvest rate C/(C+SB) (Day 365)"

          if(hrType == 'VB_SR')
            xLab  <- "Harvest rate C/VB (SR)"

          if(hrType == 'startOfYrB')
            xLab  <- "Harvest rate C/B (Day 1)"
        }

        if( indepVar == "E" )
        {

          X     <- E
          Xvec  <- Evec
          Yeq   <- YeqRefCurve_spe[sIdx,pIdx,]
          Beq   <- BeqRefCurve_spe[sIdx,pIdx,]

          xLab  <- "Fishing Effort"
        }

        ubX <- max(which(Yeq >= 0) )

        CXspline <- splinefun(x = X, y = C)
        BXspline <- splinefun(x = X, y = SB2)

        empXmsy  <- try(uniroot(  interval = c(0.01,maxXspline),
                              f = CXspline,
                              deriv = 1 )$root)
        if( class(empXmsy) == "try-error")
          empXmsy <- 0

        empMSY  <- CXspline(empXmsy)
        empSBmsy <- BXspline(empXmsy)

        empXmsy <- round(empXmsy,2)
        empSBmsy <- round(empSBmsy,1)
        empMSY  <- round(empMSY,1)

        plot( x = c(0,maxX), y = c(0,max(B)),
            type = "n", axes = FALSE,
            xlab = "", ylab = "" )
        axis( side = 1 )
        axis( side = 2, las = 1 )
        grid()
        box()

        lines( x = X, y = B, col = "grey50", lwd = 2, lty = 1 )
        lines( x = X, y = SB1, col = "#b2df8a", lwd = 2, lty = 1 )
        lines( x = X, y = SB2, col = "#33a02c", lwd = 2, lty = 1 )
        lines( x = X, y = C, col = "#1f78b4", lwd = 2, lty = 1 )

        legNames <- c('Biomass (Day 1)',
                      'SB (Day 1)',
                      'SB (Day 365)',
                      'VB (Seine Roe)',
                      'VB (Gillnet)',
                      'Fisheries catch and ponded fish',
                      'Predator consumption')
        colrs <- c('grey50', '#b2df8a', '#33a02c',"#6a3d9a", "#cab2d6",  '#1f78b4', '#a6cee3')

        if(plotYPRcurves)
        {
          lines( x = Xvec, y = Yeq,
                  col = "salmon", lty = 2 )
          lines( x = Xvec, y = Beq,
                  col = "black", lty = 2 )
        }

        if(is.null(predG) & length(predG) == 0)
        {

          if(plotVB)
          {
            lines( x = X, y = VB2, col = "#6a3d9a", lwd = 2, lty = 1 )
            lines( x = X, y = VB3, col = "#cab2d6", lwd = 2, lty = 1 )
            legend( "topright", bty = "n",
                    legend = legNames[1:6],
                    lty=1, lwd=2,
                    col=colrs[1:6])
          }

          if(!plotVB)
          {
            legend( "topright", bty = "n",
                    legend = legNames[c(1:3,6)],
                    lty=1, lwd=2,
                    col=colrs[c(1:3,6)])
          }
        }

        if(!is.null(predG) & length(predG) > 0)
        {

          lines( x = X, y = Cp, col = "#a6cee3", lwd = 2, lty = 1 )

          if(plotVB)
          {
            lines( x = X, y = VB2, col = "#6a3d9a", lwd = 2, lty = 1 )
            lines( x = X, y = VB3, col = "#cab2d6", lwd = 2, lty = 1 )
            legend( "topright", bty = "n",
                    legend = legNames[1:7],
                    lty=1, lwd=2,
                    col=colrs[1:7])
          }
s
          if(!plotVB)
          {
            legend( "topright", bty = "n",
                    legend = legNames[c(1:3,6,7)],
                    lty=1, lwd=2,
                    col=colrs[c(1:3,6,7)])
          }
        }

        legend( x = maxX*0.6, y = max(B)*0.7, bty = "n", cex=0.9,
                legend = c( paste("SB0 (Day 1) = ", round(SB1[1],1), sep = "" ),
                            paste("SB0 (Day 365) = ", round(SB2[1],1), sep = "" ),
                            paste("Umsy = ", empXmsy, sep = "" ),
                            paste("MSY = ", empMSY, sep = ""),
                            paste("SBmsy = ", empSBmsy, sep = "") ) )


        # Plot some guidelines
        segments(  x0 = empXmsy, x1 = empXmsy,
                    y0 = 0, y1 = empSBmsy,
                    lty = 2, col = "red" )
        segments(  x0 = 0, x1 = empXmsy,
                    y0 = empMSY, y1 = empMSY,
                    lty = 2, col = "red" )
        segments(  x0 = 0, x1 = empXmsy,
                    y0 = empSBmsy, y1 = empSBmsy,
                    lty = 2, col = "red" )

        Xmsy_sp[s,p] <- empXmsy
        SBmsy_sp[s,p] <- empSBmsy
        MSY_sp[s,p]  <- empMSY
    }

  mtext( side = 1, text = xLab, outer = T, line = 2)
  mtext( side = 2, text = "Equilibrium biomass and catch (kt)", outer = T, line = 2 )

  out <- list(  Xmsy_sp = Xmsy_sp,
                SBmsy_sp = SBmsy_sp,
                MSY_sp  = MSY_sp )

  out
} # END plotEmpYieldCurves()

plotFgridDiag <- function(nHRs=NULL,
                                folder = "expPred_conM_Fgrid_noProcErr",
                                indepVar = "U",
                                sIdx=1, pIdx=1,
                                curveType='mean',
                                hrType = 'endOfYrSB')
{

  load(here::here("Outputs",folder,"empRefCurves.RData"))

  mdC_spk  <- saveEmpRefCurves$mdC_spk
  mdCp_spk <- saveEmpRefCurves$mdCp_spk
  mdSB_spk  <- saveEmpRefCurves$mdSB_spk
  mdB_spk  <- saveEmpRefCurves$mdB_spk
  mdStartSB_k  <- saveEmpRefCurves$mdStartSB_k
  mdU1_spk  <- saveEmpRefCurves$mdU1_spk
  mdU2_spk  <- saveEmpRefCurves$mdU2_spk

  mnC_spk  <- saveEmpRefCurves$mnC_spk
  mnCp_spk <- saveEmpRefCurves$mnCp_spk
  mnSB_spk  <- saveEmpRefCurves$mnSB_spk
  mnB_spk  <- saveEmpRefCurves$mnB_spk
  mnStartSB_k  <- saveEmpRefCurves$mnStartSB_k
  mnU1_spk  <- saveEmpRefCurves$mnU1_spk
  mnU2_spk  <- saveEmpRefCurves$mnU2_spk

  F_spk  <- saveEmpRefCurves$F_spk
  E_spk  <- saveEmpRefCurves$E_spk
  E_pk   <- saveEmpRefCurves$E_pk

  inputF_k  <- saveEmpRefCurves$inputF_k
  mpNames_k <- saveEmpRefCurves$mpNames_k

  nS    <- dim(mdC_spk)[1]
  nP    <- dim(mdC_spk)[2]
  nSims <- dim(mdC_spk)[3]

  if(is.null(pIdx))
    pIdx <- 1:nP

  if(is.null(sIdx))
    sIdx <- 1:nS

  if(curveType=='median')
  {
    C   <- mdC_spk[sIdx,pIdx,] # fisheries/survey catch
    Cp  <- mdCp_spk[sIdx,pIdx,] # predator catch
    B   <- mdB_spk[sIdx,pIdx,]
    SB1 <- mdStartSB_k
    SB2 <- mdSB_spk[sIdx,pIdx,]
    F   <- F_spk[sIdx,pIdx,]
    E   <- E_spk[sIdx,pIdx,]

    if(hrType == 'startOfYrSB')
      U   <- mdU1_spk[sIdx,pIdx,]

    if(hrType == 'endOfYrSB')
      U   <- mdU2_spk[sIdx,pIdx,]

  }

  if(curveType=='mean')
  {
    C   <- mnC_spk[sIdx,pIdx,] # fisheries/survey catch
    Cp  <- mnCp_spk[sIdx,pIdx,] # predator catch
    SB1 <- mnStartSB_k
    SB2 <- mnSB_spk[sIdx,pIdx,]
    B   <- mnB_spk[sIdx,pIdx,]
    F   <- F_spk[sIdx,pIdx,]
    E   <- E_spk[sIdx,pIdx,]

    if(hrType == 'startOfYrSB')
      U   <- mnU1_spk[sIdx,pIdx,]

    if(hrType == 'endOfYrSB')
      U   <- mnU2_spk[sIdx,pIdx,]
  }

    actualOrder <- order(inputF_k)
    inputF_k <- inputF_k[actualOrder]

    if(!is.null(nHRs))
      actualOrder <- actualOrder[1:nHRs]

    C   <- C[actualOrder]
    Cp  <- Cp[actualOrder]
    B   <- B[actualOrder]
    SB1 <- SB1[actualOrder]
    SB2 <- SB2[actualOrder]
    U   <- U[actualOrder]
    F   <- F[actualOrder]
    E   <- E[actualOrder]

    par(mfrow=c(2,2), mar=c(3,3,1,1), mgp=c(1.8,0.6,0), tck=-0.02)

    plot( x = inputF_k , y = C,
          xlab='Input F', ylab = 'Catch')
    plot( x = inputF_k , y = U, xlab='Input F', ylab = 'HR', las=1)
    plot( x = U , y = C, xlab='HR', ylab = 'Catch', las=1)
    plot( x = inputF_k , y = SB2, las=1, ylim=c(0, max( SB2)),
          xlab='Input F', ylab = 'SB (Day 365)', )

}

# plotEmpYieldCurves
# Function to plot median simuated yield
# resulting from a grid of constant fishing
# mortality rates - used for checking the
# reference point calculations
plot2dimEmpYieldCurves <- function( sims = 1:100,
                                    folder = "CC_refPts",
                                    indepVar = "F",
                                    plotYPRcurves = TRUE,
                                    redoEmpRefCurves = FALSE )
{
  nSims <- length(sims)
  blobList <- vector(mode = "list", length = nSims)

  .loadSim(sims[1], folder = folder)

  nS <- blob$om$nS
  nP <- blob$om$nP
  nT <- blob$om$nT
  tMP <- blob$om$tMP

  goodReps <- blob$goodReps


  # Arrays to hold empirical eqbm yields
  C_spk     <- array( 0, dim = c(nS,nP,nSims) )
  B_spk     <- array( 0, dim = c(nS,nP,nSims) )
  initB_spk <- array( 0, dim = c(nS,nP,nSims) )
  F_spk     <- array( 0, dim = c(nS,nP,nSims) )
  E_spk     <- array( 0, dim = c(nS,nP,nSims) )
  M_spk     <- array( 0, dim = c(nS,nP,nSims) )
  E_pk      <- array( 0, dim = c(nP,nSims) )

  if(!file.exists(here::here("Outputs",folder,"empRefCurves.RData")) |
      redoEmpRefCurves )
  {
    for( x in sims )
    {
      .loadSim(x, folder = folder)
      C_spk[,,x]  <- apply(X = blob$om$C_ispt[goodReps,,,nT,drop = FALSE], FUN = median, MARGIN = c(2,3), na.rm = T )
      B_spk[,,x]  <- apply(X = blob$om$SB_ispt[goodReps,,,nT,drop = FALSE], FUN = median, MARGIN = c(2,3), na.rm = T )
      initB_spk[,,x]  <- apply(X = blob$om$SB_ispt[goodReps,,,tMP-1,drop = FALSE], FUN = median, MARGIN = c(2,3), na.rm = T )
      F_spk[,,x]  <- apply(X = blob$om$F_ispft[goodReps,,,2,nT,drop = FALSE], FUN = median, MARGIN = c(2,3), na.rm = T )
      M_spk[,,x]  <- apply(X = blob$om$M_iaxspt[goodReps,2,1,,,nT,drop = FALSE], FUN = median, MARGIN = c(2,3), na.rm = T )
      E_pk[,x]    <- apply(X = blob$om$E_ipft[goodReps,,2,nT,drop = FALSE], FUN = median, MARGIN = c(2), na.rm = T )

      # Clean up
      gc()
    }

    saveEmpRefCurves <- list( C_spk = C_spk,
                              B_spk = B_spk,
                              initB_spk = initB_spk,
                              F_spk = F_spk,
                              E_spk = E_spk,
                              M_spk = M_spk,
                              E_pk  = E_pk )

    save( saveEmpRefCurves, file = here::here("Outputs",folder,"empRefCurves.RData"))
  } else {
    # Load ref curve list
    load(here::here("Outputs",folder,"empRefCurves.RData"))



    C_spk <- saveEmpRefCurves$C_spk
    B_spk <- saveEmpRefCurves$B_spk
    initB_spk <- saveEmpRefCurves$initB_spk
    F_spk <- saveEmpRefCurves$F_spk
    E_spk <- saveEmpRefCurves$E_spk
    E_pk  <- saveEmpRefCurves$E_pk
  }



  # Pull F based ref curves from RP object
  refCurves       <- blob$rp[[1]]$refCurves
  Fvec            <- refCurves$F
  BeqRefCurve_spf <- refCurves$Beq_spf
  YeqRefCurve_spf <- refCurves$Yeq_spf

  # Pull effort based ref points
  EmsyRefPts      <- blob$rp[[1]]$EmsyRefPts
  Evec            <- refCurves$EffCurves$E

  YeqRefCurve_spe <- refCurves$EffCurves$Yeq_spe
  BeqRefCurve_spe <- refCurves$EffCurves$Beq_spe

  # There should only be a few initB (based on histFmult levels)
  initBvals_spc <- aperm(apply(X = initB_spk, FUN = unique, MARGIN = c(1,2) ),c(2,3,1))

  nCurves <- dim(initBvals_spc)[3]

  Xmsy_spc <- array(0, dim = c(nS,nP,nCurves))
  Bmsy_spc <- array(0, dim = c(nS,nP,nCurves))
  MSY_spc <- array(0, dim = c(nS,nP,nCurves))


  par(  mfcol = c(nP,nS),
        mar = c(1,1.5,1,1.5),
        oma = c(5,5,3,3) )

  for( s in 1:nS )
    for( p in 1:nP )
    {
      if( indepVar == "F" )
        maxX <- max(F_spk[s,p,])
      if( indepVar == "E" )
        maxX <- max(E_pk[p,],E_spk[s,p,])

      plot( x = c(0,maxX), y = c(0,max(B_spk[s,p,])),
            type = "n", axes = FALSE,
            xlab = "", ylab = "" )
        axis( side = 1 )
        axis( side = 2, las = 1 )
        grid()
        box()

        if( indepVar == "F" )
        {
          X     <- F
          Xvec  <- Fvec
          Yeq   <- YeqRefCurve_spf[s,p,]
          Beq   <- BeqRefCurve_spf[s,p,]

          xLab  <- "Fishing Mortality"
        }

        if( indepVar == "E" )
        {

          X     <- E
          Xvec  <- Evec
          Yeq   <- YeqRefCurve_spe[s,p,]
          Beq   <- BeqRefCurve_spe[s,p,]

          xLab  <- "Fishing Effort"
        }

        ubX <- max(which(Yeq >= 0) )

        # Loop over histFmult factors
        for( j in 1:nCurves)
        {
          initB <- initBvals_spc[s,p,j]

          simIdx <- which(initB_spk[s,p,] == initB)

          F <- F_spk[s,p,simIdx]
          C <- C_spk[s,p,simIdx]
          B <- B_spk[s,p,simIdx]
          E <- E_spk[s,p,simIdx]

          actualOrder <- order(F)

          C <- C[actualOrder]
          B <- B[actualOrder]
          F <- F[actualOrder]
          E <- E[actualOrder]

          if(indepVar == "F")
            X <- F
          if(indepVar == "E")
            X <- E

          CXspline <- splinefun(x = X, y = C)
          BXspline <- splinefun(x = X, y = B)

          empXmsy  <- try(uniroot(  interval = range(X),
                                f = CXspline,
                                deriv = 1 )$root)
          if( class(empXmsy) == "try-error")
            empXmsy <- 0

          empMSY  <- CXspline(empXmsy)
          empBmsy <- BXspline(empXmsy)

          empXmsy <- round(empXmsy,2)
          empBmsy <- round(empBmsy,2)
          empMSY  <- round(empMSY,2)




            lines( x = X, y = C,
                    col = "steelblue", lwd = 2, lty = 1 )
            lines( x = X, y = B,
                    col = "black", lwd = 2, lty = 1 )


          legend( "topright", bty = "n",
                  legend = c( paste(indepVar,"msy = ", empXmsy, sep = "" ),
                              paste(" MSY = ", empMSY, sep = ""),
                              paste("Bmsy = ", empBmsy, sep = "") ) )


          # Plot some guidelines
          segments(  x0 = empXmsy, x1 = empXmsy,
                      y0 = 0, y1 = empBmsy,
                      lty = 2, col = "red" )
          segments(  x0 = 0, x1 = empXmsy,
                      y0 = empMSY, y1 = empMSY,
                      lty = 2, col = "red" )
          segments(  x0 = 0, x1 = empXmsy,
                      y0 = empBmsy, y1 = empBmsy,
                      lty = 2, col = "red" )

          Xmsy_spc[s,p,j] <- empXmsy
          Bmsy_spc[s,p,j] <- empBmsy
          MSY_spc[s,p,j]  <- empMSY
        }
        if(plotYPRcurves)
        {
          lines( x = Xvec, y = Yeq,
                  col = "salmon", lty = 2 )
          lines( x = Xvec, y = Beq,
                  col = "black", lty = 2 )
        }
    }

  mtext( side = 1, text = xLab, outer = T, line = 2)
  mtext( side = 2, text = "Eqbm biomass and catch (kt)", outer = T, line = 2 )



  out <- list(  Xmsy_spc  = Xmsy_spc,
                Bmsy_spc  = Bmsy_spc,
                MSY_spc   = MSY_spc,
                initB_spc = initBvals_spc )

  out
} # END plotEmpYieldCurves()


# compareRefPts()
# Compares OM and conditioning fit reference points
# for a sanity check
compareYieldCurves <- function( obj )
{
  refPtsOM <- obj$rp[[1]]
  refPtsAM <- obj$ctlList$opMod$histRpt$refPts

  # Pull model dimensions
  nS <- obj$om$nS
  nP <- obj$om$nP

  speciesNames  <- obj$om$speciesNames
  stockNames    <- obj$om$stockNames

  # Now, loop and plot reference points
  par(  mfcol = c(nP,nS),
        mar = c(1,1.5,1,1.5),
        oma = c(5,5,3,3) )

  for( s in 1:nS )
    for( p in 1:nP )
    {
      refCurvesOM <- refPtsOM$refCurves
      refCurvesAM <- refPtsAM$refCurves
      plot( x = range( refCurvesOM$F ),
            y = c(0, max(refCurvesOM$Yeq_spf[s,p,], refCurvesOM$Beq_spf[s,p,]  )),
            type = "n", axes = FALSE )
        axis( side = 1 )
        axis( side = 2, las = 1 )

        lines(  x = refCurvesOM$F, y = refCurvesOM$Yeq_spf[s,p,],
                col = "black", lwd = 3 )
        lines(  x = refCurvesOM$F, y = refCurvesOM$Beq_spf[s,p,],
                col = "steelblue", lwd = 3 )

        lines(  x = refCurvesAM$F, y = refCurvesAM$Yeq_spf[s,p,],
                col = "black", lty = 2, lwd = 3 )
        lines(  x = refCurvesAM$F, y = refCurvesAM$Beq_spf[s,p,],
                col = "black", lty = 2, lwd = 3 )


        points( x = refPtsOM$FmsyRefPts$Fmsy_sp[s,p],
                y = refPtsOM$FmsyRefPts$YeqFmsy_sp[s,p],
                col = "red", pch = 16, cex = 1.5 )
        points( x = refPtsAM$FmsyRefPts$Fmsy_sp[s,p],
                y = refPtsAM$FmsyRefPts$YeqFmsy_sp[s,p],
                col = "steelblue", pch = 21, cex = 1.5 )

    }
}


# plotCvsB()
# Fishing mortality as a function of biomass,
# used for showing how well an MP meets the
# proposed HCR, and for determining the HCR
# implied by the omniscient manager sim
plotCvsB <- function( obj )
{
  goodReps <- obj$goodReps

  SB_ispt   <- obj$om$SB_ispt[goodReps,,,,drop = FALSE]
  C_ispt    <- obj$om$C_ispt[goodReps,,,,drop = FALSE]
  F_ispt    <- obj$om$F_ispft[goodReps,,,2,]
  vB_ispt   <- obj$om$vB_ispft[goodReps,,,2,]

  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP
  nT      <- obj$om$nT
  nReps   <- sum(goodReps)

  # Now, let's start plotting. We can fix it later
  speciesNames  <- (obj$om$speciesNames)
  stockNames    <- (obj$om$stockNames)


  projYrs <- tMP:nT

  dep_ispt <- SB_ispt
  Bmsy_sp  <- obj$rp[[1]]$FmsyRefPts$BeqFmsy_sp
  Fmsy_sp  <- obj$rp[[1]]$FmsyRefPts$Fmsy_sp
  MSY_sp   <- obj$rp[[1]]$FmsyRefPts$YeqFmsy_sp

  ctlList <- obj$ctlList

  par(  mfcol = c(nP,nS),
        mar = c(1,1.5,1,1.5),
        oma = c(5,5,3,3) )

  U_ispt <- C_ispt / vB_ispt

  for(s in 1:nS)
  {
    for( p in 1:nP )
    {
      # Convert biomass to Bmsy depletion
      dep_ispt[,s,p,] <- SB_ispt[,s,p,] / Bmsy_sp[s,p]

      # Depletion and F
      maxDep <- max(dep_ispt[,s,p,projYrs], na.rm = T)
      maxC   <- max(C_ispt[,s,p,projYrs]/MSY_sp[s,p], na.rm = T)

      # Plot window
      plot( x = c(0,3), y = c(0,2),
            type = "n", axes = F)
        mfg <- par("mfg")
        # axes
        axis( side = 1 )

        if( mfg[1] == 1 )
          mtext( side = 3, text = speciesNames[s], font = 2, line = 0.2, cex = 1.5 )

        if( mfg[2] == mfg[4] )
        {
          corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
          par(xpd = TRUE) #Draw outside plot area
          text( x = corners[2]+0.5,
                y = mean(corners[3:4]),
                labels = stockNames[p], srt = 270,
                font = 2, cex = 1.5 )
          par(xpd = FALSE)
        }

        axis( side = 2, las = 1 )
        box()

        ptCol <- scales::alpha("grey70", alpha = .3)

        points( x = dep_ispt[,s,p,projYrs], y = C_ispt[,s,p,projYrs]/MSY_sp[s,p],
                col = ptCol, pch = 1 )

        for( i in 1:nReps )
        {
          # Plot a smoother for each replicate
          lineCol <- scales::alpha("red", alpha = .3)
          smoother <- loess.smooth( x = dep_ispt[i,s,p,projYrs],
                                    y = C_ispt[i,s,p,projYrs]/MSY_sp[s,p] )
          lines( smoother, col = "grey30", lwd = .8 )
          flB <- dep_ispt[i,s,p,projYrs[c(1,length(projYrs))]]
          flC <- C_ispt[i,s,p,projYrs[c(1,length(projYrs))]]/MSY_sp[s,p]
          points( x = flB,
                  y = flC,
                  col = c("blue","red"), cex = .5,
                  pch = 16 )
        }


        abline( v = 1, lty = 2, lwd = .8)
        abline( h = 1, lty = 2, lwd = .8)

    }
  }

  mtext( side = 1, text = expression(B/B[MSY]), outer = TRUE, line = 2, font = 2)
  mtext( side = 2, text = expression(C/MSY), outer = TRUE, line = 2, font = 2)

}


# plotFvsB()
# Fishing mortality as a function of biomass,
# used for showing how well an MP meets the
# proposed HCR, and for determining the HCR
# implied by the omniscient manager sim
plotFvsB <- function( obj )
{
  goodReps <- obj$goodReps

  SB_ispt   <- obj$om$SB_ispt[goodReps,,,,drop = FALSE]
  C_ispt    <- obj$om$C_ispt[goodReps,,,,drop = FALSE]
  F_ispt    <- obj$om$F_ispft[goodReps,,,2,]
  vB_ispt   <- obj$om$vB_ispft[goodReps,,,2,]

  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP
  nT      <- obj$om$nT
  nReps   <- sum(goodReps)

  # Now, let's start plotting. We can fix it later
  speciesNames  <- (obj$om$speciesNames)
  stockNames    <- (obj$om$stockNames)


  projYrs <- tMP:nT

  dep_ispt <- SB_ispt
  Bmsy_sp  <- obj$rp[[1]]$FmsyRefPts$BeqFmsy_sp
  Fmsy_sp  <- obj$rp[[1]]$FmsyRefPts$Fmsy_sp

  ctlList <- obj$ctlList

  par(  mfcol = c(nP,nS),
        mar = c(1,1.5,1,2.5),
        oma = c(5,5,3,3) )

  U_ispt <- C_ispt / vB_ispt

  for(s in 1:nS)
  {
    for( p in 1:nP )
    {
      # Convert biomass to Bmsy depletion
      dep_ispt[,s,p,] <- SB_ispt[,s,p,] / Bmsy_sp[s,p]

      # Depletion and F
      maxDep <- max(dep_ispt[,s,p,projYrs], na.rm = T)
      maxF   <- max(F_ispt[,s,p,projYrs]/Fmsy_sp[s,p], na.rm = T)

      # Plot window
      plot( x = c(0,3), y = c(0,2),
            type = "n", axes = F)
        mfg <- par("mfg")
        # axes
        axis(side =1)

        if( mfg[1] == 1 )
          mtext( side = 3, text = speciesNames[s], font = 2, line = 0.2, cex = 1.5 )

        if( mfg[2] == mfg[4] )
        {
          corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
          par(xpd = TRUE) #Draw outside plot area
          text( x = corners[2]+0.5,
                y = mean(corners[3:4]),
                labels = stockNames[p], srt = 270,
                font = 2, cex = 1.5 )
          par(xpd = FALSE)
        }

        axis( side = 2, las = 1 )
        box()

        ptCol <- scales::alpha("grey70", alpha = .3)

        points( x = dep_ispt[,s,p,projYrs], y = F_ispt[,s,p,projYrs]/Fmsy_sp[s,p],
                col = ptCol, pch = 1 )

        for( i in 1:nReps )
        {
          # Plot a smoother for each replicate
          lineCol <- scales::alpha("red", alpha = .3)
          smoother <- loess.smooth( x = dep_ispt[i,s,p,projYrs],
                                    y = F_ispt[i,s,p,projYrs]/Fmsy_sp[s,p] )
          lines( smoother, col = "grey30", lwd = .8 )
          flB <- dep_ispt[i,s,p,projYrs[c(1,length(projYrs))]]
          flF <- F_ispt[i,s,p,projYrs[c(1,length(projYrs))]]/Fmsy_sp[s,p]
          points( x = flB,
                  y = flF,
                  col = c("blue","red"), cex = .5,
                  pch = 16 )
        }


        abline( v = 1, lty = 2, lwd = .8)
        abline( h = 1, lty = 2, lwd = .8)

    }
  }

  mtext( side = 1, text = expression(B/B[MSY]), outer = TRUE, line = 2, font = 2)
  mtext( side = 2, text = expression(F/F[MSY]), outer = TRUE, line = 2, font = 2)

}

# plotBatchPerf_sp()
# Plots multipanels (faceted by species/stock)
# of performance statistics on the y axis, with respect
# to the OM grid on the x axis
# CAUTION: currently only works for numeric xAxis and yAxis
plotBatchPerf_sp <- function( batchFolder = "fourthBatch",
                              xAxis = "projObsErrMult",
                              yAxis = "PBtGt.8Bmsy",
                              yRangeIn = NULL )
{
  # First load full stats table from the batch folder
  statsFile <- here::here(  "Outputs",batchFolder,"statistics",
                            "fullStatTable.csv")
  statTable <- read.csv(statsFile, header = TRUE, stringsAsFactors = FALSE )

  # Now, let's start plotting. We can fix it later
  speciesNames  <- unique(statTable$species)
  stockNames    <- unique(statTable$stock)

  scenarios <- unique(statTable$scenario)
  mps       <- unique(statTable$mp)

  nS <- length(speciesNames)
  nP <- length(stockNames)

  nScen <- length(scenarios)
  nMPs  <- length(mps)



  xLabs <- unique(statTable[,xAxis])
  xLabs <- xLabs[order(xLabs)]
  xMax  <- max(xLabs)
  xMin  <- min(xLabs)

  cols <- RColorBrewer::brewer.pal(n = nMPs, "Dark2")


  par(  mfcol = c( nP, nS ),
        mar = c( 0,2,0,2 ),
        oma = c(5,7,3,4) )

  yRange <- yRangeIn

  mpJitter <- seq( from = -.3, to = .3, length = nMPs )

  xDiff <- mean(diff(xLabs))
  mpJitter <- mpJitter * xDiff


  for( s in 1:nS )
    for( p in 1:nP )
    {
      subTable <- statTable %>%
                  filter( species == speciesNames[s],
                          stock == stockNames[p] )

      if(is.null(yRangeIn))
        yRange <- range(subTable[,yAxis])


      plot( x = c(-.5 + xMin,xMax + .5), y = yRange,
            type = "n", axes = FALSE )
        mfg <- par("mfg")
        # axes
        if( mfg[1] == mfg[3])
          axis( side = 1, at = xLabs )

        if( mfg[1] == 1 )
          mtext( side = 3, text = speciesNames[s], font = 2, line = 0 )

        axis( side = 2, las = 1 )

        box()
        # Add vertical grid lines to split OMs
        grid()

        # Now, plot the stuff
        for( m in 1:nMPs )
        {
          mpTable <- subTable %>% filter( mp == mps[m] )
          mpTable <- mpTable[order(mpTable[,xAxis]),]

          points( x = mpTable[,xAxis] + mpJitter[m],
                  y = mpTable[,yAxis],
                  col = cols[m], pch = 14 + m )
          lines(  x = mpTable[,xAxis] + mpJitter[m],
                  y = mpTable[,yAxis],
                  col = cols[m], lwd = .8 )
        }

        if( mfg[2] == mfg[4] )
        {
          corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
          par(xpd = TRUE) #Draw outside plot area
          text( x = corners[2]+0.2,
                y = mean(corners[3:4]),
                labels = stockNames[p], srt = 270,
                font = 2, cex = 1.5 )
          par(xpd = FALSE)
        }



    }
  legend( "topleft",
          legend = mps,
          col = cols,
          pch = 14 + 1:nMPs,
          lwd = 1 )



  mtext( side = 1, text = xAxis,
         outer = TRUE, cex = 1.5,
         font = 2, line = 3 )

  mtext( side = 2, text = yAxis,
         outer = TRUE, cex = 1.5, font = 2, line = 2 )



}



# plotHCR()
# Plots a generic hockey-stick
# harvest control rule, with given
# lower, upper control points, and
# low/high Fs. Stock status is calculated
# as a proportion of Bmsy, and fishing
# mortality rate as a proportion of
# Fmsy.
plotHCR <- function(  LCP = .4,
                      UCP = .6,
                      lowF = .1,
                      highF = 1 )
{
  x <- seq(0,1.3, length.out = 100 )
  y <- rep(lowF, length(x))


  par( mar = c(4,5,1,1), oma = c(2,2,1,1) )

  plot( x = range(x), y = c(0,1.2), type = "n",
        las = 1,
        xlab = expression(B/B[MSY]),
        ylab = expression(F/F[MSY]),
        cex.lab = 1.5,
        cex.axis = 1.5 )
    segments( x0 = 0, x1 = LCP,
              y0 = lowF, y1 = lowF,
              col = "grey50", lwd = 3 )
    segments( x0 = LCP, x1 = UCP,
              y0 = lowF, y1 = highF,
              col = "grey50", lwd = 3 )
    segments( x0 = UCP, x1 = max(x),
              y0 = highF, y1 = highF,
              col = "grey50", lwd = 3 )
    abline( v = c(LCP, UCP),
            col = c("red","darkgreen"),
            lty = 2, lwd = 2 )
    abline( h = highF, lty = 2, col = "grey70" )


} # END plotHCR()

plotUtilityFunction <- function( )
{
  x <- seq(1e-1,1, length = 1000)

  y <- (5*x - 1)/4/x

  plot( x = range(x), y = range(y),
        xlab = expression(C[spt]/TAC[spt]),
        ylab = "Utility", type = "n", las = 1 )
    mtext( side = 3, text = "TAC utilisation utility", font = 2)
    lines( x = x, y = y, lwd = 3, col = "grey40" )
    abline(h = c(0,1), lty = 2, col = "black" )
    abline( v = .2, lty = 2, col = "red", lwd = .8 )

}


# plotConvStats()
plotConvStats <- function( obj = blob )
{
  goodReps      <- obj$goodReps

  # Pull max gradient value and hessian indicator
  maxGrad_itsp  <- obj$mp$assess$maxGrad_itsp[goodReps,,,,drop = FALSE]
  pdHess_itsp   <- obj$mp$assess$pdHess_itsp[goodReps,,,,drop = FALSE]

  nReps <- dim(maxGrad_itsp)[1]

  # Now we want to get the mean and SD
  # of these values over the replicates
  quantsMaxGrad_qtsp <- apply(  X = maxGrad_itsp, FUN = quantile,
                                MARGIN = 2:4, probs = c(0.05, 0.5, 0.95),
                                na.rm = T )

  propPDHess_tsp  <- apply( X = pdHess_itsp, FUN = mean,
                            MARGIN = 2:4, na.rm = T )

  pT <- obj$ctlList$opMod$pT
  nS <- obj$om$nS
  nP <- obj$om$nP

  speciesNames  <- obj$om$speciesNames
  stockNames    <- obj$om$stockNames

  xJitter     <- seq( from = -.3, to = .3, length.out = nP )
  stockCols   <- RColorBrewer::brewer.pal( nP, "Dark2" )
  stockPts    <- seq( from = 21, by = 1, length.out = nP )
  rectWidth   <- .6 / nP

  par( mfcol = c(nS,2), mar = c(0,2,0,1), oma = c(3,3,3,2) )

  # First, plot the maxGrads, using
  # different colours for each species, and
  # pch for each stock... Jitter!
  for( s in 1:nS )
  {
    plot( x = c(1,pT),
          y = range(quantsMaxGrad_qtsp, na.rm = T),
          type = "n",
          axes = FALSE )
      mfg <- par("mfg")
      if( mfg[1] == mfg[3])
        axis( side = 1 )
      if(mfg[1] == 1 )
        mtext( side = 3, text = "Max Gradient Component")

      axis( side = 2, las = 1 )
      box()
      grid()

      for( p in 1:nP )
      {
        points( x = 1:pT + xJitter[p], y = quantsMaxGrad_qtsp[2,,s,p],
                col = stockCols[p], pch = stockPts[p], bg = stockCols[p] )
        segments( x0 = 1:pT + xJitter[p], x1 = 1:pT + xJitter[p],
                  y0 = quantsMaxGrad_qtsp[1,,s,p],
                  y1 = quantsMaxGrad_qtsp[3,,s,p],
                  col = stockCols[p], lty = 1 )
      }
      if( s == 1 & !is.null(stockNames) )
        legend( "topright", bty = "n",
                legend = stockNames,
                col = stockCols,
                pch = stockPts,
                pt.bg = stockCols )
  }

  # now do proportion of PD hessians
  for( s in 1:nS )
  {
    plot( x = c(0,pT+1), y = c(0,1.3),
          axes = FALSE, type = "n" )
      # Axes
      if( mfg[1] == mfg[3])
        axis( side = 1 )
      axis( side = 2, las = 1 )
      mtext( side = 4, text = speciesNames[s], font = 2, line = 2)
      if(mfg[1] == 1 )
        mtext( side = 3, text = "Proportion of PD Hessians")

      box()
      abline( h = 1.0, lty = 2, lwd = 2, col = "grey40" )

      # Plot rectangles of PD Hessians
      for( p in 1:nP )
      {
        rect( xleft = 1:pT + xJitter[p] - rectWidth/2,
              xright = 1:pT + xJitter[p] + rectWidth/2,
              ybottom = 0,
              ytop = propPDHess_tsp[,s,p],
              col = stockCols[p] )
      }
      if( s == 1 & (!is.null(stockNames)) )
        legend( "topright", bty = "n",
                legend = stockNames,
                col = stockCols,
                pch = 22, pt.bg = stockCols )

  }

}

plotTulipHR <- function( obj = blob, nTrace = 3 )
{
  # Model dimensions
  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP
  nT      <- obj$om$nT
  nF      <- obj$om$nF

  # Good replicates
  goodReps <- obj$goodReps

  # Pull target HR
  # Uref_p <- obj$ctlList$mp$hcr$Uref_p
  Uref_p <- rep(obj$ctlList$mp$hcr$inputF_s[1], nP)


  fleetType_f <- obj$om$fleetType_f
  sokFleets   <- which(fleetType_f %in% 2:3)
  predGears   <- obj$ctlList$opMod$predGears

  if(is.null(predGears))
    predGears <- c()

  commGears   <- (1:nF)[!1:nF %in% predGears]

  # Catch
  C_ispft  <- obj$om$C_ispft[goodReps,,,,,drop = FALSE]
  C_ispt       <- apply(X = C_ispft[,,,-sokFleets,,drop = FALSE], FUN = sum, MARGIN = c(1:3,5))
  predC_ispt   <- apply(X = C_ispft[,,,predGears[!predGears %in% sokFleets],,drop = FALSE], FUN = sum, MARGIN = c(1:3,5))
  commC_ispt   <- apply(X = C_ispft[,,,commGears[!commGears %in% sokFleets],,drop = FALSE], FUN = sum, MARGIN = c(1:3,5))

  # Dead ponded fish
  dP_ispft     <- obj$om$dP_ispft[goodReps,,,,,drop = FALSE]
  dP_ispt      <- apply(X = obj$om$dP_ispft[goodReps,,,sokFleets,,drop = FALSE], FUN=sum, MARGIN=c(1:3,5))
  predP_ispt   <- apply(dP_ispft[,,,intersect(sokFleets,predGears),,drop=FALSE], FUN=sum, MARGIN=c(1,2,3,5))
  commP_ispt   <- apply(dP_ispft[,,,intersect(sokFleets,commGears),,drop=FALSE], FUN=sum, MARGIN=c(1,2,3,5))

  # Total oonded fish
  P_ispft     <- obj$om$P_ispft[goodReps,,,,,drop = FALSE]
  P_ispt      <- apply(X = obj$om$P_ispft[goodReps,,,sokFleets,,drop = FALSE], FUN=sum, MARGIN=c(1:3,5))
  predP_ispt   <- apply(P_ispft[,,,intersect(sokFleets,predGears),,drop=FALSE], FUN=sum, MARGIN=c(1,2,3,5))
  commP_ispt   <- apply(P_ispft[,,,intersect(sokFleets,commGears),,drop=FALSE], FUN=sum, MARGIN=c(1,2,3,5))


  # SB_ispt  <- obj$om$endSB_ispt[goodReps,,,,drop = FALSE]
  SB_ispt  <- obj$om$SB_ispt[goodReps,,,,drop = FALSE]

  # Calculate aggregate harvest rate
  # commU_ispt   <- commC_ispt/(SB_ispt + commC_ispt + commP_ispt + predP_ispt + predC_ispt)
  # predU_ispt   <- predC_ispt/(SB_ispt + commC_ispt + commP_ispt + predP_ispt + predC_ispt)
  commU_ispt   <- (commC_ispt + commP_ispt)/(SB_ispt + commC_ispt + commP_ispt + predC_ispt)
  predU_ispt   <- predC_ispt/(SB_ispt + commC_ispt + commP_ispt + predC_ispt)

  # Calculate fleet specific harvest rate
  U_ispft <- array(NA, dim=c(sum(goodReps),nS,nP,nF,nT))
  for (s in 1:nS)
    for(f in 1:nF)
    {
      if(f %in% sokFleets)
        U_ispft[,s,,f,] <- P_ispft[,s,,f,]/(SB_ispt[,s,,] + C_ispt[,s,,] + commP_ispt[,s,,] + predC_ispt[,s,,])

      if(!f %in% sokFleets )
        U_ispft[,s,,f,] <- C_ispft[,s,,f,]/(SB_ispt[,s,,] + C_ispt[,s,,] + commP_ispt[,s,,] + predC_ispt[,s,,])
    }

  # Harvest rate envelopes
  U_qspft <- apply(  X = U_ispft, FUN = quantile,
                    MARGIN = c(2,3,4,5), probs = c(0.025, 0.5, 0.975),
                    na.rm = T )

  # Quantile for aggregate harvest rate
  predU_qspt <- apply(  X = predU_ispt, FUN = quantile,
                        MARGIN = c(2,3,4), probs = c(0.025, 0.5, 0.975),
                        na.rm = T )

  commU_qspt <- apply(  X = commU_ispt, FUN = quantile,
                        MARGIN = c(2,3,4), probs = c(0.025, 0.5, 0.975),
                        na.rm = T )


  # fishing fleets
  fishG  <- obj$ctlList$opMod$commGears
  fleets <- obj$ctlList$opMod$fleets
  fColrs <- brewer.pal(nF, 'Dark2')

  nReps   <- dim(C_ispft)[1]

  speciesNames  <- obj$om$speciesNames
  stockNames    <- dimnames(obj$ctlList$opMod$histRpt$I_pgt)[[1]]
  fYear         <- obj$ctlList$opMod$fYear

  yrs <- seq( from = fYear, by = 1, length.out = nT)

  traces <- sample( 1:nReps, size = min(nTrace,nReps)  )

  par(  mfcol = c(nP,nS),
        mar = c(1,1.5,1,2),
        oma = c(5,5,3,3) )

  for(s in 1:nS)
  {
    for( p in 1:nP )
    {
      plot( x = range(yrs),
            y = c(0,max(predU_qspt[,s,p,], commU_qspt[,s,p,], na.rm = T) ),
            type = "n", axes = FALSE)

      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1 )
      if( mfg[1] == 1 )
        mtext( side = 3, text = speciesNames[s], font = 2, line = 0 )
      axis( side = 2, las = 1 )
      if( mfg[2] == mfg[4] )
      {
        corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
        par(xpd = TRUE) #Draw outside plot area
        text(x = corners[2]+1, y = mean(corners[3:4]), stockNames[p], srt = 270,
              font = 2, cex = 1.5 )
        par(xpd = FALSE)
      }

      box()
      grid()

      # plot one polygon for aggregate comm and pred U
      polygon(  x = c(yrs, rev(yrs)),
                y = c(commU_qspt[1,s,p,], rev(commU_qspt[3,s,p,])),
                col = "grey65", border = NA )
      lines( x = yrs, y = commU_qspt[2,s,p,], lwd = 1, lty=2 )

      polygon(  x = c(yrs, rev(yrs)),
                y = c(predU_qspt[1,s,p,], rev(predU_qspt[3,s,p,])),
                col = scales::alpha("darkgreen",.5), border = NA )
      lines( x = yrs, y = predU_qspt[2,s,p,], lwd = 1, lty=2,
             col = "darkgreen" )

      # plot individual lines for each fleet
      for( f in fishG )
      {

        if(sum(U_qspft[2,s,p,f,], na.rm=T)==0)
          next()

        lines( x = yrs, y = U_qspft[2,s,p,f,], col=fColrs[f], lwd = 1.5)
        # for( tIdx in traces )
        #   lines( x = yrs, y = C_ispft[tIdx,s,p,f,], lwd = .8 )
      }

      # abline(h=0, lwd=1.5)
      abline( v = yrs[tMP], col = "grey30", lty = 3 )
      abline( h = Uref_p[p], lty = 2, lwd = 1, col = "red")

      legend('topright', bty='n', cex=0.8,
              legend=c('predator HR', fleets[fishG],'Target HR'),
              lwd=c(1, rep(1,length(fishG)),1),
              lty=c(2, rep(1,length(fishG)),3),
              col=c('darkgreen',fColrs[fishG],'red')
              )

    }
  }

  mtext( side = 2, outer = TRUE, text = "Harvest Rate",
          line = 2, font = 2)

  mtext( side = 1, outer = TRUE, text = "Year",
          line = 2, font = 2)
}


plotTulipF <- function( obj = blob, nTrace = 3 )
{
  # Model dimensions
  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP
  nT      <- obj$om$nT
  nF      <- obj$om$nF

  # Good replicates
  goodReps <- obj$goodReps

  # Pull reference points
  Fmsy_sp <- obj$rp[[1]]$FmsyRefPts$Fmsy_sp

  # Fishing mortality series
  F_ispft <- obj$om$F_ispft[goodReps,,,,,drop = FALSE]

  # Fishing mortality envelopes
  F_qspft <- apply(  X = F_ispft, FUN = quantile,
                    MARGIN = c(2,3,4,5), probs = c(0.025, 0.5, 0.975),
                    na.rm = T )

  # Aggregate fishing mortality
  F_ispt <- apply(  X = F_ispft, FUN = sum,
                    MARGIN = c(1,2,3,5),
                    na.rm = T )
  F_qspt <- apply(  X = F_ispt, FUN = quantile,
                    MARGIN = c(2,3,4), probs = c(0.025, 0.5, 0.975),
                    na.rm = T )

  # fishing fleets
  fishG  <- obj$ctlList$opMod$commGears
  fleets <- obj$ctlList$opMod$fleets
  fColrs <- brewer.pal(nF, 'Dark2')


  nReps   <- dim(F_ispft)[1]

  speciesNames  <- obj$om$speciesNames
  stockNames    <- dimnames(obj$ctlList$opMod$histRpt$I_pgt)[[1]]
  fYear         <- obj$ctlList$opMod$fYear

  yrs <- seq( from = fYear, by = 1, length.out = nT)

  traces <- sample( 1:nReps, size = min(nTrace,nReps)  )

  par(  mfcol = c(nP,nS),
        mar = c(1,1.5,1,2),
        oma = c(5,5,3,3) )


  for(s in 1:nS)
  {
    for( p in 1:nP )
    {
      plot( x = range(yrs),
            y = c(0,max(F_qspt[,s,p,], na.rm = T) ),
            type = "n", axes = F, ylim=c(0,max(F_qspt[,s,p,], na.rm = T) ))

      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1 )
      if( mfg[1] == 1 )
        mtext( side = 3, text = speciesNames[s], font = 2, line = 0 )
      axis( side = 2, las = 1 )
      if( mfg[2] == mfg[4] )
      {
        corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
        par(xpd = TRUE) #Draw outside plot area
        text(x = corners[2]+2, y = mean(corners[3:4]), stockNames[p], srt = 270,
              font = 2, cex = 1.5 )
        par(xpd = FALSE)
      }

      box()
      grid()

      # plot one polygon for aggregate F
      polygon(  x = c(yrs, rev(yrs)),
                y = c(F_qspt[1,s,p,], rev(F_qspt[3,s,p,])),
                col = "grey65", border = NA )
      # lines( x = yrs, y = F_qspt[2,s,p,], lwd = 2 )

      # plot individual lines for each fleet
      for( f in fishG )
      {

        if(sum(F_qspft[2,s,p,f,], na.rm=T)==0)
          next()

        lines( x = yrs, y = F_qspft[2,s,p,f,], col=fColrs[f], lwd = 1.5)
        # for( tIdx in traces )
        #   lines( x = yrs, y = F_ispft[tIdx,s,p,f,], lwd = .8 )
      }

      # abline(h=0, lwd=1.5)
      abline( v = yrs[tMP], col = "grey30", lty = 3 )
      abline( h = Fmsy_sp[s,p], lty = 2, lwd = 1, col = "red")

      legend('topright', bty='n', cex=0.8,
              legend=c('total F', fleets[fishG],'Fmsy'),
              lwd=c(2, rep(1,length(fishG)),1),
              lty=c(1, rep(1,length(fishG)),3),
              col=c('black',fColrs[fishG],'red')
              )

    }
  }

  mtext( side = 2, outer = TRUE, text = "Fishing mortality (/yr)",
          line = 2, font = 2)

  mtext( side = 1, outer = TRUE, text = "Year",
          line = 2, font = 2)
}

# TAC utilisation envelopes
plotTulipTACu <- function( obj = blob, nTrace = 3, fIdx=1 )
{

  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP
  nT      <- obj$om$nT

  goodReps <- obj$goodReps


  # Get catch for trawl fleet, in projections only
  C_ispft     <- obj$om$C_ispft[goodReps,,,2,tMP:nT,drop = FALSE]
  TAC_ispft   <- obj$mp$hcr$TAC_ispft[goodReps,,,2,tMP:nT,drop = FALSE]

  nReps   <- dim(TAC_ispft)[1]

  TACu_ispft <- C_ispft / TAC_ispft


  TACu_qspt <- apply( X = TACu_ispft, FUN = quantile,
                      MARGIN = c(2,3,5), probs = c(0.025, 0.5, 0.975),
                      na.rm = T )

  TACu_ispt <- array(0, dim = c(nReps,nS,nP,(nT - tMP + 1)))
  TACu_ispt[1:nReps,,,] <- TACu_ispft[1:nReps,,,fIdx,]

  speciesNames  <- obj$om$speciesNames
  stockNames    <- obj$om$stockNames
  fYear         <- obj$ctlList$opMod$fYear

  yrs <- seq( from = fYear, by = 1, length.out = nT)
  yrs <- yrs[tMP:nT]


  traces <- sample( 1:nReps, size = min(nTrace,nReps)  )

  par(  mfcol = c(nP,nS),
        mar = c(1,1.5,1,2),
        oma = c(5,5,3,3) )


  for(s in 1:nS)
  {
    for( p in 1:nP )
    {
      plot( x = range(yrs),
            y = c(0,max(TACu_qspt, na.rm = T) ),
            type = "n", axes = F )

      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1 )
      if( mfg[1] == 1 )
        mtext( side = 3, text = speciesNames[s], font = 2, line = 0 )
      axis( side = 2, las = 1 )
      if( mfg[2] == mfg[4] )
      {
        corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
        par(xpd = TRUE) #Draw outside plot area
        text(x = corners[2]+2, y = mean(corners[3:4]), stockNames[p], srt = 270,
              font = 2, cex = 1.5 )
        par(xpd = FALSE)
      }
      box()
      grid()
      polygon(  x = c(yrs, rev(yrs)),
                y = c(TACu_qspt[1,s,p,], rev(TACu_qspt[3,s,p,])),
                col = "grey65", border = NA )
      lines( x = yrs, y = TACu_qspt[2,s,p,], lwd = 3 )
      for( tIdx in traces )
        lines( x = yrs, y = TACu_ispt[tIdx,s,p,], lwd = .8 )

      abline( v = yrs[tMP], col = "grey30", lty = 3 )
    }
  }
  mtext( side = 2, outer = TRUE, text = expression(C[spt]/TAC[spt]),
          line = 2, font = 2)

  mtext( side = 1, outer = TRUE, text = "Year",
          line = 2, font = 2)

} # END plotTulipTACu


# plotDistN_SOK()
# Distribution of N_SOK values across
# all simulation replicates.
plotDistN_SOK <- function( obj = blob )
{
  goodReps <- obj$goodReps

  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP
  nT      <- obj$om$nT
  nReps   <- length(goodReps)

  stockNames    <- dimnames(obj$ctlList$opMod$histRpt$I_pgt)[[1]]

  stockNames <- c(stockNames,"Agg")

  fleetType_f <- obj$om$fleetType_f
  sokFleets <- which(fleetType_f == 2)

  P_ispft   <- obj$om$P_ispft[goodReps,,,sokFleets,tMP:nT,drop = FALSE]
  P_ispt    <- apply(X = P_ispft[,1,,,,drop = FALSE], FUN = sum, MARGIN = c(1,2,3,5) )

  P_it      <- apply(X = P_ispt, FUN = sum, MARGIN = c(1,4))

  P_ispt[P_ispt > 0] <- 1
  P_it[P_it > 0] <- 1

  stamp <- paste(obj$ctlList$ctl$scenarioName,":",obj$ctlList$ctl$mpName,sep = "")


  # Get granular dist
  Nsok_isp <- apply( X = P_ispt, FUN = sum, MARGIN = c(1,2,3), na.rm = T)

  Nsok_qsp <- apply(X = Nsok_isp, FUN = quantile, probs = c(0.025,.25, .5, .75, .975),
                      MARGIN = c(2,3))

  # Get aggregate distribution as well
  Nsok_i <- apply(X = P_it, FUN = sum, MARGIN = 1)
  Nsok_q <- quantile(Nsok_i, probs = c(.025, .25, .5, .75, .975) )

  # Make a bar plot

  plot( x = c(.5,4.5), y = c(0,max(Nsok_i, Nsok_isp)),
        type = "n", axes = FALSE, xlab = "Population", ylab = "Number of years with SOK fishing" )
    axis( side = 1, at = 1:4, labels = stockNames )
    axis( side = 2, las = 1 )
    box()
    grid()
    for( p in 1:nP )
    {
      # rect( xleft = p - .3, xright= p + .3,
      #       ybottom = 0, ytop = Nsok_qsp[2,1,p], border = NA,
      #       col = "grey60" )

      segments( x0 = p, x1 = p,
                y0 = Nsok_qsp[1,1,p],
                y1 = Nsok_qsp[5,1,p], col = "black", lwd = 2 )
      rect( xleft = p-.1, xright = p+.1,
            ybottom = Nsok_qsp[2,1,p], border = "black",
            ytop = Nsok_qsp[4,1,p], col = "white", lwd = 2 )
      segments( x0 = p-.1, x1 = p + .1, y0 = Nsok_qsp[3,1,p], col = "black", lwd = 2 )
    }
    segments( x0 = nP + 1, x1 = nP + 1,
              y0 = Nsok_q[1],
              y1 = Nsok_q[5], col = "black", lwd = 2 )
    rect( xleft = nP + 1 - .1, xright= nP + 1 + .1,
          ybottom = Nsok_q[2], ytop = Nsok_q[4], border = "black",
          col = "white", lwd = 2 )
    segments( x0 = nP + 1-.1, x1 = nP + 1 + .1, y0 = Nsok_q[3], col = "black", lwd = 2 )


    mtext(  side = 1, adj = .9, line = 2, cex = .6,
            text = stamp, col = "grey60" )

} # END plotDistN_SOK

plotTulipPondedFish <- function(  obj = blob,
                                  nTrace = 3,
                                  traces = NULL,
                                  leg = TRUE,
                                  proj = FALSE )
{
  goodReps <- obj$goodReps

  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP
  nT      <- obj$om$nT
  nReps   <- length(goodReps)

  fleetType_f <- obj$om$fleetType_f
  sokFleets <- which(fleetType_f == 2)

  P_ispft   <- obj$om$P_ispft[goodReps,,,sokFleets,,drop = FALSE]
  P_ispt    <- apply(X = P_ispft[,1,,,,drop = FALSE], FUN = sum, MARGIN = c(1,2,3,5) )

  P_qspt    <- apply( X = P_ispt, FUN = quantile, probs = c(0.025, 0.5, 0.975),
                      MARGIN = c(2,3,4) )

  if( is.null(traces))
    traces <- sample( 1:nReps, size = min(nTrace,nReps)  )

  stamp <- paste(obj$ctlList$ctl$scenarioName,":",obj$ctlList$ctl$mpName,sep = "")

  speciesNames  <- obj$om$speciesNames
  stockNames    <- dimnames(obj$ctlList$opMod$histRpt$I_pgt)[[1]]
  fYear         <- obj$ctlList$opMod$fYear

  yrs <- seq( from = fYear, by = 1, length.out = nT)
  tMin <-1

  if( proj )
    tMin <- tMP - 1

  par(  mfcol = c(nP,nS),
        mar = c(1,1.5,1,1.5),
        oma = c(4,3,3,3) )


  for(s in 1:nS)
  {
    for( p in 1:nP )
    {
      plot( x = range(yrs[tMin:nT]),
            y = c(0,max(P_qspt[,s,p,tMin:nT], na.rm = T) ),
            type = "n", axes = F)

      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1 )
      if( mfg[1] == 1 )
        mtext( side = 3, text = speciesNames[s], font = 2, line = 0 )
      axis( side = 2, las = 1 )
      if( mfg[2] == mfg[4] )
        rmtext( outer = TRUE, cex = 1.5, txt = stockNames[p],
                font = 2, line = 1)
      box()
      grid()
      polygon(  x = c(yrs, rev(yrs)),
                y = c(P_qspt[1,s,p,], rev(P_qspt[3,s,p,])),
                col = "grey65", border = NA )
      lines( x = yrs, y = P_qspt[2,s,p,], lwd = 3 )

      for( tIdx in traces )
        lines( x = yrs, y = P_ispt[tIdx,s,p,], lwd = .8 )

      abline( v = yrs[tMP], col = "grey30", lty = 3 )


      if( mfg[1] == 1 & mfg[2] == 1 & leg )
        legend( x = "topleft", bty = "n",
                legend = c( "Median Ponded Fish",
                            "Central 95%",
                            "Replicate Traces"),
                col = c(  "black", "grey65", "black"),
                pch = c(NA,15, NA ),
                lty = c(1, NA, 1),
                lwd = c(3, NA, .8 ) )
    }
  }
  mtext( side = 2, outer = TRUE, text = "Ponded Fish (kt)",
          line = 1.5, font = 2)

  mtext( side = 1, outer = TRUE, text = "Year",
          line = 2, font = 2)

  mtext(  outer = TRUE, side = 1, adj = .8, line = 3, cex = .6,
            text = stamp, col = "grey60" )

} # END plotTulipPondedFish

# Biomass envelopes - with catch
# if chosen
plotTulipBt <- function(  obj = blob, nTrace = 3,
                          traces = NULL,
                          dep = FALSE,
                          ref = "B0",
                          var = "SB_ispt",
                          Ct  = FALSE,
                          leg = TRUE,
                          proj = FALSE,
                          tMin = NULL,
                          B0_p = NULL )
{
  goodReps <- obj$goodReps

  SB_ispt   <- obj$om[[var]][goodReps,,,,drop = FALSE]
  C_ispt    <- obj$om$C_ispt[goodReps,,,,drop = FALSE]

  SB_ispt[SB_ispt == 0] <- NA


  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP
  nT      <- obj$om$nT
  nReps   <- dim(SB_ispt)[1]
  pT      <- obj$ctlList$opMod$pT

  # Get reference points
  B0_sp     <- array(NA, dim = c(nS,nP))
  BmsySS_sp <- array(NA, dim = c(nS,nP))
  BmsyMS_sp <- array(NA, dim = c(nS,nP))

  tInit_sp  <- obj$om$tInit_sp

  B0_sp[1:nS,]     <- obj$rp[[1]]$B0_sp
  BmsySS_sp[1:nS,] <- obj$rp[[1]]$FmsyRefPts$BeqFmsy_sp

  speciesNames  <- obj$om$speciesNames
  stockNames    <- dimnames(obj$ctlList$opMod$histRpt$I_pgt)[[1]]
  fYear         <- obj$ctlList$opMod$fYear

  yrs <- seq( from = fYear, by = 1, length.out = nT)

  if(is.null(tMin))
    tMin <-1

  if( proj )
    tMin <- tMP - 1

  if( dep )
  {

    if( ref == "B0" )
    {
      for( s in 1:nS )
        for( p in 1:nP )
        {

          SB_ispt[,s,p,] <- SB_ispt[,s,p,] / B0_sp[s,p]
          C_ispt[,s,p,]  <- C_ispt[,s,p,] / B0_sp[s,p]

        }
      # LCP_p     <- LCP_p / B0_sp
      BmsySS_sp <- BmsySS_sp / B0_sp
      B0_sp     <- B0_sp / B0_sp

    }

    if( ref == "Bmsy" )
    {
      for( s in 1:nS )
        for( p in 1:nP )
        {
          SB_ispt[,s,p,] <- SB_ispt[,s,p,] / Bmsy_sp[s,p]
          C_ispt[,s,p,]  <- C_ispt[,s,p,] / Bmsy_sp[s,p]
        }

      # LCP_p     <- LCP_p / Bmsy_sp[1,]
      B0_sp     <- B0_sp / Bmsy_sp
      BmsySS_sp <- BmsySS_sp / Bmsy_sp

    }
  }

  if( !dep )
    yAxisLab <- "Biomass (kt)"

  if( dep )
  {
    if( ref == "Bmsy" )
      yAxisLab <- expression(B[t]/B[MSY])

    if( ref == "B0" )
      yAxisLab <- expression(B[t]/B[0])
  }

  # Now take quantiles
  SB_qspt <- apply( X = SB_ispt, FUN = quantile,
                    MARGIN = c(2,3,4), probs = c(0.025, 0.5, 0.975),
                    na.rm = T )

  C_qspt <- apply( X = C_ispt, FUN = quantile,
                    MARGIN = c(2,3,4), probs = c(0.025, 0.5, 0.975),
                    na.rm = T )

  # Estimate SB0 estimate for years 100-500
  if(pT>100)
  {
      t1 <- tMP + 50
      t2 <- tMP+pT - 1
      simB0_sp <- apply( X = SB_ispt[,,,t1:t2,drop=F], FUN = mean,
                        MARGIN = c(2,3), na.rm = T )
  }


  SB_qspt[SB_qspt == 0] <- NA

  # grab B0 from history ref pts, unless supplied
  if (is.null(B0_p))
    B0_p <- obj$ctlList$opMod$repList$repOpt$refPts$refCurves$SBeq_pf[,1]

  if( is.null(traces))
    traces <- sample( 1:nReps, size = min(nTrace,nReps)  )

  stamp <- paste(obj$ctlList$ctl$scenarioName,":",obj$ctlList$ctl$mpName,sep = "")

  par(  mfcol = c(nP,nS),
        mar = c(1,1.5,1,1.5),
        oma = c(4,3,3,3) )


  for(s in 1:nS)
  {
    for( p in 1:nP )
    {
      plot( x = range(yrs[tMin:nT]),
            y = c(0,max(SB_qspt[,s,p,tMin:nT], na.rm = T) ),
            type = "n", axes = F)

      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1 )
      if( mfg[1] == 1 )
        mtext( side = 3, text = speciesNames[s], font = 2, line = 0 )
      axis( side = 2, las = 1 )
      if( mfg[2] == mfg[4] )
        rmtext( outer = TRUE, cex = 1.5, txt = stockNames[p],
                font = 2, line = 10)
      box()
      grid()
      tPlotIdx <- tInit_sp[s,p]:nT
      polygon(  x = c(yrs[tPlotIdx], rev(yrs[tPlotIdx])),
                y = c(SB_qspt[1,s,p,tPlotIdx], rev(SB_qspt[3,s,p,tPlotIdx])),
                col = "grey65", border = NA )
      lines( x = yrs[tPlotIdx], y = SB_qspt[2,s,p,tPlotIdx], lwd = 3 )

      for( tIdx in traces )
        lines( x = yrs[tPlotIdx], y = SB_ispt[tIdx,s,p,tPlotIdx], lwd = .8 )

      if( Ct )
      {
        rect( xleft = yrs[tPlotIdx] - .3, xright = yrs[tPlotIdx] + .3,
              ybottom = 0, ytop = C_qspt[2,s,p,tPlotIdx], col = "grey65",
              border = NA )
        segments( x0 = yrs[tMP:nT], x1 = yrs[tMP:nT],
                  y0 = C_qspt[1,s,p,tMP:nT], y1 = C_qspt[3,s,p,tMP:nT],
                  col = "black" )
      }

      abline( v = yrs[tMP], col = "grey30", lty = 3 )
      abline( h = B0_p[p], lty = 3, col = "grey50", lwd = 2  )
      # abline( h = 0.3*B0_sp[s,p], lty = 2, col = "red", lwd = 3)
      # abline( h = LCP_p[p], lty = 2, col = "salmon", lwd = 3)
      # abline( h = BmsySS_sp[s,p], lty = 3, col = "darkgreen", lwd = 2)

      if( mfg[1] == 1 & mfg[2] == 1 & leg )
        legend( x = "topright", bty = "n",
                legend = c( "Median Spawning Biomass",
                            "Central 95%",
                            "Replicate Traces"),
                            "Unfished Biomass",
                            # expression(paste("LRP = 0.3",B[0]))),
                            # "Median Projection SB" ),
                            # expression(B[MSY,MS])),
                col = c(  "black", "grey65", "black",
                          "grey50", "red","salmon" ),
                pch = c(NA,15, NA, NA),#  NA,NA,NA),
                lty = c(1, NA, 1, 3), #  2, 2),
                lwd = c(3, NA, .8, 2) )#, 3, 3 ) )


      if(pT>100)
      {
        abline( h = simB0_sp[s,p], col = "blue", lty = 3, lwd=1.5 )
        legend( x = "topleft", bty = "n",
                legend = c( "B0 simulation estimate"),
                col = "blue",lty = 3, lwd=1.5 )

      }

    }
  }
  mtext( side = 2, outer = TRUE, text = yAxisLab,
          line = 1.5, font = 2)

  mtext( side = 1, outer = TRUE, text = "Year",
          line = 2, font = 2)

  mtext(  outer = TRUE, side = 1, adj = .8, line = 3, cex = .6,
            text = stamp, col = "grey60" )
}




# Catch envelopes
plotTulipCt <- function(  obj = blob, nTrace = 3,
                          ref = "B0",
                          leg = TRUE,
                          tMin = NULL,
                          plotFleets=FALSE )
{
  goodReps <- obj$goodReps

  C_ispt    <- obj$om$C_ispt[goodReps,,,,drop = FALSE]
  C_ispft   <- obj$om$C_ispft[goodReps,,,,,drop = FALSE]

  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP
  nT      <- obj$om$nT
  nF      <- obj$om$nF
  nReps   <- dim(C_ispt)[1]

  speciesNames  <- obj$om$speciesNames
  stockNames    <- dimnames(obj$ctlList$opMod$histRpt$I_pgt)[[1]]
  fYear         <- obj$ctlList$opMod$fYear

  # fishing fleets
  fishG  <- obj$ctlList$opMod$commGears
  fleets <- obj$ctlList$opMod$fleets
  fColrs <- brewer.pal(nF, 'Dark2')

  yrs <- seq( from = fYear, by = 1, length.out = nT)

  if(is.null(tMin))
    tMin <-1

  C_qspt <- apply( X = C_ispt, FUN = quantile,
                    MARGIN = c(2,3,4), probs = c(0.025, 0.5, 0.975),
                    na.rm = T )

  C_qspft <- apply( X = C_ispft, FUN = quantile,
                    MARGIN = c(2,3,4,5), probs = c(0.025, 0.5, 0.975),
                    na.rm = T )


  traces <- sample( 1:dim(C_ispt)[1], size = min(nTrace,nReps)  )

  stamp <- paste(obj$ctlList$ctl$scenarioName,":",obj$ctlList$ctl$mpName,sep = "")

  par(  mfcol = c(nP,nS),
        mar = c(1,1.5,1,1.5),
        oma = c(4,3,3,3) )

  for(s in 1:nS)
  {
    for( p in 1:nP )
    {
      plot( x = range(yrs[tMin:nT]),
            y = c(0,max(C_qspt[,s,p,tMin:nT], na.rm = T) ),
            type = "n", axes = F)

      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1 )
      if( mfg[1] == 1 )
        mtext( side = 3, text = speciesNames[s], font = 2, line = 0 )
      axis( side = 2, las = 1 )
      if( mfg[2] == mfg[4] )
      {
        corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
        par(xpd = TRUE) #Draw outside plot area
        text(x = corners[2]+0.1, y = mean(corners[3:4]), stockNames[p], srt = 270,
              font = 2, cex = 1.5 )
        par(xpd = FALSE)
      }
      box()
      grid()
      polygon(  x = c(yrs, rev(yrs)),
                y = c(C_qspt[1,s,p,], rev(C_qspt[3,s,p,])),
                col = "grey65", border = NA )
      lines( x = yrs, y = C_qspt[2,s,p,], lwd = 3 )


      for( tIdx in traces )
        lines( x = yrs, y = C_ispt[tIdx,s,p,], lwd = .8 )

      abline( v = yrs[tMP], col = "grey30", lty = 3 )

      if( mfg[1] == 1 & mfg[2] == 1 & leg )
        legend( x = "topleft", bty = "n",
                legend = c( "Median Catch (kt)",
                            "Central 95%",
                            "Replicate Traces"),
                col = c(  "black", "grey65", "black"),
                pch = c(NA,15, NA),
                lty = c(1, NA, 1),
                lwd = c(3, NA, .8) )


      # plot individual lines for each fleet
      if(plotFleets)
      {
        for( f in fishG )
        {

          if(sum(C_qspft[2,s,p,f,],na.rm=T)==0)
            next()

          lines( x = yrs, y = C_qspft[2,s,p,f,], col=fColrs[f], lwd = 1.5)
          # for( tIdx in traces )
          #   lines( x = yrs, y = F_ispft[tIdx,s,p,f,], lwd = .8 )
        }

        if( mfg[1] == 1 & mfg[2] == 1 & leg & !is.null(fleets) )
        legend('bottomleft', bty='n', cex=0.8,
                legend=c(fleets[fishG]),
                lwd=c(rep(1,length(fishG))),
                lty=c(rep(1,length(fishG))),
                col=c(fColrs[fishG])
                )
      }
    }
  }
  mtext( side = 2, outer = TRUE, text = 'Catch (1000 t)',
          line = 1.5, font = 2)

  mtext( side = 1, outer = TRUE, text = "Year",
          line = 2, font = 2)

  mtext(  outer = TRUE, side = 1, adj = .8, line = 3, cex = .6,
            text = stamp, col = "grey60" )
} # END plotTulipCt()

# Catch envelopes
plotTulipPredC_ft <- function( obj = blob, nTrace = 3,
                            leg = TRUE,
                            tMin = NULL )
{
  goodReps <- obj$goodReps

  C_ispft   <- obj$om$C_ispft[goodReps,,,,,drop = FALSE]

  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP
  nT      <- obj$om$nT
  nF      <- obj$om$nF
  nReps   <- dim(C_ispft)[1]

  speciesNames  <- obj$om$speciesNames
  stockNames    <- dimnames(obj$ctlList$opMod$histRpt$I_pgt)[[1]]
  fYear         <- obj$ctlList$opMod$fYear

  # fishing fleets
  fishG  <- obj$ctlList$opMod$commGears
  predG  <- obj$ctlList$opMod$predGears

  nPreds <- length(predG)

  fleets <- obj$ctlList$opMod$fleets
  fColrs <- brewer.pal(nPreds, 'Dark2')

  yrs <- seq( from = fYear, by = 1, length.out = nT)

  if(is.null(tMin))
    tMin <-1

  C_qspft <- apply( X = C_ispft, FUN = quantile,
                    MARGIN = c(2,3,4,5), probs = c(0.025, 0.5, 0.975),
                    na.rm = T )

  traces <- sample( 1:nReps, size = min(nTrace,nReps)  )

  stamp <- paste(obj$ctlList$ctl$scenarioName,":",obj$ctlList$ctl$mpName,sep = "")



  par(  mfcol = c(nPreds,nP),
        mar = c(.1,1.5,.1,1.5),
        oma = c(4,3,3,3) )

  for( p in 1:nP )
  {
    for( fIdx in 1:nPreds )
    {
      f <- predG[fIdx]
      plot( x = range(yrs[tMin:nT]),
            y = c(0,max(C_qspft[,1,p,f,tMin:nT], na.rm = T) ),
            type = "n", axes = F)

      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1 )

      if(mfg[2] == 1)
        axis( side = 2, las = 1 )
      box()
      grid()
      polygon(  x = c(yrs, rev(yrs)),
                y = c(C_qspft[1,1,p,f,], rev(C_qspft[3,1,p,f,])),
                col = fColrs[fIdx], border = NA )
      lines( x = yrs, y = C_qspft[2,1,p,f,], lwd = 3 )


      for( tIdx in traces )
        lines( x = yrs, y = C_ispft[tIdx,1,p,f,], lwd = .8 )

      abline( v = yrs[tMP], col = "grey30", lty = 3 )

      legend(x = "topleft", legend = fleets[f], bty = "n")


    }

  }

  mtext( side = 2, outer = TRUE, text = 'Consumption (kt)',
          line = 1.5, font = 2)

  mtext( side = 1, outer = TRUE, text = "Year",
          line = 2, font = 2)

  mtext(  outer = TRUE, side = 1, adj = .8, line = 3, cex = .6,
            text = stamp, col = "grey60" )
} # END plotTulipCt()

# Catch envelopes
plotTulipPredVulnB_ft <- function(  obj = blob, nTrace = 3,
                                    leg = TRUE,
                                    tMin = NULL )
{
  goodReps <- obj$goodReps

  vB_ispft   <- obj$om$vB_ispft[goodReps,,,,,drop = FALSE]

  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP
  nT      <- obj$om$nT
  nF      <- obj$om$nF
  nReps   <- dim(vB_ispft)[1]

  speciesNames  <- obj$om$speciesNames
  stockNames    <- dimnames(obj$ctlList$opMod$histRpt$I_pgt)[[1]]
  fYear         <- obj$ctlList$opMod$fYear

  # fishing fleets
  fishG  <- obj$ctlList$opMod$commGears
  predG  <- obj$ctlList$opMod$predGears

  nPreds <- length(predG)

  fleets <- obj$ctlList$opMod$fleets
  fColrs <- brewer.pal(nPreds, 'Dark2')

  yrs <- seq( from = fYear, by = 1, length.out = nT)

  if(is.null(tMin))
    tMin <-1

  vB_qspft <- apply( X = vB_ispft, FUN = quantile,
                    MARGIN = c(2,3,4,5), probs = c(0.025, 0.5, 0.975),
                    na.rm = T )

  traces <- sample( 1:nReps, size = min(nTrace,nReps)  )

  stamp <- paste(obj$ctlList$ctl$scenarioName,":",obj$ctlList$ctl$mpName,sep = "")



  par(  mfcol = c(nPreds,nP),
        mar = c(.1,1.5,.1,1.5),
        oma = c(4,3,3,3) )

  for( p in 1:nP )
  {
    for( fIdx in 1:nPreds )
    {
      f <- predG[fIdx]
      plot( x = range(yrs[tMin:nT]),
            y = c(0,max(vB_qspft[,1,p,f,tMin:nT], na.rm = T) ),
            type = "n", axes = F)

      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1 )

      if(mfg[2] == 1)
        axis( side = 2, las = 1 )
      box()
      grid()
      polygon(  x = c(yrs, rev(yrs)),
                y = c(vB_qspft[1,1,p,f,], rev(vB_qspft[3,1,p,f,])),
                col = fColrs[fIdx], border = NA )
      lines( x = yrs, y = vB_qspft[2,1,p,f,], lwd = 3 )


      for( tIdx in traces )
        lines( x = yrs, y = vB_ispft[tIdx,1,p,f,], lwd = .8 )

      abline( v = yrs[tMP], col = "grey30", lty = 3 )

      legend(x = "topleft", legend = fleets[f], bty = "n")


    }

  }

  mtext( side = 2, outer = TRUE, text = 'Vulnerable Herring Biomass (kt)',
          line = 1.5, font = 2)

  mtext( side = 1, outer = TRUE, text = "Year",
          line = 2, font = 2)

  mtext(  outer = TRUE, side = 1, adj = .8, line = 3, cex = .6,
            text = stamp, col = "grey60" )
} # END plotTulipCt()

# Catch envelopes
plotTulipPredF_ft <- function(  obj = blob, nTrace = 3,
                                leg = TRUE,
                                tMin = NULL )
{
  goodReps <- obj$goodReps

  F_ispft   <- obj$om$F_ispft[goodReps,,,,,drop = FALSE]

  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP
  nT      <- obj$om$nT
  nF      <- obj$om$nF
  nReps   <- dim(F_ispft)[1]

  speciesNames  <- obj$om$speciesNames
  stockNames    <- dimnames(obj$ctlList$opMod$histRpt$I_pgt)[[1]]
  fYear         <- obj$ctlList$opMod$fYear

  # fishing fleets
  fishG  <- obj$ctlList$opMod$commGears
  predG  <- obj$ctlList$opMod$predGears
  nPreds <- length(predG)

  fleets <- obj$ctlList$opMod$fleets
  fColrs <- brewer.pal(nPreds, 'Dark2')

  yrs <- seq( from = fYear, by = 1, length.out = nT)

  if(is.null(tMin))
    tMin <-1

  F_qspft <- apply( X = F_ispft, FUN = quantile,
                    MARGIN = c(2,3,4,5), probs = c(0.25, 0.5, 0.75),
                    na.rm = T )

  traces <- sample( 1:nReps, size = min(nTrace,nReps)  )

  stamp <- paste(obj$ctlList$ctl$scenarioName,":",obj$ctlList$ctl$mpName,sep = "")


  par(  mfcol = c(nPreds,nP),
        mar = c(.1,1.5,.1,1.5),
        oma = c(4,3,3,3) )

  for( p in 1:nP )
  {
    for( fIdx in 1:nPreds )
    {
      f <- predG[fIdx]
      plot( x = range(yrs[tMin:nT]),
            y = c(0,max(F_qspft[,1,p,f,tMin:nT], na.rm = T) ),
            type = "n", axes = F)

      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1 )

      if(mfg[2] == 1)
        axis( side = 2, las = 1 )
      box()
      grid()
      polygon(  x = c(yrs, rev(yrs)),
                y = c(F_qspft[1,1,p,f,], rev(F_qspft[3,1,p,f,])),
                col = fColrs[fIdx], border = NA )
      lines( x = yrs, y = F_qspft[2,1,p,f,], lwd = 3 )


      for( tIdx in traces )
        lines( x = yrs, y = F_ispft[tIdx,1,p,f,], lwd = .8 )

      abline( v = yrs[tMP], col = "grey30", lty = 3 )

      legend(x = "topleft", legend = fleets[f], bty = "n")


    }

  }

  mtext( side = 2, outer = TRUE, text = 'Predation Mortality (/yr)',
          line = 1.5, font = 2)

  mtext( side = 1, outer = TRUE, text = "Year",
          line = 2, font = 2)

  mtext(  outer = TRUE, side = 1, adj = .8, line = 3, cex = .6,
            text = stamp, col = "grey60" )
} # END plotTulipCt()


# Catch envelopes
plotTulipMt <- function(  obj = blob, nTrace = 3,
                          leg = TRUE, tMin = NULL)
{
  goodReps <- obj$goodReps

  M_iaxspt <- obj$om$M_iaxspt[goodReps,,,,,,drop = FALSE]


  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP
  nT      <- obj$om$nT
  nF      <- obj$om$nF
  nReps   <- dim(M_iaxspt)[1]

  speciesNames  <- obj$om$speciesNames
  stockNames    <- dimnames(obj$ctlList$opMod$histRpt$I_pgt)[[1]]
  fYear         <- obj$ctlList$opMod$fYear

  yrs <- seq( from = fYear, by = 1, length.out = nT)

  if(is.null(tMin))
    tMin <-1

  M_qspt <- apply( X = M_iaxspt[,5,1,,,,drop=FALSE], FUN = quantile,
                    MARGIN = c(4,5,6), probs = c(0.025, 0.5, 0.975),
                    na.rm = T )

  traces <- sample( which(goodReps), size = min(nTrace,nReps)  )

  stamp <- paste(obj$ctlList$ctl$scenarioName,":",obj$ctlList$ctl$mpName,sep = "")

  par(  mfcol = c(nP,nS),
        mar = c(1,1.5,1,1.5),
        oma = c(4,3,3,3) )

  for(s in 1:nS)
  {
    for( p in 1:nP )
    {
      plot( x = range(yrs[tMin:nT]),
            y = c(0,max(M_iaxspt[,5,,s,p,tMin:nT], na.rm = T) ),
            type = "n", axes = F)

      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1 )
      if( mfg[1] == 1 )
        mtext( side = 3, text = speciesNames[s], font = 2, line = 0 )
      axis( side = 2, las = 1 )
      if( mfg[2] == mfg[4] )
      {
        corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
        par(xpd = TRUE) #Draw outside plot area
        text(x = corners[2]+1, y = mean(corners[3:4]), stockNames[p], srt = 270,
              font = 2, cex = 1.5 )
        par(xpd = FALSE)
      }
      box()
      grid()
      polygon(  x = c(yrs, rev(yrs)),
                y = c(M_qspt[1,s,p,], rev(M_qspt[3,s,p,])),
                col = "grey65", border = NA )
      lines( x = yrs, y = M_qspt[2,s,p,], lwd = 3 )


      for( tIdx in traces )
        lines( x = yrs, y = M_iaxspt[tIdx,5,1,s,p,], lwd = .8 )

      abline( v = yrs[tMP], col = "grey30", lty = 3 )


        legend( x = "topleft", bty = "n",
                legend = c( "Median M",
                            "Central 95%",
                            "Replicate Traces"),
                col = c(  "black", "grey65", "black"),
                pch = c(NA,15, NA),
                lty = c(1, NA, 1),
                lwd = c(3, NA, .8) )

    }
  }
  mtext( side = 2, outer = TRUE, text = 'Natural mortality',
          line = 1.5, font = 2)

  mtext( side = 1, outer = TRUE, text = "",
          line = 2, font = 2)

  mtext(  outer = TRUE, side = 1, adj = .8, line = 3, cex = .6,
            text = stamp, col = "grey60" )
} # END plotTulipMt()

# plotTulipEffort_p()
# Effort over time gridded
# by stock area and fleet (right now predators) - envelopes
plotTulipEffort_p <- function(  obj = blob,
                                nTrace = 3,
                                pIdx = 1,
                                fIdx = 7:14,
                                cols = NULL,
                                tdx = NULL )
{


  goodReps <- obj$goodReps

  E_ipft <- obj$om$E_ipft[goodReps,pIdx,,,drop = FALSE ]

  E_ipft[E_ipft == 0] <- NA

  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP
  nT      <- obj$om$nT
  nF      <- obj$om$nF
  nReps   <- dim(E_ipft)[1]

  if(is.null(tdx))
    tdx <- 1:nT

  E_qpft <- apply(  X = E_ipft, FUN = quantile,
                    MARGIN = c(2,3,4), probs = c(0.025, 0.5, 0.975),
                    na.rm = T )

  E_qpft[E_qpft == 0] <- NA

  traces <- sample( 1:nReps, size = min(nTrace,nReps)  )

  speciesNames  <- obj$om$speciesNames
  stockNames    <- obj$om$stockNames
  fleetNames    <- obj$om$fleetNames
  fYear         <- obj$ctlList$opMod$fYear

  yrs <- seq( from = fYear, by = 1, length.out = nT)

  if(is.null(cols))
    fleetCols <- RColorBrewer::brewer.pal( length(fIdx), "Dark2" )
  else fleetCols <- cols

  par(  mfcol = c(length(fIdx),nP),
        mar = c(.1,1.5,.1,1.5),
        oma = c(5,4,3,3) )

  for(p in pIdx)
    for( i in 1:length(fIdx) )
    {
      f <- fIdx[i]
      plot( x = range(yrs[tdx]),
            y = c(0,max(E_qpft[,,f,],na.rm = T) ),
            type = "n", axes = F )
        mfg <- par("mfg")
        # Plot axes and facet labels
        if( mfg[1] == mfg[3] )
          axis( side = 1 )
        axis( side = 2, las = 1 )
        if( mfg[2] == mfg[4] )
        {
          rmtext( txt =  fleetNames[f], font = 2, line = .05,
                  cex = 1.5, outer = TRUE )
        }
        box()
        grid()

        polygon(  x = c(yrs[tdx],rev(yrs[tdx])), y = c(E_qpft[1,p,f,tdx],rev(E_qpft[3,p,f,tdx])),
                  col = scales::alpha(fleetCols[i], alpha = .3), border = NA )

        for( tIdx in traces )
          lines( x = yrs[tdx], y = E_ipft[tIdx,p,f,tdx], col = fleetCols[i], lwd = .8 )

        lines( x = yrs[tdx], y = E_qpft[2,p,f,tdx], col = fleetCols[i], lwd = 3)
        abline( v = yrs[tMP] - 0.5, lty = 2, lwd = 0.5 )
    }
  mtext( side = 2, outer = TRUE, text = "Predator Effort (Biomass/Abundance)",
          line = 2 )

  mtext( side = 1, outer = TRUE, text = "Year",
          line = 2 )
} # END plotTulipEffort_p()


# plotTulipCatch_pft()
# Catch over time gridded
# by stock area and fleet (right now predators) - envelopes
plotTulipCatch_pft <- function( obj = blob,
                                nTrace = 3,
                                pIdx = 1,
                                fIdx = 7:13 )
{
  goodReps <- obj$goodReps

  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP
  nT      <- obj$om$nT
  nF      <- obj$om$nF
  nReps   <- dim(obj$om$C_ispft)[1]


  C_ipft <- array(0,dim = c(length(goodReps),nP,nF,nT) )

  C_ipft[,pIdx,fIdx,] <- obj$om$C_ispft[goodReps,1,pIdx,fIdx,]


  C_ipft[C_ipft == 0] <- NA



  C_qpft <- apply(  X = C_ipft, FUN = quantile,
                    MARGIN = c(2,3,4), probs = c(0.025, 0.5, 0.975),
                    na.rm = T )

  C_qpft[C_qpft == 0] <- NA

  traces <- sample( 1:nReps, size = min(nTrace,nReps)  )

  speciesNames  <- obj$om$speciesNames
  stockNames    <- obj$om$stockNames
  fleetNames    <- obj$om$fleetNames
  fYear         <- obj$ctlList$opMod$fYear

  yrs <- seq( from = fYear, by = 1, length.out = nT)

  fleetCols <- RColorBrewer::brewer.pal( length(fIdx), "Dark2" )

  par(  mfcol = c(length(fIdx),nP),
        mar = c(1,1.5,1,1.5),
        oma = c(3,4,3,3) )

  for(p in pIdx)
    for( i in 1:length(fIdx) )
    {
      f <- fIdx[i]
      plot( x = range(yrs),
            y = c(0,max(C_qpft[,,f,],na.rm = T) ),
            type = "n", axes = F )
        mfg <- par("mfg")
        # Plot axes and facet labels
        if( mfg[1] == mfg[3] )
          axis( side = 1 )
        axis( side = 2, las = 1 )
        if( mfg[2] == mfg[4] )
        {
          corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
          par(xpd = TRUE) #Draw outside plot area
          text(x = corners[2]+2, y = mean(corners[3:4]), stockNames[p], srt = 270,
                font = 2, cex = 1.5 )
          par(xpd = FALSE)
        }
        box()
        grid()

        polygon(  x = c(yrs,rev(yrs)), y = c(C_qpft[1,p,f,],rev(C_qpft[3,p,f,])),
                  col = scales::alpha(fleetCols[i], alpha = .3), border = NA )

        for( tIdx in traces )
          lines( x = yrs, y = C_ipft[tIdx,p,f,], col = fleetCols[i], lwd = .8 )

        lines( x = yrs, y = C_qpft[2,p,f,], col = fleetCols[i], lwd = 3)
        abline( v = yrs[tMP] - 0.5, lty = 2, lwd = 0.5 )
    }
  mtext( side = 2, outer = TRUE, text = "Trawl Effort (fishing hours?)",
          line = 2, font = 2)
} # END plotTulipCatch_pft()

# plotTulipPondedBio_pft()
# Catch over time gridded
# by stock area and fleet (right now predators) - envelopes
plotTulipPondedBio_pft <- function( obj = blob,
                                    nTrace = 3,
                                    pIdx = 1 )
{
  goodReps <- obj$goodReps

  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP
  nT      <- obj$om$nT
  nF      <- obj$om$nF
  nReps   <- dim(obj$om$C_ispft)[1]

  # Get SOK fleets
  browser()
  sokFleets <- c(6,13)

  P_ipft <- array(0,dim = c(length(goodReps),nP,nF,nT) )

  P_ipft[,pIdx,sokFleets,] <- obj$om$P_ispft[goodReps,1,pIdx,sokFleets,]


  P_ipft[P_ipft == 0] <- NA



  P_qpft <- apply(  X = P_ipft, FUN = quantile,
                    MARGIN = c(2,3,4), probs = c(0.025, 0.5, 0.975),
                    na.rm = T )

  P_qpft[P_qpft == 0] <- NA

  traces <- sample( 1:nReps, size = min(nTrace,nReps)  )

  speciesNames  <- obj$om$speciesNames
  stockNames    <- obj$om$stockNames
  fleetNames    <- obj$om$fleetNames
  fYear         <- obj$ctlList$opMod$fYear

  yrs <- seq( from = fYear, by = 1, length.out = nT)

  fleetCols <- RColorBrewer::brewer.pal( length(sokFleets), "Dark2" )

  par(  mfcol = c(length(sokFleets),nP),
        mar = c(1,1.5,1,1.5),
        oma = c(3,4,3,3) )

  for(p in pIdx)
    for( i in 1:length(sokFleets) )
    {
      f <- sokFleets[i]
      plot( x = range(yrs),
            y = c(0,max(P_qpft[,,f,],na.rm = T) ),
            type = "n", axes = F )
        mfg <- par("mfg")
        # Plot axes and facet labels
        if( mfg[1] == mfg[3] )
          axis( side = 1 )
        axis( side = 2, las = 1 )
        if( mfg[2] == mfg[4] )
        {
          corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
          par(xpd = TRUE) #Draw outside plot area
          text(x = corners[2]+2, y = mean(corners[3:4]), stockNames[p], srt = 270,
                font = 2, cex = 1.5 )
          par(xpd = FALSE)
        }
        box()
        grid()

        polygon(  x = c(yrs,rev(yrs)), y = c(P_qpft[1,p,f,],rev(P_qpft[3,p,f,])),
                  col = scales::alpha(fleetCols[i], alpha = .3), border = NA )

        for( tIdx in traces )
          lines( x = yrs, y = P_ipft[tIdx,p,f,], col = fleetCols[i], lwd = .8 )

        lines( x = yrs, y = P_qpft[2,p,f,], col = fleetCols[i], lwd = 3)
        abline( v = yrs[tMP] - 0.5, lty = 2, lwd = 0.5 )
    }
  mtext( side = 2, outer = TRUE, text = "Trawl Effort (fishing hours?)",
          line = 2, font = 2)
} # END plotTulipPondedBio_pft()



# plotTulipSOKEffort_p()
# Effort over time gridded
# by stock area - envelopes
plotTulipSOKEffort_p <- function( obj = blob,
                                  nTrace = 3,
                                  proj = TRUE,
                                  traces = NULL )
{
  goodReps <- obj$goodReps

  E_ipft <- obj$mp$hcr$sokEff_ispft[goodReps,1,,,]

  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP
  nT      <- obj$om$nT
  nF      <- obj$om$nF
  nReps   <- dim(E_ipft)[1]

  Eagg_ift <- apply(X = E_ipft, FUN = sum, MARGIN = c(1,3,4), na.rm = T)

  Etot_ipft <- array(0, dim = dim(E_ipft) + c(0,1,0,0) )
  Etot_ipft[,1:nP,,] <- E_ipft
  Etot_ipft[,nP+1,,] <- Eagg_ift


  E_ipft <- Etot_ipft


  # Pull scenario and MP labels
  scenarioName <- obj$ctlList$ctl$scenarioName
  mpName <- obj$ctlList$ctl$mpName

  stamp <- paste0(scenarioName,":",mpName)

  E_qpft <- apply(  X = E_ipft, FUN = quantile,
                    MARGIN = c(2,3,4), probs = c(0.025, 0.5, 0.975),
                    na.rm = TRUE )

  if( is.null(traces))
    traces <- sample( 1:nReps, size = min(nTrace,nReps)  )


  speciesNames  <- obj$om$speciesNames
  stockNames    <- c("C/S","JP/S","Lou","Aggregate")
  fleetNames    <- obj$om$fleetNames
  fYear         <- obj$ctlList$opMod$fYear

  yrs <- seq( from = fYear, by = 1, length.out = nT)

  if(proj)
    tdx <- (tMP - 1):nT
  else tdx <- 1:nT

  fleetCols <- RColorBrewer::brewer.pal( nF, "Dark2" )

  par(  mfcol = c(nP+1,1),
        mar = c(1,1.5,1,1.5),
        oma = c(3,4,3,3) )

  for( p in 1:(nP+1) )
  {
    plot( x = range(yrs[tdx]),
          y = c(0,max(E_qpft[,p,,tdx],na.rm = T) ),
          type = "n", axes = F )
      mfg <- par("mfg")
      # Plot axes and facet labels
      if( mfg[1] == mfg[3] )
        axis( side = 1 )
      axis( side = 2, las = 1 )
      box()
      grid()
      for( f in 6:7 )
      {
        # Skip if the fleet doesn't fish
        if( all(E_qpft[,,f,] == 0) )
          next

        polygon( x = c(yrs,rev(yrs)), y = c(E_qpft[1,p,f,],rev(E_qpft[3,p,f,])),
                  col = scales::alpha(fleetCols[f], alpha = .3), border = NA )
          for( tIdx in traces )
            lines( x = yrs, y = E_ipft[tIdx,p,f,], col = fleetCols[f], lwd = .8 )

        lines( x = yrs, y = E_qpft[2,p,f,], col = fleetCols[f], lwd = 3)

      }
      if( mfg[2] == mfg[4] )
      {
        rmtext( line = 0.25, txt = stockNames[p], font = 2, outer = FALSE, cex = 1.5 )
      }
      abline( v = yrs[tMP] - 0.5, lty = 2, lwd = 0.5 )
  }
  mtext( side = 2, outer = TRUE, text = "Active SOK licenses",
          line = 2, font = 2)
  mtext( side = 1, text = stamp, col = "grey60", outer = TRUE,
          line = 1, adj = .75, cex = .5 )
} # END plotTulipSOKEffort_p()

# rmtext()
# Refactored procedure to plot right hand inner
# margin mtext with the bottom towards the middle
# of the plot
rmtext <- function( line = 1,
                    txt = "Sample",
                    font = 1,
                    cex = 1,
                    outer = FALSE,
                    yadj = .5)
{
  corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
  if( outer )
    par(xpd = NA) #Draw outside the figure region
  if( !outer )
    par(xpd = TRUE)

  xRange <- corners[2] - corners[1]

  text( x = corners[2]+line,
        y = yadj * sum(corners[3:4]),
        labels = txt, srt = 270,
        font = font, cex = cex )
  par(xpd = FALSE)
} # END rmtext()


# plotEffort_p()
# Effort over time gridded
# by stock area
plotEffort_p <- function( obj = blob,
                          iRep = 1 )
{
  E_pft <- obj$om$E_ipft[iRep,,,]

  E_pft[E_pft == 0] <- NA

  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP
  nT      <- obj$om$nT
  nF      <- obj$om$nF

  speciesNames  <- obj$om$speciesNames
  stockNames    <- obj$om$stockNames
  fleetNames    <- obj$om$fleetNames
  fYear         <- obj$ctlList$opMod$fYear

  yrs <- seq( from = fYear, by = 1, length.out = nT)

  fleetCols <- RColorBrewer::brewer.pal( nF, "Dark2" )

  par(  mfcol = c(nP,1),
        mar = c(1,1.5,1,1.5),
        oma = c(3,4,3,3) )

  for( p in 1:nP )
  {
    plot( x = range(yrs),
          y = c(0,max(E_pft[p,,],na.rm = T) ),
          type = "n", axes = F )
      mfg <- par("mfg")
      # Plot axes and facet labels
      if( mfg[1] == mfg[3] )
        axis( side = 1 )
      axis( side = 2, las = 1 )
      if( mfg[2] == mfg[4] )
      {
        corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
        par(xpd = TRUE) #Draw outside plot area
        text(x = corners[2]+2, y = mean(corners[3:4]), stockNames[p], srt = 270,
              font = 2, cex = 1.5 )
        par(xpd = FALSE)
      }
      box()
      grid()
      for( f in 1:nF )
        lines( x = yrs, y = E_pft[p,f,], col = fleetCols[f], lwd = 3 )

      abline( v = yrs[tMP] - 0.5, lty = 2, lwd = 0.5 )
  }
  mtext( side = 2, outer = TRUE, text = "Trawl Effort (fishing hours?)",
          line = 2, font = 2)
} # END plotEffort_p()


# plotCtTACt_sp()
# Comparison plot of TAC and realized catch for each
# species/stock
plotCtTACt_sp <- function(  obj = blob,
                            iRep = 1,
                            fleets = 1:2 )
{
  C_spft   <- obj$om$C_ispft[iRep,,,,]
  TAC_spft <- obj$mp$hcr$TAC_ispft[iRep,,,,]

  C_spft[C_spft == 0] <- NA

  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP
  nT      <- obj$om$nT
  nF      <- obj$om$nF

  speciesNames  <- obj$om$speciesNames
  stockNames    <- obj$om$stockNames
  fleetNames    <- obj$om$fleetNames
  fYear         <- obj$ctlList$opMod$fYear

  yrs <- seq( from = fYear, by = 1, length.out = nT)

  fleetCols <- RColorBrewer::brewer.pal( nF, "Dark2" )


  par(  mfcol = c(nP,nS),
        mar = c(1,1.5,1,1.5),
        oma = c(3,3,3,3) )

  for(s in 1:nS)
    for( p in 1:nP )
    {
      plot( x = range(yrs),
            y = c(0,max(C_spft[s,p,,], TAC_spft[s,p,,], na.rm = T) ),
            type = "n", axes = F )

      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1 )
      if( mfg[1] == 1 )
        mtext( side = 3, text = speciesNames[s], font = 2, line = 0 )
      axis( side = 2, las = 1 )
      if( mfg[2] == mfg[4] )
      {
        corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
        par(xpd = TRUE) #Draw outside plot area
        text(x = corners[2]+2, y = mean(corners[3:4]), stockNames[p], srt = 270,
              font = 2, cex = 1.5 )
        par(xpd = FALSE)
      }
      box()
      grid()
      for( f in fleets )
      {
        rect( xleft = yrs - .3,
              xright = yrs + .3,
              ybottom = 0,
              ytop = TAC_spft[s,p,f,],
              border = NA, col = fleetCols[f] )

        lines(  x = yrs[], y = C_spft[s,p,f,],
                lty = 2, lwd = 2, col = "grey50" )
      }
      abline( v = yrs[tMP] - 0.5, lty = 2, lwd = 0.5 )
    }
    mtext( outer = TRUE, side = 2, text = "Catch and TAC (kt)" )
} # END plotCtTACt_sp()




# plotBtCtRt_p()
# Biomass and catch plots
# by stock for the historical
# and projection period
plotBtCtRt_p <- function( obj = blob, iRep = 1, sIdx=1)
{
  SB_ispt  <- obj$om$SB_ispt
  B_ispt   <- obj$om$B_ispt
  C_ispt   <- obj$om$C_ispt
  R_ispt   <- obj$om$R_ispt
  TAC_ispt <- obj$mp$hcr$TAC_ispt

  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP
  nT      <- obj$om$nT


  # stockNames    <- obj$om$stockNames
  stockNames    <- dimnames(obj$ctlList$opMod$histRpt$I_pgt)[[1]]
  fYear         <- obj$ctlList$opMod$fYear

  yrs <- seq( from = fYear, by = 1, length.out = nT)

  par(  mfcol = c(3,nP),
        mar = c(1,1.5,1,1.5), mgp=c(1,0.4,0),
        oma = c(2,3,0,0), tck=-.01 )

    for( p in 1:nP )
    {

      # biomass plots
      plot( x = range(yrs),
            y = c(0,max(SB_ispt[iRep,sIdx,p,]) ),
            type = "n", axes = T, las=1, xlab='', ylab='' )

      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1 )
      if( mfg[1] == 1 )
        mtext( side = 3, text = stockNames[p], font = 2, line = 0 )
      # axis( side = 2, las = 1 )

      box()
      grid()
      lines( x = yrs, y = SB_ispt[iRep,sIdx,p,], col = "red", lwd = 1 )
      abline( v = yrs[tMP] - 0.5, lty = 2, lwd = 0.5 )

      # catch plots
      plot( x = range(yrs),
            y = c(0,max(C_ispt[iRep,sIdx,p,],na.rm=T) ),
            type = "n", axes = T, las=1, xlab='', ylab='' )
      # axis( side = 2, las = 1 )
      box()
      grid()
      lines( x = yrs, y = C_ispt[iRep,sIdx,p,], col = "black", lwd = 1 )
      abline( v = yrs[tMP] - 0.5, lty = 2, lwd = 0.5 )

      # recruitment plots
      plot( x = range(yrs),
            y = c(0,max(R_ispt[iRep,sIdx,p,]) ),
            type = "n", axes = T, las=1, xlab='', ylab='' )
      # axis( side = 2, las = 1 )
      box()
      grid()
      lines( x = yrs, y = R_ispt[iRep,sIdx,p,], col = "orange", lwd = 1 )
      abline( v = yrs[tMP] - 0.5, lty = 2, lwd = 0.5 )

    }

  mtext( side =2, outer = TRUE, text = "Stock SB (kt), catch (kt), and age-1 recruitment (1000s)", line = 1.5)

}

# plotBtCtRt_p()
# Biomass and catch plots
# by stock for the historical
# and projection period
plotBtCtRtMt_p <- function( obj = blob, iRep = 1, sIdx=1)
{
  SB_ispt  <- obj$om$endSB_ispt
  B_ispt   <- obj$om$B_ispt
  C_ispt   <- obj$om$C_ispt
  C_ispft  <- obj$om$C_ispft
  R_ispt   <- obj$om$R_ispt
  M_iaxspt <- obj$om$M_iaxspt

  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP
  nT      <- obj$om$nT

  # replace zeros with NAs for Mortality object
  M_iaxspt[M_iaxspt==0] <- NA

  # Get predator fleets
  predGears <- obj$ctlList$opMod$predGears
  predC_ispt <- array(NA,dim = dim(C_ispt))

  if(length(predGears) > 0 )
  {
    C_ispt <- apply( X = C_ispft[,,,-predGears,,drop = FALSE],
                     FUN = sum, MARGIN = c(1:3,5), na.rm = T )
    predC_ispt <- apply(  X = C_ispft[,,,predGears,,drop = FALSE],
                          FUN = sum, MARGIN = c(1:3,5), na.rm = T )

    # Need to add "dead" ponded fish here
  }

  # stockNames    <- obj$om$stockNames
  stockNames    <- dimnames(obj$ctlList$opMod$histRpt$I_pgt)[[1]]
  fYear         <- obj$ctlList$opMod$fYear

  yrs <- seq( from = fYear, by = 1, length.out = nT)

  par(  mfcol = c(4,nP),
        mar = c(1,1.5,1,1.5), mgp=c(1,0.4,0),
        oma = c(2,3,0.5,0), tck=-.01 )

    for( p in 1:nP )
    {

      # biomass plots
      plot( x = range(yrs),
            y = c(0,max(SB_ispt[iRep,sIdx,p,]) ),
            type = "n", axes = T, las=1, xlab='', ylab='' )
        mfg <- par("mfg")
        if( mfg[1] == mfg[3] )
          axis( side = 1 )
        if( mfg[1] == 1 )
          mtext( side = 3, text = stockNames[p], font = 2, line = 0 )
        # axis( side = 2, las = 1 )
        if(p==1)
          mtext( side = 2, text = 'Spawning Biomass (kt)', line=2.5, cex=0.8 )

        box()
        grid()
        lines( x = yrs, y = SB_ispt[iRep,sIdx,p,], col = "red", lwd = 2 )
        abline( v = yrs[tMP] - 0.5, lty = 2, lwd = 1 )
      # legend('topright',bty='n', cex=0.8,
      #        legend='SBt',lty=1, col='red')

      # catch plots
      plot( x = range(yrs),
            y = c(0,max(C_ispt[iRep,sIdx,p,],
                        predC_ispt[iRep,sIdx,p,],na.rm=T) ),
            type = "n", axes = T, las=1, xlab='', ylab='' )
      # axis( side = 2, las = 1 )
        if(p==1)
          mtext( side = 2, text = 'Catch (kt)', line=2.5, cex=0.8 )
        box()
        grid()
        lines( x = yrs, y = C_ispt[iRep,sIdx,p,], col = "black", lwd = 1 )
        abline( v = yrs[tMP] - 0.5, lty = 2, lwd = 1 )
        if(length(predGears) > 0)
          lines(  x = yrs,
                  y = predC_ispt[iRep,sIdx,p,], col = "salmon", lwd = 2)
      # legend('topright',bty='n', cex=0.8,
      #        legend='Ct',lty=1, col='black')

      # recruitment plots
      plot( x = range(yrs),
            y = c(0,max(R_ispt[iRep,sIdx,p,]) ),
            type = "n", axes = T, las=1, xlab='', ylab='' )
      # axis( side = 2, las = 1 )
      if(p==1)
        mtext( side = 2, text = 'Recruits (1000s)', line=2.5, cex=0.8 )

      box()
      grid()
      lines( x = yrs, y = R_ispt[iRep,sIdx,p,], col = "orange", lwd = 1 )
      abline( v = yrs[tMP] - 0.5, lty = 2, lwd = 1 )
      legend('topleft',bty='n', cex=0.8,
             legend=c('Age 1 Rt'),lty=1, col='orange')


      # natural mortality plots
      plot( x = range(yrs),
            y = c(0,max(M_iaxspt[iRep,,1,sIdx,p,]) ),
            type = "n", axes = T, las=1, xlab='', ylab='' )
      # axis( side = 2, las = 1 )
      if(p==1)
        mtext( side = 2, text = 'Natural Mortality', line=2.5, cex=0.8 )

      box()
      grid()

      # age1 and age2 mortality
      M1 <- M_iaxspt[iRep,1,1,sIdx,p,]
      M2 <- M_iaxspt[iRep,2,1,sIdx,p,]
      M1[M1==0] <- NA
      M2[M2==0] <- NA

      lines( x = yrs, y = M1, col = "blue", lwd = 1, lty=3 )
      lines( x = yrs, y = M2, col = "blue", lwd = 1, lty=1 )
      abline( v = yrs[tMP] - 0.5, lty = 2, lwd = 1 )

      legend('topleft',bty='n', cex=0.8,
              legend=c('age 1 M',"age 2 M"),
              lty=c(3,1), col='blue')

    }

  # mtext( side =2, outer = TRUE, text = "Stock SB (kt), catch (kt), age-1 recruitment (1000s), and natural mortality", line = 1.5)

}


# plotBtCt_sp()
# Biomass and catch plots
# by species/stock for the historical
# and projection period
plotBtCt_sp <- function(  obj = blob,
                          iRep = 1 )
{
  SB_ispt  <- obj$om$SB_ispt
  B_ispt   <- obj$om$B_ispt
  C_ispt   <- obj$om$C_ispt
  TAC_ispt <- obj$mp$hcr$TAC_ispt

  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP
  nT      <- obj$om$nT

  speciesNames  <- obj$om$speciesNames
  stockNames    <- obj$om$stockNames
  fYear         <- obj$ctlList$opMod$fYear

  yrs <- seq( from = fYear, by = 1, length.out = nT)

  par(  mfcol = c(nP,nS),
        mar = c(1,1.5,1,1.5),
        oma = c(3,3,3,3) )

  for(s in 1:nS)
    for( p in 1:nP )
    {
      plot( x = range(yrs),
            y = c(0,max(B_ispt[iRep,s,p,]) ),
            type = "n", axes = F )

      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1 )
      if( mfg[1] == 1 )
        mtext( side = 3, text = speciesNames[s], font = 2, line = 0 )
      axis( side = 2, las = 1 )
      if( mfg[2] == mfg[4] )
      {
        corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
        par(xpd = TRUE) #Draw outside plot area
        text(x = corners[2]+2, y = mean(corners[3:4]), stockNames[p], srt = 270,
              font = 2, cex = 1.5 )
        par(xpd = FALSE)
      }
      box()
      grid()
      rect( xleft = yrs - .3, xright = yrs + .3,
            ybottom = 0, ytop = C_ispt[iRep,s,p,], col = "grey70", border = NA )
      segments( x0 = yrs, x1 = yrs,
                y0 = 0, y1 = TAC_ispt[iRep,s,p,], col = "black" )
      lines( x = yrs, y = SB_ispt[iRep,s,p,], col = "red", lwd = 1 )
      lines( x = yrs, y = B_ispt[iRep,s,p,], col = "black", lwd = 1, lty = 2 )


      abline( v = yrs[tMP] - 0.5, lty = 2, lwd = 0.5 )

    }


  mtext( side =1, outer = TRUE, text = "Year", line = 2)
  mtext( side =2, outer = TRUE, text = "Stock biomass and catch (kt)", line = 1.5)

}

# plotRefPtSeries()
# A copy-forward of the mseR plot of the same
# name. Shows a time-series of control points
# as biomass estimates, true SSB, and target
# vs effective harvest rates.
plotRefPtSeries <- function(  obj = blob,
                              iRep = 1,
                              sIdx = 1, pIdx = 1,
                              noPar = FALSE )
{
  # Get dims
  tMP     <- obj$om$tMP
  pT      <- obj$ctlList$opMod$pT
  nT      <- obj$om$nT
  nF      <- obj$om$nF
  fYear   <- obj$ctlList$opMod$fYear

  # Pull SSB and U
  SB_t    <- obj$om$SB_ispt[iRep,sIdx,pIdx,]

  # Now pull retrospective values
  retroSB_t     <- diag(obj$mp$assess$retroSB_itspt[iRep,1:pT,sIdx,pIdx,tMP:nT])
  # retroUmsy_t   <- obj$mp$assess$retroUmsy_itsp[iRep,,sIdx,pIdx]
  # retroBmsy_t   <- obj$mp$assess$retroBmsy_itsp[iRep,,sIdx,pIdx]

  simLabel <- paste(obj$ctlList$ctl$scenarioName,obj$ctlList$ctl$mpName,sep=":")

  # Calculate U
  fleetType_f <- obj$om$fleetType_f
  sokFleets   <- which(fleetType_f %in% 2:3)
  predGears   <- obj$ctlList$opMod$predGears

  if(is.null(predGears))
    predGears <- c()

  commGears   <- (1:nF)[!1:nF %in% predGears]

  # Catch
  C_ft        <- obj$om$C_ispft[iRep,sIdx,pIdx,,]
  C_t         <- apply(X = C_ft[-sokFleets,,drop = FALSE], FUN = sum, MARGIN = c(2))
  predC_t     <- apply(X = C_ft[predGears[!predGears %in% sokFleets],,drop = FALSE], FUN = sum, MARGIN = c(2))
  commC_t     <- apply(X = C_ft[commGears[!commGears %in% sokFleets],,drop = FALSE], FUN = sum, MARGIN = c(2))

  # Dead ponded fish
  dP_ft       <- obj$om$dP_ispft[iRep,sIdx,pIdx,,]
  dP_t        <- apply(X = dP_ft[sokFleets,,drop = FALSE], FUN=sum, MARGIN=c(2))
  predDP_t    <- apply(dP_ft[intersect(sokFleets,predGears),,drop=FALSE], FUN=sum, MARGIN=c(2))
  commDP_t    <- apply(dP_ft[intersect(sokFleets,commGears),,drop=FALSE], FUN=sum, MARGIN=c(2))

  # Total oonded fish
  P_ft        <- obj$om$P_ispft[iRep,sIdx,pIdx,,]
  P_t         <- apply(X = P_ft[sokFleets,,drop = FALSE], FUN=sum, MARGIN=c(2))
  predP_t     <- apply(P_ft[intersect(sokFleets,predGears),,drop=FALSE], FUN=sum, MARGIN=c(2))
  commP_t     <- apply(P_ft[intersect(sokFleets,commGears),,drop=FALSE], FUN=sum, MARGIN=c(2))



  # Calculate aggregate harvest rate
  # commU_ispt   <- commC_ispt/(SB_ispt + commC_ispt + commP_ispt + predP_ispt + predC_ispt)
  # predU_ispt   <- predC_ispt/(SB_ispt + commC_ispt + commP_ispt + predP_ispt + predC_ispt)
  commUnoPred_t   <- (commC_t + commDP_t)/(SB_t + commC_t + commDP_t)
  commUincPred_t  <- (commC_t + commDP_t)/(SB_t + commC_t + commDP_t + predC_t)
  predU_t         <- predC_t/(SB_t + commC_t + commP_t + predC_t)

  # commU_t     <- commU_spt[sIdx,pIdx,]

  # Need CPs and Fref/targetF
  LCP_t     <- obj$mp$hcr$LCP_ispt[iRep,sIdx,pIdx,]
  UCP_t     <- obj$mp$hcr$UCP_ispt[iRep,sIdx,pIdx,]
  Fref_t    <- obj$mp$hcr$Fref_ispt[iRep,sIdx,pIdx,]
  Bref_t    <- obj$mp$hcr$Bref_ispt[iRep,sIdx,pIdx,]
  targF_t   <- obj$mp$hcr$targetF_ispt[iRep,sIdx,pIdx,]


  # Full years
  yrs     <- seq( from = fYear, by = 1, length.out = nT)
  if(!noPar)
    par(  mfrow = c(2,1),
          oma = c(3,3,2,2),
          mar = c(.1,2,.1,2) )
  # First plot B series, over projection period
  maxB <- max(SB_t[tMP:nT],UCP_t[tMP:nT],retroSB_t,na.rm = T)
  plot( x = range(yrs[tMP:nT]), y = c(0,maxB),
        axes = FALSE, type = "n")
    mfg <- par("mfg")
    mtext(side = 3, text = simLabel, font = 2)
    axis( side = 2, las = 1)
    if(mfg[2] == 1)
      mtext( side = 2, text = "OM and AM Biomass (kt)", line = 3)
    grid()
    box()
    lines(  x = yrs, y = SB_t, col = "salmon",
            lwd = 2 )
    points( x = yrs, y = LCP_t, pch = 16, col = "red" )
    points( x = yrs, y = UCP_t, pch = 16, col = "darkgreen" )
    # points( x = yrs[tMP:nT], y = retroBmsy_t, pch = 16, col = "royalblue" )
    lines(  x = yrs[tMP:nT], y = retroSB_t, col = "grey40", lwd = 2 )

  maxU <- max(commUnoPred_t[tMP:nT],commUincPred_t[tMP:nT],Fref_t[tMP:nT],na.rm = T)
  plot( x = range(yrs[tMP:nT]), y = c(0,maxU),
        axes = FALSE, type = "n")
    mfg <- par("mfg")
    axis( side = 2, las = 1)
    if(mfg[2] == 1)
      mtext( side = 2, text = "Harvest Rate", line = 3)
    if(mfg[1] == mfg[3])
    {
      axis( side = 1 )
    }
    if(!noPar)
      mtext( side = 1, text = "Year", line = 3)

    grid()
    box()
    segments( x0 = yrs,
              y1 = Fref_t,
              y0 = targF_t,
              col = "grey40",
              lwd = 2 )

    lines( x = yrs, y = commUincPred_t, lwd = 2, col = "black")
    lines( x = yrs, y = commUnoPred_t, lwd = 2, col = "royalblue")
    # lines( x = yrs, y = legU_t, lwd = 1, col = "royalblue")
    points( x = yrs, y = Fref_t,
            col = "salmon", pch = 16 )
    points( x = yrs, y = targF_t,
            col = "salmon", pch = 21, bg = NA )


}


plotMultiRefPtSeries <- function(objList = list(blob),
                                  ... )
{
  nSims <- length(objList)

  par(  mfcol = c(2,nSims),
        oma = c(5,3,2,2),
        mar = c(.1,2,.1,2) )

  for(k in 1:nSims)
    plotRefPtSeries(  obj = objList[[k]],
                      noPar = TRUE, ... )

  mtext(side = 1, text = "Year", outer = TRUE, line = 2)
}



plotMultiRetroMt <- function( groupFolder = "data/siscaMPs_2023_10_28",
                              iRep = 1 )
{
  # First, we'll go over the sims in the groupFolder
  # and arrange by scenario/MP. Use the same
  # routine we used to make the weighted blobs

  # Now read the infoFile in each sim folder
  dirList   <- list.dirs(groupFolder)
  dirList   <- dirList[grepl(x = dirList, pattern = "sim_")]

  infoList  <- file.path(dirList,"infoFile.txt")
  infoList  <- lapply(X = infoList, FUN = lisread)
  infoList  <- lapply(X = infoList, FUN = as.data.frame)
  info.df   <- do.call(rbind, infoList)

  # Now pull MPs and scenarios
  scenarios <- unique(info.df$scenario)
  mps       <- unique(info.df$mp)
  mps       <- mps[order(mps)]

  nScen <- length(scenarios)
  nMP   <- length(mps)

  par(  mfrow = c(nScen, nMP),
        oma   = c(4,5,3,4),
        mar   = c(.1,.1,.1,.1))

  for(sIdx in 1:nScen)
    for(mIdx in 1:nMP )
    {

      blobRow <- info.df |>
                  filter( scenario == scenarios[sIdx],
                          mp == mps[mIdx] )

      folderName <- blobRow$simLabel



      load(file.path(groupFolder,folderName,paste0(folderName,".Rdata")))

      plotRetroMt(obj = blob, iRep = iRep, noPar = TRUE)

      if(sIdx == 1)
        mtext(side = 3, text = mps[mIdx], font = 2)

      if(sIdx == 1 & mIdx == 1)
        legend( x = "topright", bty = "n",
                cex = 1.1,
                legend = c("OM Mt", "SISCAH estimates"),
                lwd = c(3,1),
                lty = c(1,1),
                pch = c(NA,NA),
                col = c("salmon","grey60"))

      legend(x = "topleft", legend = scenarios[sIdx], bty = "n")

    }

  mtext(side = 1, text = "Year", outer = TRUE, line = 3)
  mtext(side = 2, text = "Age-2 Natural Mortality (kt)", outer = TRUE, line = 3)
}


# plotRetroMt()
# Retrospective plots of AM fits for a given replicate
plotRetroMt <- function(  obj  = blob,
                          iRep = 1,
                          noPar = FALSE )
{
  nS  <- obj$om$nS
  nP  <- obj$om$nP
  nT  <- obj$om$nT
  tMP <- obj$om$tMP
  pT  <- nT - tMP + 1
  nF  <- obj$om$nF

  # Need to adjust this for blended index
  # blendIdx      <- obj$ctlList$opMod$blendIdx
  # SISCA AM is only provided blended index, i.e. obj$mp$data$I_spft[1,1:nP,5,histdx] <- repObj$combI_pt

  # sum( matAge_asp[,s,p] * vB_axspft[,nX,s,p,f,t])

  # Get OM mortality (Age-2) arrays
  M_spt         <- ClimProjDiags::Subset( x= obj$om$M_iaxspt[iRep,2,1,,,,drop=FALSE],
                                          along = 4:6, indices = list( 1:nS, 1:nP, 1:nT),
                                          drop = "non-selected" )
  F_spft        <- ClimProjDiags::Subset( x = obj$om$F_ispft[iRep,,,,,drop = FALSE],
                                          along = 2:5, indices = list(1:nS,1:nP,1:nF,1:nT),
                                          drop = "non-selected" )

  # Need to put in the age-2 predation mortality
  if(grepl("predM",obj$ctlList$ctl$scenarioName))
  {
    # browser()
    for(f in 7:12)
    {
      # Pull age-2selectivity
      s_t <- obj$ctlList$opMod$histRpt$sel_apgt[2,1,f,]
      # Pull F


      # Make quantile, add to M

      for(t in 1:nT)
      {
        if(t < tMP)
          M_spt[,,t] <- M_spt[,,t] + s_t[t] * F_spft[1:nS,1:nP,f,t]
        if(t >= tMP)
          M_spt[,,t] <- M_spt[,,t] + s_t[tMP-1] * F_spft[1:nS,1:nP,f,t]

      }


    }

  }


  retroM_tspt  <- ClimProjDiags::Subset( x = obj$mp$assess$retroM_itaxspt[iRep,,2,1,,,,drop = FALSE],
                                          along = c(2,5:7),
                                          indices = list(1:pT,1:nS,1:nP,1:nT),
                                          drop = "non-selected" )

  retroM_tspt[retroM_tspt < 0] <- NA


  speciesNames  <- obj$ctlList$opMod$species
  stockNames    <- obj$ctlList$opMod$stock
  fYear         <- obj$ctlList$opMod$fYear

  yrs <- seq( from = fYear, by = 1, length.out = nT)

  if(nS == 2)
    nS <- 1

  if(!noPar)
    par(  mfcol = c(nP,nS),
          mar = c(1,1.5,1,1.5),
          oma = c(3,3,3,3) )
  for(s in 1:nS)
    for( p in 1:nP )
    {
      maxM <- 1.2*max(M_spt[s,p,], na.rm = T)
      plot( x = range(yrs),
            y = c(0,maxM ),
            type = "n", axes = F )

      mfg <- par("mfg")
      if(noPar)
      {
        if( mfg[1] == mfg[3] )
          axis( side = 1 )
        if(mfg[2] == 1)
          axis( side = 2, las = 1 )

        if( mfg[2] == mfg[4] )
          axis(side = 4, las = 1 )
      }

      if(!noPar)
      {
        if( mfg[1] == mfg[3] )
          axis( side = 1 )
        if( mfg[1] == 1 )
          mtext( side = 3, text = speciesNames[s], font = 2, line = 0 )
        if(mfg[2] == 1)
          axis( side = 2, las = 1 )
        if( mfg[2] == mfg[4] )
          mtext( side = 4, text= stockNames[p], line = 1.5)
      }

      box()
      grid()

      lines( x = yrs, y = M_spt[s,p,], col = "salmon", lwd = 3 )

      # Plot retro fits
      for( tt in 1:pT )
      {
        # propTAC <- propTAC_spt[s,p,tMP + tt - 1]
        lines( x = yrs, y = retroM_tspt[tt,s,p,], col = "grey60", lwd = 1 )
      }

      abline( v = yrs[tMP] - 0.5, lty = 2, lwd = 0.5 )
    }

} # plotRetroMt()

plotMultiRetroSB <- function( groupFolder = "data/siscaMPs_2023_10_28",
                              iRep = 1,
                              plotData = TRUE )
{
  # First, we'll go over the sims in the groupFolder
  # and arrange by scenario/MP. Use the same
  # routine we used to make the weighted blobs

  # Now read the infoFile in each sim folder
  dirList   <- list.dirs(groupFolder)
  dirList   <- dirList[grepl(x = dirList, pattern = "sim_")]

  infoList  <- file.path(dirList,"infoFile.txt")
  infoList  <- lapply(X = infoList, FUN = lisread)
  infoList  <- lapply(X = infoList, FUN = as.data.frame)
  info.df   <- do.call(rbind, infoList)

  # Now pull MPs and scenarios
  scenarios <- unique(info.df$scenario)
  mps       <- unique(info.df$mp)
  mps       <- mps[order(mps)]

  nScen <- length(scenarios)
  nMP   <- length(mps)

  par(  mfrow = c(nScen, nMP),
        oma   = c(4,5,3,4),
        mar   = c(.1,.1,.1,.1))

  for(sIdx in 1:nScen)
    for(mIdx in 1:nMP )
    {

      blobRow <- info.df |>
                  filter( scenario == scenarios[sIdx],
                          mp == mps[mIdx] )

      folderName <- blobRow$simLabel



      load(file.path(groupFolder,folderName,paste0(folderName,".Rdata")))

      plotRetroSB(obj = blob, iRep = iRep, noPar = TRUE)

      if(sIdx == 1)
        mtext(side = 3, text = mps[mIdx], font = 2)

      if(sIdx == 1 & mIdx == 1)
        legend( x = "topright", bty = "n",
                cex = 1.1,
                legend = c("OM SB", "SISCA estimates","Spawn Index"),
                lwd = c(3,1,NA),
                lty = c(1,1,NA),
                pch = c(NA,NA,16),
                col = c("red","grey60","grey30"))

      legend(x = "topleft", legend = scenarios[sIdx], bty = "n")

    }

  mtext(side = 1, text = "Year", outer = TRUE, line = 3)
  mtext(side = 2, text = "Spawning Biomass (kt)", outer = TRUE, line = 3)
}


# plotRetroSB()
# Retrospective plots of AM fits for a given replicate
plotRetroSB <- function(  obj  = blob,
                          iRep = 1,
                          vB_f = 5,
                          plotTB = FALSE,
                          plotVB = FALSE,
                          plotData = TRUE,
                          noPar = FALSE  )
{
  nS  <- obj$om$nS
  nP  <- obj$om$nP
  nT  <- obj$om$nT
  tMP <- obj$om$tMP
  pT  <- nT - tMP + 1
  nF  <- obj$om$nF

  # Need to adjust this for blended index
  # blendIdx      <- obj$ctlList$opMod$blendIdx
  # SISCA AM is only provided blended index, i.e. obj$mp$data$I_spft[1,1:nP,5,histdx] <- repObj$combI_pt

  # sum( matAge_asp[,s,p] * vB_axspft[,nX,s,p,f,t])

  # Get biomass arrays
  SB_spt       <- ClimProjDiags::Subset(  x= obj$om$SB_ispt[iRep,,,,drop=FALSE],
                                          along = 2:4, indices = list( 1:nS, 1:nP, 1:nT),
                                          drop = "non-selected" )

  VB_spft      <- ClimProjDiags::Subset(  x = obj$om$vB_ispft[iRep,,,,,drop=FALSE],
                                          along = 2:5, indices = list(1:nS,1:nP,vB_f,1:nT),
                                          drop = "non-selected")

  totB_spt     <- ClimProjDiags::Subset(  x = obj$om$B_ispt[iRep,,,,drop=F],
                                          along = 2:4,
                                          indices = list(1:nS,1:nP,1:nT),
                                          drop = "non-selected")

  q_spft      <- ClimProjDiags::Subset( x = obj$om$q_ispft[iRep,,,,,drop = FALSE],
                                        along = 2:5,
                                        indices = list(1:nS,1:nP,1:nF,1:(nT)),
                                        drop = "non-selected" )


  I_spft       <- ClimProjDiags::Subset(  x = obj$mp$data$I_ispft[iRep,,,,,drop = FALSE],
                                          along = 2:5,
                                          indices = list(1:nS,1:nP,1:nF,1:(nT)),
                                          drop = "non-selected" )

  retroSB_tspt  <- ClimProjDiags::Subset( x = obj$mp$assess$retroSB_itspt[iRep,,,,,drop = FALSE],
                                          along = c(2:5),
                                          indices = list(1:pT,1:nS,1:nP,1:nT),
                                          drop = "non-selected" )

  retroSB_tspt[retroSB_tspt < 0] <- NA

  # Pull blended index stuff
  rI_spft       <- ClimProjDiags::Subset( x = obj$om$rI_ispft[iRep,,,,,drop = FALSE],
                                          along = 2:5,
                                          indices = list(1:nS,1:nP,1:nF,1:(nT)),
                                          drop = "non-selected" )


  omqComb_spt <- array(NA, dim=c(nS,nP,nT))

  # histFile <- obj$ctlList$opMod$histFile
  # load(here('history',histFile,paste0(histFile,'.Rdata')))
  # rI_pgt <- reports$data$rI_pgt
  # for(tIdx in 1:(tMP-1))
  #   omqComb_spt[1:nS,1,tIdx] <- sum(rI_pgt[1,,tIdx] * q_spft[1,1,,tIdx])

  # omqComb_spt[1:nS,1,tMP:nT] <- 1

  if(!is.null(obj$mp$assess$retroqComb_itspt))
    qComb_tspt <- ClimProjDiags::Subset( x = obj$mp$assess$retroqComb_itspt[iRep,,,,,drop = FALSE],
                                          along = 2:5,
                                          indices = list(1:pT,1:nS,1:nP,1:(nT)),
                                          drop = "non-selected" )

  # Get proportion of TACs for splitting aggregated biomass
  propTAC_spt   <- obj$mp$hcr$propTAC_ispt[iRep,,,]

  I_spft[I_spft < 0] <- NA


  speciesNames  <- obj$ctlList$opMod$species
  stockNames    <- obj$ctlList$opMod$stock
  fYear         <- obj$ctlList$opMod$fYear

  yrs <- seq( from = fYear, by = 1, length.out = nT)

  if(!noPar)
    par(  mfcol = c(nP,nS),
          mar = c(1,1.5,1,1.5),
          oma = c(3,3,3,3) )
  for(s in 1:nS)
    for( p in 1:nP )
    {
      maxSB <- 1.2*max(SB_spt[s,p,], plotTB * totB_spt[s,p,], plotVB*VB_spft[s,p,,], na.rm = T)
      plot( x = range(yrs),
            y = c(0,maxSB ),
            type = "n", axes = F )

      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1 )
      if( mfg[1] == 1 )
        mtext( side = 3, text = speciesNames[s], font = 2, line = 0 )
      axis( side = 2, las = 1 )
      if( mfg[2] == mfg[4] )
      {
        corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
        par(xpd = TRUE) #Draw outside plot area
        text(x = corners[2]+1.5, y = mean(corners[3:4]), stockNames[p], srt = 270,
              font = 2, cex = 1.5 )
        par(xpd = FALSE)
      }
      box()
      grid()
      lines( x = yrs, y = SB_spt[s,p,], col = "red", lwd = 3 )
      if( plotData )
      {
          for( f in 1:nF )
            points(x = yrs, y = I_spft[s,p,f, ]/q_spft[s,p,f,] )

        # if( blendIdx )
        #   points( x = yrs, y = I_spft[s,p,5,] / omqComb_spt[s,p,])
      }

      if( plotVB )
        for( f in vB_f )
          lines( x = yrs, y = VB_spft[s,p,f,], col = "grey40", lwd = 2, lty = 3 )
      if( plotTB )
        lines( x = yrs, y = totB_spt[s,p,], col = "black", lwd = 2 )

      # Plot retro fits
      for( tt in 1:pT )
      {
        # propTAC <- propTAC_spt[s,p,tMP + tt - 1]
        lines( x = yrs, y = retroSB_tspt[tt,s,p,], col = "grey60", lwd = 1 )
      }

      abline( v = yrs[tMP] - 0.5, lty = 2, lwd = 0.5 )
    }

} # plotRetroSB()

# plotBtIt_p()
# Plot Biomass and index
plotBtIt_p <- function( obj = blob, iRep = 1, f=5,
                      plotVB=FALSE, addCatch=TRUE,
                      parArg=TRUE, YlabOn= TRUE,
                      legdOn=TRUE)
{
  # Get biomass arrays
  SB_ispt     <- obj$om$SB_ispt
  VB_ispft    <- obj$om$vB_ispft
  I_ispt      <- obj$mp$data$I_ispft[,,,f,]
  # nS+1 and nP+1 for I_ispt represent aggregates

  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP
  nT      <- obj$om$nT

  # Hack to plot empirical assessment method for Herring projections
  if(addCatch)
  {
    C_ispt <- obj$om$C_ispt

    for(s in 1:nS)
      for(p in 1:nP)
        I_ispt[,s,p,]<- I_ispt[,s,p,]+ C_ispt[,s,p,]

  }

  speciesNames  <- obj$ctlList$opMod$species
  stockNames    <- dimnames(obj$ctlList$opMod$histRpt$I_pgt)[[1]]
  fYear         <- obj$ctlList$opMod$fYear
  pT            <- obj$ctlList$opMod$pT

  yrs <- seq( from = fYear, by = 1, length.out = nT)

  if(parArg)
  par(  mfcol = c(nP,nS),
        mar = c(1,1.5,1,1.5),
        oma = c(3,3,3,3) )
  for(s in 1:nS)
    for( p in 1:nP )
    {
      plot( x = range(yrs),
            y = c(0,max(VB_ispft[iRep,s,p,f,],
                        SB_ispt[iRep,s,p,],
                        I_ispt[iRep,s,p,], na.rm = T) ),
            type = "n", axes = F )

      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1 )
      if( mfg[1] == 1 )
        mtext( side = 3, text = speciesNames[s], font = 2, line = 0 )
      axis( side = 2, las = 1 )

      if( mfg[2] == mfg[4] )
      {
        corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
        par(xpd = TRUE) #Draw outside plot area
        text(x = corners[2]+1.5, y = mean(corners[3:4]), stockNames[p], srt = 270,font = 2, cex = 1.5 )
        par(xpd = FALSE)
      }
      box()
      grid()
      lines( x = yrs, y = SB_ispt[iRep,s,p,], col = "red", lwd = 2 )

      # plot vulnerable biomass
      if(plotVB)
      lines( x = yrs, y = VB_ispft[,s,p,,], col = "grey40", lwd = 2, lty = 3 )

      # plot index
      points( x = yrs, y = I_ispt[iRep,s,p,],
        bg = "green", pch=21, cex=1 )

      abline( v = yrs[tMP], lty = 2, lwd = 0.5 )

      if(s==1 & p==1 & !addCatch & legdOn)
      legend('topright',bty='n',
             legend=c('Spawning Biomass',
                      'Dive Survey Index'),
             lwd=c(2,NA), pch=c(NA,21),
             col=c('red','black'),
             pt.bg=c(NA, 'green'))

      if(s==1 & p==1 & addCatch & legdOn)
      legend('topright',bty='n',
             legend=c('Spawning Biomass',
                      'It + Ct'),
             lwd=c(2,NA), pch=c(NA,21),
             col=c('red','black'),
             pt.bg=c(NA, 'green'))
    }

  if(YlabOn)
  mtext( side =2, outer = TRUE, text = "Spawning biomass kt)", line = 1.5)

} # plotBtIt_p()



# plotRetroSBagg()
# Retrospective plots of aggregate SISCA AM fits (twoStep HCR) for a
# given replicate. Spatial pooling is always applied: OM SSB and catch
# are summed across PFMAs. One SISCA biomass trajectory is drawn per
# projection year, using blob$mp$assess$retroSB_agg_itt[iRep,,].
# If scaledIdx = TRUE, overlays the aggregated blended spawn index
# scaled by estimated combined q (I_agg / qComb) from the final
# projection year's assessment, as open circles.
plotRetroSBagg <- function( obj = blob, iRep = 1, Ct = TRUE,
                             scaledIdx = TRUE )
{
  tMP   <- obj$om$tMP
  nT    <- obj$om$nT
  nP    <- obj$om$nP
  nF    <- obj$om$nF
  pT    <- obj$ctlList$opMod$pT
  fYear <- obj$ctlList$opMod$fYear

  # Aggregate SISCA retro estimates: dim [pT, nT]
  retroSB_agg_tt        <- obj$mp$assess$retroSB_agg_itt[iRep,,]
  retroSB_agg_tt[retroSB_agg_tt <= 0] <- NA

  # OM aggregate SSB for Herring (s=1), summed over PFMAs: dim [nT]
  omSB_pt <- obj$om$SB_ispt[iRep, 1, , , drop = FALSE]
  dim(omSB_pt) <- c(nP, nT)
  omSB_t  <- colSums(omSB_pt)
  omSB_t[omSB_t == 0] <- NA

  # Coastwide catch and TAC for Herring (s=1):
  # fishing fleets only (exclude predators)
  fishGears <- 1:3
  sokGear   <- which(obj$om$fleetType_f == 2)
  allFishF  <- c(fishGears, sokGear)

  C_all     <- obj$om$C_ispft[iRep, 1, , , , drop = FALSE]
  dim(C_all) <- c(nP, nF, nT)
  C_fish    <- C_all[, allFishF, , drop = FALSE]
  C_t       <- apply(X = C_fish, MARGIN = 3,
                     FUN = sum, na.rm = TRUE)

  TAC_all     <- obj$mp$hcr$TAC_ispft[iRep, 1, , , ,
                                       drop = FALSE]
  dim(TAC_all) <- c(nP, nF, nT)
  TAC_fish    <- TAC_all[, allFishF, , drop = FALSE]
  TAC_t       <- apply(X = TAC_fish, MARGIN = 3,
                       FUN = sum, na.rm = TRUE)

  # Scaled indices: aggregated blended spawn index / combined q
  # Uses the index fleet and q from the final projection year's fit
  scaledI_t  <- NULL
  hasScaledIdx <- scaledIdx &&
                  !is.null(obj$mp$assess$retroqComb_agg_itt) &&
                  !is.null(obj$ctlList$mp$assess$idxFleets)
  if( hasScaledIdx )
  {
    idxFleet      <- obj$ctlList$mp$assess$idxFleets[1]
    # Aggregated combined index: sum blended index across PFMAs [nT]
    I_raw         <- obj$mp$data$I_ispft[
      iRep, 1, 1:nP, idxFleet, , drop = FALSE]
    dim(I_raw)    <- c(nP, nT)
    I_agg_t       <- colSums(I_raw)
    I_agg_t[I_agg_t <= 0] <- NA
    # Combined q from the final projection year's aggregate SISCA fit
    qComb_final_t <- obj$mp$assess$retroqComb_agg_itt[iRep, pT, ]
    qComb_final_t[qComb_final_t <= 0] <- NA
    scaledI_t     <- I_agg_t / qComb_final_t
  }

  yrs   <- seq(from = fYear, by = 1, length.out = nT)
  stamp <- paste(obj$ctlList$ctl$scenarioName, ":",
                 obj$ctlList$ctl$mpName, ":", iRep, sep = "")

  yMax <- max(omSB_t, retroSB_agg_tt, scaledI_t, na.rm = TRUE) * 1.2

  par( mar = c(4, 5, 2, 1) )
  plot( x = range(yrs), y = c(0, yMax),
        type = "n", axes = FALSE,
        xlab = "", ylab = "" )
  axis(side = 1)
  axis(side = 2, las = 1)
  box()
  grid()

  # Catch and TAC bars (projection period only)
  if( Ct )
  {
    rect( xleft  = yrs - 0.3, xright = yrs + 0.3,
          ytop   = C_t,       ybottom = 0,
          col    = "grey75",  border  = NA )
    rect( xleft  = yrs[tMP:nT] - 0.3, xright = yrs[tMP:nT] + 0.3,
          ytop   = TAC_t[tMP:nT],      ybottom = 0,
          col    = NA,                  border  = "black" )
  }

  # One SISCA trajectory per projection year
  for( pt in 1:pT )
    lines( x = yrs, y = retroSB_agg_tt[pt, ],
           col = "grey50", lwd = 1 )

  # Highlight most recent estimate
  lines( x = yrs, y = retroSB_agg_tt[pT, ],
         col = "dodgerblue", lwd = 2 )

  # OM truth
  lines( x = yrs, y = omSB_t, col = "red", lwd = 2 )

  # Scaled indices (I_agg / qComb from final year's fit)
  if( !is.null(scaledI_t) )
    points( x = yrs, y = scaledI_t,
            pch = 21, bg = "white", col = "darkgreen", cex = 0.9 )

  abline( v = yrs[tMP] - 0.5, lty = 2, lwd = 0.5 )

  # Reference lines
  B0_agg  <- sum(
    obj$ctlList$opMod$histRpt$refPts$refCurves$SBeq_pf[
      , 1])
  Bref    <- 17.35
  LRP_agg <- 0.3 * B0_agg
  abline(h = B0_agg,  lty = 2, col = "darkgreen", lwd = 1.5)
  abline(h = Bref,    lty = 2, col = "blue",     lwd = 1.5)
  abline(h = LRP_agg, lty = 2, col = "red",      lwd = 1.5)

  mtext(side = 1, text = "Year",                            line = 2.5, font = 2)
  mtext(side = 2, text = "Spawning Biomass (kt)", line = 3.5, font = 2)
  mtext(side = 3,
        text = "Retrospective estimates of WCVI Herring Spawning Biomass (SISCAH-DIM)",
        line = 0.5, font = 2)
  mtext(outer = FALSE, side = 1, adj = 1, line = 3.5, cex = 0.6,
        text = stamp, col = "grey60")

  refLeg <- list(
    legend = c(expression(B[0]), expression(B[ref]), "LRP"),
    col    = c("darkgreen", "blue", "red"),
    lwd    = c(1.5, 1.5, 1.5),
    lty    = c(2, 2, 2),
    pch    = c(NA, NA, NA))

  idxLegend <- if( !is.null(scaledI_t) )
                 list( legend = c("OM truth", "SISCA (final)", "SISCA (earlier)",
                                  "Catch", "Index/q",
                                  refLeg$legend),
                       col    = c("red", "dodgerblue", "grey50",
                                  "grey75", "darkgreen",
                                  refLeg$col),
                       lwd    = c(2, 2, 1, NA, NA,
                                  refLeg$lwd),
                       lty    = c(1, 1, 1, NA, NA,
                                  refLeg$lty),
                       pch    = c(NA, NA, NA, 15, 21,
                                  refLeg$pch),
                       pt.cex = c(NA, NA, NA, 1.5, 0.9,
                                  NA, NA, NA) )
               else
                 list( legend = c("OM truth", "SISCA (final)", "SISCA (earlier)",
                                  "Catch",
                                  refLeg$legend),
                       col    = c("red", "dodgerblue", "grey50", "grey75",
                                  refLeg$col),
                       lwd    = c(2, 2, 1, NA,
                                  refLeg$lwd),
                       lty    = c(1, 1, 1, NA,
                                  refLeg$lty),
                       pch    = c(NA, NA, NA, 15,
                                  refLeg$pch),
                       pt.cex = c(NA, NA, NA, 1.5,
                                  NA, NA, NA) )
  legend( "topleft", bty = "n",
          legend = idxLegend$legend,
          col    = idxLegend$col,
          lwd    = idxLegend$lwd,
          lty    = idxLegend$lty,
          pch    = idxLegend$pch,
          pt.cex = idxLegend$pt.cex )
}

# plotEBSBratio()
# Plots of ratio between vulnerable
# biomass and exploitable biomass
# for given depletion levels.
plotEBSBratio <- function( obj = blob )
{
  # Pull model dims
  nS <- obj$om$nS
  nP <- obj$om$nP
  nF <- obj$om$nF

  fleetNames <- obj$om$fleetNames

  # Pull reference points
  rp <- obj$rp[[1]]

  nu_spfk <- array(NA, dim = c(nS+1,nP+1,nF,3))
  P1_spf  <- array(NA, dim = c(nS+1,nP+1,nF))
  P2_spf  <- array(NA, dim = c(nS+1,nP+1,nF))

  # Fill!
  # nu pars
  nu_spfk[1:nS,1:nP,1:nF,] <- rp$EBSBpars$stockSpec$nu_spfk
  nu_spfk[nS+1,1:nP,1:nF,] <- rp$EBSBpars$coastWide$nu_spfk
  nu_spfk[1:nS,1+nP,1:nF,] <- rp$EBSBpars$dataPooled$nu_spfk
  nu_spfk[1+nS,1+nP,1:nF,] <- rp$EBSBpars$totalAgg$nu_spfk

  # P1
  P1_spf[1:nS,1:nP,1:nF] <- rp$EBSBpars$stockSpec$P1_spf
  P1_spf[nS+1,1:nP,1:nF] <- rp$EBSBpars$coastWide$P1_spf
  P1_spf[1:nS,1+nP,1:nF] <- rp$EBSBpars$dataPooled$P1_spf
  P1_spf[1+nS,1+nP,1:nF] <- rp$EBSBpars$totalAgg$P1_spf

  # P2
  P2_spf[1:nS,1:nP,1:nF] <- rp$EBSBpars$stockSpec$P2_spf
  P2_spf[nS+1,1:nP,1:nF] <- rp$EBSBpars$coastWide$P2_spf
  P2_spf[1:nS,1+nP,1:nF] <- rp$EBSBpars$dataPooled$P2_spf
  P2_spf[1+nS,1+nP,1:nF] <- rp$EBSBpars$totalAgg$P2_spf

  depSeq <- seq(0,1,length.out = 100 )

  par( mfcol = c(nS+1,nP+1),
        mar = c(.5,.5,.5,.5),
        oma = c(3,3,2,2) )


  fleetCols <- RColorBrewer::brewer.pal(nF,"Set1")

  for( s in 1:(nS+1) )
    for( p in 1:(nP+1) )
    {
      plot( x = c(0,1), y = c(0,2), axes = FALSE,
            type = "n" )
        mfg <- par( "mfg" )

        if( mfg[1] == mfg[3] )
          axis( side = 1 )
        if( mfg[2] == 1 )
          axis( side = 2, las = 1 )

        box()

        abline( h = 1, lwd = .8, lty = 2 )

        for( f in 1:nF )
        {
          numerator   <- (1 - exp( -nu_spfk[s,p,f,3]*( depSeq - P1_spf[s,p,f])))
          denominator <- (1 - exp( -nu_spfk[s,p,f,3]*( P2_spf[s,p,f] - P1_spf[s,p,f])))

          ratioSeq    <- nu_spfk[s,p,f,1] + (nu_spfk[s,p,f,2] - nu_spfk[s,p,f,1]) * numerator/denominator


          lines( x = depSeq, y = ratioSeq, lwd = 2, col = fleetCols[f] )
        }

    }

  legend( x = "topleft",
          legend = fleetNames,
          col = fleetCols,
          lwd = 2, bty = "n" )

} # END plotEBSBratio()

# plotRetroCatchability()
# Plot of retrospective catchability estimates for each
# fleet compared to the mean taken from the OM
plotRetroCatchability <- function(  obj = blob,
                                    iRep = 1 )
{
  # Get model dims
  tMP <- obj$om$tMP
  nT  <- obj$om$nT
  nS  <- obj$om$nS
  nP  <- obj$om$nP
  nF  <- obj$om$nF

  # Get model state arrays and AM fits
  SB_spt        <- obj$om$SB_ispt[iRep,,,1:t]
  VB_spt        <- obj$om$vB_ispft[iRep,,,2,1:t]
  totB_spt      <- obj$om$B_ispt[iRep,,,1:t]
  fitSB_spt     <- obj$mp$assess$retroSB_itspt[iRep,projt,,,1:t]
  fitVB_spft    <- obj$mp$assess$retroVB_itspft[iRep,projt,,,,1:t]
  fitq_spf      <- obj$mp$assess$retroq_itspf[iRep,,,,]
  fitq_spft     <- obj$mp$assess$retroq_itspft[iRep,,,,,]

  # Get ctlList
  ctlList <- obj$ctlList
  pT      <- ctlList$opMod$pT

  # Get mean catchability (OM fits)
  mq_spf      <- ctlList$opMod$histRpt$q_spf
  mq_sf       <- array(1, dim = c(nS,nF))
  mq_sf[,3:4] <- ctlList$opMod$histRpt$qSurv_sf
  sdlnq_f     <- ctlList$mp$assess$spsdlnq_f


}



# plotScaledIndices()
# Scaled indices with model fits for a given
# replicate and year.
plotScaledIndices <- function(  obj = blob,
                                iRep = 1,
                                Ct = TRUE,
                                t = blob$om$tMP )
{
  # Get model dims
  tMP <- obj$om$tMP
  nT  <- obj$om$nT
  nS  <- obj$om$nS
  nP  <- obj$om$nP
  nF  <- obj$om$nF

  # Calculate the projection time step
  projt <- t - tMP + 1

  # Blended index?
  blendIdx      <- obj$ctlList$opMod$blendIdx

  SB_spt        <- array(NA, dim = c(nS,nP,t))
  VB_spt        <- array(NA, dim = c(nS,nP,t))
  fitSB_spt     <- array(NA, dim = c(nS,nP,t))
  fitVB_spft    <- array(NA, dim = c(nS,nP,nF,t))
  fitq_spf      <- array(NA, dim = c(nS,nP,nF))
  fitq_spft     <- array(NA, dim = c(nS,nP,nF,t))
  fitqComb_spt  <- array(NA, dim = c(nS,nP,t))
  C_spft        <- array(NA, dim = c(nS,nP,nF,t))
  C_spt         <- array(NA, dim = c(nS,nP,t))

  # Get biomass arrays and catchabilities
  SB_spt[1:nS,,]        <- obj$om$SB_ispt[iRep,1:nS,,1:t]
  VB_spt[1:nS,,]        <- obj$om$vB_ispft[iRep,,,2,1:t]
  fitSB_spt[1:nS,,]     <- obj$mp$assess$retroSB_itspt[iRep,projt,,,1:t]
  fitVB_spft[1:nS,,,]   <- obj$mp$assess$retroVB_itspft[iRep,projt,,,,1:t]
  fitq_spf[1:nS,,]      <- obj$mp$assess$retroq_itspf[iRep,projt,,,]
  fitq_spft[1:nS,,,]    <- obj$mp$assess$retroq_itspft[iRep,projt,,,,1:t]

  if(is.null(obj$mp$assess$retroqComb_itspt[iRep,projt,,,1:t]))
  {
    cat('retroqComb_itspt not in mp$assess, uqing qComb=1 for projT...\n')
    histFile <- obj$ctlList$opMod$histFile
    load(here('history',histFile,paste0(histFile,'.Rdata')))
    rI_pgt <- reports$data$rI_pgt
    # q_pg   <- reports$repOpt$qComb_pg
    q_g   <- fitq_spf[1:nS,1,]

    for(tIdx in 1:(tMP-1))
      fitqComb_spt[1:nS,1,t] <- sum(rI_pgt[1,,tIdx] * q_g)

    fitqComb_spt[1:nS,1,tMP:t] <- 1

  } else
    fitqComb_spt[1:nS,,]  <- obj$mp$assess$retroqComb_itspt[iRep,projt,,,1:t]

  ctlList <- obj$ctlList
  spFleets <- ctlList$mp$assess$spFleets


  fitSB_spt[fitSB_spt < 0] <- NA

  C_spft[1:nS,,,]   <- obj$om$C_ispft[iRep,,,,1:t]
  C_spt[1:nS,,]     <- apply(X = C_spft, FUN = sum, MARGIN = c(1,2,4))

  # Now pull indices
  I_spft <- obj$mp$data$I_ispft[iRep,,,,1:t]
  I_spft[I_spft < 0] <- NA

  # Get ratios for blending index
  rI_spft <- obj$om$rI_ispft[iRep,,,,1:t]

  fleetPCH  <- 20 + 1:nF
  fleetBG   <- RColorBrewer::brewer.pal(nF, "Set1")
  stockCol  <- RColorBrewer::brewer.pal(nP, "Spectral")

  spTVqFleets <- ctlList$mp$assess$spTVqFleets

  SB_spt[SB_spt == 0]     <- NA
  VB_spt[VB_spt == 0]     <- NA

  nSS     <- dim( SB_spt)[1]
  nPP     <- dim( SB_spt)[2]



  scaledIdx_spft <- array( NA, dim = c(nSS, nP, nF, t ) )
  scaledIdx_spt  <- array( NA, dim = c(nSS, nP, t ) )
  for( s in 1:nSS )
    for( p in 1:nPP )
    {
      for( f in 1:nF )
      {
        if( f %in% spTVqFleets)
          scaledIdx_spft[s,p,f,] <- I_spft[s,p,f,] / fitq_spft[s,p,f,]
        else
          scaledIdx_spft[s,p,f,] <- I_spft[s,p,f,] / fitq_spf[s,p,f]

        # Scale by ratio of SB and VB
        scaledIdx_spft[s,p,f,] <- scaledIdx_spft[s,p,f,] * fitSB_spt[s,p,] / fitVB_spft[s,p,f,]
      }
      if( blendIdx ) # blendIdx not used since I_spft[s,p,5,] over
      {
        scaledIdx_spt[s,p,]   <- I_spft[s,p,5,] / fitqComb_spt[s,p,]
      }
    }


  speciesNames  <- obj$ctlList$opMod$species
  stockNames    <- obj$ctlList$opMod$stock
  fleetNames    <- obj$om$fleetNames
  fYear         <- obj$ctlList$opMod$fYear
  pT            <- obj$ctlList$opMod$pT

  yrs <- seq( from = fYear, by = 1, length.out = t)
  ppJitter <- seq(from = -.3, to = .3, length.out = nP )

  par(  mfcol = c(nPP,nSS),
        mar = c(1,1.5,1,1.5),
        oma = c(3,3,3,3) )
  for(s in 1:nSS)
    for( p in 1:nPP )
    {
      plot( x = range(yrs),
            y = c(0,max(VB_spt[s,p,],SB_spt[s,p,],na.rm = T) ),
            type = "n", axes = F )

      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1 )
      if( mfg[1] == 1 )
        mtext( side = 3, text = speciesNames[s], font = 2, line = 0 )
      axis( side = 2, las = 1 )
      if( mfg[2] == mfg[4] )
      {
        corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
        par(xpd = TRUE) #Draw outside plot area
        text(x = corners[2]+1.5, y = mean(corners[3:4]), stockNames[p], srt = 270,
              font = 2, cex = 1.5 )
        par(xpd = FALSE)
      }
      box()
      grid()

      if( Ct )
      {
        # Plot actual catch
        rect( xleft = yrs - .3, xright = yrs + .3,
              ytop = C_spt[s,p,], ybottom = 0, col = "grey40",
              border = NA )

      }

      lines( x = yrs, y = SB_spt[s,p,], col = "red", lwd = 3 )
      lines( x = yrs, y = fitSB_spt[s,p,], col = "black", lwd = 1 )

      if( !blendIdx )
        for( f in 1:nF )
          points( x = yrs, y = scaledIdx_spft[s,p,f,],
                  pch = fleetPCH[f], bg = fleetBG[f],
                  cex = 1.3, col = stockCol[p] )

      if( blendIdx )
      {
        points( x = yrs, y = scaledIdx_spt[s,p,],
                pch = 16, col = "grey40" )
        legend('topleft',bty='n',
                legend=paste0('AM qComb for last year = ',round(fitqComb_spt[s,p,t-1],2)))
      }

      abline( v = yrs[tMP] - 0.5, lty = 2, lwd = 0.5 )
    }
    # legend( x= "topright",
    #         legend = fleetNames,
    #         pch = fleetPCH,
    #         pt.bg = fleetBG )

} # END plotScaledIndices()

# plotAMIdxResids()
# Plots of standardised AM residuals
plotAMIdxResids <- function(  obj = blob,
                              iRep = 1,
                              Ct = TRUE,
                              t = blob$om$tMP )
{
  # Get model dims
  tMP <- obj$om$tMP
  nT  <- obj$om$nT
  nS  <- obj$om$nS
  nP  <- obj$om$nP
  nF  <- obj$om$nF

  blendIdx      <- obj$ctlList$opMod$blendIdx

  # Calculate the projection time step
  projt <- t - tMP + 1

  SB_spt        <- array(NA, dim = c(nS,nP,t))
  VB_spt        <- array(NA, dim = c(nS,nP,t))
  fitSB_spt     <- array(NA, dim = c(nS,nP,t))
  fitVB_spft    <- array(NA, dim = c(nS,nP,nF,t))
  fitq_spf      <- array(NA, dim = c(nS,nP,nF))
  fitq_spft     <- array(NA, dim = c(nS,nP,nF,t))
  I_spft        <- array(NA, dim = c(nS,nP,nF,t))
  tauObs_spf    <- array(NA, dim = c(nS,nP,nF))

  # Get biomass arrays
  SB_spt[1:nS,,1:t]     <- obj$om$SB_ispt[iRep,,,1:t]
  VB_spt[1:nS,,1:t]     <- obj$om$vB_ispft[iRep,,,2,1:t]
  fitSB_spt[1:nS,,1:t]  <- obj$mp$assess$retroSB_itspt[iRep,projt,,,1:t]
  fitVB_spft[1:nS,,,1:t]<- obj$mp$assess$retroVB_itspft[iRep,projt,,,,1:t]
  fitq_spf[1:nS,,]      <- obj$mp$assess$retroq_itspf[iRep,projt,,,]
  fitq_spft[1:nS,,,1:t] <- obj$mp$assess$retroq_itspft[iRep,projt,,,,1:t]

  tauObs_spf[1:nS,,]    <- obj$mp$assess$retrotauObs_itspf[iRep,projt,,,]

  if( all( is.na(tauObs_spf)) )
    tauObs_spf[1:nS,,] <- 1

  tauObs_spf[is.na(tauObs_spf)] <- 0

  ctlList <- obj$ctlList

  fitSB_spt[fitSB_spt < 0] <- NA

  # idxFleets


  spTVqFleets <- ctlList$mp$assess$spTVqFleets
  spFleets    <- ctlList$mp$assess$spFleets

  # Model dims
  tMP     <- obj$om$tMP
  nS      <- obj$om$nS
  nP      <- obj$om$nP
  nT      <- obj$om$nT

  # Now pull indices
  I_spft[1:nS,1:nP,,1:t] <- obj$mp$data$I_ispft[iRep,1:nS,1:nP,,1:t]
  I_spft[I_spft <= 0] <- NA

  nSS <- nS
  nPP <- nP

  stdResids_spft  <- array( NA, dim = c(nSS, nPP, nF, t) )
  stdResids_spt   <- array( NA, dim = c(nSS, nPP, t) )

  for( s in 1:nSS )
    for( p in 1:nPP )
    {
      for( f in spFleets )
      {
        if( tauObs_spf[s,p,f] > 0 )
        {

          if( f %in% spTVqFleets )
            stdResids_spft[s,p,f,] <- (-log( I_spft[s,p,f,]/fitq_spft[s,p,f,] ) + log(fitVB_spft[s,p,f,]))/tauObs_spf[s,p,f]
          if( !f %in% spTVqFleets )
            stdResids_spft[s,p,f,] <- (-log( I_spft[s,p,f,]/fitq_spf[s,p,f] ) + log(fitVB_spft[s,p,f,]))/tauObs_spf[s,p,f]


        }
      }
    }


  fleetPCH <- 20 + 1:nF
  fleetBG <- RColorBrewer::brewer.pal(nF, "Set1")



  speciesNames  <- obj$ctlList$opMod$species
  stockNames    <- obj$ctlList$opMod$stock
  fleetNames    <- obj$om$fleetNames
  fYear         <- obj$ctlList$opMod$fYear
  pT            <- obj$ctlList$opMod$pT

  # HACK - need to import fleetnames
  fleetNames <- c("reduction","seineRoe","gillnet","surf","dive","SOK")

  if( nSS == 1 )
    speciesNames <- "Data Pooled"

  yrs <- seq( from = fYear, by = 1, length.out = t)

  ppJitter <- seq( from = -.3, to = .3, length.out = nP )

  par(  mfcol = c(nPP,nSS),
        mar = c(1,1.5,1,1.5),
        oma = c(3,3,3,3) )
  for(s in 1:nSS)
    for( p in 1:nPP )
    {
      if( any(!is.na(stdResids_spft[s,p,,])))
        maxResid <- max(abs(stdResids_spft[s,p,,]),na.rm = T)
      else maxResid <- 1

      plot( x = range(yrs),
            y = range(-maxResid,maxResid),
            type = "n", axes = F )

      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1 )
      if( mfg[1] == 1 )
        mtext( side = 3, text = speciesNames[s], font = 2, line = 0 )
      axis( side = 2, las = 1 )
      if( mfg[2] == mfg[4] )
      {
        corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
        par(xpd = TRUE) #Draw outside plot area
        text(x = corners[2]+1.5, y = mean(corners[3:4]), stockNames[p], srt = 270,
              font = 2, cex = 1.5 )
        par(xpd = FALSE)
      }
      box()
      grid()

      for( f in 1:nF )
      {
        points( x = yrs, y = stdResids_spft[s,p,f,],
                pch = fleetPCH[f], bg = fleetBG[f],
                cex = 1.3 )
        nonNA <- which(!is.na(stdResids_spft[s,p,f,]))
        if( length(nonNA) > 0 )
        {
          yVal <- stdResids_spft[s,p,f,nonNA]
          xVal <- yrs[nonNA]



          dat <- data.frame(x = xVal, y = yVal )
          regLine <- lm( y~x, data = dat )

          pVal <- round(summary(regLine)$coefficients[2,4],3)

          pLabel <- paste("p = ", pVal, sep = "")

          dat <- dat %>%
                  mutate( regLine = predict.lm(regLine, newdata = dat) )


          lines( x = dat$x, y = dat$regLine, col = fleetBG[f], lwd = 2 )
          text( x = dat$x[1], y = 1.2*dat$y[1], label = pLabel, col = fleetBG[f], font = 2 )

        }
      }
      abline( h = 0, lty = 2)


    }

    legend( x= "bottomleft", bty = "n",
            legend = c(fleetNames),
            pch = fleetPCH,
            pt.bg = fleetBG,
            cex = 1.3 )



} # END plotAMIdxResids()


.plotDiagCondition <- function( obj = blob, iRep = 1)
{
  repObj <- obj$ctlList$opMod$histRpt
  posts <- obj$ctlList$opMod$posts

  iPost <- ifelse( obj$ctlList$opMod$posteriorSamples, obj$ctlList$opMod$postDraws_i[iRep], 0 )

  diagCondition( repObj, obj, posts, iRep, iPost )

}

# diagCondition()
# Plots that help diagnose issues with conditioning
# between the fitted OM report, and the conditioned
# ms3R operating model.
diagCondition <- function(  repObj  = blob$ctlList$opMod$histRpt,
                            ms3Obj  = blob,
                            iRep    = 1 )
{

  nP <- repObj$nP
  nT <- repObj$nT
  nF <- repObj$nG
  stockNames    <- dimnames(ms3Obj$ctlList$opMod$histRpt$I_pgt)[[1]]

  par(mfrow = c(3,2), mar = c(1,2,1.5,1), oma = c(3,5,3,3) )

  # Biomass RE
  plotRE_spt( repObj = repObj, omObj = ms3Obj$om,
              AMseries = "SB_pt",
              OMseries = "SB_ispt",
              iRep = iRep )
  mtext( side = 3, text = "SB_spt", line = 0, font = 2)

  # Recruitment
  plotRE_spt( repObj = repObj,
              omObj = ms3Obj$om,
              AMseries = "R_pt",
              OMseries = "R_ispt",
              iRep = iRep )
  mtext( side = 3, text = "R_spt", line = 0, font = 2)

  # Recruitment errors
  plotRE_spt( repObj = repObj,
              omObj = ms3Obj$om$errors,
              nS = 1,
              AMseries = "SRdevs_pt",
              OMseries = "omegaR_ispt",
              iRep = iRep )
  mtext( side = 3, text = expression(omega[R]), line = 0, font = 2)

  # Catch
  plotRE_spft(  repObj = repObj,
                omObj = ms3Obj$om,
                AMseries = "C_pgt",
                OMseries = "C_ispft",
                iRep = iRep )
  mtext( side = 3, text = expression(C[spft]), line = 0, font = 2)

  # F
  plotRE_spft(  repObj = repObj,
                omObj = ms3Obj$om,
                AMseries = "U_pgt",
                OMseries = "F_ispft",
                iRep = iRep,
                legendON=TRUE, stocks=stockNames )
  mtext( side = 3, text = expression(F[spft]), line = 0, font = 2)

  mtext( side = 1, outer = TRUE, text = "Time Step", font = 2, line = 2)
  mtext( side = 2, outer = TRUE, text = "Relative error in conditioning",
          font = 2, line = 3 )



  # # Numbers at age
  # plotRE_axspt( repObj = repObj, omObj = ms3Obj$om, series = "N_iaxspt" )
  # mtext( side = 3, text = expression(N[axspt]), line = 0, font = 2)

  # # Biomass at age
  # plotRE_axspt( repObj = repObj, omObj = ms3Obj$om, series = "N_axspt" )
  # mtext( side = 3, text = expression(B[axspt]), line = 0, font = 2)
}
calcREdist <- function( est, true, marg,
                        qProbs = c(0.025, 0.5, 0.975) )
{

  # Calculate REs
  re <- (est - true)/true

  # calculate distribution over margin
  reQuants <-  apply( X = re, FUN = quantile,
                      MARGIN = marg, probs = qProbs,
                      na.rm = TRUE )

  reQuants[!is.finite(reQuants)] <- 0

  return(reQuants)

}

plotRE_spt <- function( repObj, omObj, posts=NULL, nS = 1,
                        AMseries = "SB_spt",
                        OMseries = "SB_ispt", iRep = 1,
                        iPost = 0,
                        yRange = NULL )
{

  nP <- repObj$nP
  nT <- repObj$nT

  true_spt  <- array(NA, dim = c(nS,nP,nT))
  est_spt   <- array(NA, dim = c(nS,nP,nT))

  if( iPost )
    true_spt[1:nS,,] <- posts[[AMseries]][iPost, ,1:nT,drop = FALSE]
  else
    true_spt[1:nS,,] <- repObj[[AMseries]][,1:nT,drop = FALSE]

  est_spt[1:nS,1:nP,1:nT] <- omObj[[OMseries]][iRep,1:nS,1:nP,1:nT]

  if( OMseries == "omegaR_ispt" )
  {

    if(is.null(repObj$avgRcode_p))
      repObj$avgRcode_p <- rep(0,nP)

    for (p in 1:nP)
    {

#      if(repObj$avgRcode_p[p]==1)
#        true_spt[1:nS,p,]   <- repObj$SRdevs_pt[p,1:nT]
#      else
#        true_spt[1:nS,p,]   <- repObj$omegaR_pt[p,1:nT]

    }
    sigmaR <- repObj$sigmaR
    est_spt[1:nS,1:nP,1:nT] <- est_spt[1:nS,1:nP,1:nT] - 0.5 * sigmaR
  }

  re_qspt  <- calcREdist( true = true_spt,
                          est  = est_spt,
                          marg = c(1,2,3) )

  stockCols <- RColorBrewer::brewer.pal(n = nP, "Dark2")
  stockLty <- rep(1,nP)

  if( is.null(yRange) )
  {
    yRange <- range(re_qspt, na.rm = TRUE)
    yRange[2] <- max(.01,yRange[2])
    yRange[1] <- min(-.01,yRange[1])
  }

  plot( x = c(1,nT), y = yRange,
        axes = FALSE, type = "n" )
    mfg <- par("mfg")
    if(mfg[1] == mfg[3])
      axis( side = 1 )
    axis( side = 2, las = 1 )
    box()
    grid()
    for( s in 1:nS )
      for( p in 1:nP )
      {
        # polyY <- c(re_qspt[1,s,p,],rev(re_qspt[3,s,p,]))
        # polygon( x = c(1:nT,nT:1), y = polyY,
        #           col = scales::alpha(stockCols[p],alpha = .3),
        #           lty = stockLty[p] )
        lines( x = 1:nT, y = re_qspt[2,s,p,],
                col = stockCols[p], lty = stockLty[p] )

      }
    abline( h = 0, lty = 3, lwd = .8 )
}

plotRE_spft <- function(  repObj, omObj, posts=NULL,
                          OMseries = "C_spft",
                          AMseries = "C_pgt",
                          iRep = 1, iPost=0, nS = 1,
                          legendON=FALSE, stocks=NULL )
{
  nP <- repObj$nP
  nT <- repObj$nT
  nF <- repObj$nG

  # if(OMseries == "F_ispft")
  #   browser()

  fleetType_f <- omObj$fleetType_f
  sokFleets <- which(fleetType_f >= 2)

  true_spft  <- array(NA, dim = c(nS,nP,nF,nT))
  est_spft   <- array(NA, dim = c(nS,nP,nF,nT))

  if( iPost )
    true_spft[1:nS,,,] <- posts[[AMseries]][iPost,,,1:nT,drop = FALSE]
  else
    true_spft[1:nS,,,] <- repObj[[AMseries]][,,1:nT,drop = FALSE]


  est_spft[1:nS,1:nP,1:nF,1:nT] <- omObj[[OMseries]][iRep,1:nS,1:nP,1:nF,1:nT]


  # repObj C_pgt contains SOK product for SOK fleets
  # omObj C_spft contains ponded fish for SOK fleets

  # # convert SOK product which is what is reported in AMseries
  # if( OMseries == "C_ispft" & length(sokFleets) > 0 )
  # {
  #   psi_pgt <- repObj$psi_pgt
  #   psi_gt  <- repObj$psi_gt

  #   for (p in 1:nP)
  #     for( f in sokFleets)
  #     {
  #       if(any(psi_pgt[p,f,] != 0))
  #         true_spft[1:nS,p,f,] <- true_spft[1:nS,p,f,]/psi_pgt[p,f,]

  #       if(any(psi_gt[f,] != 0))
  #         true_spft[1:nS,p,f,] <- true_spft[1:nS,p,f,]/psi_gt[f,]
  #     }
  # }

  re_qspt  <- calcREdist( true = true_spft,
                            est  = est_spft,
                            marg = c(1,2,4) )

  re_qft  <- calcREdist(  true = true_spft,
                          est  = est_spft,
                          marg = c(3:4) )



  stockCols <- RColorBrewer::brewer.pal(n = nP, "Dark2")
  stockLty <- rep(1,nP)

  yRange <- range(re_qspt, na.rm = TRUE)
  yRange[2] <- max(.01,yRange[2])
  yRange[1] <- min(-.01,yRange[1])

  col <- RColorBrewer::brewer.pal(n = nF, "Paired")

  plot( x = c(1,nT), y = yRange,
        axes = FALSE, type = "n" )
    mfg <- par("mfg")
    if(mfg[1] == mfg[3])
      axis( side = 1 )
    axis( side = 2, las = 1 )
    box()
    grid()

    for( s in 1:nS )
      for( p in 1:nP )
      {
        polyY <- c(re_qspt[1,s,p,],rev(re_qspt[3,s,p,]))
        polygon( x = c(1:nT,nT:1), y = polyY,
                  col = scales::alpha(stockCols[p],alpha = .3),
                  lty = stockLty[p] )
        lines( x = 1:nT, y = re_qspt[2,s,p,],
                col = stockCols[p], lty = stockLty[p] )

      }
    abline( h = 0, lty = 3, lwd = .8 )

    for(f in 1:nF )
      lines(x = 1:nT, y=re_qft[2,f,], col = col[f])

  if(legendON)
    legend('topright', bty='n',
           legend=1:nF, lty=1, col=col)

}

plotRE_at <- function(  repObj, omObj,
                          OMseries = "M_iaxspt",
                          AMseries = "M_apt",
                          posts = NULL,
                          iRep = 1, iPost=0, nS = 1,
                          legendON=FALSE, stocks=NULL )
{
  nS <- repObj$nS
  nP <- repObj$nP
  nT <- repObj$nT

  if( iPost )
    true_at <- posts[[AMseries]][iPost,,1,1:nT]
  else
    true_at <- repObj[[AMseries]][,1,1:nT]

  est_at  <- omObj[[OMseries]][iRep,,1,1,1,1:nT]

  re_qt  <- calcREdist( true = true_at[,1:nT],
                        est  = est_at[,1:nT],
                        marg = c(2) )

  stockLty <- rep(1,nP)

  stockCols <- RColorBrewer::brewer.pal(n = nP, "Dark2")

  yRange <- range(re_qt, na.rm = TRUE)
  yRange[2] <- max(.01,yRange[2])
  yRange[1] <- min(-.01,yRange[1])

  plot( x = c(1,nT), y = yRange,
        axes = FALSE, type = "n" )
    mfg <- par("mfg")
    if(mfg[1] == mfg[3])
      axis( side = 1 )
    axis( side = 2, las = 1 )
    box()
    grid()
    polyY <- c(re_qt[1,],rev(re_qt[3,]))
    polygon( x = c(1:nT,nT:1), y = polyY,
              col = scales::alpha(stockCols[1],alpha = .3),
              lty = stockLty[1] )
    lines( x = 1:nT, y = re_qt[2,],
            col = stockCols[1], lty = stockLty[1] )


    abline( h = 0, lty = 3, lwd = .8 )
}

# So, x-axis is m1, colour is Mb, and pch is sigmaM
assignPCH <- function( Mb )
{
  pchList <- c(21,23,24)
  names(pchList) <- c("0.2","0.4","0.6")

  pch <- pchList[as.character(Mb)]
  pch
}

assignCol <- function( sigmaM, nameVec = c("0","0.1","0.3") )
{
  cols <- brewer.pal("Dark2", n = length(nameVec))[1:length(nameVec)]
  names(cols) <- nameVec

  col <- cols[as.character(sigmaM)]

  col
}

addJitter <- function( Mb, jitt = c(-.05,0,.05), nameVec = c("0.2","0.4","0.6") )
{
  names(jitt) <- nameVec

  xJitt <- jitt[as.character(Mb)]

  xJitt
}


plotTable <- function(  medREsTable = HG_MREs, ylab = "MRE",
                        columns = c("m1","Mb","h","B0","totB0"))
{

  medREsTable <- medREsTable |>
                  mutate( m1sim = as.numeric(m1sim),
                          Mbsim = as.numeric(Mbsim),
                          sigmaMsim = as.numeric(sigmaMsim),
                          pch = sapply( X = Mbsim, FUN = assignPCH),
                          col = sapply( X = sigmaMsim, FUN = assignCol),
                          MbJitter = sapply(X = Mbsim, FUN = addJitter),
                          sigmaMJitter = sapply(  X = sigmaMsim, FUN = addJitter,
                                                  nameVec = c("0","0.1","0.3"), jitt = c(-.2,0.,.2)),
                          xJitter = m1sim + MbJitter + sigmaMJitter  )

  nReps <- as.integer(medREsTable$nReps)
  medREsTable$bg <- medREsTable$col
  medREsTable$bg[nReps < 95] <- "white"
  medREsTable$pch[nReps < 75] <- NA


  nPlots <- length(columns)
  par(mfrow = c(nPlots,1),
      mar = c(.1,2,.1,2),
      oma = c(4,4,2,2) )

  m1Range <- c(0.5,7)

  m1s <- unique(medREsTable$m1sim)

  for( cIdx in 1:nPlots )
  {

    yRange <- abs(range(as.numeric(medREsTable[nReps >= 75,columns[cIdx]])))

    yRange <- c(-max(yRange),max(yRange))
    plot( x = m1Range, y = yRange,
          axes = FALSE, type = "n" )
      mfg <- par("mfg")
      axis( side = 2, las = 1)
      if( mfg[1] == mfg[3] )
        axis( side = 1, labels = m1s, at = m1s )
      grid()
      box()

      mtext(side = 2, text = paste0(ylab," (", columns[cIdx] ,")"),
            line = 3)

      # Plot points now, see about a jitter later
      points( x = medREsTable$xJitter, y = medREsTable[,columns[cIdx]],
              col = medREsTable$col, pch = medREsTable$pch,
              bg = medREsTable$bg,
              cex = 1.5 )
      abline(h = 0, lty = 2, col = "grey60")

      if(mfg[1] == 1 )
        legend( x = "topright", bty = "n",
                legend = c( "sigmaM = 0",
                            "sigmaM = 0.1",
                            "sigmaM = 0.3",
                            "Mb = 0.2",
                            "Mb = 0.4",
                            "Mb = 0.6"),
                col = c(brewer.pal("Dark2", n = 3),
                        rep("grey40",3)),
                bg = c(brewer.pal("Dark2", n = 3),
                        rep("white",3)),
                pch = c(15,15,15,21,23,24),
                cex = 1.2 )

  }
}

# plotSimEstDepM()
# Plots the true OM depM model,
# and the tulip of estimated depM
# models from a sim-est experiment
plotSimEstDepM <- function( obj = blob, nB = 100,
                            nTrace = 3,
                            sIdx = 1, pIdx = 1, pt = 1)
{
  # First, pull out OM values
  totB0     <- obj$rp[[1]]$totB0_sp[sIdx,pIdx]
  Mb        <- obj$rp[[1]]$Mb_sp[sIdx,pIdx]
  m1        <- obj$rp[[1]]$m1_sp[sIdx,pIdx]

  # then pull retros
  retroMb_i <- obj$mp$assess$retroM_itsp[,pt,sIdx,pIdx]
  retrom1_i <- obj$mp$assess$retrom1_itsp[,pt,sIdx,pIdx]

  nReps <- sum(obj$goodReps)
  goodReps <- which(obj$goodReps)

  Dseq <- seq(from = 0, to = 3, length.out = nB)

  depMcurve_id <- array(NA, dim = c(nReps,nB))
  for( i in 1:length(goodReps) )
  {
    idx <- goodReps[i]

    depMcurve_id[i,] <- retroMb_i[idx] + exp( -retrom1_i[idx] * Dseq )
  }

  traceIdx <- sample(1:nReps,nTrace)

  trueDepM <- Mb + exp(-m1 * Dseq )

  depMcurve_qd <- apply(  X = depMcurve_id, FUN  = quantile,
                          MARGIN = 2, probs = c(0.025, 0.5, 0.975) )


  maxM <- max(c(Mb + 1, retroMb_i + 1), na.rm = T)

  plot( x = c(0,3), y = c(0,maxM), type = "n",
        las = 1 )
    polygon(  x = c(Dseq,rev(Dseq)),
              y = c(depMcurve_qd[1,],rev(depMcurve_qd[3,])),
              border = NA, col = "grey70" )
    lines( x = Dseq, y = depMcurve_qd[2,], lwd = 3, col = "black", lty = 2 )
    lines( x = Dseq, y = trueDepM, col = "red", lwd = 3 )
    for( tIdx in traceIdx )
      lines( x = Dseq, y = depMcurve_id[tIdx,], lwd = .8)

}

plotMPSummary <- function( obj        = blob,
                           commGears  = 1:3,
                           sokGear    = 6,
                           aggregate  = FALSE,
                           showBref   = TRUE,
                           nTrace     = 3,
                           traces     = NULL,
                           tMin       = NULL,
                           areaColors = NULL,
                           gearColors = NULL )
{
  goodReps   <- obj$goodReps
  tMP        <- obj$om$tMP
  nP         <- obj$om$nP
  nT         <- obj$om$nT
  fYear      <- obj$ctlList$opMod$fYear
  yrs        <- seq(from = fYear, by = 1, length.out = nT)
  stockNames  <- obj$om$stockNames
  fleetNames  <- obj$ctlList$opMod$fleets

  if(aggregate | nP == 1L)
    stockNames <- "WCVI Herring"

  # Aggregating a single spatial stock is a no-op; treat as per-area
  if( nP == 1L )
    aggregate <- FALSE

  if( is.null(tMin) )
    tMin <- 1

  tPlotIdx <- tMin:nT

  # Herring species index fixed at s = 1
  sIdx      <- 1L
  nGoodReps <- sum(goodReps)

  # Always use drop = FALSE when subsetting arrays so that singleton
  # dimensions (e.g. nP = 1) are never silently dropped.  Strip the
  # fixed s dimension (always size 1 after sIdx selection) explicitly
  # via dim<-, leaving [nReps, nP, nT].
  SB_ipt      <- obj$om$endSB_ispt[goodReps, sIdx, , , drop = FALSE]
  dim(SB_ipt) <- c(nGoodReps, nP, nT)
  SB_ipt[SB_ipt == 0] <- NA

  # Catch: [nReps, nP, nF, nT]
  nF          <- dim(obj$om$C_ispft)[4L]
  C_ipft      <- obj$om$C_ispft[goodReps, sIdx, , , , drop = FALSE]
  dim(C_ipft) <- c(nGoodReps, nP, nF, nT)

  nReps <- nGoodReps

  # Reference points — one value per spatial patch (p)
  B0_p <- obj$ctlList$opMod$histRpt$refPts$refCurves$SBeq_pf[, 1]
  Bref <- obj$ctlList$mp$hcr$Bref_s[1]

  if( is.null(traces) )
    traces <- sample(seq_len(nReps), size = min(nTrace, nReps))

  # Helper: get SB as [nReps, nT] for patch col or coastwide aggregate
  getSB <- function( col )
  {
    if( aggregate )
      apply(SB_ipt, c(1, 3), sum, na.rm = TRUE)
    else
      SB_ipt[, col, ]
  }

  # Helper: catch summed over gears → [nReps, nT] respecting aggregate flag
  getC <- function( col, gears )
  {
    if( aggregate )
    {
      dat <- C_ipft[, , gears, , drop = FALSE]
      apply(X = dat, MARGIN = c(1, 4), FUN = sum, na.rm = TRUE)
    } else {
      dat <- C_ipft[, col, gears, , drop = FALSE]
      apply(X = dat, MARGIN = c(1, 4), FUN = sum, na.rm = TRUE)
    }
  }

  # Helper: catch for a specific area (ignoring aggregate flag) → [nReps, nT]
  getC_area <- function( area, gears )
  {
    dat <- C_ipft[, area, gears, , drop = FALSE]
    apply(X = dat, MARGIN = c(1, 4), FUN = sum, na.rm = TRUE)
  }

  # Helper: draw stacked bars from a [nGroups, nT] matrix of medians.
  # Uses yrs and tPlotIdx from the enclosing scope.
  plotStackedBars <- function( meds_gt, gColors )
  {
    nGroups <- nrow(meds_gt)
    bottom  <- rep(0, length(tPlotIdx))
    for( g in seq_len(nGroups) )
    {
      tops <- bottom + pmax(meds_gt[g, tPlotIdx], 0, na.rm = TRUE)
      rect( xleft   = yrs[tPlotIdx] - 0.3,
            xright  = yrs[tPlotIdx] + 0.3,
            ybottom = bottom,
            ytop    = tops,
            col     = gColors[g], border = NA )
      bottom <- tops
    }
  }

  # Safe per-area labels: stockNames may be NULL or shorter than nP.
  # Default to PFMA names (25 = index 1, 24 = index 2, 23 = index 3).
  defAreaLabels <- paste0("PFMA ", c(25, 24, 23))[1:nP]
  areaLabels    <- if( !is.null(stockNames) && length(stockNames) >= nP )
                     stockNames[1:nP]
                   else
                     defAreaLabels

  nCol     <- if( aggregate ) 1L else nP
  colNames <- if( aggregate ) "WCVI Aggregate" else areaLabels

  # Default colours — distinct palettes for areas vs fleets so the two
  # groupings are always visually separable across all plot rows.
  defAreaColors <- c("steelblue3", "darkorange2", "forestgreen",
                     "mediumpurple3", "goldenrod2")
  defGearColors <- c("grey35", "coral3", "darkseagreen4",
                     "sienna3", "orchid4")
  if( is.null(areaColors) )
    areaColors <- defAreaColors[1:nP]
  if( is.null(gearColors) )
    gearColors <- defGearColors[seq_along(commGears)]

  stamp <- paste(obj$ctlList$ctl$scenarioName, ":",
                 obj$ctlList$ctl$mpName, sep = "")

  par( mfrow = c(3L, nCol),
       oma  = c(4, 5, 3, 1),
       mar  = c(0.3, 0.3, 0.8, 0.3),
       mgp  = c(2.5, 0.6, 0) )

  commLabel  <- paste0("Gears ", paste(commGears, collapse = ", "))
  showAggRef <- aggregate || nP == 1L

  # Pre-compute shared y-axis limits across all columns for each row
  # so that panels within a row are on the same scale.
  yMaxSB   <- 0
  yMaxComm <- 0
  yMaxSOK  <- 0
  for( col in 1:nCol )
  {
    B0_col  <- if( aggregate ) sum(B0_p) else B0_p[col]
    LRP_col <- 0.3 * B0_col
    refVals <- if( showAggRef ) c(B0_col, LRP_col, Bref)
               else             c(B0_col, LRP_col)

    SB_it    <- getSB(col)
    SB_qt    <- apply(X     = SB_it, MARGIN = 2, FUN = quantile,
                      probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
    yMaxSB   <- max(yMaxSB, SB_qt[3, tPlotIdx], refVals, na.rm = TRUE)

    Cc_it    <- getC(col = col, gears = commGears)
    Cc_qt    <- apply(X     = Cc_it, MARGIN = 2, FUN = quantile,
                      probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
    yMaxComm <- max(yMaxComm, Cc_qt[3, tPlotIdx], na.rm = TRUE)

    # SOK catch is small; use 95th percentile so open-fishery years
    # are visible even when most replicates have zero catch.
    Cs_it   <- getC(col = col, gears = sokGear)
    Cs_q95  <- apply(X     = Cs_it,
                     MARGIN = 2,
                     FUN   = quantile, probs = 0.95, na.rm = TRUE)
    yMaxSOK <- max(yMaxSOK, Cs_q95[tPlotIdx], na.rm = TRUE)
  }
  yMaxSB   <- ifelse(is.finite(yMaxSB), yMaxSB * 1.05, 1)
  yMaxComm <- ifelse(is.finite(yMaxComm), yMaxComm * 1.05, 1)
  yMaxSOK  <- ifelse(is.finite(yMaxSOK),yMaxSOK * 1.05, 0.1 )

  # ----------------------------------------------------------------
  # Row 1: Spawning Biomass
  # ----------------------------------------------------------------
  for( col in 1:nCol )
  {
    SB_it  <- getSB(col)
    SB_qt  <- apply(X     = SB_it, MARGIN = 2, FUN = quantile,
                    probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
    B0_ref <- if( aggregate ) sum(B0_p) else B0_p[col]
    LRP    <- 0.3 * B0_ref

    plot( x    = range(yrs[tPlotIdx]),
          y    = c(0, yMaxSB),
          type = "n", las = 1, xaxt = "n", yaxt = "n",
          xlab = "", ylab = "" )

    mfg <- par("mfg")
    if( mfg[2] == 1 )
    {
      axis( side = 2, las = 1 )
      mtext( side = 2, text = "SSB (kt)",
              line = 3, las = 0 )

    }

    mtext( side = 3, text = colNames[col], font = 2, line = 0.2 )
    grid()
    box()

    polygon( x      = c(yrs[tPlotIdx], rev(yrs[tPlotIdx])),
             y      = c(SB_qt[1, tPlotIdx], rev(SB_qt[3, tPlotIdx])),
             col    = "grey65", border = NA )
    lines( x = yrs[tPlotIdx], y = SB_qt[2, tPlotIdx], lwd = 2.5 )

    for( tIdx in traces )
      lines( x = yrs[tPlotIdx], y = SB_it[tIdx, tPlotIdx],
             lwd = 0.6 )

    abline( v = yrs[tMP], col = "grey30",  lty = 3 )
    abline( h = B0_ref,   col = "grey50",  lty = 3, lwd = 2 )

    if( showBref )
    {
      abline( h = LRP,  col = "red",       lty = 2, lwd = 2 )
      if( showAggRef )
        abline( h = Bref, col = "darkgreen", lty = 2, lwd = 2 )
    }

    # When aggregating, overlay per-area median SB lines in area colours
    if( aggregate )
    {
      for( p in 1:nP )
      {
        sb_med <- apply(X      = SB_ipt[, p, ],
                        MARGIN = 2,
                        FUN    = median, na.rm = TRUE)
        lines( x   = yrs[tPlotIdx],
               y   = sb_med[tPlotIdx],
               col = areaColors[p], lwd = 1.5, lty = 1 )
      }
    }

    # Legend in the rightmost column: SB trajectories + reference lines
    if( col == nCol )
    {
      leg_lab <- character(0)
      leg_col <- character(0)
      leg_lty <- numeric(0)
      leg_lwd <- numeric(0)

      if( aggregate )
      {
        leg_lab <- c(leg_lab, "WCVI Aggregate", areaLabels)
        leg_col <- c(leg_col, "black",     areaColors)
        leg_lty <- c(leg_lty, 1,           rep(1, nP))
        leg_lwd <- c(leg_lwd, 2.5,         rep(1.5, nP))
      }

      leg_lab <- c(leg_lab, "B0")
      leg_col <- c(leg_col, "grey50")
      leg_lty <- c(leg_lty, 3)
      leg_lwd <- c(leg_lwd, 2)

      if( showBref )
      {
        leg_lab <- c(leg_lab, "LRP (0.3B0)")
        leg_col <- c(leg_col, "red")
        leg_lty <- c(leg_lty, 2)
        leg_lwd <- c(leg_lwd, 2)

        if( showAggRef )
        {
          leg_lab <- c(leg_lab, "Bref")
          leg_col <- c(leg_col, "darkgreen")
          leg_lty <- c(leg_lty, 2)
          leg_lwd <- c(leg_lwd, 2)
        }
      }

      legend( "topright",
              legend = leg_lab, col = leg_col,
              lty    = leg_lty, lwd = leg_lwd,
              bty    = "n", cex = 1.0 )
    }
  }

  # ----------------------------------------------------------------
  # Row 2: Commercial Catch — stacked by gear (!aggregate) or area (aggregate)
  # ----------------------------------------------------------------
  for( col in 1:nCol )
  {
    # Total catch for 95% whiskers on projection years
    C_tot_it <- getC(col = col, gears = commGears)
    C_tot_qt <- apply(X     = C_tot_it, MARGIN = 2, FUN = quantile,
                      probs = c(0.025, 0.5, 0.975), na.rm = TRUE)

    plot( x    = range(yrs[tPlotIdx]),
          y    = c(0, yMaxComm), 
          type = "n", las = 1, xaxt = "n", yaxt = "n",
          xlab = "", ylab = "" )

    mfg <- par("mfg")
    if( mfg[2] == 1 )
    {
      axis( side = 2, las = 1 )
      mtext(  side = 2, 
              text = paste0("Commercial catch (kt)"),
              line = 3, las = 0 )

    }

    grid()
    box()

    if( aggregate )
    {
      # Stacked by area: one colour per spatial patch
      meds_gt <- do.call(rbind, lapply(
        X   = 1:nP,
        FUN = function(p)
          apply(X      = getC_area(area = p, gears = commGears),
                MARGIN = 2, FUN = median, na.rm = TRUE) ))
      plotStackedBars(meds_gt = meds_gt, gColors = areaColors)
    } else {
      # Stacked by gear: one colour per commercial gear
      meds_gt <- do.call(rbind, lapply(
        X   = commGears,
        FUN = function(g)
          apply(X      = getC_area(area = col, gears = g),
                MARGIN = 2, FUN = median, na.rm = TRUE) ))
      plotStackedBars(meds_gt = meds_gt, gColors = gearColors)
    }

    segments( x0  = yrs[tMP:nT], x1 = yrs[tMP:nT],
              y0  = C_tot_qt[1, tMP:nT],
              y1  = C_tot_qt[3, tMP:nT],
              col = "black", lwd = 0.8 )

    abline( v = yrs[tMP], col = "grey30", lty = 3 )

    # Legend in the rightmost column panel
    if( col == nCol )
    {
      if( aggregate )
        legend( "topright",
                legend = areaLabels,
                fill   = areaColors,
                bty    = "n", cex = 1.0 )
      else
        legend( "topright",
                legend = fleetNames[commGears],
                fill   = gearColors,
                bty    = "n", cex = 1.0 )
    }
  }

  
  # ----------------------------------------------------------------
  # Row 3: SOK / FSC — stacked by area when aggregate, single colour otherwise
  # ----------------------------------------------------------------
  for( col in 1:nCol )
  {
    # Total SOK for 95% whiskers
    C_tot_it <- getC(col = col, gears = sokGear)
    C_tot_qt <- apply(X     = C_tot_it, MARGIN = 2, FUN = quantile,
                      probs = c(0.025, 0.5, 0.975), na.rm = TRUE)

    plot( x    = range(yrs[tPlotIdx]),
          y    = c(0, yMaxSOK),
          type = "n", las = 1, yaxt = "n",
          xlab = "", ylab = "" )

    mfg <- par("mfg")
    if( mfg[2] == 1 )
    {
      axis( side = 2, las = 1 )
      mtext( side = 2, text = "SOK/FSC (kt)",
         line = 3, las = 0 )

    }

    grid()
    box()

    if( aggregate )
    {
      # Stacked by area in area colours
      meds_gt <- do.call(rbind, lapply(
        X   = 1:nP,
        FUN = function(p)
          apply(X      = getC_area(area = p, gears = sokGear),
                MARGIN = 2, FUN = median, na.rm = TRUE) ))
      plotStackedBars(meds_gt = meds_gt, gColors = areaColors)
    } else {
      # Single area, single colour
      meds_gt <- matrix(
        apply(X      = getC_area(area = col, gears = sokGear),
              MARGIN = 2, FUN = median, na.rm = TRUE),
        nrow = 1)
      plotStackedBars(
        meds_gt = meds_gt,
        gColors = adjustcolor("steelblue", alpha.f = 0.7))
    }

    segments( x0  = yrs[tMP:nT], x1 = yrs[tMP:nT],
              y0  = C_tot_qt[1, tMP:nT],
              y1  = C_tot_qt[3, tMP:nT],
              col = "black", lwd = 0.8 )

    abline( v = yrs[tMP], col = "grey30", lty = 3 )
  }


  mtext( side = 1, outer = TRUE, text = "Year", line = 2.5 )

  mtext( outer = TRUE, side = 1, adj = 0.9,
         line = 3.5, cex = 0.6,
         text = stamp, col = "grey60" )

} # END plotMPSummary()
