# --------------------------------------------------------------------------
# refPts.R
#
# Takes parameters from a report object and calculates
# reference points and curves under an age structured/YPR formulation.
# Called from within runHierSCAL() after runnings an AM.
#
#
# Author: Samuel Johnson
# Date: March 7, 2019
#
# --------------------------------------------------------------------------

# calcRefPts()
# Calculates biological reference points based on given biological parameters
# assuming a delay difference biological model.
# inputs:   obj = list of biological parameters
# ouputs:   refPts = list() of reference points
calcRefPts <- function( obj )
{
  # Calculate selectivity
  # obj <- .calcSel_xsp(obj, fleetIdx = 2)
  nS <- obj$om$nS
  nP <- obj$om$nP
  nF <- obj$om$nF

  nSpec   <- nS

  # temporarily use calcPerRecruit to recalc R0
  obj$relF_spf    <- array(0,dim = c(nS,nP,nF))
  tmp             <- .calcPerRecruit( f = 0, obj = obj )
  yprList         <- tmp$yprList

  # Calculate R0_sp and totB
  h_sp      <- obj$om$h_sp
  B0_sp     <- obj$om$B0_sp

  obj$R0_sp     <- array(0,dim = c(nS,nP))
  obj$totB0_sp  <- array(0,dim = c(nS,nP))

  obj$M0_xsp    <- yprList$Meq_xsp

  # Identify valid stocks for each species based on allocation
  # A stock is valid if it has any non-zero allocation across all fleets
  validStock_sp <- apply(X = obj$om$alloc_spf, FUN = function(x) sum(x) > 0,
                         MARGIN = c(1,2))

  obj$R0_sp <- B0_sp / yprList$ssbpr_sp
  obj$totB0_sp <- obj$R0_sp * yprList$totbpr_sp

  # Set invalid stocks to NA to prevent them from being used in calculations
  obj$R0_sp[!validStock_sp] <- NA
  obj$totB0_sp[!validStock_sp] <- NA

  lastB_sp <<- obj$totB0_sp

  # Beverton-Holt a/b parameters
  obj$rec.a_sp  <- 4.*h_sp*obj$R0_sp/(B0_sp*(1.-h_sp))
  obj$rec.b_sp  <- (5.*h_sp-1.)/(B0_sp*(1.-h_sp))

  # Need to calculate relF_spf
  # Calculate relative F

  alloc_spf       <- obj$om$alloc_spf
  commIdx_spf     <- alloc_spf > 0

  # In case there's only one fishery...
  commIdx_spf[alloc_spf == 1] <- FALSE
  obj$relF_spf[alloc_spf == 1] <- 1


  if(any(commIdx_spf))
  {
    nComm           <- sum(commIdx_spf)

    # Test the objective function with initial parameters
    testVal <- .getRelFs(lnfg = rep(-2, nComm),
                         f = mean(obj$om$M_xsp),
                         obj = obj,
                         commIdx_spf = commIdx_spf,
                         allocVar = obj$allocVar)

    if(!is.finite(testVal)) {
      cat("Warning: Initial value in .getRelFs is not finite:\n")
      cat("  testVal =", testVal, "\n")
      cat("  nComm =", nComm, "\n")
      cat("  mean(M_xsp) =", mean(obj$om$M_xsp), "\n")
      cat("  Trying alternate initialization...\n")

      # Try different initial values
      initVals <- c(-2, -1, -0.5, 0, 0.5)
      for(initVal in initVals) {
        testVal <- .getRelFs(lnfg = rep(initVal, nComm),
                             f = mean(obj$om$M_xsp),
                             obj = obj,
                             commIdx_spf = commIdx_spf,
                             allocVar = obj$allocVar)
        if(is.finite(testVal)) {
          cat("  Found finite initial value with par =", initVal, "\n")
          break
        }
      }

      if(!is.finite(testVal)) {
        stop("Unable to find finite initial values for optimization. Check obj$totB0_sp, obj$om$alloc_spf, and per-recruit calculations.")
      }
    } else {
      initVal <- -2
    }

    optRelF         <- optim( par = rep(initVal, nComm), fn = .getRelFs,
                              method = "BFGS", control=list(maxit = 50),
                              f = mean(obj$om$M_xsp), obj = obj,
                              commIdx_spf = commIdx_spf,
                              allocVar = obj$allocVar )

    # Overwrite relF_pf
    idx <- 0
    for( s in 1:nS )
    {
      for(p in 1:nP)
      {
        fIdx <- which(commIdx_spf[s,p,])
        nF <- length(fIdx)
        obj$relF_spf[s,p,fIdx] <- exp(optRelF$par[idx + 1:nF])
      }

      idx <- idx + nF
    }

    # Normalise
    for(s in 1:nS)
      for( p in 1:nP)
          obj$relF_spf[s,p,] <- obj$relF_spf[s,p,]/sum(obj$relF_spf[s,p,])
  }

  # Calculate generation time
  obj$genTime_sp <- .calcGenTime(obj)

  lastB_sp <<- obj$totB0_sp

  # Calculate reference curves
  refCurves <- .calcRefCurves( obj )

  # First, let's just do Fmsy reference points
  FmsyRefPts <- .getFmsy_sp(  obj = obj,
                              refCurves = refCurves )




  obj$refPts <- list()
  obj$refPts$refCurves    <- refCurves
  obj$refPts$FmsyRefPts   <- FmsyRefPts

  # message("Fmsy = ", FmsyRefPts$Fmsy_sp, "\n")
  # browser()
  # message("Bmsy = ", FmsyRefPts$S_sp, "\n")

  # totB0_sp[1:nS,1:nP] <- refCurves$totBeq_spf[,,1]


  # Get survivorship
  obj$refPts$Surv_axsp    <- tmp$Surv_axsp
  obj$refPts$ssbpr_sp     <- yprList$ssbpr_sp
  obj$refPts$R0_sp        <- obj$R0_sp
  obj$refPts$totB0_sp     <- obj$totB0_sp
  obj$refPts$totbpr_sp    <- yprList$totbpr_sp
  obj$refPts$rec.a_sp     <- obj$rec.a_sp
  obj$refPts$rec.b_sp     <- obj$rec.b_sp
  obj$refPts$B0_sp        <- refCurves$Beq_spf[,,1,drop = FALSE]
  obj$refPts$M_xsp        <- obj$om$M_xsp
  obj$refPts$m1_sp        <- obj$om$m1_sp
  obj$refPts$Mb_sp        <- obj$om$Mb_sp
  obj$refPts$h_sp         <- h_sp
  obj$refPts$relF_spf     <- obj$relF_spf



  return(obj$refPts)

} # END calcRefPts()

# .getRelFs()
# Function used for determining relative fishing
# mortality rates that meet a nominated allocation
# among gears.
.getRelFs <- function( lnfg, obj, f=0, commIdx_spf, allocVar = "Biomass" )
{
  nS <- obj$om$nS
  nP <- obj$om$nP
  nF <- obj$om$nF

  # Overwrite relF_pf
  idx <- 0
  for( s in 1:nS )
  {
    for(p in 1:nP)
    {
      fIdx <- which(commIdx_spf[s,p,])
      nF <- length(fIdx)
      obj$relF_spf[s,p,fIdx] <- exp(lnfg[idx + 1:nF])
    }

    idx <- idx + nF
  }


  # Normalise
  for(s in 1:nS)
    for( p in 1:nP)
        obj$relF_spf[s,p,] <- obj$relF_spf[s,p,]/sum(obj$relF_spf[s,p,])

  lastB_sp <<- obj$totB0_sp

  obj$relF_spf[!is.finite(obj$relF_spf)] <- 0
  # Calculate legal YPR
  tmp         <- .calcPerRecruit( f=f, obj )$yprList

  funcVal <- 0
  for(s in 1:nS)
    for(p in 1:nP)
    {
      commGears <- which(obj$om$alloc_spf[s,p,] > 0)

      # Skip if no commercial gears for this stock/species
      if(length(commGears) == 0)
        next

      # Skip if only one fleet (100% allocation) - already handled in lines 67-69
      if(any(obj$om$alloc_spf[s,p,] == 1))
        next

      # Switch here between egg yield and bio/ponded fish yield
      if(allocVar == "Eggs")
        prop <- tmp$epr_spf[s,p,commGears]/sum(tmp$epr_spf[s,p,commGears])
      if( allocVar == "Biomass")
        prop <- tmp$ypr_spf[s,p,commGears]/sum(tmp$ypr_spf[s,p,commGears])

      # Add small epsilon to prevent log(0) and handle NaN/Inf
      epsilon <- 1e-10
      prop <- pmax(prop, epsilon)
      alloc <- pmax(obj$om$alloc_spf[s,p,commGears], epsilon)

      # Check for non-finite values before calculating
      if(any(!is.finite(prop)) || any(!is.finite(alloc))) {
        cat("Warning in .getRelFs: non-finite values detected\n")
        cat("  s =", s, "p =", p, "\n")
        cat("  prop =", prop, "\n")
        cat("  alloc =", alloc, "\n")
        return(Inf)
      }

      # cat("s = ",s,"; propYield = ",prop,"\n")
      funcVal <- funcVal + sum((log(prop)-log(alloc))^2.)
      # cat("f = ", funcVal,"\n")
    }


  funcVal
}


# .calcRefCurves()
# Calculates equilibrium curves of equilbrium biomass, numbers,
# yield and recruitment as a function of input fishing mortality rates
# inputs:   obj = list of biological parameters
# ouputs:   refCurves = list() of reference curves (vectors)
.calcRefCurves <- function( obj, nFs = 200 )
{
  # First, compute max F (tolerance of 1e-5)
  nT   <- dim(obj$om$qF_spft)[4]
  maxF <- 10 * max(obj$om$M_xsp,na.rm = T)

  if( obj$condModel == "hierSCAL")
    maxE <- max( maxF / obj$om$qF_spft[,,2,nT])

  # We're going to need to fill each species' ref curves,
  # so labeling and dimensions are needed
  nS          <- obj$om$nS
  nP          <- obj$om$nP
  nA          <- obj$om$nA
  nX          <- obj$om$nX

  specNames   <- c("Herring","Hake")[1:nS]
  stockNames  <- dimnames(obj$M_apt)[[2]]

  f <- seq( from = 0.0, to = maxF, length = nFs )

  if( obj$condModel == "hierSCAL")
    e <- seq( from = 0.0, to = maxE, length = nFs )



  # Create matrices to hold Recruitment reference curve, name rows and copy
  # for each of Beq, Neq and Yeq
  Req_spf      <- array( NA,  dim = c(nS, nP, nFs),
                              dimnames = list(  species = specNames,
                                                stock = stockNames,
                                                F = f ) )

  surv_axspf  <- array( NA, dim = c(nA,nX,nS,nP,nFs) )
  Z_axspf     <- array( NA, dim = c(nA,nX,nS,nP,nFs) )
  Meq_xspf    <- array( NA, dim = c(nX,nS,nP,nFs) )

  Beq_spf       <- Req_spf
  totBeq_spf    <- Req_spf
  Yeq_spf       <- Req_spf
  ypr_spf       <- Req_spf
  ssbpr_spf     <- Req_spf
  Ueq_spf       <- Req_spf
  Yeq_spf       <- Req_spf
  EYeq_spf      <- Req_spf
  ypr_spf       <- Req_spf
  epr_spf       <- Req_spf
  ssbpr_spf     <- Req_spf
  matbpr_spf    <- Req_spf
  Ueq_spf       <- Req_spf
  Udead_spf     <- Req_spf


  # Loop and fill
  for( i in 1:length(f) )
  {
    tmp               <- .calcEquil( f = f[i], obj = obj )
    Req_spf[,,i]      <- tmp$Req_sp
    Beq_spf[,,i]      <- tmp$Beq_sp
    # expBeq_spf[,,i]   <- tmp$expBeq_sp
    totBeq_spf[,,i]   <- tmp$totBeq_sp
    Yeq_spf[,,i]      <- tmp$Yeq_sp
    EYeq_spf[,,i]     <- tmp$EYeq_sp
    ypr_spf[,,i]      <- tmp$ypr_sp
    epr_spf[,,i]      <- tmp$epr_sp
    ssbpr_spf[,,i]    <- tmp$ssbpr_sp
    Ueq_spf[,,i]      <- tmp$Ueq_sp
    surv_axspf[,,,,i] <- tmp$surv_axsp
    Meq_xspf[,,,i]    <- tmp$Meq_xsp
    Z_axspf[,,,,i]    <- tmp$Z_axsp
  }


  # Save F based ref points
  refCurves <- list()
    refCurves$F           <- f
    refCurves$ypr_spf     <- ypr_spf
    refCurves$Req_spf     <- Req_spf
    refCurves$Beq_spf     <- Beq_spf
    # refCurves$expBeq_spf  <- expBeq_spf
    refCurves$totBeq_spf  <- totBeq_spf
    refCurves$Yeq_spf     <- Yeq_spf
    refCurves$EYeq_spf    <- EYeq_spf
    refCurves$Ueq_spf     <- Yeq_spf
    refCurves$Meq_xspf    <- Meq_xspf
    refCurves$surv_axspf  <- surv_axspf
    refCurves$Z_axspf     <- Z_axspf


  return( refCurves )
} # END .calcRefCurves

# .calcEquil()
# Calculates equilibrium Biomass, Yield and Recruitment
# with respect to given fishing mortality rates F
# inputs:   f = input fishing mortality rate
#           obj = list of biological parameters
# ouputs:   equil = list() of equilibrium biomass, yield and recruitment
.calcEquil <- function( f = 0, obj, type = "fmort" )
{
  nS    <- obj$om$nS
  nMICE <- obj$om$nMICE
  nP    <- obj$om$nP
  A_s   <- obj$om$A_s


  # Now calculate eqbm recruitment at given f value
  tmp <- .calcPerRecruit( f = f, obj = obj, type = type )
  yprList <- tmp$yprList



  recruits_sp <- ( obj$rec.a_sp * yprList$ssbpr_sp - 1) / (obj$rec.b_sp * yprList$ssbpr_sp)


  equil <- list()
    equil$Req_sp      <- recruits_sp
    equil$Beq_sp      <- recruits_sp * yprList$ssbpr_sp
    # equil$expBeq_sp  <- recruits_sp * yprList$expbpr_sp
    equil$totBeq_sp   <- recruits_sp * yprList$totbpr_sp
    equil$Yeq_sp      <- recruits_sp * yprList$ypr_sp
    equil$EYeq_sp     <- recruits_sp * yprList$epr_sp
    equil$ypr_sp      <- yprList$ypr_sp
    equil$epr_sp      <- yprList$epr_sp
    equil$ssbpr_sp    <- yprList$ssbpr_sp
    equil$Ueq_sp      <- equil$Yeq_sp / (equil$Beq_sp + equil$Yeq_sp)
    equil$surv_axsp   <- tmp$Surv_axsp
    equil$Meq_xsp     <- yprList$Meq_xsp
    equil$Z_axsp      <- yprList$Z_axsp

  return(equil)
}


.calcSel_xsp <- function( obj, fleetIdx = 2)
{
  # Model dimensions
  nS    <- obj$nS
  nP    <- obj$nP
  nX    <- obj$nX
  A_s   <- obj$A_s
  nA    <- obj$nA
  nL    <- obj$nL

  # Probability of being a given length-at-age
  probLenAge_laspx <- obj$probLenAge_laspx

  # Selectivity - this is mean over each fleet's
  # time period, not time-varying
  # Might be useful to take the group mean for the comm fleet...
  xSel50_sp   <- obj$xSel50_spf[ , , fleetIdx, drop = FALSE ]
  xSel95_sp   <- obj$xSel95_spf[ , , fleetIdx, drop = FALSE ]
  xSelStep_sp <- obj$xSelStep_spf[ , , fleetIdx, drop = FALSE ]

  # Harcode for comm.mod for now,
  # and length based selectivity
  selLen_lsp  <- array(NA, dim = c(nL,nS,nP) )
  selAge_axsp <- array(NA, dim = c(nA,nX,nS,nP) )

  # Loop over species and stocks, calculate
  # selLen so we can get selAge
  for( s in 1:nS )
    for(p in 1:nP )
    {
      selLen_lsp[,s,p] <- (1 + exp(-(1:nL - xSel50_sp[s,p,])/xSelStep_sp[s,p,]/log(19)))^(-1)
      for( x in 1:nX)
      {
        for( a in 1:nA )
        {
          selAge_axsp[a,x,s,p] <- sum(probLenAge_laspx[,a,s,p,x] * selLen_lsp[,s,p])
        }
        selAge_axsp[,x,s,p] <- selAge_axsp[,x,s,p] / max(selAge_axsp[,,s,p],na.rm = T)
      }
    }

  obj$selLen_lsp <- selLen_lsp
  obj$selAge_axsp <- selAge_axsp

  obj
}

# .calcPerRecruit
# Purpose:     Calculate all equilibrium per-recruit quantities of interest for an
#              input fishing mortality.
# Parameters:  f=scalar input fishing mortality rate; obj=list of all operating
#              model parameters.
# Returns:     a list with equilibrium quantities - (i) spawning stock biomass-per-recruit
#              and (ii) yield-per-recruit (ypr)
# Source:      S.P. Cox, modified for hierSCAL by SDNJ
.calcPerRecruit <- function( f, obj, type = "fmort" )
{

  # Compute eqbm spawning biomass per recruit for
  # given f and species/stock pars
  nS      <- obj$om$nS
  nP      <- obj$om$nP
  nF      <- obj$om$nF
  A_s     <- obj$om$A_s
  nA      <- obj$om$nA
  nL      <- obj$om$nL
  nX      <- obj$om$nX
  M_xsp   <- obj$om$M_xsp
  Mb_sp   <- obj$om$Mb_sp
  m1_sp   <- obj$om$m1_sp
  nT      <- obj$om$nT
  tMP     <- obj$om$tMP
  qF_spf  <- obj$om$qF_spft[,,,nT]

  nSpec   <- nS

  # Life history schedules
  matAge_asp        <- obj$om$matAge_asp
  wtAge_axsp        <- obj$om$meanWtAge_axsp

  # Spawn timing
  spawnTiming_s     <- obj$om$spawnTiming_s

  # Recruitment pars
  if(f > 0)
  {
    rec.a_sp <- obj$rec.a_sp
    rec.b_sp <- obj$rec.b_sp
    totB0_sp <- obj$totB0_sp
  }



  # Selectivity
  selAge_axspf           <- array(NA, dim = c(nA,nX,nS,nP,nF))
  selAge_axspf[1:nA,1:nX,1:nS,1:nP,1:nF]  <- obj$om$sel_axspft[1:nA,1:nX,1:nS,1:nP,1:nF,tMP-1]

  # Fishing mortality
  relF_spf  <- obj$relF_spf
  relF_spf[relF_spf < 1e-5] <- 0
  f_spf <- f * relF_spf

  Mb_sp         <- array(0,dim = c(nS,nP))
  totB0_sp      <- array(0,dim = c(nS,nP))
  m1_sp         <- obj$om$m1_sp

  for(s in 1:nS )
  {
    totB0_sp[s,] <- obj$totB0_p
    for(p in 1:nP )
      Mb_sp     <- M_xsp[1,s,p] - exp(-m1_sp[s,p])
  }

  if(obj$densityDepM == 1)
  {
    if(f == 0 )
    {
      Meq_sp    <- Mb_sp + exp(-m1_sp)
      lastB_sp  <<- totB0_sp
    }
    if( f > 0 )
    {
      lnB_sp <- log(lastB_sp)
      if( nS > 1 | nP > 1)
      {
        opt <- optim( par = lnB_sp,
                      f = solveForMeq,
                      ff = f, obj = obj,
                      fit = TRUE )
        lnB_sp[1:nS,1:nP] <- opt$par
      }

      if(nS == 1 & nP == 1)
      {
        opt <- optimize(  f = solveForMeq,
                          interval = c(log(0.5*lastB_sp), lnB_sp),
                          ff = f, obj = obj, fit = TRUE )
        lnB_sp[1:nS,1:nP] <- opt$minimum
      }


      lastB_sp[1:nS,1:nP] <<- exp(lnB_sp)



      Meq_sp <- solveForMeq( lnB_sp = lnB_sp, obj = obj,ff = f, fit = FALSE)
    }

    # cat("lastB_sp = ", lastB_sp, "\n", sep = "")

    M_xsp[1,,] <- Meq_sp
  }

  # message("lastB = ", lastB_sp, "\n")

  # Zero indexing
  juveMage <- obj$juveMage


  # Compute Z_asp
  Z_axsp    <- array( NA, dim = c(nA,nX,nS,nP))
  Surv_axsp <- array( NA, dim = c(nA,nX,nS,nP))
  Surv_axsp[1,,,] <- 1/nX


  for( x in 1:nX )
    for( s in 1:nS )
      for( p in 1:nP )
      {

        if( juveMage > 0 )
          Z_axsp[1:(juveMage),x,s,p] <- obj$Mjuve_p[p]

        Z_axsp[(juveMage+1):A_s[s],x,s,p] <- M_xsp[x,s,p]

        for( a in 1:A_s[s])
        {
          for( f in 1:nF)
            Z_axsp[a,x,s,p] <- Z_axsp[a,x,s,p] + selAge_axspf[a,x,s,p,f] * f_spf[s,p,f]

          if( a > 1 )
            Surv_axsp[a,x,s,p] <- Surv_axsp[a-1,x,s,p] * exp( -Z_axsp[a-1,x,s,p])
          if( a == A_s[s])
            Surv_axsp[a,x,s,p] <- Surv_axsp[a,x,s,p] / (1 - exp(-Z_axsp[a,x,s,p]))
        }
      }


  # Calculate yield-per-recruit, ssb per recruit, and exp biomass per recruit
  # by using survival array
  ssbpr_asp     <- array( NA, dim = c(nA,nS,nP) )
  matbpr_asp    <- array( NA, dim = c(nA,nS,nP) )
  expbpr_axspf  <- array( NA, dim = c(nA,nX,nS,nP,nF) )
  totbpr_axsp   <- array( NA, dim = c(nA,nX,nS,nP) )
  C_axspf       <- array(0,   dim = c(nA,nX,nS,nP,nF))
  P_axspf       <- array(0,   dim = c(nA,nX,nS,nP,nF))
  E_axspf       <- array(0,   dim = c(nA,nX,nS,nP,nF))


  for( s in 1:nS )
    for( p in 1:nP)
    {
      ssbpr_asp[,s,p]  <- Surv_axsp[,nX,s,p] *
                          wtAge_axsp[,nX,s,p] *
                          matAge_asp[,s,p]

      # Keep track of separate ssbpr and mature biomass. SSB is
      # fish that are used in the SR relationship, while mature biomass
      # includes spawners whose eggs are harvested by SOK (or predators)
      ssbpr_asp[,s,p] <- ssbpr_asp[,s,p] * exp(-spawnTiming_s[s] * Z_axsp[,nX,s,p])
      matbpr_asp[,s,p] <- ssbpr_asp[,s,p]


      for( x in 1:nX )
      {
        for( f in 1:nF)
        {
          # Catch - total ponded fish for SOK fleets
          C_axspf[,x,s,p,f] <-  Surv_axsp[,x,s,p] * wtAge_axsp[,x,s,p] *
                                selAge_axspf[,x,s,p,f] * f_spf[s,p,f] *
                                (1 - exp(-Z_axsp[,x,s,p]))/Z_axsp[,x,s,p]

          # Egg yield
          E_axspf[,x,s,p,f]   <- C_axspf[,x,s,p,f] * matAge_asp[,s,p] * 0.5 * 200 * 0.35
          # Convert eggs to SOK and overwite catch,
          # record "catch" based on F as ponded fish
          if( obj$om$fleetType_f[f] == 2 )
          {
            P_axspf[,x,s,p,f] <-  C_axspf[,x,s,p,f]

            C_axspf[,x,s,p,f] <-  E_axspf[,x,s,p,f] *
                                  obj$om$gamma_f[f]

            # Return ponded fish to mature biomass, but
            # they don't contribute to spawners
            matbpr_asp[,s,p]  <- matbpr_asp[,s,p] + exp(-0.315)* P_axspf[,x,s,p,f]
          }

          expbpr_axspf[,x,s,p,f] <- Surv_axsp[,x,s,p] *
                                    selAge_axspf[,x,s,p,f] *
                                    wtAge_axsp[,x,s,p]


        }


        totbpr_axsp[,x,s,p] <- Surv_axsp[,x,s,p] * wtAge_axsp[,x,s,p]
      }

    }

  # Replace NAs with 0 (unmodeled ages)
  # Z_axsp[is.na(Z_axsp)] <- 0
  # Surv_axsp[is.na(Surv_axsp)] <- 0
  # C_axsp[is.na(C_axsp)] <- 0

  # Compute ratios
  ypr_sp        <- apply( X = C_axspf, FUN = sum, MARGIN = c(3,4),na.rm = T)
  ypr_spf       <- apply( X = C_axspf, FUN = sum, MARGIN = c(3:5),na.rm = T)
  epr_sp        <- apply( X = E_axspf, FUN = sum, MARGIN = c(3,4),na.rm = T)
  epr_spf       <- apply( X = E_axspf, FUN = sum, MARGIN = c(3:5),na.rm = T)
  ssbpr_sp      <- apply( X = ssbpr_asp, FUN = sum, MARGIN = c(2,3), na.rm = T )
  matbpr_sp     <- apply( X = matbpr_asp, FUN = sum, MARGIN = c(2,3), na.rm = T )
  expbpr_spf    <- apply( X = expbpr_axspf, FUN = sum, MARGIN = 3:5, na.rm = T )
  totbpr_sp     <- apply( X = totbpr_axsp[2:nA,,,,drop = FALSE], FUN = sum, MARGIN = c(3,4), na.rm = T )


  # if(any(expbpr_sp - ypr_sp < 0))

  # compile output list
  yprList <- list(  ssbpr_sp    = ssbpr_sp,
                    ypr_sp      = ypr_sp,
                    ypr_spf     = ypr_spf,
                    epr_sp      = epr_sp,
                    epr_spf     = epr_spf,
                    expbpr_spf  = expbpr_spf,
                    totbpr_sp   = totbpr_sp,
                    matbpr_sp   = matbpr_sp,
                    Meq_xsp     = M_xsp,
                    Z_axsp      = Z_axsp )

  obj$yprList <- yprList
  obj$Surv_axsp <- Surv_axsp

  return(obj)
}

# solve for density dependent M
solveForMeq <- function( lnB_sp = log(totB0_sp), obj, ff, fit = TRUE )
{
  # Compute eqbm spawning biomass per recruit for
  # given f and species/stock pars
  # Compute eqbm spawning biomass per recruit for
  # given f and species/stock pars
  nS      <- obj$om$nS
  nP      <- obj$om$nP
  nF      <- obj$om$nF
  A_s     <- obj$om$A_s
  nA      <- obj$om$nA
  nL      <- obj$om$nL
  nX      <- obj$om$nX
  M_xsp   <- obj$om$M_xsp
  Mb_sp   <- obj$om$Mb_sp
  m1_sp   <- obj$om$m1_sp
  nT      <- obj$om$nT
  tMP     <- obj$om$tMP
  qF_spf  <- obj$om$qF_spft[,,,nT]


  totB0_sp <- obj$totB0_sp

  # Life history schedules
  matAge_asp        <- obj$om$matAge_asp
  wtAge_axsp        <- obj$om$meanWtAge_axsp

  # Spawn timing
  spawnTiming_s     <- obj$om$spawnTiming_s

  # Selectivity
  selAge_axspf           <- array(NA, dim = c(nA,nX,nS,nP,nF))
  selAge_axspf[1:nA,1:nX,1:nS,1:nP,1:nF]  <- obj$om$sel_axspft[1:nA,1:nX,1:nS,1:nP,1:nF,tMP-1]


  # Fishing mortality
  relF_spf  <- obj$relF_spf
  relF_spf[relF_spf < 1e-3] <- 0
  f_spf <- ff * relF_spf

  # DepM model pars
  Mb_sp         <- obj$om$Mb_sp
  m1_sp         <- obj$om$m1_sp

  if(is.null(Mb_sp))
    for(s in 1:nS )
      for(p in 1:nP )
        Mb_sp     <- M_xsp[1,s,p] - exp(-m1_sp[s,p])

  juveMage      <- obj$juveMage
  Mjuve_p       <- obj$Mjuve_p

  # Recruitment pars
  if(ff > 0)
  {
    rec.a_sp <- obj$rec.a_sp
    rec.b_sp <- obj$rec.b_sp
  }

  # message("totB0 = ", totB0_sp, "\n")
  initTotBeq_sp <- exp(lnB_sp)
  Meq_sp        <- Mb_sp + exp(-m1_sp * initTotBeq_sp/totB0_sp)

  if(ff == 0)
    return(Meq_sp)

  ssbpr_sp <- array(0,dim = c(nS,nP))
  totbpr_sp <- array(0,dim = c(nS,nP))

  # Compute Z_asp
  Z_axsp    <- array( NA, dim = c(nA,nX,nS,nP))
  Surv_axsp <- array( NA, dim = c(nA,nX,nS,nP))
  Surv_axsp[1,,,] <- 1/nX


  for( x in 1:nX )
    for( s in 1:nS )
      for( p in 1:nP )
      {

        if( juveMage > 0 )
          Z_axsp[1:(juveMage),x,s,p] <- Mjuve_p[p]

        Z_axsp[(juveMage+1):A_s[s],x,s,p] <- Meq_sp[s,p]

        for( a in 1:A_s[s])
        {
          for( fIdx in 1:nF)
            Z_axsp[a,x,s,p] <- Z_axsp[a,x,s,p] + selAge_axspf[a,x,s,p,fIdx] * f_spf[s,p,fIdx]

          if( a > 1 )
            Surv_axsp[a,x,s,p] <- Surv_axsp[a-1,x,s,p] * exp( -Z_axsp[a-1,x,s,p])
          if( a == A_s[s])
            Surv_axsp[a,x,s,p] <- Surv_axsp[a,x,s,p] / (1 - exp(-Z_axsp[a,x,s,p]))
        }
      }


  # Calculate ssbpr and totbpr
  for( p in 1:nP )
  {
    for( s in 1:nS )
    {
      ssbpr_sp[s,p]  <- sum(Surv_axsp[,,s,p] * wtAge_axsp[,,s,p] * matAge_asp[,s,p] * exp(-spawnTiming_s[s] * Z_axsp[,,s,p]))
      totbpr_sp[s,p] <- sum(Surv_axsp[(juveMage+1):nA,,s,p] * wtAge_axsp[(juveMage+1):nA,,s,p])
    }

  }

  guessR_sp     <- initTotBeq_sp/totbpr_sp
  guessSB_sp    <- guessR_sp * ssbpr_sp
  guessR2_sp    <- (rec.a_sp * ssbpr_sp - 1) / (rec.b_sp * ssbpr_sp)

  guessTotB_sp  <- guessR2_sp * totbpr_sp

  # guessR22_p <- rec.a_sp[1,] * guessSB_p / (1 + rec.b_sp[1,] * guessSB_p)

  # message( "totbpr = ", totbpr_sp, "\n")
  # message( "ssbpr = ", ssbpr_sp, "\n")
  # message( "guessB = ", initTotBeq_sp, "\n")
  # message( "guessB2 = ", guessTotB_sp, "\n")
  # message( "guessB/tbpr = ", guessR_sp, "\n")
  # message( "ssbpr Rec = ", guessR2_sp, "\n")
  # message( "ssbpr rec * totbpr = ", guessTotB_sp, "\n")

  if( fit )
  {
    resid <- sum((log(guessR2_sp) - log(guessR_sp))^2)
    # message("optB f = ", resid, "\n")
    return(resid)
  }

  return(Meq_sp)

}


# # Calculates recruitment parameters, and equilibrium unfished
# # numbers and recruitment.
# .calcRecPars <- function( obj )
# {
#   # Calculate eqbm parameters
#   # Survival


#   # Beverton-Holt a parameters
#   rec.a <- 4.*obj$rSteepness*R0/(obj$B0*(1.-obj$rSteepness))
#   rec.b <- (5.*obj$rSteepness-1.)/(obj$B0*(1.-obj$rSteepness))

#   # Now return everything in a list
#   recList <- list(  S0 = S0,
#                     wbar0 = wbar0,
#                     N0 = N0,
#                     R0 = R0,
#                     rec.a = rec.a,
#                     rec.b = rec.b  )

#   return(recList)
# }

# .getEmsy_p()
# Calculates effort based reference curves
# for each stock and species, then computes
# multispecies optimum yield and effort
# values for each stock area. Returns
# arrays for plotting of yield curves.
.getEmsy_p <- function( obj, refCurves, FmsyRefPts )
{
  nS  <- obj$nS
  nP  <- obj$nP

  Eseq <- refCurves$E

  .getEmsy <- function( yieldCurve, E = Eseq )
  {
    minE <- 0
    maxE <- max(E)

    # Now check if yieldCurve goes negative anywhere
    if( any(yieldCurve < 0) )
    {
      minE <- E[min(which(yieldCurve >= 0))]
      maxE <- E[max(which(yieldCurve >= 0))]
    }

    eySplineFun <- splinefun( x=E, y=yieldCurve )

    # Find stat point for Fmsy
    Emsy <- try( uniroot( f = eySplineFun, interval = c(minE, maxE),
                          deriv = 1 )$root )
    if(class(Emsy) == "try-error")
    {
      browser(cat("uniroot failed\n\n"))
      Emsy <- 0
    }

    Emsy <- min(Emsy, maxE)

    return(Emsy)
  }

  # calculate Fmsy for each species/stock
  Emsy_sp <-  apply(  X = refCurves$Yeq_spe, FUN = .getEmsy,
                      MARGIN = c(1,2) )

  Yeq_spe <- refCurves$Yeq_spe
  Yeq_spe[Yeq_spe < 0] <- NA

  # Now we want the total stock area yield curve
  Yeq_pe  <- apply( X  = refCurves$Yeq_spe, FUN = sum, MARGIN = c(2,3),
                    na.rm = T )
  Yeq_pe[Yeq_pe == 0] <- NA
  Yeq_pe[,1] <- 0
  Emsy_p  <- apply( X = Yeq_pe, FUN = .getEmsy, MARGIN = 1)

  # Now create a matrix to hold the species/stock values
  # on each curve
  spMat       <- matrix(0, nrow = nS, ncol = nP)
  dimnames(spMat) <- dimnames(Emsy_sp)

  EmsyRefPts  <- list(  yprEmsy_sp    = spMat,
                        YeqEmsy_sp    = spMat,
                        BeqEmsy_sp    = spMat,
                        ReqEmsy_sp    = spMat,
                        expBeqEmsy_sp = spMat,
                        Umsy_sp       = spMat )

  EmsyMSRefPts <- list( yprEmsy_sp    = spMat,
                        YeqEmsy_sp    = spMat,
                        BeqEmsy_sp    = spMat,
                        ReqEmsy_sp    = spMat,
                        expBeqEmsy_sp = spMat,
                        Umsy_sp       = spMat,
                        YeqEmsy_p     = rep(0,nP) )


  # Calculate ref points
  EmsyRefPts$Emsy_sp    <- Emsy_sp
  EmsyMSRefPts$EmsyMS_p   <- Emsy_p

  # Loop and get each stock/species ref pt
  for( p in 1:nP )
  {
    for( s in 1:nS )
    {
      tmp <- .calcEquil( f = Emsy_sp[s,p], obj = obj, type = "effort" )
      EmsyRefPts$yprEmsy_sp[s,p]    <- tmp$ypr_sp[s,p]
      EmsyRefPts$YeqEmsy_sp[s,p]    <- tmp$Yeq_sp[s,p]
      EmsyRefPts$BeqEmsy_sp[s,p]    <- tmp$Beq_sp[s,p]
      EmsyRefPts$expBeqEmsy_sp[s,p] <- tmp$expBeq_sp[s,p]
      EmsyRefPts$NeqEmsy_sp[s,p]    <- tmp$Neq_sp[s,p]
      EmsyRefPts$ReqEmsy_sp[s,p]    <- tmp$Req_sp[s,p]
      EmsyRefPts$Umsy_sp[s,p]       <- tmp$Ueq_sp[s,p]

      tmpMS <- .calcEquil( f = Emsy_p[p], obj = obj, type = "effort" )
      EmsyMSRefPts$yprEmsy_sp[s,p]    <- tmpMS$ypr_sp[s,p]
      EmsyMSRefPts$YeqEmsy_sp[s,p]    <- tmpMS$Yeq_sp[s,p]
      EmsyMSRefPts$BeqEmsy_sp[s,p]    <- tmpMS$Beq_sp[s,p]
      EmsyMSRefPts$expBeqEmsy_sp[s,p] <- tmpMS$expBeq_sp[s,p]
      EmsyMSRefPts$NeqEmsy_sp[s,p]    <- tmpMS$Neq_sp[s,p]
      EmsyMSRefPts$ReqEmsy_sp[s,p]    <- tmpMS$Req_sp[s,p]
      EmsyMSRefPts$Umsy_sp[s,p]       <- tmpMS$Ueq_sp[s,p]

    }
    # Now sum species eq yield for complex opt yield
    EmsyMSRefPts$YeqEmsy_p[p]         <- sum( EmsyMSRefPts$YeqEmsy_sp[,p] )
  }

  outlist <- list(  EmsyRefPts    = EmsyRefPts,
                    EmsyMSRefPts  = EmsyMSRefPts )
}

# .getFmsy     ()
# Purpose:     fit a spline function to f vs yield, then use a root finder to get Fmsy.
# Parameters:  obj=list of all operating model parameters, schedules, equilibrium functions
# Returns:     a list with all equilibrium quantities for Fmsy
# Source:      S.P. Cox, modified for hierSCAL by SDNJ
.getFmsy_sp <- function( obj, refCurves )
{
  Fseq  <- refCurves$F
  nS    <- obj$om$nS
  nP    <- obj$om$nP
  nMICE <- obj$om$nMICE

  .getFmsy <- function( yieldCurve, F = Fseq )
  {
    minF <- 0
    maxF <- max(F)

    # Check if yield curve is all zeros or near-zero (e.g., stocks with zero selectivity)
    if(max(yieldCurve, na.rm = TRUE) < 1e-10) {
      return(0)
    }

    # Now check if yieldCurve goes negative anywhere
    if( any(yieldCurve < 0) )
    {
      minF <- F[min(which(yieldCurve >= 0))]
      maxF <- F[max(which(yieldCurve >= 0))]
    }

    fySplineFun <- splinefun( x=F, y=yieldCurve )
    # Find stat point for Fmsy
    Fmsy <- try( uniroot( f = fySplineFun, interval = c(minF, maxF),
                          deriv = 1 )$root, silent = TRUE )

    if(class(Fmsy) == "try-error")
    {
      # browser(cat("uniroot failed\n\n"))
      Fmsy <- 0
    }

    Fmsy <- min(Fmsy, maxF)

    return(Fmsy)
  }

  # calculate Fmsy for each species/stock
  Fmsy_sp <-  apply(  X = refCurves$Yeq_spf, FUN = .getFmsy,
                      MARGIN = c(1,2) )
  # Now create a matrix to hold the species/stock values
  # on each curve
  spMat       <- matrix(0, nrow = nS, ncol = nP)
  dimnames(spMat) <- dimnames(Fmsy_sp)

  FmsyRefPts  <- list(  yprFmsy_sp    = spMat,
                        YeqFmsy_sp    = spMat,
                        EYeqFmsy_sp   = spMat,
                        BeqFmsy_sp    = spMat,
                        ReqFmsy_sp    = spMat,
                        totBeqFmsy_sp = spMat,
                        Umsy_sp       = spMat )


  # Calculate ref points
  FmsyRefPts$Fmsy_sp    <- Fmsy_sp
  # # Reset lastB_p
  nothing <- .calcEquil( f = 0, obj = obj )


  # Loop and get each stock/species ref pt
  for( s in 1:nS )
    for( p in 1:nP )
    {
      tmp <- .calcEquil( f = Fmsy_sp[s,p], obj = obj )
      FmsyRefPts$yprFmsy_sp[s,p]    <- tmp$ypr_sp[s,p]
      FmsyRefPts$YeqFmsy_sp[s,p]    <- tmp$Yeq_sp[s,p]
      FmsyRefPts$EYeqFmsy_sp[s,p]   <- tmp$EYeq_sp[s,p]
      FmsyRefPts$BeqFmsy_sp[s,p]    <- tmp$Beq_sp[s,p]
      FmsyRefPts$totBeqFmsy_sp[s,p] <- tmp$totBeq_sp[s,p]
      FmsyRefPts$NeqFmsy_sp[s,p]    <- tmp$Neq_sp[s,p]
      FmsyRefPts$ReqFmsy_sp[s,p]    <- tmp$Req_sp[s,p]
      FmsyRefPts$Umsy_sp[s,p]       <- tmp$Ueq_sp[s,p]
    }

  FmsyRefPts
}     # END function .getFmsy


# .calcGenTime
# Purpose:     Calculate generation time as average age of mature stock
# Parameters:  natural mortality and maturity
# Returns:     generation time
# Source:      S.P. Cox
.calcGenTime <- function( obj )
{
  # Compute eqbm spawning biomass per recruit for
  # given f and species/stock pars
  nS    <- obj$om$nS
  nP    <- obj$om$nP
  A_s   <- obj$om$A_s
  nA    <- obj$om$nA
  nL    <- obj$om$nL
  nX    <- obj$om$nX
  M_xsp <- obj$om$M_xsp

  # Life history schedules
  matAge_asp        <- obj$om$matAge_asp

  genTime_sp <- array(0, dim = c(nS,nP))

  for(s in 1:nS )
    for( p in 1:nP )
    {
      surv <- rep(1, length.out = A_s[s])
      a <- c(1:A_s[s])
      surv[a] <- exp( -M_xsp[nX,s,p] * (a - 1) )
      surv[A_s[s]] <- surv[A_s[s]] / (1 - exp(-M_xsp[nX,s,p]))
      genTime_sp[s,p] <- sum( a * matAge_asp[a,s,p] * surv) / sum( surv * matAge_asp[a,s,p])
    }

  genTime_sp


}
