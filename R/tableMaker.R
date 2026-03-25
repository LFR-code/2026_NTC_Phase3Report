# tableMaker.R
#
# Table-generation functions for the NTC Phase 3 report.
# Returns kableExtra tables formatted for PDF output.
#
# Functions:
#   makeAllocTable()   -- Allocation proportions by Stat Area (slides 14-15)
#   makePerfTables()   -- Performance metrics tables     (slides 20-21)


# makeAllocTable()
# Creates a table of TAC allocation proportions by Stat Area for all three
# allocation methods, at the start of the projection (2024) and as
# long-run medians over the full projection period.
# Inputs:
#   mp1  -- sim3S_MP1 blob (SI allocation)
#   mp2  -- sim3S_MP2 blob (Catch allocation)
#   mp3  -- sim3S_MP3 blob (All-Stat Area-23 allocation)
makeAllocTable <- function( mp1 = sim3S_MP1,
                            mp2 = sim3S_MP2,
                            mp3 = sim3S_MP3 )
{
  stocks <- c("Stat Area 25", "Stat Area 24", "Stat Area 23")

  # Helper: extract median propTAC at tMP and over projection
  extractProp <- function(obj) {
    tMP      <- obj$om$tMP
    nT       <- obj$om$nT
    prop_isp <- obj$mp$hcr$propTAC_ispt[, 1, , ]   # i x nP x nT

    # 2024 allocation (tMP year)
    prop_2024 <- apply(prop_isp[, , tMP], MARGIN = 2, FUN = median)

    # Long-run median (tMP to nT)
    prop_lr   <- apply(prop_isp[, , tMP:nT], MARGIN = 2, FUN = median)

    list(start = prop_2024, longrun = prop_lr)
  }

  si    <- extractProp(mp1)
  catch <- extractProp(mp2)
  alloc <- extractProp(mp3)

  # Format as percentages
  fmt <- function(x) paste0(round(x * 100, 0), "\\%")

  df <- data.frame(
    "Stat Area"  = stocks,
    SI_2024      = fmt(si$start),
    SI_LR        = fmt(si$longrun),
    Catch_2024   = fmt(catch$start),
    Catch_LR     = fmt(catch$longrun),
    All23_2024   = fmt(alloc$start),
    All23_LR     = fmt(alloc$longrun),
    stringsAsFactors = FALSE,
    check.names  = FALSE
  )

  kableExtra::kbl(
    x       = df,
    format  = "latex",
    booktabs = TRUE,
    escape  = FALSE,
    col.names = c(
      "Stat Area",
      "2024", "Long-run",
      "2024", "Long-run",
      "2024", "Long-run"
    ),
    caption = paste0(
      "Median TAC allocation proportions by Stat Area for the three ",
      "allocation methods. Values shown for the first projection ",
      "year (2024) and as long-run medians over the full 25-year ",
      "projection period."
    )
  ) |>
    kableExtra::add_header_above(
      header = c(
        " " = 1,
        "SI allocation" = 2,
        "Catch allocation" = 2,
        "All Stat Area 23" = 2
      )
    ) |>
    kableExtra::kable_styling(
      latex_options = c("hold_position"),
      font_size     = 10
    )
} # END makeAllocTable()


# makePerfTables()
# Reads the pre-computed fullStatTable.csv and returns two kableExtra
# tables: probability of low biomass and probability of low catch,
# each pivoted to Stat Area x allocation method layout.
# Inputs:
#   statFile  -- path to fullStatTable.csv
makePerfTables <- function(
    statFile = here::here(
      'Models/Sims/3S_allocTest/statistics/fullStatTable.csv'
    )
)
{
  df <- read.csv(file = statFile, stringsAsFactors = FALSE)

  # Map mp values to readable labels
  mpMap <- c(
    "twoStep-minE-SI"    = "SI",
    "twoStep-minE-catch" = "Catch",
    "twoStep-minE-area3" = "All Stat Area 23"
  )
  df$mpLab <- mpMap[df$mp]

  # Map stock to readable labels
  stockMap <- c(pfma25 = "Stat Area 25", pfma24 = "Stat Area 24", pfma23 = "Stat Area 23")
  df$stockLab <- stockMap[df$stock]

  # ---- Table 1: probability of low biomass ----
  pBio <- reshape(
    data      = df[, c("stockLab", "mpLab", "pBltMinHistB")],
    idvar     = "stockLab",
    timevar   = "mpLab",
    direction = "wide"
  )
  names(pBio) <- gsub("^pBltMinHistB\\.", "", names(pBio))
  rownames(pBio) <- NULL
  # Reorder columns
  pBio <- pBio[order(match(pBio$stockLab, c("Stat Area 25","Stat Area 24","Stat Area 23"))),
               c("stockLab", "SI", "Catch", "All Stat Area 23")]

  tbl1 <- kableExtra::kbl(
    x        = pBio,
    format   = "latex",
    booktabs = TRUE,
    digits   = 3,
    label    = "perfBiomass",
    col.names = c("Stat Area", "SI", "Catch", "All Stat Area 23"),
    caption  = paste0(
      "Probability of low biomass (biomass below the historical ",
      "minimum) by Stat Area and allocation method under the ",
      "\\emph{SISCAH-3S-Pred} operating model."
    )
  ) |>
    kableExtra::kable_styling(
      latex_options = c("hold_position"),
      font_size     = 10
    )

  # ---- Table 2: probability of low catch ----
  pCat <- reshape(
    data      = df[, c("stockLab", "mpLab", "pCltMinHistC")],
    idvar     = "stockLab",
    timevar   = "mpLab",
    direction = "wide"
  )
  names(pCat) <- gsub("^pCltMinHistC\\.", "", names(pCat))
  rownames(pCat) <- NULL
  pCat <- pCat[order(match(pCat$stockLab, c("Stat Area 25","Stat Area 24","Stat Area 23"))),
               c("stockLab", "SI", "Catch", "All Stat Area 23")]

  tbl2 <- kableExtra::kbl(
    x        = pCat,
    format   = "latex",
    booktabs = TRUE,
    digits   = 3,
    label    = "perfCatch",
    col.names = c("Stat Area", "SI", "Catch", "All Stat Area 23"),
    caption  = paste0(
      "Probability of low catch (catch below the historical ",
      "minimum) by Stat Area and allocation method under the ",
      "\\emph{SISCAH-3S-Pred} operating model."
    )
  ) |>
    kableExtra::kable_styling(
      latex_options = c("hold_position"),
      font_size     = 10
    )

  list(biomass = tbl1, catch = tbl2)
} # END makePerfTables()
