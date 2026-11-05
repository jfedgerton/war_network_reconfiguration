###############################################
# Code to prepare initial data sets for:      #
# "A Modern Peace? Schumpeter, the Decline of #    
# Conﬂict, and the Investment-War Tradeoff"   #
#                                             #
# J. Tyson Chatagnier and Emanuele Castelli   #
#                                             #
# Forthcoming in Political Research Quarterly #
#                                             #
# File created 29 July 2016                   #
###############################################

rm(list = ls())

# port: data-build step for P8_049 (authors' modern_peace_prepare_data.r, reshape->reshape2 edit from the sweep).
# Reads raw author inputs from SWEEP_DATA/work_p8/P8_049; writes the built .dta files to SWEEP_OUT/work_p8/P8_049.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
WD <- file.path(ROOT, "work_p8/P8_049"); WO <- file.path(OUT, "work_p8/P8_049")
dir.create(WO, recursive = TRUE, showWarnings = FALSE)

library(foreign)
library(reshape2)  # P8 edit: reshape not installed; reshape2::melt equivalent
library(countrycode)

# setwd("path")
                        # Use to change default directory

source(file.path(CODE, "work_p8/P8_049/modern_peace_auxiliary_functions.r"))  # port

################
# MONADIC DATA #
################

#################################
# Base data sets being employed #
#################################

gdp <- read.csv(file.path(WD, "maddison_ad1_2010.csv"),
                check.names = FALSE
                )
                        # Maddison data set gives historical GDPpc for
                        # countries in constant 1990 USD (GK$), from 
                        # AD 1 through 2010

polity <- read.csv(file.path(WD, "p4v2014.csv"),
                   check.names = FALSE
                   )
                        # Polity IV data, which run from 1800--2015

m <- read.csv(file.path(WD, "m_scores.csv"))
                        # Arena's M Score data, 1821--2007

cinc <- read.csv(file.path(WD, "nmc_v4.csv"))
                        # CINC Scores, 1816--2007

mid <- read.csv(file.path(WD, "MIDB_401.csv"))
                        # Disputant-level MID data, 1816--2010

trade <- read.csv(file.path(WD, "national_trade_3_0.csv"))
                        # COW national trade data, 1870--2009

mp <- read.csv(file.path(WD, "Majors.csv"))
                        # Info on when states became major powers

nukes <- read.csv(file.path(WD, "nukes.csv"))
                        # Info on when states successfully tested nuclear
                        # weapons; includes Israel from 1961

industry <- read.csv(file.path(WD, "industry_value_added.csv"),
                     check.names = FALSE
                     )

service <- read.csv(file.path(WD, "service_value_added.csv"),
                    check.names = FALSE
                    )
                        # World Bank data on % value added to GDP for
                        # industrial and service sectors, 1960--2014

gartzke <- read.dta(file.path(WD, "IMFcapdata.dta"))
                        # Replication data for Gartzke's capital
                        # openness variable, 1977--2000

mcdonald <- read.dta(file.path(WD, "mcdonald_replication.dta"))
                        # Replication for McDonald's public property
                        # variable, 1969--2001

mousseau <- read.dta(file.path(WD, "Mousseau_Contracts.dta"))
                        # Replication for Mousseau's contract 
                        # intensiveness variable, 1960--2005

################################
# Prepare the data for merging #
################################

#######
# GDP #
#######

gdp <- melt(gdp,
            id.vars = "year"
            )

gdp[, 2] <- countrycode(gdp[, 2],
                        "country.name",
                        "cown"
                        )

colnames(gdp) <- c("year",
                   "ccode",
                   "gdp"
                   )

########################
# Industry value added #
########################

industry <- melt(industry,
                 id.vars = "Country Name"
                 )

industry[, 1] <- countrycode(industry[, 1],
                             "country.name",
                             "cown"
                             )

colnames(industry) <- c("ccode",
                        "year",
                        "industry"
                        )

industry$year <- as.numeric(as.character(industry$year))
                        # The year variable is a factor here,
                        # so make sure it's numeric.

#######################
# Service value added #
#######################

service <- melt(service,
                id.vars = "Country Name"
                )

service[, 1] <- countrycode(service[, 1],
                            "country.name",
                            "cown"
                            )

colnames(service) <- c("ccode",
                       "year",
                       "service"
                       )

service$year <- as.numeric(as.character(service$year))
                        # The year variable is a factor here,
                        # so make sure it's numeric.

##########
# Polity #
##########

polity <- polity[, c("ccode", 
                     "year", 
                     "polity2"
                     )
                 ]
                        # Polity data set has a lot of variables.
                        # Clean it up a bit.

#########
# Trade #
#########

trade[trade == -9] <- NA
                        # Get rid of NAs

trade <- data.frame(ccode = trade$ccode,
                    year  = trade$year,
                    trade = trade$imports + trade$exports
                    )
                        # We just want total trade

########
# CINC #
########

cinc[cinc == -9] <- NA
                        # Get rid of NAs

cinc <- cinc[order(cinc$year, 
                   cinc$ccode
                   ), 
             ]
                        # Reorder by year, and then by ccode,
                        # which will be necessary for the 
                        # by() arguments below
irst.vals <- by(cinc$irst,
                list(cinc$year),
                function(x) {
                  x / sum(x, na.rm= T)
                }
                )

pec.vals <- by(cinc$pec,
               list(cinc$year),
               function(x) {
                 x / sum(x, na.rm= T)
               }
               )
                        # Get the proportion of world's iron and
                        # steel production and energy consumption
                        # for each country in each year; this is the
                        # value used when computing CINC scores

cinc <- data.frame(year   = cinc$year,
                   ccode  = cinc$ccode,
                   irst   = unlist(irst.vals),
                   pec    = unlist(pec.vals),
                   urbpct = cinc$upop / cinc$tpop
                   )

########
# MIDs #
########

mid[mid == -9] <- NA
                        # Get rid of NAs

################
# Gartzke data #
################

ccode.gartzke <- unique(c(gartzke$statea,
                          gartzke$stateb
                          )
                        )

year.gartzke <- unique(gartzke$year)

gartzke.grid <- expand.grid(ccode.gartzke,
                            year.gartzke
                            )

open <- vector(length = nrow(gartzke.grid))

for (i in 1 : nrow(gartzke.grid)) {
  rows.val <- which((gartzke$statea == gartzke.grid[i, 1] |
                     gartzke$stateb == gartzke.grid[i, 1]
                     ) &
                    gartzke$year == gartzke.grid[i, 2]
                    )
  cap.vals <- c(gartzke$capopenh[rows.val], gartzke$capopenl[rows.val])
  open[i] <- ifelse(length(table(cap.vals)) > 0,
                    as.numeric(names(which.max(table(cap.vals)))),
                    NA
                    )
                        # Get all values of capital openness associated with
                        # a given country-year, and sort by frequency. The 
                        # most frequent should be the state's true value.
                        # If there are no values, this is NA.
}

gartzke <- data.frame(ccode = gartzke.grid[, 1],
                      year  = gartzke.grid[, 2],
                      open  = open
                      )

#################
# McDonald data #
#################

mcdonald <- mcdonald[, c("ccode",
                         "year",
                         "public"
                         )
                     ]
                        # Drop other variables

#################
# Mousseau data #
#################

mousseau <- data.frame(ccode   = mousseau$ccode,
                       year    = mousseau$year,
                       cie     = mousseau$CIE,
                       cie.bin = mousseau$CIEbinary
                       )
                        # Keep both the continuous and binary
                        # versions of the CIE variable

###################
# Create data set #
###################

monadic <- Reduce(function(...) {
                  merge(..., 
                        all = TRUE
                        )
                }, 
                list(gdp, 
                     polity, 
                     m, 
                     cinc,
                     trade, 
                     industry, 
                     service,
                     gartzke,
                     mcdonald,
                     mousseau
                     )
                )
                        # Merge the country-year data together, using function
                        # taken from Stack Overflow:
                        # http://stackoverflow.com/questions/14096814/r-merging-a-lot-of-data-frames

monadic$year <- as.numeric(monadic$year)

monadic$majpow <- 0
monadic$nuclear <- 0

for (i in 1 : nrow(mp)) {
  mp.obs <- which(monadic$ccode == mp$ccode[i] &
                  monadic$year >= mp$styear[i] &
                  monadic$year < mp$endyear[i]
                  )
                        # Get the relevant observations, which include
                        # the year a state becomes a major power and 
                        # do not include the year it loses major power
                        # status
  monadic$majpow[mp.obs] <- 1
}

for (i in 1 : nrow(nukes)) {
  nukes.obs <- which(monadic$ccode == nukes$ccode[i] &
                     monadic$year >= nukes$styear[i]
                     )
                        # Get the relevant observations, which include
                        # the year a state conducts a nuclear test
  monadic$nuclear[nukes.obs] <- 1
}

#############
# MIDs Data #
#############

monadic.mids <- data.frame(start   = mid$StYear,
                           end     = mid$EndYear,
                           ccode   = mid$ccode,
                           fatal   = mid$Fatality,
                           hihost  = mid$HostLev
                           )

mids.one <- monadic.mids[monadic.mids$start == monadic.mids$end, ]
                        # Pull out those conflicts that lasted less than
                        # a year.

mids.one$init.mid <- 1
mids.one$ongo.mid <- 0
mids.one$init.hostile <- 1 * (mids.one$hihost >= 4)
mids.one$ongo.hostile <- 0
mids.one$init.fatal <- 1 * (mids.one$fatal >= 1)
mids.one$ongo.fatal <- 0

                        # Create an initiation variable, which is necessarily
                        # one for all of these observations, and an ongoing
                        # variable, which is necessarily zero.

mids.one$year <- mids.one$start

mids.one <- mids.one[, !names(mids.one) %in% c("start", 
                                               "end"
                                               )
                     ]
                        # Replace the start and end variables with a 
                        # year variable

mids.mult <- monadic.mids[monadic.mids$start < monadic.mids$end, ]
mult.list <- vector("list",
                    length = nrow(mids.mult)
                    )
                        # Pull out those conflicts that lasted more than
                        # one year, and initialize a list to turn them
                        # into matrices.

for (i in 1 : nrow(mids.mult)) {
  start <- mids.mult[i, "start"]
  end <- mids.mult[i, "end"]
  mids.sub <- mids.mult[i, ]
  mid.mat <- mids.sub[rep(row.names(mids.sub), 
                          end - start + 1
                          ), 
                      ]
                        # Create data frame with rows for each year that the
                        # conflict lasted.
  mid.mat <- cbind(mid.mat,
                   c(1,
                     rep(0,
                         times = end - start
                         )
                     ),
                   c(0, 
                     rep(1,
                         times = end - start
                         )
                     ),
                   seq(start,
                       end
                       )
                   )
                        # Combine the data with initiation and ongoing
                        # vectors, as well as a year vector.
  colnames(mid.mat) <- c(names(mids.mult), 
                         "init.mid",
                         "ongo.mid",
                         "year"
                         )
  mid.mat <- mid.mat[, !colnames(mid.mat) %in% c("start", 
                                                 "end"
                                                 )
                     ]
  mid.mat$init.hostile <- mid.mat$init.mid * (mid.mat$hihost >= 4)
  mid.mat$ongo.hostile <- mid.mat$ongo.mid * (mid.mat$hihost >= 4)
  mid.mat$init.fatal <- mid.mat$init.mid * (mid.mat$fatal >= 1)
  mid.mat$ongo.fatal <- mid.mat$ongo.mid * (mid.mat$fatal >= 1)
  mult.list[[i]] <- mid.mat
}

monadic.mids <- rbind(mids.one,
                      do.call("rbind",
                              mult.list
                              )
                      )

monadic.mids <- monadic.mids[c("ccode",
                               "year",
                               "init.mid",
                               "init.hostile",
                               "init.fatal",
                               "ongo.mid",
                               "ongo.hostile",
                               "ongo.fatal"
                               )
                             ]

monadic.id <- paste(monadic.mids$ccode,
                    monadic.mids$year,
                    sep = "."
                    )

mids.combine <- by(monadic.mids,
                   list(monadic.id),
                   function(x) {
                     apply(x,
                           2,
                           max,
                           na.rm = TRUE
                           )
                   }
                   )
                        # We will have some overlap (i.e., states with
                        # multiple MIDs initiated and/or ongoing for a
                        # given year), so take the largest values here.

mids.combine <- do.call("rbind",
                        mids.combine
                        )

mids.combine[is.infinite(mids.combine)] <- NA
                        # Some of our fatal values were NAs and the maxima
                        # became -Inf above. So turn them back into NAs.

mids.combine <- data.frame(mids.combine)
                        # Turn this into a data frame to make it easier 
                        # to work with.

cy <- paste(monadic$ccode,
            monadic$year,
            sep = "."
            )

cy.combine <- paste(mids.combine$ccode,
                    mids.combine$year,
                    sep = "."
                    )

monadic$mid.onset <- ifelse(monadic$year < 1816 | 
                            monadic$year > 2010,
                            NA,
                            ifelse(cy %in% unique(cy.combine[mids.combine$init.mid == 1]),
                                   1,
                                   0
                                   )
                            )

monadic$hostile.onset <- ifelse(monadic$year < 1816 | 
                                monadic$year > 2010,
                                NA,
                                ifelse(cy %in% unique(cy.combine[mids.combine$init.hostile == 1]),
                                       1,
                                       0
                                       )
                                )

monadic$fatal.onset <- ifelse(monadic$year < 1816 | 
                              monadic$year > 2010,
                              NA,
                              ifelse(cy %in% unique(cy.combine[mids.combine$init.fatal == 1]),
                                     1,
                                     0
                                     )
                              )

monadic$mid.ongo <- ifelse(monadic$year < 1816 | 
                            monadic$year > 2010,
                            NA,
                            ifelse(cy %in% unique(cy.combine[mids.combine$ongo.mid == 1]),
                                   1,
                                   0
                                   )
                            )

monadic$hostile.ongo <- ifelse(monadic$year < 1816 | 
                            monadic$year > 2010,
                            NA,
                            ifelse(cy %in% unique(cy.combine[mids.combine$ongo.hostile == 1]),
                                   1,
                                   0
                                   )
                            )

monadic$fatal.ongo <- ifelse(monadic$year < 1816 | 
                            monadic$year > 2010,
                            NA,
                            ifelse(cy %in% unique(cy.combine[mids.combine$ongo.fatal == 1]),
                                   1,
                                   0
                                   )
                            )

######################################
# Calculate time since last conflict #
######################################

t.mid <- GetTime(monadic$mid.onset,
                 monadic$mid.ongo,
                 monadic$ccode,
                 monadic$year
                 )

t.hostile <- GetTime(monadic$hostile.onset,
                     monadic$hostile.ongo,
                     monadic$ccode,
                     monadic$year
                     )

t.fatal <- GetTime(monadic$fatal.onset,
                   monadic$fatal.ongo,
                   monadic$ccode,
                   monadic$year
                   )
                        # Create measure of time since the last 
                        # MID of each type

t.df <- data.frame(year      = t.mid[, 2],
                   ccode     = t.mid[, 1],
                   t.mid     = t.mid[, 3],
                   t.hostile = t.hostile[, 3],
                   t.fatal   = t.fatal[, 3]
                   )
                        # Put these together in a data frame

monadic <- merge(monadic,
                 t.df,
                 by = c("ccode",
                        "year"
                        )
                 )
                        # Merge into monadic data set

###############################
# Write data set for analysis #
###############################

monadic.subset <- monadic[monadic$year >= 1960 & monadic$year <= 2007, ]
                        # The data sets that we use are limited to this
                        # range of years, so subset the file in order
                        # to remove necessarily incomplete observations

write.dta(monadic.subset,
          file.path(WO, "monadic_data_set.dta")
          )

write.dta(monadic,
          file.path(WO, "monadic_data_full.dta")
          )

###############
# DYADIC DATA #
###############

trade <- read.csv(file.path(WD, "dyadic_trade_3_0.csv"))
                        # COW bilateral trade data, 1870--2009

alliance <- read.csv(file.path(WD, "alliance_v4_1_by_dyad_yearly.csv"))
                        # COW formal alliance data, 1816--2012

affinity <- read.dta(file.path(WD, "affinity_03102006.dta"))
                        # Gartzke's affinity data, 1946--2002

contig <- read.csv(file.path(WD, "contdir.csv"))
                        # COW Direct Contiguity data, 1816--2006

################################
# Prepare the data for merging #
################################

#########
# Trade #
#########

trade[trade == -9] <- NA
                        # Get rid of NAs

trade <- data.frame(ccode1  = trade$ccode1,
                    ccode2  = trade$ccode2,
                    year    = trade$year,
                    trade   = trade$flow1 + trade$flow2,
                    lntrade = log(1 + trade$flow1 + trade$flow2)
                    )

############
# Affinity #
############

names(affinity)[4 : 5] <- c("ccode1",
                            "ccode2"
                            )
                        # Rename the country code columns so that
                        # it will be compatible with other data sets

affinity <- affinity[, -c(1 : 3)]
                        # Drop unnecessary columns

#############
# Alliances #
#############

allies.error <- which(alliance$year == 0)
n.egypt <- nrow(alliance[alliance$year == 0 & alliance$ccode1 == 651, ])
alliance$year[allies.error] <- alliance$dyad_st_year[allies.error]
alliance$dyad_st_year[allies.error] <- rep(c(1979,
                                             1981
                                             ),
                                           times = c(n.egypt,
                                                     length(allies.error) - 
                                                       n.egypt
                                                     )
                                           )
                        # Correct an error with the coding of Israeli
                        # alliances with the US and Egypt

alliance$id <- paste(alliance$ccode1,
                     alliance$ccode2,
                     alliance$year,
                     sep = "."
                     )

alliance.agg <- aggregate(alliance[,14 : 17], 
                          list(alliance$ccode1, 
                               alliance$ccode2, 
                               alliance$year
                               ), 
                          sum
                          )

names(alliance.agg)[1 : 3] <- c("ccode1",
                                "ccode2",
                                "year"
                                )

alliance.agg[, 4 : 7] <- apply(alliance.agg[, 4 : 7],
                               2,
                               function(x) {
                                 ifelse(x >= 1,
                                        1,
                                        0
                                        )
                               }
                               )
                        # Data are in alliance-year format. Aggregate up to
                        # dyad-year format, and use dummies for categories

alliance.agg$alliance <- 1
                        # These are all instances in which an alliance of some
                        # kind existed, so create a variable to note that

alliance.agg[is.na(alliance.agg)] <- -99
                        # Recode this to -99 for now, and then code
                        # back to NA later.

##############
# Contiguity #
##############

contig$contig <- ifelse(contig$conttype <= 3,
                        1,
                        0
                        )

contig$start <- as.numeric(substr(contig$begin,
                                  1,
                                  4
                                  )
                           )

contig$end <- as.numeric(substr(contig$end,
                                1,
                                4
                                )
                         )
                        # Get start and end years

contig <- data.frame(ccode1 = contig$statelno,
                     ccode2 = contig$statehno,
                     start  = as.numeric(substr(contig$begin,
                                                1,
                                                4
                                                )
                                         ),
                     end    = as.numeric(substr(contig$end,
                                                1,
                                                4
                                                )
                                         ),
                     contig = ifelse(contig$conttype <= 3,
                                     1,
                                     0
                                     )
                     )
                        # Create a data set with just the variables we want,
                        # recode the start and end variables to years, and 
                        # get contiguity for separation by less than 25 miles

contig.one <- contig[contig$start == contig$end, ]

contig.one$year <- contig.one$start

contig.one <- contig.one[, !names(contig.one) %in% c("start", 
                                                     "end"
                                                     )
                         ]
                        # Replace the start and end variables with a 
                        # year variable

contig.mult <- contig[contig$start < contig$end, ]

mult.list <- vector("list",
                    length = nrow(contig.mult)
                    )
                        # Pull out those pairs that lasted more than
                        # one year, and initialize a list to turn them
                        # into matrices.

for (i in 1 : nrow(contig.mult)) {
  start <- contig.mult[i, "start"]
  end <- contig.mult[i, "end"]
  contig.sub <- contig.mult[i, ]
  contig.mat <- contig.sub[rep(row.names(contig.sub), 
                        end - start + 1
                        ), 
                        ]
                        # Create data frame with rows for each year that the
                        # conflict lasted.

  contig.mat <- cbind(contig.mat,
                      seq(start,
                          end
                          )
                   )
                        # Combine the data with  a year vector.
  colnames(contig.mat) <- c(names(contig.mult), 
                            "year"
                            )
  contig.mat <- contig.mat[, !colnames(contig.mat) %in% c("start", 
                                                          "end"
                                                          )
                     ]
  mult.list[[i]] <- contig.mat
}

dyadic.contig <- rbind(contig.one,
                       do.call("rbind",
                               mult.list
                               )
                       )

##########################
# Create dyadic MID data #
##########################

unique.mids <- unique(mid$DispNum3)

dyads.list <- vector("list",
                     length = length(unique.mids)
                     )

for (i in 1 : length(dyads.list)) {
  mid.sub <- mid[mid$DispNum3 == unique.mids[i], ]
                        # Look at one dispute at a time
  mid.pairs <- expand.grid(mid.sub$ccode[mid.sub$SideA == 1],
                           mid.sub$ccode[mid.sub$SideA == 0]
                           )
                        # Get all pairs of Side A vs. Side B
  names(mid.pairs) <- c("ccode1",
                        "ccode2"
                        )
  for (j in 1 : nrow(mid.pairs)) {
    .pair <- c(mid.pairs$ccode1[j],
               mid.pairs$ccode2[j]
               )
    mid.pairs$start[j] <- max(mid.sub$StYear[mid.sub$ccode %in% .pair])
    mid.pairs$end[j] <- min(mid.sub$EndYear[mid.sub$ccode %in% .pair])
                        # Get the latest start date and earliest end date
                        # for the members of each pair of combatants
    mid.pairs$hostile[j] <- min(mid.sub$HostLev[mid.sub$ccode %in% .pair])
    mid.pairs$fatal[j] <- min(mid.sub$Fatality[mid.sub$ccode %in% .pair])
                        # Get the lowest hostility level among the pair,
                        # and the lowest fatality level among the pair
    mid.pairs$origin[j] <- prod(mid.sub$Orig[mid.sub$ccode %in% .pair])
                        # Were both of the actors originators?
  }
  dyads.list[[i]] <- mid.pairs
}

dyadic.mids <- do.call("rbind",
                       dyads.list
                       )
                        # Create a data frame the list of data frames we have
                        # above. This will give us all of our conflict dyads,
                        # with the relevant information.

dyadic.mids <- dyadic.mids[dyadic.mids$start <= dyadic.mids$end, ]
                        # Remove those instances in which one of the pair
                        # exited the conflict before the other entered.

mids.one <- dyadic.mids[dyadic.mids$start == dyadic.mids$end, ]
                        # Pull out those conflicts that lasted less than
                        # a year.

mids.one$mid.onset <- 1
mids.one$mid.ongo <- 0
mids.one$hostile.onset <- 1 * (mids.one$hostile >= 4)
mids.one$hostile.ongo <- 0
mids.one$fatal.onset <- 1 * (mids.one$fatal >= 1)
mids.one$fatal.ongo <- 0
                        # Create an initiation variable, which is necessarily
                        # one for all of these observations, and an ongoing
                        # variable, which is necessarily zero.

mids.one$year <- mids.one$start

mids.one <- mids.one[, !names(mids.one) %in% c("start", 
                                               "end",
                                               "fatal",
                                               "hostile"
                                               )
                     ]
                        # Remove columns that are no longer needed.

mids.mult <- dyadic.mids[dyadic.mids$start < dyadic.mids$end, ]
mult.list <- vector("list",
                    length = nrow(mids.mult)
                    )
                        # Pull out those conflicts that lasted more than
                        # one year, and initialize a list to turn them
                        # into matrices.

for (i in 1 : nrow(mids.mult)) {
  start <- mids.mult[i, "start"]
  end <- mids.mult[i, "end"]
  mids.sub <- mids.mult[i, ]
  mid.mat <- mids.sub[rep(row.names(mids.sub), 
                          end - start + 1
                          ), 
                      ]
                        # Create data frame with rows for each year that the
                        # conflict lasted.
  mid.mat <- cbind(mid.mat,
                   c(1,
                     rep(0,
                         times = end - start
                         )
                     ),
                   c(0, 
                     rep(1,
                         times = end - start
                         )
                     )
                   )
                        # Combine the data with initiation and ongoing
                        # vectors.
  colnames(mid.mat) <- c(names(mids.mult), 
                         "mid.onset",
                         "mid.ongo"
                         )
  mid.mat$hostile.onset <- mid.mat$mid.onset * (mid.mat$hostile >= 4)
  mid.mat$hostile.ongo <- mid.mat$mid.ongo * (mid.mat$hostile >= 4)
  mid.mat$fatal.onset <- mid.mat$mid.onset * (mid.mat$fatal >= 1)
  mid.mat$fatal.ongo <- mid.mat$mid.ongo * (mid.mat$fatal >= 1)
                        # Create columns for hostile and fatal MIDs
  mid.mat <- mid.mat[, !colnames(mid.mat) %in% c("start", 
                                                 "end",
                                                 "fatal",
                                                 "hostile"
                                                 )
                     ]
                        # Remove columns that are no longer needed.
  mid.mat$year <- seq(start, 
                      end
                      )
                        # Create a column for the year variable
  mult.list[[i]] <- mid.mat
}

dyadic.mids <- rbind(mids.one,
                     do.call("rbind",
                             mult.list
                             )
                     )

dyadic.mids <- dyadic.mids[c("ccode1",
                             "ccode2",
                             "year",
                             names(dyadic.mids)[3 : (ncol(dyadic.mids) - 1)]
                             )
                           ]

dyadic.mids[is.na(dyadic.mids)] <- -99
                        # Recode this to -99 for now, and then code
                        # back to NA later.

################
# Combine data #
################

dyadic <- merge(trade,
                affinity,
                by  = c("ccode1",
                        "ccode2",
                        "year"
                        ),
                all = TRUE
                )

names.ally <- names(alliance.agg[, 4 : ncol(alliance.agg)])
names.mids <- names(dyadic.mids[, 4 : ncol(dyadic.mids)])

dyads.id <- paste(dyadic$ccode1,
                  dyadic$ccode2,
                  dyadic$year,
                  sep = "."
                  )

allies.id <- paste(alliance.agg$ccode1,
                   alliance.agg$ccode2,
                   alliance.agg$year,
                   sep = "."
                   )

mids.id <- paste(dyadic.mids$ccode1,
                 dyadic.mids$ccode2,
                 dyadic.mids$year,
                 sep = "."
                 )

contig.id <- paste(dyadic.contig$ccode1,
                   dyadic.contig$ccode2,
                   dyadic.contig$year,
                   sep = "."
                   )

allies.obs <- ifelse(dyads.id %in% allies.id,
                     match(dyads.id,
                           allies.id
                           ),
                     NA
                     )

mids.obs <- ifelse(dyads.id %in% mids.id,
                   match(dyads.id,
                         mids.id
                         ),
                   NA
                   )

contig.obs <- ifelse(dyads.id %in% contig.id,
                     match(dyads.id,
                           contig.id
                           ),
                     NA
                     )

dyadic[, names.ally] <- alliance.agg[allies.obs, names.ally]
dyadic[, names.mids] <- dyadic.mids[mids.obs, names.mids]
dyadic[, "contig"] <- dyadic.contig[contig.obs, "contig"]

dyadic[, c(names.ally,
           names.mids,
           "contig"
           )
       ][is.na(dyadic[, c(names.ally,
                          names.mids,
                          "contig"
                          )
                      ]
               )
         ] <- 0
                        # Those observations that should be zero
                        # will be NA, so recode them

dyadic[dyadic == -99] <- NA
                        # Recode true NAs to NA

##################################
# Add in monadic variables that  #
# depend on dyad characteristics #
##################################

dyads.sep <- cbind(paste(dyadic$ccode1,
                         dyadic$year,
                         sep = "."
                         ),
                   paste(dyadic$ccode2,
                         dyadic$year,
                         sep = "."
                         )
                   )

###############
# CINC values #
###############

cinc.id <- paste(cinc$ccode,
                 cinc$year,
                 sep = "."
                 )

irst.vals <- cbind(cinc$irst[match(dyads.sep[, 1],
                                   cinc.id
                                   )
                             ],
                   cinc$irst[match(dyads.sep[, 2],
                                   cinc.id
                                   )
                             ]
                   )

dyadic$irst.hi <- pmax(irst.vals[, 1], 
                       irst.vals[, 2]
                       )

pec.vals <- cbind(cinc$pec[match(dyads.sep[, 1],
                                  cinc.id
                                  )
                             ],
                  cinc$pec[match(dyads.sep[, 2],
                                 cinc.id
                                 )
                           ]
                  )

dyadic$pec.hi <- pmax(pec.vals[, 1], 
                       pec.vals[, 2]
                       )
                        # Get the higher value of iron and steel production
                        # and energy consumption for robustness check

#####################
# Industrialization #
#####################

ind.id <- paste(industry$ccode,
                industry$year,
                sep = "."
                )

ind.vals <- cbind(industry$industry[match(dyads.sep[, 1],
                                          ind.id
                                          )
                                    ],
                  industry$industry[match(dyads.sep[, 2],
                                          ind.id
                                          )
                                    ]
                  )

dyadic$ind.hi <- pmax(ind.vals[, 1], 
                      ind.vals[, 2]
                      )

##################
# Service sector #
##################

serv.id <- paste(service$ccode,
                 service$year,
                 sep = "."
                 )

serv.vals <- cbind(service$service[match(dyads.sep[, 1],
                                         serv.id
                                         )
                                   ],
                  service$service[match(dyads.sep[, 2],
                                        serv.id
                                        )
                                  ]
                  )

dyadic$serv.ind.hi <- ifelse(ind.vals[, 1] > ind.vals[, 2],
                             serv.vals[, 1],
                             serv.vals[, 2]
                             )
                        # Get the size of the service sector
                        # for the more industrialized country

################################
# Reverse industry and service #
################################

dyadic$serv.hi <- pmax(serv.vals[, 1], 
                       serv.vals[, 2]
                       )

dyadic$ind.serv.hi <- ifelse(serv.vals[, 1] > serv.vals[, 2],
                             ind.vals[, 1],
                             ind.vals[, 2]
                             )
                        # For robustness, get the size of the
                        # service sector for the more service-heavy
                        # country, and the size of its industrial 
                        # sector


#######
# GDP #
#######

gdp.id <- paste(gdp$ccode,
                gdp$year,
                sep = "."
                )

gdp.vals <- cbind(gdp$gdp[match(dyads.sep[, 1],
                                gdp.id
                                )
                          ],
                  gdp$gdp[match(dyads.sep[, 2],
                                gdp.id
                                )
                          ]
                  )

dyadic$gdp.lo <- pmin(gdp.vals[, 1], 
                      gdp.vals[, 2]
                      )

################
# Capabilities #
################

m.id <- paste(m$ccode,
              m$year,
              sep = "."
              )

M.vals <- cbind(m$m[match(dyads.sep[, 1],
                          m.id
                          )
                    ],
                m$m[match(dyads.sep[, 2],
                          m.id
                          )
                    ]
                )

m.large <- pmax(M.vals[, 1],
                M.vals[, 2]
                )

dyadic$Mrat <- m.large / (rowSums(M.vals))

######################
# Major Power Status #
######################

dyadic$majpow <- 0

for (i in 1 : nrow(mp)) {
  mp.obs <- which((dyadic$ccode1 == mp$ccode[i] |
                   dyadic$ccode2 == mp$ccode[i]
                   ) &
                  dyadic$year >= mp$styear[i] &
                  dyadic$year < mp$endyear[i]
                  )
                        # Get the relevant observations, which include
                        # the year a state becomes a major power and 
                        # do not include the year it loses major power
                        # status
  dyadic$majpow[mp.obs] <- 1
}

###################
# Nuclear Weapons #
###################

dyadic$nuclear <- 0

for (i in 1 : nrow(nukes)) {
  nuclear.obs <- which((dyadic$ccode1 == mp$ccode[i] |
                        dyadic$ccode2 == mp$ccode[i]
                        ) &
                       dyadic$year >= mp$styear[i] &
                       dyadic$year < mp$endyear[i]
                       )
                        # Get the relevant observations, which include
                        # the year a state becomes a major power and 
                        # do not include the year it loses major power
                        # status
  dyadic$nuclear[nuclear.obs] <- 1
}

#############
# Democracy #
#############

polity.id <- paste(polity$ccode,
                   polity$year,
                   sep = "."
                   )

polity.vals <- cbind(polity$polity[match(dyads.sep[, 1],
                                   polity.id
                                   )
                             ],
                     polity$polity[match(dyads.sep[, 2],
                                   polity.id
                                   )
                             ]
                     )

dyadic$lodem <- pmin(polity.vals[, 1], 
                     polity.vals[, 2]
                     )

dyadic$hidem <- pmax(polity.vals[, 1], 
                     polity.vals[, 2]
                     )

###########
# Gartzke #
###########

gartzke.id <- paste(gartzke$ccode,
                    gartzke$year,
                    sep = "."
                    )

gartzke.vals <- cbind(gartzke$open[match(dyads.sep[, 1],
                                         gartzke.id
                                         )
                                   ],
                      gartzke$open[match(dyads.sep[, 2],
                                         gartzke.id
                                         )
                                   ]
                      )

dyadic$openl <- pmin(gartzke.vals[, 1], 
                     gartzke.vals[, 2]
                     )

############
# McDonald #
############

mcdonald.id <- paste(mcdonald$ccode,
                     mcdonald$year,
                     sep = "."
                     )

mcdonald.vals <- cbind(mcdonald$public[match(dyads.sep[, 1],
                                             mcdonald.id
                                             )
                                       ],
                       mcdonald$public[match(dyads.sep[, 2],
                                             mcdonald.id
                                             )
                                       ]
                       )

dyadic$publich <- pmax(mcdonald.vals[, 1], 
                       mcdonald.vals[, 2]
                       )

############
# Mousseau #
############

mousseau.id <- paste(mousseau$ccode,
                     mousseau$year,
                     sep = "."
                     )

cie.vals <- cbind(mousseau$cie[match(dyads.sep[, 1],
                                     mousseau.id
                                     )
                               ],
                  mousseau$cie[match(dyads.sep[, 2],
                                     mousseau.id
                                     )
                               ]
                  )

ciebin.vals <- cbind(mousseau$cie.bin[match(dyads.sep[, 1],
                                            mousseau.id
                                            )
                                      ],
                     mousseau$cie.bin[match(dyads.sep[, 2],
                                            mousseau.id
                                            )
                                      ]
                     )

dyadic$cielo <- pmin(cie.vals[, 1], 
                     cie.vals[, 2]
                     )
dyadic$cie.binlo <- pmin(ciebin.vals[, 1], 
                         ciebin.vals[, 2]
                         )

######################################
# Calculate time since last conflict #
######################################

dyadic$dyadid <- paste(dyadic$ccode1,
                       dyadic$ccode2,
                       sep = "."
                       )

t.mid <- GetTime(dyadic$mid.onset,
                 dyadic$mid.ongo,
                 dyadic$dyadid,
                 dyadic$year
                 )

t.hostile <- GetTime(dyadic$hostile.onset,
                     dyadic$hostile.ongo,
                     dyadic$dyadid,
                     dyadic$year
                     )

t.fatal <- GetTime(dyadic$fatal.onset,
                   dyadic$fatal.ongo,
                   dyadic$dyadid,
                   dyadic$year
                   )
                        # Create measure of time since the last 
                        # MID of each type

t.df <- data.frame(year      = t.mid[, 2],
                   dyadid    = t.mid[, 1],
                   t.mid     = as.numeric(t.mid[, 3]),
                   t.hostile = as.numeric(t.hostile[, 3]),
                   t.fatal   = as.numeric(t.fatal[, 3])
                   )
                        # Put these together in a data frame

dyadic <- merge(dyadic,
                t.df,
                by = c("dyadid",
                       "year"
                       )
                )
                        # Merge into monadic data set

###############################
# Write data set for analysis #
###############################

dyadic.subset <- dyadic[dyadic$year >= 1960 & dyadic$year <= 2002, ]
                        # The data sets that we use are limited to this
                        # range of years, so subset the file in order
                        # to remove necessarily incomplete observations

write.dta(dyadic.subset,
          file.path(WO, "dyadic_data_set.dta")
          )

write.dta(dyadic,
          file.path(WO, "dyadic_data_full.dta")
          )
                          