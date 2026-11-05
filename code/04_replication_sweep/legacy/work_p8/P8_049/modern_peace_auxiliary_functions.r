###############################################
# Auxiliary functions used in replication of: #
# "A Modern Peace? Schumpeter, the Decline of #    
# Conﬂict, and the Investment-War Tradeoff"   #
#                                             #
# J. Tyson Chatagnier and Emanuele Castelli   #
#                                             #
# Forthcoming in Political Research Quarterly #
#                                             #
# File created 8 August 2016                  #
###############################################

#############################
# Function to generate time #
# since last event          #
#############################

GetTime <- function(event,
                    ongo,
                    group,
                    time
                    ) {
  require(data.table)
  if (length(group) != length(time)) {
    stop("'group' and 'time' vectors must be equal lengths.")
  }
  noinit.last <- (sequence(rle(event)$lengths) * event == 1) +
                 ((sequence(rle(event)$lengths) * !event) > 1)
                        # This tells us whether we had an onset last year
  peace.year <- (noinit.last == 1) * (ongo == 0)
                        # We want to count actual peace years. These are years
                        # in which there are no ongoing conflicts and there
                        # was no conflict the previous year.
  time.since <- by(peace.year,
                   group,
                   function(x) {
                     sequence(rle(x)$lengths) * x
                   }
                   )
  
  group.time <- na.omit(cbind(group,
                              time
                              )
                        )
                        # Combine group and time and eliminate those
                        # rows in which the group is NA

  val.mat <- cbind(group.time[order(group.time[, 1]), ],
                   unlist(time.since)
                   )
                        # Sort the group.time object and combine it 
                        # with time since
  
  return(val.mat)
}

####################################
# Code to calculate average        #
# effects for substantive results. #
####################################

AvgPreds <- function(mod,
                     xmat       = NULL,
                     varcov     = NULL,
                     var.change,
                     vars       = NULL,
                     change.lo  = NULL,
                     change.hi  = NULL,
                     change.vec = NULL,
                     pol.rel    = FALSE,
                     n.out      = 1000
                     ) {
    if (is.null(xmat) & !is.numeric(mod$x)) {
        cat("Please provide an x matrix.\n")
    } else if (is.null(xmat) & is.numeric(mod$x)) {
        xmat <- mod$x
    }
    if (pol.rel == TRUE) {
      xmat <- xmat[xmat[, "cont"] == 1 | xmat[, "majpow"] == 1, ]
    }
    k <- ncol(xmat)
    b <- mod$coefficients
    if (is.null(change.vec)) {
        change.lo <- ifelse(is.null(change.lo),
                            min(xmat[, var.change]),
                            change.lo
        )
        change.hi <- ifelse(is.null(change.hi),
                            max(xmat[, var.change]),
                            change.hi
        )
        change.vec <- seq(from       = change.lo,
                          to         = change.hi,
                          length.out = n.out
        )
    }
    xc <- xmat
    pred.mat <- matrix(nrow = length(change.vec),
                       ncol = 3
                       )
    for (i in 1 : length(change.vec)) {
      xc[, var.change] <- change.vec[i]
      if (!is.null(vars)) {
        for (j in 1 : nrow(vars)) {
          xc[, vars[j, 1]] <- vars[j, 2]
        }
                        # vars is an optional variable that takes
                        # a matrix as its argument, with the first
                        # column being variable locations to be held
                        # constant and the second being the values
                        # at which to hold those variables.
      }
      pp <- plogis(xc %*% b)    
      deriv <- dlogis(xc %*% b)
                      # dlogis() is the derivative of plogis()
      deriv.mat <- sweep(xc,
                         1,
                         deriv,
                         '*'
                   )
      if (is.null(varcov)) {
        vhat <- vcov(mod)
      } else {
        vhat <- varcov
      }
      se.pp <- sqrt(rowSums((deriv.mat %*% vhat) * deriv.mat))
                        # Just need the standard errors, so this is
                        # quicker than matrix multiplication for the
                        # entire covariance matrix, and it avoids 
                        # problems with very large data sets.
      ci.lo <- pp - (1.96 * se.pp)
      ci.lo <- ifelse(ci.lo < 0,
                      0,
                      ifelse(ci.lo > 1,
                             1,
                             ci.lo)
                      )
      ci.hi <- pp + (1.96 * se.pp)
      ci.hi <- ifelse(ci.hi < 0,
                      0,
                      ifelse(ci.hi > 1,
                             1,
                             ci.hi)
                      )
      pred.mat[i, 1] <- mean(pp)
      pred.mat[i, 2] <- mean(ci.lo)
      pred.mat[i, 3] <- mean(ci.hi)
    }
  return(list(pred = pred.mat[, 1],
              lo   = pred.mat[, 2],
              hi   = pred.mat[, 3]
              )
         )
}

#########################
# Convenience functions #
#########################

MinMax <- function(var,
                   n.out = NULL
                   ) {
  min <- min(na.omit(var))
  max <- max(na.omit(var))
  if(is.null(n.out)) {
    vec <- c(min, 
             max
             )
  } else {
    vec <- seq(min,
               max,
               length.out = n.out
               )
  }
  return(vec)
}

UnitLength <- function(n) {
  seq(0,
      1,
      length.out = n
      )
}

############################
# Functions for imputation #
############################

ImputeMonadic <- function(data) {
  monadic.df <- data.frame(industry = data$industry,
                           service  = data$service,
                           gdp      = log(data$gdp),
                           M        = data$m,
                           majpow   = data$majpow,
                           polity   = data$polity2,
                           trade    = log(data$trade),
                           nuclear  = data$nuclear,
                           t        = data$t_fatal,
                           t2       = (data$t_fatal ^ 2),
                           t3       = (data$t_fatal ^ 3),
                           public   = data$public,
                           open     = data$open,
                           contract = data$cie,
                           irst     = data$irst,
                           pec      = data$pec,
                        # Include these two variables to help with
                        # imputation of industrialization variable
                           mid      = data$fatal_onset,
                           ccode    = data$ccode,
                           year     = data$year
                           )
                        # Create data frame

  data.imp <- amelia(monadic.df,
                     ts       = "year",
                     cs       = "ccode",
                     polytime = 1,
                     p2s      = 2,
                     noms     = c("mid",
                                  "majpow",
                                  "nuclear"
                                  ),
                     ords     = "polity"
                     )

  data.imp <- transform(data.imp, 
                        t2 = t ^ 2
                        )
  data.imp <- transform(data.imp, 
                        t3 = t ^ 3
                        )
  return(data.imp)
}

MonadicMI <- function(data.imp,
                      type = c("interpolate",
                               "extrapolate"
                               )
                      ) {
  b.base <- NULL
  b.mcd <- NULL
  b.open <- NULL
  b.cont <- NULL
  se.base <- NULL
  se.mcd <- NULL
  se.open <- NULL
  se.cont <- NULL
  
  for (i in 1 : data.imp$m) {
    data.set <- data.imp$imputations[[i]]
    if (type == "interpolate") {
      data.base <- data.set[data.set$year >= 1960 & data.set$year <= 2007, ]
      data.mcd <- data.set[data.set$year >= 1971 & data.set$year <= 2001, ] 
      data.open <- data.set[data.set$year >= 1966 & data.set$year <= 1992, ] 
      data.cont <- data.set[data.set$year >= 1960 & data.set$year <= 2000, ] 
    } else if (type == "extrapolate") {
      data.base <- data.set
      data.mcd <- data.set
      data.open <- data.set
      data.cont <- data.set
    } else {
      stop("'type' must be either 'interpolate' or 'extrapolate'")
    }
    x.base <- as.matrix(data.base[, c(1 : 11)])
    x.mcd <- as.matrix(cbind(data.mcd[, c(1 : 11, 12)]))
    x.open <- as.matrix(cbind(data.open[, c(1 : 11, 13)]))
    x.cont <- as.matrix(cbind(data.cont[, c(1 : 11, 14)]))
    mod.base <- glm(data.base$mid ~ x.base,
                    family = binomial(link = "logit"),
                    x      = TRUE
                    )
  
    rob.base <- clx(mod.base,
                    1,
                    factor(data.base$ccode)
                    )
  
    mod.mcd <- glm(data.mcd$mid ~ x.mcd,
                   family = binomial(link = "logit"),
                   x      = TRUE
                   )
  
    rob.mcd <- clx(mod.mcd,
                   1,
                   factor(data.mcd$ccode)
                   )
  
    mod.open <- glm(data.open$mid ~ x.open,
                    family = binomial(link = "logit"),
                    x      = TRUE
                    )
  
    rob.open <- clx(mod.open,
                    1,
                    factor(data.open$ccode)
                    )
  
    mod.cont <- glm(data.cont$mid ~ x.cont,
                    family = binomial(link = "logit"),
                    x      = TRUE
                    )
  
    rob.cont <- clx(mod.cont,
                    1,
                    factor(data.cont$ccode)
                    )
  
    b.base <- rbind(b.base, 
                    mod.base$coef
                    )
    b.mcd <- rbind(b.mcd, 
                   mod.mcd$coef
                   )
    b.open <- rbind(b.open, 
                    mod.open$coef
                    )
    b.cont <- rbind(b.cont, 
                    mod.cont$coef
                    )
  
    se.base <- rbind(se.base, 
                     rob.base[, 2]
                     )
    se.mcd <- rbind(se.mcd, 
                    rob.mcd[, 2]
                    )
    se.open <- rbind(se.open, 
                     rob.open[, 2]
                     )
    se.cont <- rbind(se.cont, 
                     rob.cont[, 2]
                     )
  }
  
  monadic.base <- mi.meld(q  = b.base, 
                          se = se.base
                          )
  monadic.mcd <- mi.meld(q  = b.mcd, 
                         se = se.mcd
                         )
  monadic.open <- mi.meld(q  = b.open, 
                          se = se.open
                          )
  monadic.cont <- mi.meld(q  = b.cont, 
                          se = se.cont
                          )
  
  
  b.monadic <- data.frame("B.Base"     = c(as.numeric(monadic.base$q.mi), 
                                           NA
                                           ),
                        # The baseline model has no capitalism variable, so
                        # we include an NA to ensure equal lengths
                          "B.Public"   = as.numeric(monadic.mcd$q.mi),
                          "B.Open"     = as.numeric(monadic.open$q.mi),
                          "B.Contract" = as.numeric(monadic.cont$q.mi)
                          )
  
  se.monadic <- data.frame("se.Base"      = c(as.numeric(monadic.base$se.mi), 
                                              NA
                                              ),
                           "se.Public"    = as.numeric(monadic.mcd$se.mi),
                           "se.Open"      = as.numeric(monadic.open$se.mi),
                           "se.Contract"  = as.numeric(monadic.cont$se.mi)
                            )

  sig.monadic <- ifelse(abs(b.monadic / se.monadic) >= 2.33,
                        "***",
                        ifelse(abs(b.monadic / se.monadic) >= 1.96,
                               "**",
                               ifelse(abs(b.monadic / se.monadic) >= 1.645,
                                      "*",
                                      ""
                                      )
                               )
                        )
   
  varnames.monad <- c("Constant",
                      "Industrialization",
                      "Service",
                      "GDP",
                      "M Score",
                      "Major Power",
                      "Polity",
                      "Trade",
                      "Nuclear",
                      "t",
                      "t^2",
                      "t^3",
                      "Capitalism"
                      )
  
  rownames(b.monadic) <- varnames.monad
  rownames(se.monadic) <- varnames.monad
  rownames(sig.monadic) <- varnames.monad

  mods.list <- list(mod.base,
                    mod.mcd,
                    mod.open,
                    mod.cont
                    )

  return(list(b    = b.monadic,
              se   = se.monadic,
              sig  = sig.monadic,
              mods = mods.list
              )
         )
}

ImputeDyadic <- function(data) {
  dyadic.df <- data.frame(ind_hi       = data$ind_hi,
                          serv_ind_hi  = data$serv_ind_hi,
                          logdp        = data$gdp_lo,
                          Mrat         = data$Mrat,
                          majpow       = data$majpow,
                          lodem        = data$lodem,
                          hidem        = data$hidem,
                          cont         = data$contig,
                          nuclear      = data$nuclear,
                          allies       = data$alliance,
                          sun2cati     = data$s2uni,
                          dyadic_trade = data$lntrade,
                          t_fatal      = data$t_fatal,
                          t2_fatal     = data$t_fatal ^ 2,
                          t3_fatal     = data$t_fatal ^ 3,
                          public_hi    = data$publich,
                          capopenl     = data$openl,
                          cie_lo       = data$cielo,
                          mid_fatal    = data$fatal_onset,
                          pec_hi       = data$pec_hi,
                          irst_hi      = data$irst_hi,
                        # Include these two variables to help with
                        # imputation of industrialization variable
                          dyadid       = data$dyadid,
                          year         = data$year
                          )
                        # Create data frame

  data.imp <- amelia(dyadic.df,
                     ts       = "year",
                     cs       = "dyadid",
                     polytime = 1,
                     p2s      = 2,
                     noms     = c("mid_fatal",
                                  "majpow",
                                  "cont",
                                  "allies"
                                  ),
                     ords     = c("lodem",
                                  "hidem"
                                  )
                     )

  data.imp <- transform(data.imp, 
                        t2_fatal = t_fatal ^ 2
                        )
  data.imp <- transform(data.imp, 
                        t3_fatal = t_fatal ^ 3
                        )
  return(data.imp)
}

DyadicMI <- function(data.imp,
                     type = c("interpolate",
                              "extrapolate"
                              )
                     ) {
  b.base <- NULL
  b.mcd <- NULL
  b.open <- NULL
  b.cont <- NULL
  se.base <- NULL
  se.mcd <- NULL
  se.open <- NULL
  se.cont <- NULL
  b.base <- NULL
  b.mcd <- NULL
  b.open <- NULL
  b.cont <- NULL
  se.base <- NULL
  se.mcd <- NULL
  se.open <- NULL
  se.cont <- NULL
  
  for (i in 1 : data.imp$m) {
    data.set <- data.imp$imputations[[i]]
    if (type == "interpolate") {
      data.base <- data.set[data.set$year >= 1960 & data.set$year <= 2002, ]
      data.mcd <- data.set[data.set$year >= 1971 & data.set$year <= 2001, ] 
      data.open <- data.set[data.set$year >= 1966 & data.set$year <= 1992, ] 
      data.cont <- data.set[data.set$year >= 1960 & data.set$year <= 2000, ] 
    } else if (type == "extrapolate") {
      data.base <- data.set
      data.mcd <- data.set
      data.open <- data.set
      data.cont <- data.set
    } else {
      stop("'type' must be either 'interpolate' or 'extrapolate'")
    }
    x.base <- as.matrix(data.base[, c(1 : 15)])
    x.mcd <- as.matrix(cbind(data.mcd[, c(1 : 15, 16)]))
    x.open <- as.matrix(cbind(data.open[, c(1 : 15, 17)]))
    x.cont <- as.matrix(cbind(data.cont[, c(1 : 15, 18)]))
    mod.base <- glm(data.base$mid_fatal ~ x.base,
                    family = binomial(link = "logit"),
                    x      = TRUE
                    )

    rob.base <- clx(mod.base,
                    1,
                    factor(data.base$dyadid)
                    )

    mod.mcd <- glm(data.mcd$mid_fatal ~ x.mcd,
                   family = binomial(link = "logit"),
                   x      = TRUE
                   )

    rob.mcd <- clx(mod.mcd,
                    1,
                    factor(data.mcd$dyadid)
                    )

    mod.open <- glm(data.open$mid_fatal ~ x.open,
                    family = binomial(link = "logit"),
                    x      = TRUE
                    )

    rob.open <- clx(mod.open,
                    1,
                    factor(data.open$dyadid)
                    )

    mod.cont <- glm(data.cont$mid_fatal ~ x.cont,
                    family = binomial(link = "logit"),
                    x      = TRUE
                    )

    rob.cont <- clx(mod.cont,
                    1,
                    factor(data.cont$dyadid)
                    )

    b.base <- rbind(b.base, 
                    mod.base$coef
                    )
    b.mcd <- rbind(b.mcd, 
                   mod.mcd$coef
                   )
    b.open <- rbind(b.open, 
                    mod.open$coef
                    )
    b.cont <- rbind(b.cont, 
                    mod.cont$coef
                    )

    se.base <- rbind(se.base, 
                     rob.base[, 2]
                     )
    se.mcd <- rbind(se.mcd, 
                    rob.mcd[, 2]
                    )
    se.open <- rbind(se.open, 
                     rob.open[, 2]
                     )
    se.cont <- rbind(se.cont, 
                     rob.cont[, 2]
                     )
  }

  dyadic.base <- mi.meld(q  = b.base, 
                         se = se.base
                         )
  dyadic.mcd <- mi.meld(q  = b.mcd, 
                        se = se.mcd
                        )
  dyadic.open <- mi.meld(q  = b.open, 
                         se = se.open
                         )
  dyadic.cont <- mi.meld(q  = b.cont, 
                         se = se.cont
                         )


  b.dyadic <- data.frame("B.Base" = c(as.numeric(dyadic.base$q.mi), 
                                          NA
                                          ),
                          # The baseline model has no capitalism variable, so
                          # we include an NA to ensure equal lengths
                         "B.Public"   = as.numeric(dyadic.mcd$q.mi),
                         "B.Open"      = as.numeric(dyadic.open$q.mi),
                         "B.Contract" = as.numeric(dyadic.cont$q.mi)
                         )

  se.dyadic <- data.frame("se.Base" = c(as.numeric(dyadic.base$se.mi), 
                                            NA
                                            ),
                          "se.Public"    = as.numeric(dyadic.mcd$se.mi),
                          "se.Open"      = as.numeric(dyadic.open$se.mi),
                          "se.Contract"  = as.numeric(dyadic.cont$se.mi)
                           )
   
  sig.dyadic <- ifelse(abs(b.dyadic / se.dyadic) >= 2.33,
                       "***",
                       ifelse(abs(b.dyadic / se.dyadic) >= 1.96,
                              "**",
                              ifelse(abs(b.dyadic / se.dyadic) >= 1.645,
                                     "*",
                                     ""
                                     )
                              )
                       )

  varnames.dyad <- c("Constant",
                     "Industrialization",
                     "Service",
                     "GDP",
                     "M Ratio",
                     "Major Power",
                     "Low Democracy",
                     "High Democracy",
                     "Contiguity",
                     "Nuclear",
                     "Allies",
                     "S Score",
                     "Trade",
                     "t",
                     "t^2",
                     "t^3",
                     "Capitalism"
                     )

  rownames(b.dyadic) <- varnames.dyad
  rownames(se.dyadic) <- varnames.dyad
  rownames(sig.dyadic) <- varnames.dyad

  mods.list <- list(mod.base,
                    mod.mcd,
                    mod.open,
                    mod.cont
                    )

  return(list(b    = b.dyadic,
              se   = se.dyadic,
              sig  = sig.dyadic,
              mods = mods.list
              )
         )
}

####################################
# Cluster Robust Standard          #
# Errors, downloaded from          #
# http://people.su.se/~ma/clmclx.R #
####################################

ClusterVarCov <- function(mod,
                          .data,
                          dfcw = 1,
                          cluster
                          ) {
  rel.vals <- complete.cases(.data[, c(as.character(mod$formula[2]), 
                                       colnames(mod$x)[-1])
                                   ]
                             )
  cluster.mod <- cluster[rel.vals]
  varcov <- clx(fm      = mod,
                dfcw    = dfcw,
                cluster = cluster.mod,
                vcov    = TRUE
                )
  return(varcov)
}

clx <- function(fm, 
                dfcw, 
                cluster,
                vcov = FALSE
                ) {
                        # R-codes (www.r-project.org) for computing
                        # clustered-standard errors. Mahmood Arai, Jan 26, 2008.
             
                       # The arguments of the function are:
                        # fitted model, cluster1 and cluster2
                        # You need to install libraries `sandwich' and `lmtest'
             
                        # reweighting the var-cov matrix for the within model
         require(sandwich)
         require(lmtest)
         M <- length(unique(cluster))   
         N <- length(cluster)           
         K <- fm$rank                        
         dfc <- (M / (M - 1)) * ((N - 1) / (N - K))  
         uj  <- apply(estfun(fm),
                      2, 
                      function(x) {
                        tapply(x, 
                               cluster, 
                               sum
                               )
                      }
                      )
         vcovCL <- dfc*sandwich(fm, 
                                meat = crossprod(uj) / N
                                ) * dfcw
         if (vcov) {
           return(vcovCL)
         } else {
           coeftest(fm, vcovCL)
         }
}

##################################
# Calculation of hazard function #
##################################

GetHaz <- function(t,
                   b,
                   varcov,
                   t.loc
                   ) {
                        # This function calculates a predicted hazard and
                        # its associated confidence interval. It takes the
                        # following arguments:
                        # t is a sequence of time points
                        # b is a vector of coefficient estimates
                        # varcov is a covariance matrix
                        # t.loc is a vector of locations for time variables
  b.t <- b[t.loc]
  haz <- t * b.t[1] + t ^ 2 * b.t[2] + t ^ 3 * b.t[3]
  se.haz <- sqrt(VarMFX(varcov,
                        t.loc[1],
                        t.loc[2],
                        t.loc[3]
                        )
                 )
  haz.lo <- haz - 1.96 * se.haz
  haz.hi <- haz + 1.96 * se.haz
  return(list(haz = haz,
              lo  = haz.lo,
              hi  = haz.hi
              )
         )
}

##########################
# Function to calculate  #
# descriptive statistics #
##########################

GetDesc <- function(vars) {
  funs <- c("min",
            "max",
            "mean",
            "sd"
            )
  desc.df <- data.frame(Variable = names(vars),
                        Minimum  = NA,
                        Maximum  = NA,
                        Mean     = NA,
                        SD       = NA
                        )
  for (i in 1 : length(funs)) {
    desc.df[, i + 1] <- apply(vars,
                               2,
                               funs[i],
                               na.rm = TRUE
                               )
  }
  return(desc.df)
}

##########################################
# Function to get covariance matrix      #
# for marginal effects on an interaction #
##########################################

VarMFX <- function(varcov,
                   num.x1,
                   num.x2,
                   num.x3
                   ) {
  var.x1 <- varcov[num.x1, num.x1]
  var.x2 <- varcov[num.x2, num.x2]
  var.x3 <- varcov[num.x3, num.x3]
  cov.x1.x2 <- varcov[num.x1, num.x2]
  cov.x1.x3 <- varcov[num.x1, num.x3]
  cov.x2.x3 <- varcov[num.x2, num.x3]
  var.sum <- var.x1 + var.x2 + var.x3 + 2 * (cov.x1.x2 + cov.x1.x3 + cov.x2.x3)
  return(var.sum)
}

####################################
# Function to print to tables from #
# multiple MI models, for appendix #
####################################

MakeTable <- function(.list,
                      filename = NULL
                      ) {
  .list$mods[[1]]$coefficients <- .list$b[, 1][-length(.list$b[, 1])]
                        # There will be an empty space in our lists
                        # for the first column. So we need to exclude it.
  .list$mods[[2]]$coefficients <- .list$b[, 2]
  .list$mods[[3]]$coefficients <- .list$b[, 3]
  .list$mods[[4]]$coefficients <- .list$b[, 4]
  .list$mods[[1]]$se <- .list$se[, 1][-length(.list$b[,1])]
  .list$mods[[2]]$se <- .list$se[, 2]
  .list$mods[[3]]$se <- .list$se[, 3]
  .list$mods[[4]]$se <- .list$se[, 4]
  names(.list$mods[[1]]$coefficients) <- rownames(.list$b)[-nrow(.list$b)]
  names(.list$mods[[2]]$coefficients) <- rownames(.list$b)
  names(.list$mods[[3]]$coefficients) <- rownames(.list$b)
  names(.list$mods[[4]]$coefficients) <- rownames(.list$b)
  if (!is.null(filename)) {
    write(c("\\documentclass[12pt]{article} \n\\usepackage{amsmath, amssymb, natbib, setspace, fullpage, dcolumn} \n\\begin{document} \n \n",
            apsrtable(.list$mods[[1]],
                      .list$mods[[2]],
                      .list$mods[[3]],
                      .list$mods[[4]],
                      stars = "default",
                      se    = "robust"
                      ),
            "\n\\end{document}"),
            file = filename
           )
  } else {
    apsrtable(.list$mods[[1]],
              .list$mods[[2]],
              .list$mods[[3]],
              .list$mods[[4]],
              stars = "default",
              se    = "robust"
              )
  }
}
