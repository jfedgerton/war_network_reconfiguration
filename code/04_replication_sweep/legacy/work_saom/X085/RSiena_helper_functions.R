######################################
##### HELPER FUNCTIONS FOR RSIENA ####
######################################

##############################################
### Generate lagged dyadic values          ###
##############################################
# Input:                                     #
#   data:   dyad-year data frame             #
#   lag:    maximum length of lag in years   #
#   ccodes: vector of two column names       #
#           containing country codes         #
# Output:                                    #
#   Dyad-year data frame with lagged entries #
##############################################
lag_dyad <- function(data, lag, ccodes) {
  results <- data
  for (i in c(1:lag)) {
    new_data <- data
    new_data$year <- data$year + i
    results <- rbind(results, new_data)
  }
  results <- results[!duplicated(results[,c(ccodes, 'year')]),]
}

########################################
### Convert undirected dyads to      ###
### symmetric directed dyads         ###
########################################
# Input:                               #
#   data:   dyad-year data frame       #
#           (one entry per dyad-year)  #
#   ccodes: vector of two column names #
#           containing country codes   #
# Output:                              #
#      directed dyad-year data frame   #
#      (two entries per dyad-year)     #
########################################
undirect <- function(data, ccodes) {
  cols <- c(ccodes, names(data)[!names(data) %in% ccodes])
  cols_rev <- c(rev(ccodes), names(data)[!names(data) %in% ccodes])
  data <- data[, cols]
  data_rev <- data[, cols_rev]
  colnames(data_rev) <- cols
  data <- rbind(data, data_rev)
}

########################################################
### Convert from edge list to adjacency matrix       ###
########################################################
# Input:                                               #
#  elist: edge list with node ids in first two columns #
#         and edge values in third column              #
#  n:     total number of nodes                        #
#  zeros: option to replace missing edge values with 0 #
########################################################
# Output:                                              #
#  n x n adjacency matrix                              #
########################################################
net2mat <- function(elist, n, zeros=TRUE) {
  if (zeros==FALSE) {
    adj_mat <- matrix(NA, n, n)
  } else {
    adj_mat <- matrix(0, n, n)
  }  
  adj_mat[elist[, 1:2]] <- elist[, 3]
  adj_mat
}

###############################################
### Convert from state-year format to       ###
### wide matrix format                      ###
###############################################
# Input:                                      #
#   data:    dyad-year data frame             #
#   val_col: column name with edge values     #
#   sts:     numeric vector of state ids      #
#   yrs:     numeric vector of years          #
#   ccode:   column name with country codes   #
# Output:                                     #
#   Wide format matrix (years across columns) #
#   formatted for sienaNet() or varCovar()    #
###############################################
dat2mat <- function(data, val_col, sts, yrs, ccode) {
  ns <- length(sts)
  ny <- length(yrs)
  mat <- matrix(NA, ns, ny)
  for (row in c(1:ns)) {
    for (col in c(1:ny)) {
      val <- data[[val_col]][data[[ccode]] == sts[row] & data$year == yrs[col]]
      if (length(val) > 0) {
        mat[row, col] <- val
      }
    }
  }
  mat
}

##################################################
### Convert dyadic data to adjacency matrix    ###
##################################################
# Input:                                         #
#   data:    dyadic data frame                   #
#   val_col: column name with edge values        #
#   sts:     numeric vector of state ids         #
#   ccodes:  vector of two column names          #
#            containing country codes            #
#   zeros:   option to replace missing           #
#            edge values with 0                  #
# Output:                                        #
#   Adjacency matrix formatted for coDyadCovar() #
##################################################
constantNet <- function(data, val_col, sts, ccodes, zeros=FALSE) {
  n <- length(sts)
  data <- data[ data[[ ccodes[1] ]] %in% sts &
                data[[ ccodes[2] ]] %in% sts , ]
  data <- data[, c(ccodes, val_col)]
  data <- as.matrix(data, rownames.force=FALSE)
  for (i in c(1:n)) {
    data[data == sts[i]] <- i
  }  
  if (zeros==FALSE) {
    adj_mat <- matrix(NA, n, n)
  } else {
    adj_mat <- matrix(0, n, n)
  }  
  adj_mat[data[, 1:2]] <- data[, 3]
  adj_mat
}

################################################
### Convert dyad-year data to                ###
### multi-period adjacency matrices          ###
################################################
# Input:                                       #
#   data:    dyad-year data frame              #
#   val_col: column name with edge values      #
#   sts:     numeric vector state ids          #
#   yrs:     numeric vector of years           #
#   ccodes:  vector of two column names        #
#            containing country codes          #
#   ee_tab:  table of entry-exit times         #
#            (created by FormatEntryExit())    #
#   absent:  option to add missing values      #
#            for absent nodes                  #
# Output:                                      #
#   Multi-period adjacency matrix              #
#   formatted for sienaNet() or varDyadCovar() #
################################################
changingNet <- function(data, val_col, sts, yrs, ccodes, 
                        ee_tab, absent = 'missing') {
  ns <- length(sts)
  ny <- length(yrs)
  cols <- c(ccodes, 'year', val_col)
  data <- data[ data[[ ccodes[1] ]] %in% sts &
                data[[ ccodes[2] ]] %in% sts , ]
  data <- data[, cols]
  out <- array(NA, dim=c(ns, ns, ny))
  for (j in c(1:length(yrs))) {
    y <- yrs[j]
    dat <- data[data$year == y,]
    dat$year <- NULL
    net <- constantNet(dat, val_col = val_col, sts = sts, ccodes = ccodes, zeros = TRUE)
    sts_in <- unique(ee_tab$ccode[ee_tab$styear <= y & ee_tab$endyear >= y])
    cond <- (!(sts %in% sts_in))
    if (absent=='structural') {
      net[cond, ] <- 10
      net[, cond] <- 10
    }
    if (absent=='missing') {
      net[cond, ] <- NA
      net[, cond] <- NA
    }
    if (absent=='zeros') {
      net[cond, ] <- 0
      net[, cond] <- 0
    }
    out[, , j] <- net
  }
  out  
}


################################################
### Format Entry-Exit Data                   ###
################################################
# Input:                                       #
#   ee:  entry-exit data                       #
#        (e.g. COW system membership table)    #
#   yrs: vector of years starting each period  #
#   sts: vector of state ids                   #
# Output:                                      #
#   List with entry-exit information formatted #
#   for sienaCompositionChange()               #
################################################
formatEntryExit <- function(ee, yrs, sts) {
  num_years <- length(yrs)
  start_yr <- min(yrs)
  end_yr <- max(yrs)
  incr_yr <- abs(yrs[1] - yrs[2])
  ee$styear <- year(as.Date(ee$start))
  ee$endyear <- year(as.Date(ee$end))
  ee <- ee[ee$styear < end_yr, ]
  ee <- ee[ee$endyear > start_yr, ]
  ee$styear[ee$styear < start_yr] <- start_yr
  ee$endyear[ee$endyear > end_yr] <- end_yr
  ee <- ee[((ee$endyear - ee$styear) > 1), ]
  ee <- ee[ee$ccode %in% sts,]
  num_periods <- dim(ee)[1]
  nodes <- sort(unique(ee$ccode))
  num_states <- length(nodes)
  ee_list <- vector("list", num_states)
  for (p in seq(1:num_periods)) {
    p_start <- round(((ee$styear[p] - start_yr) / incr_yr) + 1, 2)
    p_end <- round(((ee$endyear[p] - start_yr) / incr_yr) + 1, 2)
    siena_id <- which(nodes == ee$ccode[p])
    if (length(siena_id) > 0) {
      ee_list[[siena_id]] <- c(ee_list[[siena_id]], p_start, p_end)
    }
  }
  out <- vector('list')
  out$nodes <- nodes
  out$ee_list <- ee_list
  out$ee_tab <- ee
  return(out)
}


####################################################################
### Create Data Object for Co-Evolutionary RSiena Model          ###
####################################################################
# Input:                                                           #
#   ee:     entry-exit data                                        #
#           (e.g. COW system membership table)                     #
#   yrs:    vector of years starting each period                   #
#   ndv:    vector of names of network dependent variables         #
#   bdv:    vector of names of behavioral dependent variables      #
#   iv:     vector of names of exogenous individual variables      #
#   dyadv:  vector of names of exogenous dynamic dyadic variables  #
#   dyadc:  vector of names of exogenous constant dyadic variables #
#   ccodes: vector of two column names with dyadic country codes   #
#   ccode:  single column name with country codes                  #
# Output:                                                          #
#   RSiena formatted data object from sienaDataCreate()            #
####################################################################
formatData <- function(ee, yrs, ndv, bdv, iv, dyadv, dyadc, ccodes, ccode) {
  
  #Define Siena NodeSet
  years_est <- yrs
  ustates <- vector('character')
  for (varname in ndv) {
    us <- unique(c(get(paste0(varname, '_data'))[[ccodes[1]]], 
                   get(paste0(varname, '_data'))[[ccodes[2]]]))
    if (length(ustates > 0)) {
      ustates <- intersect(ustates, us)
    } else {
      ustates <- unique(us)
    }
  }
  ustates <- sort(ustates)
  ee_est <- formatEntryExit(ee, yrs = years_est, sts = ustates)
  nodes_siena <- sienaNodeSet(length(ee_est$nodes), nodeSetName="States")
  
  #Create composition change object
  compchange <- sienaCompositionChange(
    ee_est$ee_list, nodeSet = "States", option = 3)
  
  #################
  ## Network DVs ##
  #################
  
  for (varname in ndv) {
    ndv_est <- sienaNet(changingNet(
      get(paste0(varname, '_data')), 'val', ccodes=ccodes,
      sts=ee_est$nodes, yrs=years_est, 
      ee_tab=ee_est$ee_tab), 
      type="oneMode", nodeSet="States")
    assign(paste0(varname, '_est'), ndv_est)
  }
  
  ####################
  ## Behavioral DVs ##
  ####################
  
  for (varname in bdv) {
    bdv_est <- sienaNet(dat2mat(
      get(paste0(varname, '_data')), 'val', ccode=ccode, 
      sts=ee_est$nodes, yrs=years_est), 
      type="behavior", nodeSet="States")
    assign(paste0(varname, '_est'), bdv_est)
  }
  
  ##############################################
  ## Exogenous Dyanamic Individual Covariates ##
  ##############################################
  
  for (varname in iv) {
    iv_est <- varCovar(dat2mat(
      get(paste0(varname, '_data')), 'val', ccode=ccode, 
      sts=ee_est$nodes, yrs=years_est), 
      nodeSet="States")
    assign(paste0(varname, '_est'), iv_est)
  }
  
  ##########################################
  ## Exogenous Dyanamic Dyadic Covariates ##
  ##########################################
  
  for (varname in dyadv) {
    dyadv_est <- varDyadCovar(changingNet(
      get(paste0(varname, '_data')), 'val', sts=ee_est$nodes, 
      yrs=years_est[-length(years_est)], ee_tab=ee_est$ee_tab, 
      ccodes=ccodes), 
      type="oneMode", nodeSets=c("States","States"))
    assign(paste0(varname, '_est'), dyadv_est)
  }
  
  ##########################################
  ## Exogenous Constant Dyadic Covariates ##
  ##########################################
  
  for (varname in dyadc) {
    z <- FALSE
    if (length(unique(get(paste0(varname, '_data'))$val)) == 1) {
      z <- TRUE
    }
    dyadc_est <- coDyadCovar(constantNet(
      get(paste0(varname, '_data')), 'val', 
      sts=ee_est$nodes, ccodes=ccodes, zeros=z), 
      type="oneMode", nodeSets=c("States","States"))
    assign(paste0(varname, '_est'), dyadc_est)
  }
  
  ##############################################
  ## COMBINE DATA SETS INTO SIENA DATA OBJECT ##
  ##############################################
  
  data_list <- c(ndv, bdv, iv, dyadv, dyadc)
  
  cmd_txt <- paste0(
    'sienaDataCreate(compchange, ', 
    paste(paste0(data_list, '_est'), collapse=', '), 
    ', nodeSets = list(nodes_siena))')
  
  data_est <- eval(parse(text=cmd_txt))
  
  return(data_est)
}

###################################################################
### Create Full Effects Object for Co-Evolutionary RSiena Model ###
###################################################################
# Input:                                                          #
#   data:  RSiena data object generated by formatData()           #
#   ndv:   vector of names of network dependent variables         #
#   bdv:   vector of names of behavioral dependent variables      #
#   iv:    vector of names of exogenous individual variables      #
#   dyadv: vector of names of exogenous dynamic dyadic variables  #
#   dyadc: vector of names of exogenous constant dyadic variables #
# Output:                                                         #
#   RSiena model effects object with fully specified effects      #
###################################################################
createFullEffects <- function(data, ndv, bdv, iv, dyadv, dyadc, verbose = FALSE) {
  
  if (!verbose) {
    if (.Platform$OS.type == 'windows') {
      sink('NUL')
    } else {
      sink('/dev/null')
    }
  }
  
  effects_siena <- getEffects(data)
  
  # Already (automatically) included are:
  #     1. constant rates, for each dependent variable, for each period
  #     2. density for each dependent network
  #     3. linear behaviorial effect
  # In addition, include the following effects:
  
  #################################
  ## [Network DV --> Network DV] ##
  #################################
  
  for (varname in ndv) {
    vname_est <- paste0(varname, '_est')
    
    #Popularity
    effects_siena <- includeEffects(
      effects_siena, inPop, name = vname_est)
    
    #Remove transitive triads
    effects_siena <- includeEffects(
      effects_siena, transTriads, name = vname_est, include=FALSE)
    
    #Add transitive ties
    effects_siena <- includeEffects(
      effects_siena, transTies, name = vname_est)
  }
  
  ###################################
  ## [Network DV1 --> Network DV2] ##  
  ###################################
  
  #Network cross-products
  if (length(ndv) > 1) {
    for (varname1 in ndv) {
      for (varname2 in ndv) {
        if (varname1 != varname2) {
          vname1_est <- paste0(varname1, '_est')
          vname2_est <- paste0(varname2, '_est')
          effects_siena <- includeEffects(
            effects_siena, crprod,
            name = vname1_est, 
            interaction1 = vname2_est)
        }
      }
    }
  }
  
  ####################################
  ## [Behavioral DV --> Network DV] ##
  ####################################
  
  for (varname1 in ndv) {
    for (varname2 in bdv) {
      vname1_est <- paste0(varname1, '_est')
      vname2_est <- paste0(varname2, '_est')
      
      #Ego = Alter 
      effects_siena <- includeEffects(
        effects_siena, sameX, 
        name = vname1_est, 
        interaction1 = vname2_est)
      
      #(Ego = Alter) * Alter
      effects_siena <- includeInteraction(
        effects_siena, sameX, altX, 
        name = vname1_est, 
        interaction1 = c(vname2_est, vname2_est))
    }
  }
  
  ##################################
  ## [Network DV --> Behavior DV] ##                                 
  ##################################
  
  for (varname1 in bdv) {
    for (varname2 in ndv) {
      vname1_est <- paste0(varname1, '_est')
      vname2_est <- paste0(varname2, '_est')
      
      #Network degree
      effects_siena <- includeEffects(
        effects_siena, outdeg, 
        name = vname1_est, 
        interaction1 = vname2_est)
      
      #Average alter value
      effects_siena <- includeEffects(
        effects_siena, avAlt, 
        name = vname1_est, 
        interaction1 = vname2_est)
    }
  }
  
  ###########################################
  ## [Individual Covariate --> Network DV] ##
  ###########################################
  
  for (varname1 in ndv) {
    for (varname2 in iv) {
      vname1_est <- paste0(varname1, '_est')
      vname2_est <- paste0(varname2, '_est')
      
      #Alter value
      effects_siena <- includeEffects(
        effects_siena, altX, 
        name = vname1_est, 
        interaction1 = vname2_est)
    }
  }
  
  ##############################################
  ## [Individual Covariate --> Behavioral DV] ##
  ##############################################
  
  for (varname1 in bdv) {
    for (varname2 in iv) {
      vname1_est <- paste0(varname1, '_est')
      vname2_est <- paste0(varname2, '_est')
      
      #Self value
      effects_siena <- includeEffects(
        effects_siena, effFrom, 
        name = vname1_est, 
        interaction1 = vname2_est)
    }
  }
  
  #######################################
  ## [Dyadic Covariate --> Network DV] ##
  #######################################
  
  for (varname1 in ndv) {
    for (varname2 in c(dyadv, dyadc)) {
      vname1_est <- paste0(varname1, '_est')
      vname2_est <- paste0(varname2, '_est')
      
      #Dyadic value
      effects_siena <- includeEffects(
        effects_siena, X, 
        name = vname1_est, 
        interaction1 = vname2_est)
    }
  }
  if (!verbose) {
    sink()
  }
  return(effects_siena)
}


#################################################################
### Create Contrained Effects Object for RSiena Model         ###
### (removes co-evolutionary effects from full model)         ###
#################################################################
# Input:                                                        #
#   efx: RSiena effects object generated by CreateFullEffects() #
#   ndv: vector of names of network dependent variables         #
#   bdv: vector of names of behavioral dependent variables      #
# Output:                                                       #
#   RSiena model effects object with constrained effects        #
#################################################################
constrainEffects <- function(efx, ndv, bdv, verbose = FALSE) {
  
  if (!verbose) {
    if (.Platform$OS.type == 'windows') {
      sink('NUL')
    } else {
      sink('/dev/null')
    }
  }
  
  #Network cross-products
  if (length(ndv) > 1) {
    for (varname1 in ndv) {
      for (varname2 in ndv) {
        if (varname1 != varname2) {
          vname1_est <- paste0(varname1, '_est')
          vname2_est <- paste0(varname2, '_est')
          efx <- includeEffects(
            efx, crprod, include = FALSE,
            name = vname1_est, 
            interaction1 = vname2_est)
        }
      }
    }
  }
  
  for (varname1 in ndv) {
    for (varname2 in bdv) {
      vname1_est <- paste0(varname1, '_est')
      vname2_est <- paste0(varname2, '_est')
      
      #Ego = Alter 
      efx <- includeEffects(
        efx, sameX, include = FALSE,
        name = vname1_est, 
        interaction1 = vname2_est)
      
      #(Ego = Alter) * Alter
      efx <- includeInteraction(
        efx, sameX, altX, include = FALSE,
        name = vname1_est, 
        interaction1 = c(vname2_est, vname2_est))
    }
  }
  
  for (varname1 in bdv) {
    for (varname2 in ndv) {
      vname1_est <- paste0(varname1, '_est')
      vname2_est <- paste0(varname2, '_est')
      
      #Network degree
      efx <- includeEffects(
        efx, outdeg, include = FALSE,
        name = vname1_est, 
        interaction1 = vname2_est)
      
      #Average alter value
      efx <- includeEffects(
        efx, avAlt, include = FALSE,
        name = vname1_est, 
        interaction1 = vname2_est)
    }
  }
  if (!verbose) {
    sink()
  }
  return(efx)
}


######################################################################
### Estimate RSiena Model                                          ###
######################################################################
### The algorithm begins searching by using sub-phases             ###
### with a relatively small number of iterations, and then         ###  
### increass the number of iterations as convergence improves,     ###
### taking the previous estimates as the starting values           ###
### for each subsequent attempt, until full convergence is reached ###
######################################################################
# Input:                                                             #
#   data:     RSiena data object generated by formatData()           #
#   efx:      RSiena effects object generated by CreateFullEffects() #
#   nreps:    maximum number of repeated attempts at convergence     #
#   n2nstart: two integers, giving starting and ending               #
#             vaues for n2start (minimum iterations in phase 2)      #
#   n3sims:   two integers, giving starting and ending               #
#             values for n3 (total iterations in phase 3)            #
#   diag:     diagonalization weight for derivative matrix           #
#   firstg:   step size for Robbins-Monro procedure                  #
#   ntconv:   tconv convergence threshold                            #
#   ntcmax:   tcmax convergence threshold                            #
#   ncpu:     number of cpu cores to use for parallel computations   #
#   stub:     model name stub for output files                       #
# Output:                                                            #
#   RSiena sienaFit object from siena07()                            #
######################################################################
est_siena <- function(data, efx, nreps, n2start, n3sims, diag, firstg, ntconv, ntcmax, ncpu, stub) {
  
  st <- as.integer(Sys.time())
  
  n3sims1 <- n3sims[1]
  n3sims2 <- n3sims[2]
  n3sims  <- n3sims1
  
  n2start1 <- n2start[1]
  n2start2 <- n2start[2]
  n2start  <- n2start1
  
  est <- NULL
  n2_switched <- FALSE
  n3_switched <- FALSE
  tcmax_last <- 20
  
  cat("\n\n", "Starting calculations with:", "\n",
      "n2start =", n2start, "\n", "n3 =", n3sims, "\n\n")
  
  for (rep in c(1:nreps)) {
    pname <- paste0("AllyConDem_EST_", stub)
    model <- sienaModelCreate(projname = pname, modelType = 3, 
                              truncation = 10, nsub = 1, 
                              n2start = n2start, n3 = n3sims, 
                              diagonalize = diag, firstg = firstg)
    new_est <- siena07(
      model, 
      data = data, 
      effects = efx, 
      useCluster = as.logical(ncpu > 1), 
      nbrNodes = ncpu, 
      prevAns = est)
    
    #save(new_est, file= paste0('EST.', version_num, '.n', rep, '.RData'))
    
    tconv <- max(new_est$tconv)
    tcmax <- new_est$tconv.max
    et <- as.integer(Sys.time()) 
    cat("Max tconv: ", tconv, "\n")
    cat("Max tconv ratio: ", tcmax, "\n")
    cat("Caculation time:", (et - st) / 3600, "hours", "\n\n")
    if (is.na(tcmax) ||  tcmax > tcmax_last * 1.5) {
      cat("Discarding estimate", "\n\n")
      next
    }
    est <- new_est
    tcmax_last <- tcmax
    
    if (tconv < ntconv*3 && tcmax < ntcmax*3 && !n2_switched) {
     cat("Switching to n2start =", n2start2, "\n\n")
     n2start <- n2start2
     n2_switched <- TRUE
     next
    }
    if (tconv < ntconv*2 && tcmax < ntcmax*2 && !n3_switched) {
      cat("Switching to n3 =", n3sims2, "\n\n")
      n3sims <- n3sims2
      n3_switched <- TRUE
      next
    }
    if (tconv < ntconv && tcmax < ntcmax) {
      cat("Estimation converged", "\n\n")
      break
    } 
  }
  return(est)
}


######################################################################
### Generate Simulated Realizations from RSiena Model              ### 
######################################################################
# Input:                                                             #
#   data:     RSiena data object generated by formatData()           #
#   efx:      RSiena effects object generated by CreateFullEffects() #
#   n3sims:   number of simulated realizations to generate           #
#   ncpu:     number of cpu cores to use for parallel computations   #
#   stub:     model name stub for output files                       #
#   est:      optional results from prior estimate                   #
# Output:                                                            #
#   RSiena sienaFit object from siena07() with simulated outcomes    #
######################################################################
sim_siena <- function(data, efx, n3sims, ncpu, stub, est = NULL) {
  
  pname <- paste0("AllyConDem_SIM_", stub)
  model <- sienaModelCreate(projname = pname, 
                            modelType = 3, nsub = 0, 
                            n3 = n3sims, simOnly = TRUE)
  
  sim <- siena07(
    model, 
    data = data, 
    effects = efx, 
    useCluster = as.logical(ncpu > 1), 
    nbrNodes = ncpu, 
    prevAns = est,
    returnDeps = TRUE)
  
  return(sim)
}


