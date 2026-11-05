
########################
##### FORMAT DATA  #####
########################

ccode <-  'cow_ccode'
ccode1 <- 'cow_ccode1'
ccode2 <- 'cow_ccode2'
ccodes <- c(ccode1, ccode2)
years <- seq(start_yr, pred_yr, by=incr_yr)

#Load table of entry-exit dates
entry_exit_data <- read.csv(file.path(ROOT, 'work_saom/X085/cow_tab.csv'), as.is=TRUE)

#Load constant dyadic data
dyad_cons <- read.csv(file.path(ROOT, "work_saom/X085/dyad_data.csv"), as.is = TRUE)

#Load dyad-year data
dyad_year <- read.csv(file.path(ROOT, "work_saom/X085/dyad_year_data.csv"), as.is = TRUE)

#Load state-year data
state_year <- read.csv(file.path(ROOT, "work_saom/X085/state_year_data.csv"), as.is = TRUE)


## Alliances ##
ally_data <- dyad_year[!is.na(dyad_year$defense) & dyad_year$defense == 1, c(ccodes, 'year')]
ally_data <- undirect(ally_data, ccodes = ccodes)
ally_data <- ally_data[!duplicated(ally_data[,c(ccodes, 'year')]),]
ally_data$val <- 1


## Conflict ##
con_data <- dyad_year[!is.na(dyad_year$minhost) & dyad_year$minhost >= con_level, c(ccodes, 'year')]
con_data <- undirect(con_data, ccodes = ccodes)
con_data <- lag_dyad(con_data, lag = con_lag, ccodes = ccodes)
con_data <- con_data[!duplicated(con_data[,c(ccodes, 'year')]),]
con_data$val <- 1


## Democracy ##
dem_data <- state_year[!is.na(state_year$polity2), c(ccode, 'year', 'polity2')]
dem_data$val <- NA
dem_data$val[dem_data$polity2 >= dem_level] <- 1 
dem_data$val[dem_data$polity2 < dem_level] <- 0  


## Power ##
## (logged and rescaled to 0-1)
power_data <- state_year[!is.na(state_year$cinc), c(ccode, 'year', 'cinc')]
power_data$cinc[power_data$cinc < 0] <- NA
power_data$lcinc <- log(power_data$cinc)
power_data$lcinc <- power_data$lcinc - min(power_data$lcinc, na.rm=TRUE)
power_data$lcinc <- power_data$lcinc / max(power_data$lcinc, na.rm=TRUE)
power_data$val <- power_data$lcinc


## Power ratio ##
## (logged and rescaled to 0-1)
powrat_data <- matrix(NA, nrow=sum(table(power_data$year)^2), ncol=5)
rownum <- 0
for (y in years) {
  rows <- which(power_data$year == y)
  for (i in rows) {
    for (j in rows) {
      rownum <- rownum + 1
      powrat_data[rownum,] <- 
        c(y, power_data[[ccode]][i], power_data[[ccode]][j], 
          power_data$cinc[i],     power_data$cinc[j])
    }
  }
}
powrat_data <- as.data.frame(powrat_data)
names(powrat_data) <- c('year', ccode1, ccode2, 'cinc1', 'cinc2')
powrat_data$cinc_high <- ifelse(powrat_data$cinc1 > powrat_data$cinc2, 
                                powrat_data$cinc1,  powrat_data$cinc2)
powrat_data$cinc_low  <- ifelse(powrat_data$cinc1 < powrat_data$cinc2, 
                                powrat_data$cinc1,  powrat_data$cinc2)
powrat_data$prat <- powrat_data$cinc_low / powrat_data$cinc_high
powrat_data$lprat <- log(powrat_data$prat)
powrat_data$lprat <- powrat_data$lprat - min(powrat_data$lprat, na.rm=TRUE)
powrat_data$lprat <- powrat_data$lprat / max(powrat_data$lprat, na.rm=TRUE)
powrat_data$val <- powrat_data$lprat


## Distance ##
capdist_data <- dyad_cons[!is.na(dyad_cons$kmdist) ,c(ccodes, 'kmdist')]
capdist_data <- undirect(capdist_data, ccodes = ccodes)
capdist_data$val <- capdist_data$kmdist / 1000


## Contiguity ##
contig_data <- dyad_cons[!is.na(dyad_cons$contig) & dyad_cons$contig == 1, c(ccodes, 'contig')]
contig_data <- undirect(contig_data, ccodes = ccodes)
contig_data$val <- 1


## Language ##
lang_data <- dyad_cons[!is.na(dyad_cons$lang) & dyad_cons$lang == 1, c(ccodes, 'lang')]
lang_data <- undirect(lang_data, ccodes = ccodes)
lang_data$val <- 1


##################################
### Assisgn data sets to roles ###
##################################

#Network dependent variables
ndv_list <- c('ally', 'con')

#Behavioral dependent variables
bdv_list <- c('dem')

#Exogenous dynamic individual covariate
iv_list <- c('power')

#Exogenous dynamic dyadic covariate
dyadv_list <- c('powrat')

#Exogenous constant dyadic covariates
dyadc_list <- c('capdist', 'contig', 'lang')


##############################
### Format data for RSiena ###
##############################

est_years <- seq(start_yr, end_yr, incr_yr)

est_data <- formatData(
  ee     = entry_exit_data, 
  yrs    = est_years, 
  ndv    = ndv_list, 
  bdv    = bdv_list, 
  iv     = iv_list, 
  dyadv  = dyadv_list,
  dyadc  = dyadc_list,
  ccodes = ccodes,
  ccode  = ccode)


####################################
