# Data preparation for cid 65 (Hegre, Oneal & Russett 2010), mirrors HOR_JPR_table3.do up to estimation.
# Returns data.frame `d` (1950-2000) with all constructed variables; spline knots passed as argument.
suppressPackageStartupMessages({library(haven); library(dplyr)})
W65 <- file.path(Sys.getenv("SWEEP_DATA"), "work/65"); stopifnot(nzchar(Sys.getenv("SWEEP_DATA")))

raw65 <- function() {
  d <- as.data.frame(read_dta(file.path(W65, "HOR_JPR_table3.dta")))
  stopifnot(nrow(d) == 652944)
  d$dyadid <- d$statea * 1000 + d$stateb
  d <- d[order(d$dyadid, d$year), ]
  # btscs fatal year dyadid, gen(py): years since last event within dyad; 0 at the dyad's first obs,
  # reset to 0 in the year after an event (event==1 at t-1). Missing event treated as non-event (Stata: .==1 false).
  ev <- !is.na(d$fatal) & d$fatal == 1
  first <- c(TRUE, d$dyadid[-1] != d$dyadid[-nrow(d)])
  py <- numeric(nrow(d))
  for (i in seq_len(nrow(d))) {
    if (first[i]) py[i] <- 0 else if (ev[i - 1]) py[i] <- 0 else py[i] <- py[i - 1] + 1
  }
  d$py <- py
  d
}

add_splines <- function(d, knots, Tmax) {
  # BKT/btscs-type cubic spline basis: s_i = (py-k_i)^3_+ - ((Tmax-k_i)/Tmax) * py^3  (form verified on btscs output in Maoz 2006 data)
  for (i in seq_along(knots)) d[[paste0("_spline", i)]] <- pmax(d$py - knots[i], 0)^3 - (Tmax - knots[i]) / Tmax * d$py^3
  d
}

build65 <- function(d) {
  d <- d[d$year > 1949 & d$year < 2001, ]
  sm <- function(a, b) ifelse(is.na(a) | is.na(b), NA, pmin(a, b))
  lg <- function(a, b) ifelse(is.na(a) | is.na(b), NA, pmax(a, b))
  d$smldem <- sm(d$polity_a, d$polity_b); d$lrgdem <- lg(d$polity_a, d$polity_b)
  d$lrgcap <- lg(d$cap_a, d$cap_b); d$pwin_cap <- d$lrgcap / (d$cap_a + d$cap_b)
  d$lnlrgcap <- log(d$lrgcap)
  d$lnsmlrgdp <- log(sm(d$rgdp_a, d$rgdp_b)); d$lnlrgrgdp <- log(lg(d$rgdp_a, d$rgdp_b))
  d$lnsmlpop <- log(sm(d$tpop_a, d$tpop_b)); d$lnlrgpop <- log(lg(d$tpop_a, d$tpop_b))
  d$logdstab <- log(d$distance)
  d$dircont <- ifelse(d$contig < 6, 1, ifelse(d$contig == 6, 0, NA))
  # Stata: jtdem=1 if smldem>6 & lrgdem>6 -- missing counts as +inf, so missing-polity dyads get jtdem=1 (authors' coding kept)
  big <- function(x) ifelse(is.na(x), Inf, x)
  d$jtdem <- as.numeric(big(d$smldem) > 6 & big(d$lrgdem) > 6)
  # lnNt: number of US dyads in year + 1 (after 1950-2000 restriction), ln, centred, 0 if contiguous
  nt <- d %>% filter(statea == 2) %>% count(year, name = "Nt") %>% mutate(Nt = Nt + 1)
  stopifnot(all(1950:2000 %in% nt$year))
  d$lnNt <- log(nt$Nt[match(d$year, nt$year)]) - 4.3175
  d$lnNt[d$dircont == 1] <- 0   # Stata: replace lnNt=0 if dircont==1 (dircont missing never ==1)
  d
}
