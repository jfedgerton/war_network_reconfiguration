# Shared helpers for the 2026-10-11 protocol-B rebuilds (work_p9 P9_066+). Sourced after work_p8/common_b3.R.
# Builds dyad-years from {peacesciencer} (package data; version printed) plus the COW alliance v4.1 dyad-year file already on disk
# (work_p8/P8_049, listed in data_sources.csv). Sequential, fail-loud.
suppressPackageStartupMessages(library(peacesciencer))
cat("peacesciencer", as.character(packageVersion("peacesciencer")), "\n")

# non-directed (ccode1 < ccode2) or directed COW dyad-years with GML MIDs and the standard covariates
ps_dyads <- function(years, directed = FALSE, igos = FALSE, trade = TRUE, fpsim = FALSE) {
  d <- create_dyadyears(system = "cow", directed = directed, subset_years = years, mry = FALSE)
  d <- suppressMessages(d %>% add_gml_mids() %>% add_contiguity() %>% add_capital_distance() %>% add_nmc() %>%
         add_cow_majors() %>% add_democracy() %>% add_atop_alliance() %>% add_sim_gdp_pop())
  if (trade) d <- suppressMessages(d %>% add_cow_trade())
  if (igos) d <- suppressMessages(d %>% add_igos())
  if (fpsim) d <- suppressMessages(d %>% add_fpsim())
  d <- as.data.frame(d)
  stopifnot(!anyDuplicated(d[, c("ccode1", "ccode2", "year")]))
  d
}

# COW formal alliance v4.1 dyad-years (non-directed; any of defense/neutrality/nonaggression/entente)
add_cow_ally <- function(d) {
  a <- read.csv(file.path(ROOT, "work_p8/P8_049/alliance_v4_1_by_dyad_yearly.csv"))
  a <- a %>% mutate(c1 = pmin(ccode1, ccode2), c2 = pmax(ccode1, ccode2)) %>% group_by(c1, c2, year) %>%
    summarise(cow_def = max(defense), cow_neu = max(neutrality), cow_nonagg = max(nonaggression), cow_ent = max(entente), .groups = "drop") %>%
    mutate(cow_any = 1L)
  n0 <- nrow(d)
  d <- d %>% mutate(c1 = pmin(ccode1, ccode2), c2 = pmax(ccode1, ccode2)) %>% left_join(a, by = c("c1", "c2", "year")) %>%
    mutate(across(c(cow_def, cow_neu, cow_nonagg, cow_ent, cow_any), ~ coalesce(as.integer(.x), 0L))) %>% select(-c1, -c2)
  stopifnot(nrow(d) == n0); as.data.frame(d)
}

# peace years since the last onset within unit (data sorted by unit, year), btscs rule; splines with quartile knots (btscs)
py_calc <- function(d, unit, ev) {
  o <- order(d[[unit]], d$year); d <- d[o, ]
  u <- d[[unit]]; e <- d[[ev]]; first <- c(TRUE, u[-1] != u[-length(u)]); gap <- c(TRUE, diff(d$year) != 1)
  py <- numeric(nrow(d)); cur <- -1; last <- FALSE
  stopifnot(!is.null(e), !anyNA(first), !anyNA(gap))
  for (i in seq_len(nrow(d))) { if (first[i] || gap[i]) { cur <- -1; last <- FALSE }
    cur <- if (cur < 0 || last) 0 else cur + 1; py[i] <- cur; last <- !is.na(e[i]) && e[i] == 1 }
  d$py <- py; d
}
btscs_splines <- function(d, py = "py") {
  kn <- quantile(d[[py]], c(.25, .5, .75), type = 2, na.rm = TRUE); Tm <- max(d[[py]], na.rm = TRUE)
  cat("btscs knots:", kn, " T:", Tm, "\n")
  for (i in 1:3) d[[paste0("spl", i)]] <- (pmax(d[[py]] - kn[i], 0)^3 - (Tm - kn[i]) / Tm * d[[py]]^3) / 1000
  d
}
# politically relevant: contiguous (land or <= 150 miles water, conttype 1-4) or either a major power
pr_flag <- function(d) as.integer(d$conttype %in% 1:4 | d$cowmaj1 == 1 | d$cowmaj2 == 1)

# Polity IV 2014 (p4v2014.csv, work_p8/P8_049) mapped to COW codes (same mapping as P9_061); returns ccode, year, polity2, durable, democ, autoc
polity_cow <- function(vars = c("polity2", "durable")) {
  p <- read.csv(file.path(ROOT, "work_p8/P8_049/p4v2014.csv"))[, c("ccode", "year", vars)]
  for (v in vars) p[[v]][p[[v]] %in% c(-66, -77, -88)] <- NA
  p <- p %>% mutate(ccode = case_when(ccode == 364 ~ 365L, ccode == 769 ~ 770L, ccode == 818 ~ 816L, ccode == 342 & year < 1921 ~ 345L,
                                      TRUE ~ as.integer(ccode))) %>%
    group_by(ccode, year) %>% summarise(across(all_of(vars), ~ if (all(is.na(.x))) NA_real_ else max(.x, na.rm = TRUE)), .groups = "drop")
  as.data.frame(p)
}

# relogit fitter returning the term list (King-Zeng bias correction; work/relogit_kz_ea1.R)
fit_relogit <- function(dv, rhs, data, terms, cluster = NULL, vcov_type = NULL) {
  o <- relogit_kz(as.formula(paste(dv, "~", paste(rhs, collapse = " + "))), data, cluster = cluster, vcov_type = vcov_type)
  cat("\n---- relogit (KZ-corrected)", dv, "~", paste(rhs, collapse = " + "), "| N =", attr(o, "n"), "\n")
  print(round(o[!grepl("__unc$", rownames(o)), ], 6))
  res <- lapply(setNames(terms, terms), function(t) o[t, c("b", "se", "p")]); attr(res, "n") <- attr(o, "n"); res
}

# Oneal-Russett/Dafoe all-dyad file (work_demcap/D1/dor.dta, 1816-2000) joined with the derived trade/GDP file (work_p4/X074/or_trade.rds)
or_dyads <- function(years) {
  x <- as.data.frame(haven::zap_labels(haven::read_dta(file.path(ROOT, "work_demcap/D1/dor.dta"))))
  x <- x[x$year %in% years, c("statea", "stateb", "year", "mzmid1_orig", "mzfatal1_orig", "mzcowwar1_orig", "mzmid1_all", "mzfatal1_all", "democa", "autoca", "democb", "autocb")]
  t <- readRDS(file.path(ROOT, "work_p4/X074/or_trade.rds")); t <- t[t$year %in% years, ]
  stopifnot(!anyDuplicated(x[, c("statea", "stateb", "year")]), !anyDuplicated(t[, c("statea", "stateb", "year")]))
  x <- x[, c("statea", "stateb", "year", setdiff(names(x), names(t)))]
  d <- merge(t, x, by = c("statea", "stateb", "year"), all.x = TRUE)
  stopifnot(nrow(d) == nrow(t), all(d$statea < d$stateb))
  d
}

# COW national trade 3.0 totals (imports + exports, current $m) by ccode-year (work_p8/P8_049)
cow_nat_trade <- function() {
  x <- read.csv(file.path(ROOT, "work_p8/P8_049/national_trade_3_0.csv"))[, c("ccode", "year", "imports", "exports")]
  x$imports[x$imports < 0] <- NA; x$exports[x$exports < 0] <- NA
  stopifnot(!anyDuplicated(x[, c("ccode", "year")]))
  x$tot <- x$imports + x$exports; x[, c("ccode", "year", "tot")]
}

# GWF All Political Regimes 1.2 (work_p9/gwf/GWF_AllPoliticalRegimes.dta): personalist incl. hybrids, prior regime
gwf_cy <- function() {
  g <- as.data.frame(haven::zap_labels(haven::read_dta(file.path(ROOT, "work_p9/gwf/GWF_AllPoliticalRegimes.dta"))))
  g <- g %>% transmute(ccode = as.integer(cowcode), year = as.integer(year), gwf_regimetype, gwf_prior, gwf_nonautocracy, gwf_casename,
    autoc = as.integer(!is.na(gwf_regimetype) & gwf_regimetype != ""), pers = as.integer(grepl("personal", gwf_regimetype)))
  stopifnot(!anyDuplicated(g[, c("ccode", "year")])); g
}
# V-Dem v6.2 (shortened file staged for P9_062): v2x_libdem by COW code-year
vdem_libdem <- function() {
  v <- as.data.frame(haven::zap_labels(haven::read_dta(file.path(ROOT, "work_p9/P9_062/V-Dem-DS-CY-v6.2_shortened.dta"),
        col_select = c("COWcode", "year", "v2x_libdem"))))
  v <- v[!is.na(v$COWcode), ] %>% group_by(ccode = as.integer(COWcode), year = as.integer(year)) %>%
    summarise(libdem = if (all(is.na(v2x_libdem))) NA_real_ else max(v2x_libdem, na.rm = TRUE), .groups = "drop")
  as.data.frame(v)
}
