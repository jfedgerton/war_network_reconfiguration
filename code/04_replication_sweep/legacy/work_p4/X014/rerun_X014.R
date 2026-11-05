# X014 Oneal & Tir 2006 ISQ, "Does the Diversionary Use of Force Threaten the Democratic Peace?"
# Key term (decided from abstract BEFORE estimating): democracy x growth interaction (demXgrow).
#   Abstract: "economic conditions do affect the likelihood that a democracy ... will initiate a fatal militarized dispute"
# Headline: Table 1 col 3 (OT_TBL1.DO): directed dyads, logit finitsda (fatal MID initiation, t+1), dyad-clustered SE.
# b_pub from authors' log ot_tbl1.log: demXgrow = -13.9047 (p = .001), N = 509,604.
# No joint-democracy/capitalist term in the headline (democracy enters monadically: democ_st, democ_tgt) -> no CONTROL rows.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p4/b6/common_b6.R"))
WO <- file.path(OUT, "work_p4/X014"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work_p4/X014"); sink(file.path(WO, "log_X014.txt"), split = TRUE)
m <- list(xid = "X014", short_cite = "Oneal & Tir 2006", journal = "ISQ", year = 2006, doi = "10.1111/j.1468-2478.2006.00424.x")

if (!file.exists(file.path(WO, "OTDirect.dta"))) unzip(file.path(W, "OT.zip"), files = "OTDirect.dta", exdir = WO)  # port: zip member unzipped into OUT
stopifnot(file.exists(file.path(WO, "OTDirect.dta")))
d <- haven::read_dta(file.path(WO, "OTDirect.dta"), col_select = c(state, target, year, democ_st, democ_tgt, depend, allies, lncaprat,
  dircont, distance, numstate, finitsda, pyfatA, `_spline1`, `_spline2`, `_spline3`, leadid_st, democ_st_cg, grow_st, majpower))
d <- as.data.frame(haven::zap_labels(d)); names(d) <- sub("^_spline", "spline", names(d))
stopifnot(!anyDuplicated(d[, c("state", "target", "year")]))

# ---- replicate OT_TBL1.DO data steps, in order ----
d$lncaprat <- d$lncaprat + 10.3; d$lncapsq <- d$lncaprat^2
d <- d %>% arrange(state, target, year) %>% group_by(state, target) %>%
  mutate(grow_stldr = ifelse(!is.na(leadid_st) & leadid_st == dplyr::lag(leadid_st, 1) & leadid_st == dplyr::lag(leadid_st, 2), grow_st, NA)) %>%
  ungroup() %>% as.data.frame()
# Stata: a==b with both missing is TRUE (missing==missing); leadid missing -> treat as not same leader (conservative); check count
d$logdist <- log(d$distance)
d$systsize <- ifelse(!(d$dircont == 1 | d$majpower == 1), -log(d$numstate), NA); d$systsize[is.na(d$systsize)] <- 0
d$dyadid <- 1000 * d$state + d$target
d <- d[!(is.na(d$depend) | d$depend > 1), ]            # Stata 'drop if depend>1' also drops missing
d <- d[is.na(d$grow_stldr) | d$grow_stldr != 0, ]      # drop exactly-zero growth
d$demXgrow <- ifelse(!is.na(d$democ_st_cg) & d$democ_st_cg == 1 & !is.na(d$democ_st), d$grow_stldr, NA)
d$demXgrow[is.na(d$demXgrow) & !is.na(d$grow_stldr)] <- 0
stopifnot(all(d$finitsda %in% c(0, 1, NA)))

f0 <- finitsda ~ democ_st + democ_tgt + depend + grow_stldr + demXgrow + lncaprat + lncapsq + allies + dircont +
  logdist + majpower + systsize + pyfatA + spline1 + spline2 + spline3
K <- "demXgrow"
r <- fit_keys(f0, d, K, "~dyadid", "logit", file.path(WO, "cache_repro.rds")); n_orig <- r$n
cat("reproduction N:", n_orig, " (published 509,604)\n")
b_pub <- -13.9047
st <- repro_status(b_pub, TRUE, r$k$demXgrow); cat("status:", st, "\n")

# History anchored on the DV year (finitsda is measured in t+1)
d$dvyear <- d$year + 1
ds <- common_sample(d, f0, "dyadid", "state", "target", "dvyear"); rm(d); gc()
M <- run_six(f0, ds, K, "~dyadid", "logit", WO)
# ---- sensitivity hook (split only): runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  # fit_keys with cache = NULL (what run_six runs, uncached): run_six would read the cache_m0/mD/mDM files
  sens_split(id = "X014", dat = ds,
             fit = function(dat, add) { r <- fit_keys(if (length(add)) update(f0, as.formula(paste(". ~ . +", paste(add, collapse = " + ")))) else f0, dat, K, "~dyadid", "logit")$k$demXgrow
               c(b = unname(r[["b"]]), se = unname(r[["se"]]), p = unname(r[["p"]])) })
}
write_p4(m, "DEMOCRACY", "KEY", "demXgrow",
  "abstract: 'economic conditions do affect the likelihood that a democracy... will initiate a fatal militarized dispute'",
  st, b_pub, unname(r$k$demXgrow["b"]), "logit, dyad-clustered SE (directed dyads, fatal MID initiation t+1)", n_orig, M,
  "Y (pyfatA peace-years + 3 splines)",
  paste0("Table 1 col 3 (OT_TBL1.DO), b_pub from ot_tbl1.log (p=.001, N=509,604). Product term = growth x consolidated democracy (democ_st_cg). ",
         "History anchored on DV year t+1; sample 1921-2001 fully inside GML range; dropped vs complete-case: ", n_orig - M$m0$n,
         " (benchmark/history missingness). Democracy enters monadically (democ_st/democ_tgt), no joint-democracy term; bD is not a placebo for a democracy-family key term."))
sink()
