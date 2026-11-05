# X068 Chatagnier & Lim 2020 JPR "Does the WTO exacerbate international conflict?" (10.1177/0022343320960203)
# Data: PRIO replication zip (wto_conflict_replication_data.dta); spec from authors' R log (wto_conflict_replication.txt).
# Headline: Table 2 col 1 "All MIDs": OLS (LPM, DV x100) mid_onset ~ joint + treat.wto + did.wto + majpow + contig + lodem + fta + alliance + Mrat + t + t2 + t3
# KEY terms (fixed a priori from abstract): did.wto ("judicialization erases these benefits": joint x post-1995 DSU) and joint
#   ("membership in a trade institution facilitates peaceful interaction"). CONTROL: lodem (lower Polity, weak-link democracy).
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p4/common_p4_b04.R"))
suppressPackageStartupMessages(library(readstata13))
WO <- file.path(OUT, "work_p4/X068"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work_p4/X068"); sink(file.path(WO, "log_X068.txt"), split = TRUE)
wto <- read.dta13(file.path(W, "wto_conflict_replication_data.dta"), convert.factors = FALSE)
cat("raw rows", nrow(wto), "\n")
# ---- authors' sample construction (verbatim logic) ----
s0 <- 1980; s1 <- 2009
wto <- wto[wto$polrel == 1, ]
jl <- ifelse((!is.na(wto$join_1) & wto$join_1 > s0 & wto$join_1 <= s1) | (!is.na(wto$join_2) & wto$join_2 > s0 & wto$join_2 <= s1), 1, 0)
wto <- wto[jl == 0, ]
cc <- unique(c(wto$ccode1, wto$ccode2))
ny <- sapply(cc, function(i) length(unique(wto$year[wto$ccode1 == i | wto$ccode2 == i])))
full <- cc[ny >= s1 - s0]
wto <- wto[wto$ccode1 %in% full & wto$ccode2 %in% full & wto$year >= s0 & wto$year <= s1, ]
wto$treat.wto <- as.numeric(wto$year >= 1995); wto$did.wto <- wto$treat.wto * wto$joint
wto$mid_onset <- wto$mid_onset * 100
wto$t <- wto$t_mid; wto$t2 <- wto$t_mid^2; wto$t3 <- wto$t_mid^3
cat("analysis rows", nrow(wto), " directed?", any(wto$ccode1 > wto$ccode2), "\n")
f0 <- mid_onset ~ joint + treat.wto + did.wto + majpow + contig + lodem + fta + alliance + Mrat + t + t2 + t3
terms <- c("did.wto", "joint", "lodem")
# ---- reproduction on full sample ----
r_iid <- ff_ols_iid(f0, wto); print(round(r_iid, 4))
b_pub <- c(did.wto = 1.78, joint = -1.68, lodem = -0.06); p_pub_sig <- c(did.wto = TRUE, joint = TRUE, lodem = TRUE)
rep <- sapply(terms, function(t) { ok <- abs(r_iid[t, 1] - b_pub[t]) / abs(b_pub[t]) <= 0.10 || abs(r_iid[t, 1] - b_pub[t]) <= 0.005   # table rounded to 2 dp
  if (ok && (r_iid[t, 3] < .05) == p_pub_sig[t]) "REPRODUCED" else "NOT_REPRODUCED" })
print(rep); n_orig <- attr(r_iid, "n")
# ---- history battery ----
ds <- common_sample(wto, f0, character(0), "ccode1", "ccode2", "year")
M <- run_six(f0, ds, terms, ff_ols_iid, WO)
# ---- sensitivity hook (split only): runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  # calls ff_ols_iid directly (what fit_keys runs): fit_keys via run_six would read the cache_m0/mD/mDM files
  sens_split(id = "X068", dat = ds,
             fit = function(dat, add) { ct <- ff_ols_iid(if (length(add)) update(f0, as.formula(paste(". ~ . +", paste(add, collapse = " + ")))) else f0, dat)
               c(b = unname(ct["did.wto", 1]), se = unname(ct["did.wto", 2]), p = unname(ct["did.wto", 3])) })
}
meta <- list(xid = "X068", short_cite = "Chatagnier & Lim 2020", journal = "JPR", year = 2020, doi = "10.1177/0022343320960203")
est <- "OLS linear probability (DV x100), conventional SEs (authors' lm; texreg 'robust' arg has no effect on lm)"
nt <- "Authors' R-log spec and sample (PR dyads, 1980-2009, no late WTO joiners, full-period states). All years inside 1822-2011; 6 obs dropped for missing lagged CINC/Polity benchmarks (common sample 13,164). Published Table 2 values rounded to 2 dp."
write_p4(meta, "did.wto", "TIE", "KEY", "\"membership ... facilitates peaceful interaction, but that judicialization erases these benefits\"", rep["did.wto"],
  b_pub["did.wto"], r_iid["did.wto", 1], est, n_orig, M, "N (peace-years cubic only)", paste(nt, "Key = Joint x WTO-era (DiD) product term."))
write_p4(meta, "joint", "TIE", "KEY", "\"membership in a trade institution facilitates peaceful interaction\"", rep["joint"],
  b_pub["joint"], r_iid["joint", 1], est, n_orig, M, "N (peace-years cubic only)", paste(nt, "Joint GATT/WTO membership constituent (pre-1995 effect)."))
write_p4(meta, "lodem", "DEMOCRACY", "CONTROL", "control: Polity (Lower) in headline model", rep["lodem"],
  b_pub["lodem"], r_iid["lodem", 1], est, n_orig, M, "N (peace-years cubic only)", paste(nt, "Weak-link democracy control."))
sink()
