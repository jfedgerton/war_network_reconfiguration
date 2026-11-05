# X077 Hwang 2010 II "Power, Preferences, and Multiple Levels of Interstate Conflict".
# Headline = Model 1.1 (authors' do file): probit cwmid lnpower_1 I_1 Ipow democracy contig peaceyrs _spline1-3, robust   (MID onset, all dyads)
#   (Model 1.2 is a heckprob for war -> selection model, INFEASIBLE under protocol; not used.)
# Key term fixed from abstract before estimating: Ipow = I_1 x lnpower_1 ("conditional effect of preferences on the relationship between power and conflict").
# CONTROL: democracy (dyadic democracy control; variable unlabelled in deposit).
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p4/P4_common.R"))
WO <- file.path(OUT, "work_p4/X077"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work_p4/X077"); sink(file.path(WO, "log_X077.txt"), split = TRUE)
d <- as.data.frame(haven::zap_labels(haven::read_dta(file.path(W, "hwang.tab"))))
names(d) <- sub("^_spline", "spline", names(d))
dk <- duplicated(d[, c("ccode1","ccode2","year")]) | duplicated(d[, c("ccode1","ccode2","year")], fromLast = TRUE)
cat("rows in duplicated (ccode1,ccode2,year) keys:", sum(dk), " exact duplicate rows:", sum(duplicated(d)), " ccode1<ccode2 share:", mean(d$ccode1 < d$ccode2), "\n")
print(head(d[dk, ], 6))
stopifnot(all(d$cwmid %in% c(0, 1, NA)))   # duplicates kept as in the authors' estimation data
cat("N:", nrow(d), " years:", range(d$year), " Ipow check:", max(abs(d$Ipow - d$I_1 * d$lnpower_1), na.rm = TRUE), "\n")
print(summary(d$democracy))
f0 <- cwmid ~ lnpower_1 + I_1 + Ipow + democracy + contig + peaceyrs + spline1 + spline2 + spline3
K <- c("Ipow", "democracy")
PR <- binomial("probit")
r <- fit_keys(f0, d, K, "hetero", file.path(WO, "cache_repro.rds"), PR)
ds <- common_sample(d, f0, character(0), "ccode1", "ccode2", "year")
M <- run_six(f0, ds, K, "hetero", WO, "", PR)
# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "X077", dat = ds, y = "cwmid", tie = "democracy",
           ctrl = c("lnpower_1", "I_1", "Ipow", "contig", "peaceyrs", "spline1", "spline2", "spline3"), fe = NULL,
           fit = function(dat, add) { r <- fit_keys(if (length(add)) update(f0, as.formula(paste(". ~ . +", paste(add, collapse = " + ")))) else f0, dat, "democracy", "hetero", NULL, PR)$k$democracy
             c(b = unname(r[["b"]]), se = unname(r[["se"]]), p = unname(r[["p"]])) })
}
sc <- "Hwang 2010 II"; est <- "probit, robust (HC) SE (fixest feglm probit)"
hc <- "Y (peace years + 3 cubic splines only)"
nt <- paste0("Model 1.1 (MID onset, all dyads). Paper table inaccessible (T&F paywall; deposit has no log) -> UNVERIFIED. N repro ", r$n,
  "; deposit contains duplicated dyad-year keys (kept, as in authors' data); common sample drops ", r$n - M$m0$n, " obs.")
write_p4("X077", sc, "TIE", "KEY", "Ipow", "abstract: 'conditional effect of preferences on the relationship between power and conflict'",
  "UNVERIFIED", NA, r$k$Ipow[["b"]], est, r$n, M, hc, paste(nt, "Ipow = I_1 (lagged preference similarity) x lnpower_1 (ln power ratio), product term."))
write_p4("X077", sc, "DEMOCRACY", "CONTROL", "democracy", "control in Model 1.1 (dyadic democracy; unlabelled in deposit)",
  "UNVERIFIED", NA, r$k$democracy[["b"]], est, r$n, M, hc, paste(nt, "bD NOT a placebo for this term."))
sink()
