# Rerun 54: Bearce & Omori 2005 JPR -- Model 2 (full institutional-mechanism model):
#   logit midon trade_interdependence democracy_low democracy_high joint_alliance capabilities_ratio hegemony contiguity joint_pta
#         economic_integration nested_military high_organs midonspl midonsp1 midonsp2 midonsp3, cluster(dyadid)
# Headline tie = high_organs (the abstract's supported mechanism: PTAs with more high-level organs reduce MID onset).
# joint_pta (commercial-institution membership) is also reported in the log for reference.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work/rerun_common_k6.R"))
W <- file.path(ROOT, "work/54")
WO <- file.path(OUT, "work/54"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)
sink(file.path(WO, "log_54.txt"), split = TRUE)

# 1. Authors' analysis data (PRIO JPR 2005/42/6 zip -> file47159_bearceomorijpr.dta)
d <- as.data.frame(zap_labels(read_dta(file.path(W, "file47159_bearceomorijpr.dta"))))
stopifnot(nrow(d) == 20990, all(c("state1","state2","year","midon","high_organs","dyadid") %in% names(d)))
BASE <- "midon ~ trade_interdependence + democracy_low + democracy_high + joint_alliance + capabilities_ratio + hegemony + contiguity + joint_pta + economic_integration + nested_military + high_organs + midonspl + midonsp1 + midonsp2 + midonsp3"
MV <- c("midon","trade_interdependence","democracy_low","democracy_high","joint_alliance","capabilities_ratio","hegemony","contiguity","joint_pta","economic_integration","nested_military","high_organs","midonspl","midonsp1","midonsp2","midonsp3","dyadid")
fitter <- function(add, dat) fit_tie(as.formula(paste0(BASE, add)), dat, link = "logit", tie = "high_organs", cluster = "dyadid")

# 2. Reproduce on the full published sample (verbatim do-file model 2). No author log in deposit; paper table not retrievable.
cat("\n==== Reproduction: Model 2 (full sample) ====\n")
d_full <- d[complete.cases(d[, MV]), ]
rp <- fitter("", d_full)
cat("\nReference: joint_pta-only Model 1 on full sample\n")
rp1 <- fit_tie(midon ~ trade_interdependence + democracy_low + democracy_high + joint_alliance + capabilities_ratio + hegemony + contiguity + joint_pta + midonspl + midonsp1 + midonsp2 + midonsp3,
               d_full, tie = "joint_pta", cluster = "dyadid")
b_pub <- NA
status <- "UNVERIFIED"
n_orig <- rp$n

# 3. Merge history (undirected state1/state2; 1951-1985 inside GML range)
d <- add_history(d, "state1", "state2", "year")
d <- restrict_sample(d, MV)

# 4. Six models on one sample
R <- run_six(fitter, d)

# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "54", dat = d, y = "midon", tie = "high_organs",
           ctrl = c("trade_interdependence", "democracy_low", "democracy_high", "joint_alliance", "capabilities_ratio", "hegemony", "contiguity",
                    "joint_pta", "economic_integration", "nested_military", "midonspl", "midonsp1", "midonsp2", "midonsp3"), fe = NULL,
           fit = function(dat, add) { r <- fitter(if (length(add)) paste0(" + ", paste(add, collapse = " + ")) else "", dat); c(b = r$b, se = r$se, p = r$p) })
}

# 5. Results row
write_row(54, "Bearce & Omori 2005", status, b_pub, rp$b, "logit, SE clustered on dyadid (glm + vcovCL HC0*G/(G-1))", "high_organs",
  n_orig, R$m0$n, R, ahc = "Y",
  notes = sprintf("Model 2 (econ integration + nested military + high organs). Tie = high_organs (headline mechanism); Model 1 joint_pta b=%.3f (p=%.3f) in log. b_pub not obtainable (paper closed, no log); do-file run verbatim. Author controls peace-years + 3 splines. Dropped %d obs vs orig N=%d.", rp1$b, rp1$p, n_orig - R$m0$n, n_orig))
sink()
