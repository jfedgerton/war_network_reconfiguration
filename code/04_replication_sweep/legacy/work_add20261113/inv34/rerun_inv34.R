# Rerun inv34: Suhas 2023 ISQ (10.1093/isq/sqad070), "How Alliances Shape Rivalries" -- only the weaker rival allied ->
#   hazard that the weaker rival initiates a MID against the stronger rival, asymmetric strategic rivalries 1816-2014
# Headline: Table 5.0 model 1, author do-file 02-19-2023-main-analysis_do.do:
#   stset dur, id(rcode) failure(weakini)
#   stcox i.onlyweakallied onlystrongallied wsallied multino_mutual otherboth spatial positional demdem majmaj conttype coldwar,
#         cluster(rcode)
#   author log: HR 3.880095 (robust SE 1.639346), 236 subjects, 66 failures, N 5,239 -> b = log(3.880095) = 1.355859
# Key term = onlyweakallied. Single additive term -> sens_run (scope full).
# Estimator: survival::coxph, Breslow ties (Stata default), robust SE clustered on rcode. Stata's stset is rebuilt by hand:
#   interval (dur[_n-1], dur] within rcode (first interval starts at 0); rows after a rivalry's first failure are dropped.
# Data: author's ISQ supplemental files (sqad070_supplemental_files.zip, OUP), in _incoming/Suhas2023.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
suppressPackageStartupMessages({library(survival); library(dplyr); library(haven)})
set.seed(123)
WO <- file.path(OUT, "work_add20261113/inv34"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)
sink(file.path(WO, "log_inv34.txt"), split = TRUE)

# 1. Load the author's analysis data and rebuild stset
TIE  <- "onlyweakallied"
CTRL <- c("onlystrongallied", "wsallied", "multino_mutual", "otherboth", "spatial", "positional", "demdem", "majmaj",
          "conttype", "coldwar")
d <- as.data.frame(zap_labels(read_dta(file.path(ROOT, "_incoming/Suhas2023/02-19-2023-ISQ-Dataset.dta"))))
cat("rows in author file:", nrow(d), "\n")
stopifnot(nrow(d) == 7179, all(c("weakini", "dur", "rcode", "ccode1", "ccode2", "year", TIE, CTRL) %in% names(d)))
d <- d %>% arrange(rcode, dur) %>% group_by(rcode) %>%
  mutate(tstart = lag(dur, default = 0), fail = as.integer(!is.na(weakini) & weakini != 0),
         prior_fail = cumsum(lag(fail, default = 0L))) %>%
  ungroup() %>% filter(prior_fail == 0) %>% as.data.frame()
stopifnot(all(d$dur > d$tstart))
cat("rows after stset (exit on first failure):", nrow(d), "| subjects", n_distinct(d$rcode), "| failures", sum(d$fail), "\n")
stopifnot(nrow(d) == 5239, n_distinct(d$rcode) == 236, sum(d$fail) == 66)   # author log

f_rhs <- paste(c(TIE, CTRL), collapse = " + ")
fit <- function(rhs_add, dat) {
  m <- coxph(as.formula(paste("Surv(tstart, dur, fail) ~", f_rhs, rhs_add)), data = dat, ties = "breslow", cluster = rcode)
  ct <- summary(m)$coefficients
  print(round(ct[TIE, , drop = FALSE], 5))
  list(b = ct[TIE, "coef"], se = ct[TIE, "robust se"], p = ct[TIE, "Pr(>|z|)"], n = m$n, ev = m$nevent)
}

# 2. Reproduce Table 5.0 model 1 on the author's full sample
cat("\n==== Reproduction (author sample) ====\n")
r_pub <- fit("", d)
b_pub <- log(3.880095); se_pub <- 1.639346 / 3.880095
cat(sprintf("REPRO %s: %.5f (se %.5f) N %d events %d | published %.5f (%.5f) N 5239 events 66\n",
            TIE, r_pub$b, r_pub$se, r_pub$n, r_pub$ev, b_pub, se_pub))
stopifnot(r_pub$n == 5239, r_pub$ev == 66, abs(r_pub$b - b_pub) < 1e-3, abs(r_pub$se - se_pub) / se_pub < 0.01)   # fail loudly
n_orig <- r_pub$n

# 3. Merge history (undirected rival dyad, outcome year = year); restrict to 1822-2011
source(file.path(CODE, "hist/add_history.R"))
d <- add_history(d, "ccode1", "ccode2", "year")
d <- filter(d, year >= 1822, year <= 2011)
d <- d[complete.cases(d[, c("fail", "tstart", "dur", "rcode", TIE, CTRL, "hist_dyad_w15", "hist_mono_w15", "hist_dyad_w26",
                            "hist_mono_w26", "bm_cinc_w15", "bm_democ_w15")]), ]
cat("\nRows after year restriction + complete cases:", nrow(d), "(author model N", n_orig, ") | failures", sum(d$fail), "\n")
# Common-sample guard: PI exception 2026-11-13 at 0.75 (as P9_092). The rivalry panel runs to 2014 (history ends 2011) and the
# capability-ratio benchmark is missing for 784 nineteenth-century rows; the common sample keeps 77.4% of rows and 63 of 66 failures.
stopifnot(nrow(d) > 0.75 * n_orig)

# 4. Models on the same sample
cat("\n==== m0 ====\n");  r0  <- fit("", d)
cat("\n==== mD ====\n");  rD  <- fit("+ hist_dyad_w15", d)
cat("\n==== mDM ====\n"); rDM <- fit("+ hist_dyad_w15 + hist_mono_w15", d)
cat("\n==== m26 ====\n"); r26 <- fit("+ hist_dyad_w26 + hist_mono_w26", d)
cat("\n==== bC ====\n");  rC  <- fit("+ bm_cinc_w15", d)
cat("\n==== bD ====\n");  rBD <- fit("+ bm_democ_w15", d)
stopifnot(length(unique(c(r0$n, rD$n, rDM$n, r26$n, rC$n, rBD$n))) == 1, r0$n == nrow(d))
# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "inv34", dat = d, y = "fail", tie = TIE, ctrl = CTRL, fe = NULL,
           fit = function(dat, add) { r <- fit(if (length(add)) paste("+", paste(add, collapse = " + ")) else "", dat)
             c(b = unname(r$b), se = unname(r$se), p = unname(r$p)) })
}
pct <- function(x) 100 * (x$b - r0$b) / abs(r0$b)

out <- data.frame(cid = "inv34", short_cite = "Suhas 2023", key_term = TIE, reproduced = "REPRODUCED", b_pub = b_pub, b_repro = r_pub$b,
  estimator = "Cox PH (Breslow), robust SE clustered by rivalry (survival::coxph)", tie_var = TIE,
  n_orig = n_orig, n_used = r0$n, b0 = r0$b, se0 = r0$se, p0 = r0$p, b_dm = rDM$b, se_dm = rDM$se, p_dm = rDM$p,
  pct_dm = pct(rDM), b_d = rD$b, pct_d = pct(rD), b_26 = r26$b, se_26 = r26$se, p_26 = r26$p, pct_26 = pct(r26),
  pct_bC = pct(rC), pct_bD = pct(rBD), author_has_history_control = "N",
  notes = sprintf("Table 5.0 model 1 (HR 3.880); dropped %d obs vs author sample (1822-2011 + history/benchmark completeness).", n_orig - r0$n))
stopifnot(all(is.finite(unlist(out[, c("b0", "se0", "p0", "b_dm", "p_dm", "b_26", "p_26", "pct_bC", "pct_bD")]))))
print(Filter(is.numeric, out))
dir.create(file.path(OUT, "results"), recursive = TRUE, showWarnings = FALSE)
write.csv(out, file.path(OUT, "results/inv34.csv"), row.names = FALSE)
sink()
