# Rerun ADD_T03: Leeds & Mattes 2007 CMPS (10.1080/07388940701473054) -- ATOP alliance commitment types in a replication of
#   Russett & Oneal (2001); politically relevant dyads, MID onset (dispute1), 1885-1991 (covariate years)
# Headline: Table 6, full sample, commitment-type model (author do-file LeedsMattesCMPS2007rep.do, third xtgee):
#   xtgee dispute1 smldmat smldep smigoabi lcaprat nonaggonly notnonagg noncontg logdstab minrpwrs,
#         family(binomial) link(logit) corr(ar1) force robust      (iis dyad, tis year; `lcaprat` abbreviates lcaprat2)
# Key term = notnonagg (shared ATOP alliance with commitments beyond non-aggression); nonaggonly (non-aggression-only pact) is
#   carried as a non-headline row. Single additive key term in an additive model -> sens_run (scope full).
# Published Table 6 coefficients not obtainable here (T&F paywall); b_pub = our reproduction of the authors' own model on the
#   authors' data (published_check = unverified). The Cold War / pre-Cold War split models of Table 6 are printed for reference only.
# Estimator port: Stata xtgee corr(ar1) force takes lags from data ORDER within dyad; geeglm with waves = within-dyad order and
#   single-observation dyads dropped, exactly as in legacy/work/20/rerun_est_20.R (Oneal & Russett 2001, same model family).
# History is attached to the outcome year dvyear = year + 1, as in legacy/work/20 (Russett-Oneal files: dispute1 refers to t+1).
# Duplication note: same Russett-Oneal base model and politically relevant dyad sample as registry id 20 (Oneal & Russett 2001,
#   key term smldep, rebuild); here the key term is the ATOP alliance type, which id 20 does not estimate.
# Data: authors' replication package (Harvard Dataverse doi:10.7910/DVN/L4STW5), in _incoming/LeedsMattes2007.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
suppressPackageStartupMessages({library(geepack); library(dplyr); library(haven)})
set.seed(123)
WO <- file.path(OUT, "work_add20261110/ADD_T03"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)
sink(file.path(WO, "log_ADD_T03.txt"), split = TRUE)

TIE <- "notnonagg"; TIE2 <- "nonaggonly"
CTRL <- c("smldmat", "smldep", "smigoabi", "lcaprat2", "nonaggonly", "noncontg", "logdstab", "minrpwrs")
f0 <- as.formula(paste("dispute1 ~", paste(c("smldmat", "smldep", "smigoabi", "lcaprat2", "nonaggonly", "notnonagg",
                                             "noncontg", "logdstab", "minrpwrs"), collapse = " + ")))
rd <- function(f) { x <- as.data.frame(zap_labels(read_dta(file.path(ROOT, "_incoming/LeedsMattes2007", f)))); x[order(x$dyad, x$year), ] }
d <- rd("LeedsMattesROrep.dta")
cat("rows in author file:", nrow(d), " dyads:", length(unique(d$dyad)), "\n")
stopifnot(nrow(d) == 39996, all(all.vars(f0) %in% names(d)))

fit_gee <- function(f, data, keys = c(TIE, TIE2)) {
  data <- data[complete.cases(data[, c(all.vars(f), "dyad", "year")]), ]; data <- data[order(data$dyad, data$year), ]
  tb <- table(data$dyad); ns <- sum(tb == 1); data <- data[data$dyad %in% as.numeric(names(tb)[tb > 1]), ]
  data$w <- ave(data$year, data$dyad, FUN = seq_along)                     # order-based AR(1) (Stata force)
  m <- geeglm(f, data = data, id = dyad, waves = w, family = binomial, corstr = "ar1")
  ct <- summary(m)$coefficients
  cat("\n----", deparse(f, width.cutoff = 500), "| N =", nrow(data), "| singletons dropped", ns, "| alpha", round(unlist(m$geese$alpha), 4), "\n")
  print(round(ct, 5))
  setNames(lapply(keys, function(k) list(b = ct[k, 1], se = ct[k, 2], p = ct[k, 4], n = nrow(data))), keys)
}

# 1. Reproduce Table 6 (full sample, commitment-type model) with the authors' data; split-sample models for reference
cat("\n==== Reproduction (author sample) ====\n")
r_pub <- fit_gee(f0, d)
cw  <- fit_gee(f0, rd("LeedsMattesROrepCW.dta"))
pcw <- fit_gee(f0, rd("LeedsMattesROrepPCW.dta"))
cat(sprintf("REPRO full: %s %.4f (se %.4f, p %.4f) %s %.4f N %d | Cold War %s %.4f (p %.4f) | pre-CW %s %.4f (p %.4f)\n",
            TIE, r_pub[[TIE]]$b, r_pub[[TIE]]$se, r_pub[[TIE]]$p, TIE2, r_pub[[TIE2]]$b, r_pub[[TIE]]$n,
            TIE, cw[[TIE]]$b, cw[[TIE]]$p, TIE, pcw[[TIE]]$b, pcw[[TIE]]$p))
b_pub <- r_pub[[TIE]]$b; b_pub2 <- r_pub[[TIE2]]$b          # unverified (published table not obtained)
stopifnot(is.finite(b_pub), is.finite(b_pub2))
n_orig <- r_pub[[TIE]]$n

# 2. Merge history on the outcome year; restrict to 1822-2011; complete cases on every model column
source(file.path(CODE, "hist/add_history.R"))
d$dvyear <- d$year + 1
d <- add_history(d, "statea", "stateb", "dvyear")
s <- subset(d, dvyear >= 1822 & dvyear <= 2011)
s <- s[complete.cases(s[, c(all.vars(f0), "dyad", "hist_dyad_w15", "hist_mono_w15", "hist_dyad_w26", "hist_mono_w26",
                            "bm_cinc_w15", "bm_democ_w15")]), ]
s <- s[ave(s$year, s$dyad, FUN = length) > 1, ]                           # model sample (fit_gee drops singleton dyads)
cat("\nRows on common sample:", nrow(s), "(author model N", n_orig, ")\n")
stopifnot(nrow(s) > 0.80 * n_orig)

# 3. Models on the same sample
R <- list(m0 = fit_gee(f0, s), mD = fit_gee(update(f0, . ~ . + hist_dyad_w15), s),
          mDM = fit_gee(update(f0, . ~ . + hist_dyad_w15 + hist_mono_w15), s), m26 = fit_gee(update(f0, . ~ . + hist_dyad_w26 + hist_mono_w26), s),
          bC = fit_gee(update(f0, . ~ . + bm_cinc_w15), s), bD = fit_gee(update(f0, . ~ . + bm_democ_w15), s))
stopifnot(length(unique(sapply(R, function(r) r[[TIE]]$n))) == 1, R$m0[[TIE]]$n == nrow(s))
# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "ADD_T03", dat = s, y = "dispute1", tie = TIE, ctrl = CTRL, fe = NULL,
           fit = function(dat, add) { r <- fit_gee(if (length(add)) update(f0, as.formula(paste(". ~ . +", paste(add, collapse = " + ")))) else f0, dat)[[TIE]]
             c(b = unname(r$b), se = unname(r$se), p = unname(r$p)) })
}

# 4. Results, one row per key term, same columns as other studies
row_for <- function(key, bp) {
  g <- function(m) R[[m]][[key]]; b0 <- g("m0")$b; pct <- function(m) 100 * (g(m)$b - b0) / abs(b0)
  data.frame(cid = "ADD_T03", short_cite = "Leeds & Mattes 2007", key_term = key, reproduced = "REPRODUCED", b_pub = bp, b_repro = bp,
    estimator = "GEE logit, AR(1) (order-based, Stata force), robust SE by dyad (geeglm)", tie_var = key,
    n_orig = n_orig, n_used = g("m0")$n, b0 = b0, se0 = g("m0")$se, p0 = g("m0")$p,
    b_dm = g("mDM")$b, se_dm = g("mDM")$se, p_dm = g("mDM")$p, pct_dm = pct("mDM"), b_d = g("mD")$b, pct_d = pct("mD"),
    b_26 = g("m26")$b, se_26 = g("m26")$se, p_26 = g("m26")$p, pct_26 = pct("m26"), pct_bC = pct("bC"), pct_bD = pct("bD"),
    author_has_history_control = "N",
    notes = sprintf("Table 6 full-sample commitment-type model; published table not obtained (b_pub = reproduction, unverified); no peace years (C_none); dropped %d obs vs author sample.", n_orig - g("m0")$n))
}
out <- rbind(row_for(TIE, b_pub), row_for(TIE2, b_pub2))
stopifnot(all(is.finite(as.matrix(out[, c("b0", "se0", "p0", "b_dm", "p_dm", "b_26", "p_26", "pct_bC", "pct_bD")]))))
print(Filter(is.numeric, out))
dir.create(file.path(OUT, "results"), recursive = TRUE, showWarnings = FALSE)
write.csv(out, file.path(OUT, "results/ADD_T03.csv"), row.names = FALSE)
sink()
