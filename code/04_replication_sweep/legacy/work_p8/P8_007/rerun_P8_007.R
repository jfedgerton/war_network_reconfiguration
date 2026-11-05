# P8_007 Boehmke, Chyzh & Thies 2016 PSRM "Addressing Endogeneity in Actor-Specific Network Measures" (10.1017/psrm.2015.34)
# Data: authors' Dataverse (doi:10.7910/DVN/BKPIEV) table1.dta = output of their table1.do (SEq instruments already built, 5 imputations).
# Paper (author PDF, olgachyzh.com): Table 1, logit, all dyads 1962-96, MID (dichmid20), IVs lagged one year, cubic time (t1..t3) + peace years,
#   plain logit SEs (Stata -logit-, no clustering), multiple-imputation combination over 5 draws (authors' micalc = Rubin's rules).
#   Model 3 "instrumented network": Ltrade_yhat_seq1..5  b = -4.372 (0.496); Model 4 "instrumented SEq": Ltrade_seq_hat1..5  b = -7.616 (0.642);
#   (Model 2 naive Ltradesecorr -0.762 (0.228)); N = 154,570.
# KEY TERMS (a priori): the paper's substantive claim is the IV-corrected trade-SEq effect; both IV models (3, 4) are presented as co-equal
#   corrections ("effect larger after IV correction"), so both are key terms. Logit ERGM columns are not rerun (ERGM = DEFERRED per user rule).
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/p8_common.R"))
suppressMessages(library(haven))
WD <- file.path(ROOT, "work_p8/P8_007"); W <- file.path(OUT, "work_p8/P8_007"); dir.create(W, recursive = TRUE, showWarnings = FALSE); sink(file.path(W, "log_P8_007.txt"), split = TRUE)
iv1 <- paste0("Ltrade_yhat_seq", 1:5); iv2 <- paste0("Ltrade_seq_hat", 1:5)
ctl <- c("Lminreg302", "Lcaprat301", "Ldistance", "Latopsecorr", "Ligosecorr", "peaceyrs", "t1", "t2", "t3")
d <- read_dta(file.path(WD, "table1.dta"), col_select = any_of(c("statea", "stateb", "year", "dyad", "dichmid20", "Ltradesecorr", iv1, iv2, ctl))) %>%
  zap_labels() %>% as.data.frame()
cat("rows:", nrow(d), " years:", range(d$year), " self-dyads:", sum(d$statea == d$stateb), " directed (a>b rows):", sum(d$statea > d$stateb), "\n")
stopifnot(!anyDuplicated(d[, c("statea", "stateb", "year")]))
d <- d[!is.na(d$Ltrade_seq_hat1), ]      # authors: Model 3 'if !missing(Ltrade_seq_hat1)'; Model 4 sample identical (N 154,570 both)

# MI logit: fit on each of 5 imputations, combine with Rubin's rules (authors' micalc: V = W + (1+1/m) B)
fit_mi <- function(ivs, keyname) function(r, data) {
  r0 <- setdiff(r, keyname); bs <- c(); ws <- c(); n <- NULL
  for (v in ivs) {
    data$KEY <- data[[v]]
    f <- fit_bin("dichmid20", c("KEY", r0), data, "logit", "iid", "KEY")
    bs <- c(bs, f$KEY["b"]); ws <- c(ws, f$KEY["se"]^2); n <- c(n, attr(f, "n"))
  }
  stopifnot(length(unique(n)) == 1)
  m <- length(bs); b <- mean(bs); V <- mean(ws) + (1 + 1 / m) * var(bs); se <- sqrt(V)
  out <- setNames(list(c(b = b, se = se, p = 2 * pnorm(-abs(b / se)))), keyname); attr(out, "n") <- n[1]
  cat("MI combined", keyname, ": b =", round(b, 4), " se =", round(se, 4), "\n"); out
}
meta <- list(pid = "P8_007", short_cite = "Boehmke, Chyzh & Thies 2016", journal = "PSRM", year = 2016, doi = "10.1017/psrm.2015.34",
  family = "TIE", headline_table_model = NA, article_estimator = "logit (Stata default SE), IVs lagged 1 yr, cubic time + peace years; 5-draw MI combined by Rubin's rules",
  est_used_all_steps = "logit iid SE on each of 5 imputations + Rubin combination, all six steps", author_has_history_control = "N (peace years only)")
spec <- list(
  list(key = "trade_seq_ivnetwork", ivs = iv1, model = "Table 1, Model 3 (logit, instrumented network)", pub = c(-4.372, 0.496),
       ev = "Abstract/Table 1: IV-corrected trade SEq 'many times' larger than naive; instrumented network model"),
  list(key = "trade_seq_ivscore", ivs = iv2, model = "Table 1, Model 4 (logit, instrumented SEq score)", pub = c(-7.616, 0.642),
       ev = "Table 1: instrumented-score correction of trade SEq (Maoz et al. 2006) -> MIDs"))
for (s in spec) {
  fitter <- fit_mi(s$ivs, s$key)
  rhs <- c(s$key, ctl); vars <- c("dichmid20", s$ivs, ctl)
  d0 <- d[complete.cases(d[, vars]), ]; cat("\n==== ", s$key, " complete-case N:", nrow(d0), "(published 154,570)\n")
  full <- fitter(rhs, d0)
  ds <- common_sample(d, vars, "statea", "stateb", "year")
  M <- run_six(rhs, ds, fitter)
  meta$headline_table_model <- s$model
  write_p8(meta, s$key, s$ev, s$pub[1], s$pub[2], 154570, TRUE, FALSE, M,
    paste0("Authors' constructed data (table1.dta). Full-sample MI b = ", round(full[[1]]["b"], 3), " (SE ", round(full[[1]]["se"], 3), "). ",
           "Key IV enters as its 5 imputed versions; history/benchmarks added to every imputation. ERGM columns (Models 5-7) not rerun (DEFERRED)."))
}
sink()
