# Key-audit rerun 42: Lee & Mitchell 2012. Key tie claim (abstract): "Higher bilateral FDI flows between two disputants
# significantly reduce the chances for escalation to high levels of violence". Headline onset model's key is GLOBAL FDI
# (not a tie), bilateral FDI explicitly null there. Key model for the tie: Table 2 M3, logit midfatyr ~ dflow_avegdplag + ...,
# cluster(dyad). b_pub_key = -5.072 (authors' do-file/paper). Six models on one common sample.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work/rerun_common_b7.R"))
W <- file.path(ROOT, "work/42"); WO <- file.path(OUT, "work/42"); dir.create(WO, recursive = TRUE, showWarnings = FALSE); sink(file.path(WO, "log_key_42.txt"), split = TRUE)
d <- as.data.frame(haven::zap_labels(haven::read_dta(file.path(W, "FDI_ICOW_LEEMITCHELL.dta"))))
stopifnot(nrow(d) == 39311)
LG <- binomial("logit"); T <- "dflow_avegdplag"
f3 <- midfatyr ~ dflow_avegdplag + fdi_gdpLOWlag + lwfdilag + capLOW + jointdem + lgdppcLOW + cwpceyrs + alliance +
  icowsal + otherMID1 + otherMID2
r3 <- fit_tie(f3, d, LG, T, "dyad"); cat("REPRO b_pub -5.072 b_repro", r3$b, "N", r3$n, "\n")
stopifnot(abs(r3$b - (-5.072)) / 5.072 <= 0.10)
d <- add_history(d, "ccode1", "ccode2", "year")
d <- d[d$year >= 1822 & d$year <= 2011, ]
ds <- d[complete.cases(d[, c(all.vars(f3), "dyad", "hist_dyad_w15", "hist_mono_w15", "hist_dyad_w26", "hist_mono_w26", "bm_cinc_w15", "bm_democ_w15")]), ]
cat("N_orig", r3$n, "N common", nrow(ds), "\n")
add <- list(m0 = "", mD = "+ hist_dyad_w15", mDM = "+ hist_dyad_w15 + hist_mono_w15", m26 = "+ hist_dyad_w26 + hist_mono_w26",
            bC = "+ bm_cinc_w15", bD = "+ bm_democ_w15")
out <- do.call(rbind, lapply(names(add), function(k) {
  r <- fit_tie(update(f3, as.formula(paste(". ~ .", add[[k]]))), ds, LG, T, "dyad")
  stopifnot(r$n == nrow(ds)); data.frame(model = k, b = r$b, se = r$se, p = r$p, n = r$n) }))
out$b_repro <- r3$b; out$n_orig <- r3$n
print(out); write.csv(out, file.path(WO, "key_42_six.csv"), row.names = FALSE)

# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "42", dat = ds, y = "midfatyr", tie = T,
           ctrl = c("fdi_gdpLOWlag", "lwfdilag", "capLOW", "jointdem", "lgdppcLOW", "cwpceyrs", "alliance", "icowsal", "otherMID1", "otherMID2"),
           fe = NULL,
           fit = function(dat, add) { r <- fit_tie(update(f3, as.formula(paste(". ~ .", if (length(add)) paste("+", paste(add, collapse = " + ")) else ""))), dat, LG, T, "dyad")
             c(b = r$b, se = r$se, p = r$p) })
}

# port (extraction only, no estimation): one wide row with the same numbers as key_42_six.csv
# (the sweep assembled results_key/42.csv from key_42_six.csv in write_key_mine.py; pct = 100*(b - b0)/|b0| as there).
g <- function(m, v) { x <- out[[v]][out$model == m]; stopifnot(length(x) == 1); x }
b0 <- g("m0", "b")
ex <- data.frame(id = "42", key_term = "dflow_avegdplag", b0 = b0, se0 = g("m0", "se"), p0 = g("m0", "p"), b_dm = g("mDM", "b"), p_dm = g("mDM", "p"),
  b_26 = g("m26", "b"), p_26 = g("m26", "p"), pct_bC = 100 * (g("bC", "b") - b0) / abs(b0), pct_bD = 100 * (g("bD", "b") - b0) / abs(b0),
  n_used = g("m0", "n"))
stopifnot(all(is.finite(unlist(ex[, -(1:2)]))))
write.csv(ex, file.path(WO, "key_42_extract.csv"), row.names = FALSE)
sink()
