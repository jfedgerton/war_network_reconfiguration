# Key-audit rerun 89: Wright & Rider 2014. Key claim = does target-defense-pact deterrence hold when a territorial issue
# is at stake -> INTERACTION; key term = terrXdef (terriss x ptargdef). Same six models/sample as rerun_89.R.
# Also reports the reparametrized conditional effect of ptargdef in territorial dyads (ptargdef*terriss), same model.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work/rerun_common_k6.R"))
W <- file.path(ROOT, "work/89"); WO <- file.path(OUT, "work/89"); dir.create(WO, recursive = TRUE, showWarnings = FALSE); sink(file.path(WO, "log_key_89.txt"), split = TRUE)
d <- as.data.frame(zap_labels(read_dta(file.path(W, "WrightRider2014CMPS_Replication/icow_sample.dta"))))
stopifnot(nrow(d) == 174278)
ok <- complete.cases(d[, c("terriss","ptargdef","terrXdef")])
stopifnot(max(abs(d$terrXdef[ok] - d$terriss[ok] * d$ptargdef[ok])) < 1e-9)
d$ptd_terr <- d$ptargdef * d$terriss; d$ptd_non <- d$ptargdef * (1 - d$terriss)
BASE <- "dispute ~ terriss + ptargdef + terrXdef + pchaloff + pchalneu + jdem + ln_distance + capprop + defense_ally + s_un_glo + peaceyrs1 + py2 + py3"
BASE2 <- "dispute ~ terriss + ptd_non + ptd_terr + pchaloff + pchalneu + jdem + ln_distance + capprop + defense_ally + s_un_glo + peaceyrs1 + py2 + py3"
MV <- c("dispute","terriss","ptargdef","terrXdef","pchaloff","pchalneu","jdem","ln_distance","capprop","defense_ally","s_un_glo","peaceyrs1","py2","py3","challenger","target")
fitter  <- function(add, dat) fit_tie(as.formula(paste0(BASE, add)), dat, link = "probit", tie = "terrXdef")
fitter2 <- function(add, dat) fit_tie(as.formula(paste0(BASE2, add)), dat, link = "probit", tie = "ptd_terr")
d_full <- d[complete.cases(d[, setdiff(MV, c("challenger","target"))]), ]
rp <- fitter("", d_full); n_orig <- rp$n; stopifnot(n_orig == 173722); rm(d_full)
d <- add_history(d, "challenger", "target", "year")
d <- restrict_sample(d, MV)
R <- run_six(fitter, d)
R2 <- run_six(fitter2, d)
# ---- sensitivity hook (split only): runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_split(id = "89", dat = d,
             fit = function(dat, add) { r <- fitter(if (length(add)) paste0(" + ", paste(add, collapse = " + ")) else "", dat); c(b = r$b, se = r$se, p = r$p) })
}
tab <- function(R, term) do.call(rbind, lapply(names(R), function(k) data.frame(term = term, model = k, b = R[[k]]$b, se = R[[k]]$se, p = R[[k]]$p, n = R[[k]]$n)))
out <- rbind(tab(R, "terrXdef"), tab(R2, "ptargdef|terriss=1"))
out$b_repro <- rp$b; out$n_orig <- n_orig
print(out); write.csv(out, file.path(WO, "key_89_six.csv"), row.names = FALSE)
# port (extraction only, no estimation): one wide row for the key term with the same numbers as key_89_six.csv
# (the sweep assembled results_key/89.csv from key_89_six.csv in write_key_mine.py; pct = 100*(b - b0)/|b0| as there).
k <- out[out$term == "terrXdef", ]; g <- function(m, v) { x <- k[[v]][k$model == m]; stopifnot(length(x) == 1); x }
b0 <- g("m0", "b")
ex <- data.frame(id = "89", key_term = "terrXdef", b0 = b0, se0 = g("m0", "se"), p0 = g("m0", "p"), b_dm = g("mDM", "b"), p_dm = g("mDM", "p"),
  b_26 = g("m26", "b"), p_26 = g("m26", "p"), pct_bC = 100 * (g("bC", "b") - b0) / abs(b0), pct_bD = 100 * (g("bD", "b") - b0) / abs(b0),
  n_used = g("m0", "n"))
stopifnot(all(is.finite(unlist(ex[, -(1:2)]))))
write.csv(ex, file.path(WO, "key_89_extract.csv"), row.names = FALSE)
sink()
