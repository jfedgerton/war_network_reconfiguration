# KEY rerun 61: key claim = trade x challenger-democracy interaction (abstract: democratic leaders more averse to initiating
# conflict with trading partners). Key model = authors' do-file "NOW WITH DEMOCRACY & TRADE INTERACTION" (O&R data):
#   estsimp logit initiate ordepa ordepb orjntdep ordm1dpa dmat21_1 jntdemoc relcap sdiff majdyd border lndist peacyr1-4, robust mi(5)
# Key term = ordm1dpa (ordepa x dmat21_1, challenger democracy 1-21). Same data/imputations/sample rules as rerun_61.R.
# Results saved to key_61.rds (results/61.csv untouched).
# Rerun 61: Gelpi & Grieco 2008 JPR, trade dependence -> MID initiation (directed dyads 1950-1992)
# Headline: first model in authors' do-file (Oneal & Russett trade data, joint-dependence interaction):
#   estsimp logit initiate ordepa ordepb orjntdep jntdemoc relcap sdiff majdyd border lndist peacyr1-4,
#   robust mi(5 imputed files)  -> Rubin's rules over 5 imputations, Stata robust (HC1-type) SE.
# Tie = ordepa (challenger's trade dependence), constituent term of ordepa x ordepb (orjntdep).
# Published table not obtainable (paywalled, no OA copy found) -> UNVERIFIED.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work/rerun_helpers_s68.R"))
W <- file.path(ROOT, "work/61"); WO <- file.path(OUT, "work/61"); dir.create(WO, recursive = TRUE, showWarnings = FALSE); sink(file.path(WO, "log_key_61.txt"), split = TRUE)
library(haven)

ctrl <- c("ordepa","ordepb","orjntdep","dmat21_1","jntdemoc","relcap","sdiff","majdyd","border","lndist","peacyr1","peacyr2","peacyr3","peacyr4")
TIE <- "ordm1dpa"
f0 <- reformulate(c(TIE, ctrl), "initiate")

# 1. Load the 5 imputations
D <- lapply(1:5, function(i) as.data.frame(zap_labels(read_dta(file.path(W, sprintf("intltrade.5.18.04.impute%d.dta", i))))))
stopifnot(all(sapply(D, nrow) == 40100))
stopifnot(all(sapply(D[-1], function(x) identical(x[, c("ccode1","ccode2","year","initiate")], D[[1]][, c("ccode1","ccode2","year","initiate")]))))

# Rubin's rules combining over imputations (as Clarify estsimp mi())
fit_mi <- function(f, Dl, label) {
  r <- lapply(seq_along(Dl), function(i) fit_glm(f, Dl[[i]], TIE, cluster = "robust", label = paste(label, "imp", i)))
  b <- sapply(r, `[[`, "b"); v <- sapply(r, `[[`, "se")^2; m <- length(b)
  qbar <- mean(b); T <- mean(v) + (1 + 1/m) * var(b)
  dfr <- (m - 1) * (1 + mean(v) / ((1 + 1/m) * var(b)))^2
  out <- list(b = qbar, se = sqrt(T), p = 2 * pt(-abs(qbar / sqrt(T)), dfr), n = r[[1]]$n)
  cat(sprintf("== MI combined [%s]: b=%.5f se=%.5f p=%.5g n=%d (per-imp b: %s)\n", label, out$b, out$se, out$p, out$n,
              paste(round(b, 4), collapse = ",")))
  out
}

# 2. Reproduce (full original sample)
rep <- fit_mi(f0, D, "repro")
n_orig <- rep$n
b_pub <- NA
status <- "UNVERIFIED"

# 3. Merge history (directed dyads; add_history handles order)
S <- lapply(D, function(x) {
  x <- add_history(x, "ccode1", "ccode2", "year")
  x <- x[x$year >= 1822 & x$year <= 2011, ]
  x[complete.cases(x[, c(all.vars(f0), "hist_dyad_w15","hist_mono_w15","hist_dyad_w26","hist_mono_w26")]), ]
})
cat("Estimation sample", nrow(S[[1]]), "; dropped vs original", n_orig - nrow(S[[1]]), "\n")

# 4. Models
m0  <- fit_mi(f0, S, "m0")
mD  <- fit_mi(update(f0, . ~ . + hist_dyad_w15), S, "mD")
mDM <- fit_mi(update(f0, . ~ . + hist_dyad_w15 + hist_mono_w15), S, "mDM")
m26 <- fit_mi(update(f0, . ~ . + hist_dyad_w26 + hist_mono_w26), S, "m26")
SC <- lapply(S, function(x) x[!is.na(x$bm_cinc_w15), ]); SD <- lapply(S, function(x) x[!is.na(x$bm_democ_w15), ])
cat("bC sample", nrow(SC[[1]]), " bD sample", nrow(SD[[1]]), "\n")
b0C <- if (nrow(SC[[1]]) < nrow(S[[1]])) fit_mi(f0, SC, "m0 on bC sample")$b else NULL
bC  <- fit_mi(update(f0, . ~ . + bm_cinc_w15), SC, "bC")
b0D <- if (nrow(SD[[1]]) < nrow(S[[1]])) fit_mi(f0, SD, "m0 on bD sample")$b else NULL
bD  <- fit_mi(update(f0, . ~ . + bm_democ_w15), SD, "bD")

saveRDS(list(rep = rep, n_orig = n_orig, R = list(m0 = m0, mD = mD, mDM = mDM, m26 = m26, bC = bC, bD = bD), b0C = b0C, b0D = b0D), file.path(WO, "key_61.rds"))
print(sapply(list(m0 = m0, mD = mD, mDM = mDM, m26 = m26, bC = bC, bD = bD), unlist)); cat("b0C", if (is.null(b0C)) NA else b0C, "b0D", if (is.null(b0D)) NA else b0D, "\n")
# port (extraction only, no estimation): one CSV row with the numbers saved in key_61.rds (the sweep copied them into
# results_key/61.csv via write_key_batch.py; pct_bC / pct_bD vs m0 re-fit on the benchmark's own sample when smaller, as there).
rC <- if (is.null(b0C)) m0$b else b0C; rD <- if (is.null(b0D)) m0$b else b0D
# Estimator audit (sweep phase "estaudit"): the article's Clarify/estsimp reports normal-based inference on the MI-combined
# estimate, so p-values are 2 * pnorm(-|b / se|) from the Rubin-combined b and se (the t-reference p above is not used).
pn <- function(m) 2 * pnorm(-abs(m$b / m$se))
ex <- data.frame(id = "61", key_term = "ordm1dpa", b0 = m0$b, se0 = m0$se, p0 = pn(m0), b_dm = mDM$b, p_dm = pn(mDM), b_26 = m26$b, p_26 = pn(m26),
  pct_bC = 100 * (bC$b - rC) / abs(rC), pct_bD = 100 * (bD$b - rD) / abs(rD), n_used = m0$n)
stopifnot(all(is.finite(unlist(ex[, -(1:2)]))))
write.csv(ex, file.path(WO, "key_61_extract.csv"), row.names = FALSE)
sink()
