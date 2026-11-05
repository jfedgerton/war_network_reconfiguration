# Rerun 61: Gelpi & Grieco 2008 JPR, trade dependence -> MID initiation (directed dyads 1950-1992)
# Headline: first model in authors' do-file (Oneal & Russett trade data, joint-dependence interaction):
#   estsimp logit initiate ordepa ordepb orjntdep jntdemoc relcap sdiff majdyd border lndist peacyr1-4,
#   robust mi(5 imputed files)  -> Rubin's rules over 5 imputations, Stata robust (HC1-type) SE.
# Tie = ordepa (challenger's trade dependence), constituent term of ordepa x ordepb (orjntdep).
# Published table not obtainable (paywalled, no OA copy found) -> UNVERIFIED.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work/rerun_helpers_s68.R"))
W <- file.path(ROOT, "work/61"); WO <- file.path(OUT, "work/61"); dir.create(WO, recursive = TRUE, showWarnings = FALSE); sink(file.path(WO, "log_61.txt"), split = TRUE)
library(haven)

ctrl <- c("ordepb","orjntdep","jntdemoc","relcap","sdiff","majdyd","border","lndist","peacyr1","peacyr2","peacyr3","peacyr4")
TIE <- "ordepa"
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

# 5. Write
write_row(61, "Gelpi & Grieco 2008", status, b_pub, rep$b, "logit, robust SE, MI (5 imps, Rubin's rules)", "ordepa",
  n_orig, m0$n, m0, mD, mDM, m26, bC, bD, b0C, b0D, ahc = "Y",
  notes = sprintf("Headline = first model in do-file (O&R trade data, ordepa x ordepb interaction); tie is constituent term ordepa. Published table paywalled -> UNVERIFIED. Dropped %d obs (1950-1992 all in range). bC N=%d, bD N=%d.",
                  n_orig - m0$n, bC$n, bD$n))
sink()
