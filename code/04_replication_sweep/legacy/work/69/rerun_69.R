# Rerun 69: Mitchell & Zawahri 2015 JPR -- river treaty design -> militarization of ICOW river claims
# Headline: Table 1 Model 1: logit midissyr rbo2 infoexchange2 monitoring2 enforcement2 conflictresolution2
#   icowsal relcaps recmidwt upstreampower downstreampower if year>1944, cluster(dyad)
# Tie = infoexchange2 (abstract: information exchange + enforcement most effective at preventing militarization).
# Published Table 1 not obtainable (paywalled). Authors' appendix log for a robustness variant (minlgdpen in
# place of up/downstream power; N=489, infoexchange2=-1.647443, rbo2=3.067952) is used to verify data/estimator.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work/rerun_helpers_s68.R"))
W <- file.path(ROOT, "work/69"); WO <- file.path(OUT, "work/69"); dir.create(WO, recursive = TRUE, showWarnings = FALSE); sink(file.path(WO, "log_69.txt"), split = TRUE)
library(haven)

d <- as.data.frame(zap_labels(read_dta(file.path(W, "icowdataMZ.dta"))))
stopifnot(nrow(d) == 1000)
d <- d[!is.na(d$year) & d$year > 1944, ]
TIE <- "infoexchange2"
inst <- c("rbo2","infoexchange2","monitoring2","enforcement2","conflictresolution2")
f0 <- reformulate(c(inst, "icowsal","relcaps","recmidwt","upstreampower","downstreampower"), "midissyr")

# 2a. Verify against authors' appendix log (robustness spec)
fL <- reformulate(c(inst, "icowsal","relcaps","recmidwt","minlgdpen"), "midissyr")
lg <- fit_glm(fL, d, TIE, cluster = "dyad", label = "appendix-log check")
lgR <- fit_glm(fL, d, "rbo2", cluster = "dyad", label = "appendix-log check rbo2")
cat(sprintf("Appendix log: infoexchange2 -1.647443 (se .7266) N=489 | repro %.6f (se %.4f) N=%d; rbo2 3.067952 | repro %.6f\n",
            lg$b, lg$se, lg$n, lgR$b))
log_ok <- abs(lg$b - (-1.647443)) < 1e-3 && lg$n == 489 && abs(lgR$b - 3.067952) < 1e-3

# 2b. Headline
rep <- fit_glm(f0, d, TIE, cluster = "dyad", label = "headline repro")
n_orig <- rep$n
cat("Headline N", n_orig, "(appendix descriptives N=538)\n")
b_pub <- NA
status <- "UNVERIFIED"

# 3. Merge history (chal/tgt COW codes; directed claim-dyad)
d <- add_history(d, "chal", "tgt", "year")
s <- d[d$year >= 1822 & d$year <= 2011, ]
s <- s[complete.cases(s[, c(all.vars(f0), "dyad", "hist_dyad_w15","hist_mono_w15","hist_dyad_w26","hist_mono_w26")]), ]
cat("Estimation sample", nrow(s), "; dropped vs original", n_orig - nrow(s), "\n")

# 4. Models (also track enforcement2 and rbo2)
run <- function(f, dat, lab) {
  r <- fit_glm(f, dat, TIE, cluster = "dyad", label = lab)
  r$bE <- fit_glm(f, dat, "enforcement2", cluster = "dyad", label = paste(lab, "(enf)"))$b; r
}
m0  <- run(f0, s, "m0")
mD  <- run(update(f0, . ~ . + hist_dyad_w15), s, "mD")
mDM <- run(update(f0, . ~ . + hist_dyad_w15 + hist_mono_w15), s, "mDM")
m26 <- run(update(f0, . ~ . + hist_dyad_w26 + hist_mono_w26), s, "m26")
sC <- s[!is.na(s$bm_cinc_w15), ]; sD <- s[!is.na(s$bm_democ_w15), ]
cat("bC sample", nrow(sC), " bD sample", nrow(sD), "\n")
b0C <- if (nrow(sC) < nrow(s)) fit_glm(f0, sC, TIE, cluster = "dyad", label = "m0 on bC sample")$b else NULL
bC  <- fit_glm(update(f0, . ~ . + bm_cinc_w15), sC, TIE, cluster = "dyad", label = "bC")
b0D <- if (nrow(sD) < nrow(s)) fit_glm(f0, sD, TIE, cluster = "dyad", label = "m0 on bD sample")$b else NULL
bD  <- fit_glm(update(f0, . ~ . + bm_democ_w15), sD, TIE, cluster = "dyad", label = "bD")

write_row(69, "Mitchell & Zawahri 2015", status, b_pub, rep$b, "logit, claim-dyad clustered SE", "infoexchange2",
  n_orig, m0$n, m0, mD, mDM, m26, bC, bD, b0C, b0D, ahc = "Y",
  notes = sprintf("Table 1 pub not obtainable; authors' appendix Stata log (minlgdpen variant) %s; N=%d matches appendix descriptives (538). Author controls recmidwt (claim MID history). enforcement2: m0 %.3f, mDM %.3f, m26 %.3f. Only ~16 MID events; dropped %d obs.",
                  if (log_ok) "matched exactly" else "NOT matched", n_orig, m0$bE, mDM$bE, m26$bE, n_orig - m0$n))

# port (extraction only): the co-key enforcement2 master row was copied by vb3/copy_rows.py from log_69.txt, where fit_glm prints
# round(ct, 5) for every model above. Re-fit the SAME models/samples with enforcement2 as the reported term (full precision);
# pct_bD base = m0 on the bD sample when that sample is smaller (as in copy_rows.py).
E <- function(f, dat, lab) { r <- fit_glm(f, dat, "enforcement2", cluster = "dyad", label = lab); list(b = r$b, se = r$se, p = r$p, n = r$n) }
e0  <- E(f0, s, "m0 (enf extract)")
eDM <- E(update(f0, . ~ . + hist_dyad_w15 + hist_mono_w15), s, "mDM (enf extract)")
e26 <- E(update(f0, . ~ . + hist_dyad_w26 + hist_mono_w26), s, "m26 (enf extract)")
eC  <- E(update(f0, . ~ . + bm_cinc_w15), sC, "bC (enf extract)")
eD  <- E(update(f0, . ~ . + bm_democ_w15), sD, "bD (enf extract)")
stopifnot(abs(e0$b - m0$bE) < 1e-10, abs(eDM$b - mDM$bE) < 1e-10)
rC <- if (nrow(sC) < nrow(s)) E(f0, sC, "m0 on bC sample (enf extract)")$b else e0$b
rD <- if (nrow(sD) < nrow(s)) E(f0, sD, "m0 on bD sample (enf extract)")$b else e0$b
ex <- data.frame(id = "69", key_term = "enforcement2", b0 = e0$b, se0 = e0$se, p0 = e0$p, b_dm = eDM$b, p_dm = eDM$p, b_26 = e26$b, p_26 = e26$p,
  pct_bC = 100 * (eC$b - rC) / abs(rC), pct_bD = 100 * (eD$b - rD) / abs(rD), n_used = e0$n)
stopifnot(all(is.finite(unlist(ex[, -(1:2)]))))
write.csv(ex, file.path(WO, "key_69_enforcement2_extract.csv"), row.names = FALSE)
sink()
