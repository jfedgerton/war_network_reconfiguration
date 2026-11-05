# Rerun 68: Kinne 2013 JPR, IGO network convergence -> MID onset
# Headline: Table 1 Model 4 (author's preferred, best fit; used in all Table 2 models):
#   logit mid requivlag5 + base, vce(cluster dyadid), year >= 1970.
# Published (Table 1 M4): requivlag5 = -1.760 (0.317), N = 660,366.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work/rerun_helpers_s68.R"))
W <- file.path(ROOT, "work/68"); WO <- file.path(OUT, "work/68"); dir.create(WO, recursive = TRUE, showWarnings = FALSE); sink(file.path(WO, "log_68.txt"), split = TRUE)
library(haven); options(use_fixest = TRUE)  # memory-lean estimation (shared machine; glm path was OOM-killed)

# 1. Load authors' data (as in 01.Estimate.Stata.do)
base <- c("s_wt_glo","igos","demlow","dep_ab","sum_lncap","allies","contig","lndist","systsize",
          "peaceyrs_mid","_spline1_mid","_spline2_mid","_spline3_mid")
d <- read_dta(file.path(W, "JPR_Replication/Data/zData.dta"),
              col_select = all_of(c("dyadid","ccode1","ccode2","year","mid","requivlag5", base)))
d <- as.data.frame(zap_labels(d))
names(d) <- sub("^_", "x_", names(d)); base <- sub("^_", "x_", base)
d <- d[d$year >= 1970, ]
stopifnot(nrow(d) == 739692)

# 2. Reproduce
f0 <- reformulate(c("requivlag5", base), "mid")
b_pub <- -1.760
rep <- fit_glm(f0, d, "requivlag5", cluster = "dyadid", label = "repro")
n_orig <- rep$n
cat("b_pub", b_pub, "b_repro", rep$b, "rel diff", abs(rep$b - b_pub) / abs(b_pub), "N", n_orig, "(pub 660366)\n")
status <- if (abs(rep$b - b_pub) / abs(b_pub) <= 0.10 && rep$p < .05) "REPRODUCED" else "NOT_REPRODUCED"

# 3. Merge history (undirected dyads; ccode1<ccode2 in data)
d <- d[complete.cases(d[, all.vars(f0)]), ]; invisible(gc())   # drop rows unusable in any model (memory)
d <- add_history(d, "ccode1", "ccode2", "year")
s <- d[d$year >= 1822 & d$year <= 2011, ]; rm(d); invisible(gc())
s <- s[complete.cases(s[, c(all.vars(f0), "dyadid", "hist_dyad_w15", "hist_mono_w15", "hist_dyad_w26", "hist_mono_w26")]), ]
cat("Estimation sample", nrow(s), "; dropped vs original", n_orig - nrow(s), "\n")

# 4. Models
m0  <- fit_glm(f0, s, "requivlag5", cluster = "dyadid", label = "m0")
mD  <- fit_glm(update(f0, . ~ . + hist_dyad_w15), s, "requivlag5", cluster = "dyadid", label = "mD")
mDM <- fit_glm(update(f0, . ~ . + hist_dyad_w15 + hist_mono_w15), s, "requivlag5", cluster = "dyadid", label = "mDM")
m26 <- fit_glm(update(f0, . ~ . + hist_dyad_w26 + hist_mono_w26), s, "requivlag5", cluster = "dyadid", label = "m26")
sC <- s[!is.na(s$bm_cinc_w15), ]; sD <- s[!is.na(s$bm_democ_w15), ]
cat("bC sample", nrow(sC), " bD sample", nrow(sD), "\n")
b0C <- if (nrow(sC) < nrow(s)) fit_glm(f0, sC, "requivlag5", cluster = "dyadid", label = "m0 on bC sample")$b else NULL
bC  <- fit_glm(update(f0, . ~ . + bm_cinc_w15), sC, "requivlag5", cluster = "dyadid", label = "bC")
b0D <- if (nrow(sD) < nrow(s)) fit_glm(f0, sD, "requivlag5", cluster = "dyadid", label = "m0 on bD sample")$b else NULL
bD  <- fit_glm(update(f0, . ~ . + bm_democ_w15), sD, "requivlag5", cluster = "dyadid", label = "bD")

# 5. Write
write_row(68, "Kinne 2013", status, b_pub, rep$b, "logit, dyad-clustered SE (Stata G/(G-1))", "requivlag5",
  n_orig, m0$n, m0, mD, mDM, m26, bC, bD, b0C, b0D, ahc = "Y",
  notes = sprintf("Headline = Table 1 M4 (5-yr lagged convergence, author-preferred). Years 1970-2000 all in GML range; dropped %d obs. Benchmarks vs baseline on own non-missing sample (bC N=%d, bD N=%d).",
                  n_orig - m0$n, bC$n, bD$n))
sink()
