# Rerun 109: Johnson & Leeds 2011 FPA, "Defense Pacts: A Prescription for Peace?"
# Headline: Table 1 probit dispute ~ ptargdef + ... + peaceyrs(^2,^3); default (non-robust) ML SEs.
# b_pub from authors' smcl log: ptargdef = -.0615538 (se .0154694), N = 1,077,992.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work/rerun_common_b7.R"))
WO <- file.path(OUT, "work/109"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work/109"); sink(file.path(WO, "log_109.txt"), split = TRUE)
C <- function(k) file.path(WO, paste0("cache_", k, ".rds"))

# 1. Load authors' data (directed dyad-year: challenger, target, year)
d <- haven::read_dta(file.path(W, "jl.dta"),
  col_select = c(challenger, target, year, dispute, ptargdef, pchalally, pchaloff, pchalneu,
                 ln_distance, capprop, jdem, s_un_glo, peaceyrs, peaceyrs2, peaceyrs3))
d <- as.data.frame(haven::zap_labels(d))
stopifnot(nrow(d) == 1077992, all(d$dispute %in% c(0, 1)))
n_orig <- 1077992; b_pub <- -0.0615538

ctrl <- "pchalally + pchaloff + pchalneu + ln_distance + capprop + jdem + s_un_glo + peaceyrs + peaceyrs2 + peaceyrs3"
f0 <- as.formula(paste("dispute ~ ptargdef +", ctrl))
PB <- binomial(link = "probit")

# 2. Reproduce on full sample
r <- fit_tie(f0, d, PB, "ptargdef", cache = C("repro"))
b_repro <- r$b
stopifnot(r$n == n_orig)
repro_ok <- abs(b_repro - b_pub) / abs(b_pub) <= 0.10 && (r$p < .05) == TRUE
cat("\nb_pub", b_pub, "b_repro", b_repro, "rel diff", abs(b_repro - b_pub) / abs(b_pub), "\n")

# 3. Merge history; restrict to 1822-2011 (data end 2000)
d <- add_history(d, "challenger", "target", "year")
ds <- d[d$year >= 1822 & d$year <= 2011, ]
ds <- ds[complete.cases(ds[, c(all.vars(f0), "hist_dyad_w15", "hist_mono_w15", "hist_dyad_w26", "hist_mono_w26")]), ]
cat("dropped by 1822+ restriction / missing:", nrow(d) - nrow(ds), "\n")
cat("bm_cinc missing:", sum(is.na(ds$bm_cinc_w15)), " bm_democ missing:", sum(is.na(ds$bm_democ_w15)), "\n")
rm(d); gc()

# 4. Models on the common restricted sample
m0 <- fit_tie(f0, ds, PB, "ptargdef", cache = C("m0"))
mD <- fit_tie(update(f0, . ~ . + hist_dyad_w15), ds, PB, "ptargdef", cache = C("mD"))
mDM <- fit_tie(update(f0, . ~ . + hist_dyad_w15 + hist_mono_w15), ds, PB, "ptargdef", cache = C("mDM"))
m26 <- fit_tie(update(f0, . ~ . + hist_dyad_w26 + hist_mono_w26), ds, PB, "ptargdef", cache = C("m26"))
bC <- fit_tie(update(f0, . ~ . + bm_cinc_w15), ds, PB, "ptargdef", cache = C("bC"))
bD <- fit_tie(update(f0, . ~ . + bm_democ_w15), ds, PB, "ptargdef", cache = C("bD"))

# 5. Results row
write_row(109, "Johnson & Leeds 2011", if (repro_ok) "REPRODUCED" else "NOT_REPRODUCED",
  b_pub, b_repro, "probit (ML SE)", "ptargdef", n_orig, m0$n, m0, mD, mDM, m26, bC, bD, ahc = "Y",
  notes = paste0("Table 1; repro exact vs smcl log. Sample 1822-2000: dropped ", n_orig - m0$n,
                 " obs (1816-1821). bC N=", bC$n, ", bD N=", bD$n, " (benchmark NA rows listwise-dropped)."))
sink()
