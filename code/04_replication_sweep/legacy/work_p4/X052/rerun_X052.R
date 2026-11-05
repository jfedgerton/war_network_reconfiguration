# X052 Daxecker 2007 EJIR "Perilous Polities?" (10.1177/1354066107083146). Data + do-file: Dataverse doi:10.7910/DVN/PMOP0J.
# MONADIC state-year 1950-2000. Headline = Table 2 Onset Model 2 (full controls, contains the interaction the abstract cites):
#   logit mid democratization abschg10 abs_democ demneighbor civil ln_cinc ln_cont spline1-3, robust   (abschg10 is 'abchg' in posted data)
# KEY terms (fixed a priori from abstract): abs_democ = regime change x democratization ("such exacerbating effects are reversed for
#   democratizing states"); abchg = magnitude of regime change ("regime change itself may increase conflict propensities").
# History (monadic adaptation, see notes): hist_mono_w15/w26 = log1p(state's own GML MID onsets in t-1..t-5 / t-2..t-6);
#   no dyadic term exists -> mD not estimable (b_d NA); mDM = + hist_mono_w15; m26 = + hist_mono_w26.
#   Benchmarks (monadic analogues): bC = state's mean lagged CINC (t-1..t-5); bD = share of t-1..t-5 with own polity2 >= 6.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p4/common_p4_b04.R"))
suppressPackageStartupMessages(library(haven))
WO <- file.path(OUT, "work_p4/X052"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work_p4/X052"); sink(file.path(WO, "log_X052.txt"), split = TRUE)
d <- zap_labels(read_dta(file.path(W, "data.dta"))) %>% as.data.frame()
stopifnot(anyDuplicated(d[, c("ccode", "year")]) == 0)
cat("rows", nrow(d), "years", range(d$year), "\n")
chk <- with(d, abs_democ - abchg * democratization); cat("abs_democ == abchg*democratization:", all(abs(chk) < 1e-8, na.rm = TRUE), "\n")
f0 <- mid ~ democratization + abchg + abs_democ + demneighbor + civil + ln_cinc + ln_cont + spline1 + spline2 + spline3
terms <- c("abs_democ", "abchg")
r0 <- ff_logit_hc1(f0, d); print(round(r0, 4)); n_orig <- attr(r0, "n")
# ---- monadic history ----
H <- readRDS(Sys.getenv("SWEEP_HIST"))
d <- d %>% left_join(H$sth, by = c("ccode", "year")) %>% mutate(s15 = coalesce(s15, 0), s26 = coalesce(s26, 0),
  hist_mono_w15 = log1p(s15), hist_mono_w26 = log1p(s26))
bm <- lapply(1:5, function(L) { x <- data.frame(ccode = d$ccode, yl = d$year - L) %>%
  left_join(H$nmc, by = c("ccode", "yl" = "year")) %>% left_join(H$pol, by = c("ccode", "yl" = "year"))
  data.frame(c = x$cinc, dm = as.numeric(x$polity2 >= 6)) })
d$bm_cinc_w15 <- rowMeans(sapply(bm, `[[`, "c"), na.rm = TRUE); d$bm_democ_w15 <- rowMeans(sapply(bm, `[[`, "dm"), na.rm = TRUE)
ds <- d[complete.cases(d[, all.vars(f0)]), ]; n_cc <- nrow(ds)
ds <- ds[ds$year >= 1822 & ds$year <= 2011 & is.finite(ds$bm_cinc_w15) & is.finite(ds$bm_democ_w15), ]
cat("complete-case N:", n_cc, " common sample N:", nrow(ds), "\n"); stopifnot(nrow(ds) > 0.8 * n_cc)
C <- function(k) file.path(WO, paste0("cache_", k, ".rds"))
M <- list(m0  = fit_keys(f0, ds, terms, ff_logit_hc1, C("m0")), mD = NULL,
          mDM = fit_keys(update(f0, . ~ . + hist_mono_w15), ds, terms, ff_logit_hc1, C("mDM")),
          m26 = fit_keys(update(f0, . ~ . + hist_mono_w26), ds, terms, ff_logit_hc1, C("m26")),
          bC  = fit_keys(update(f0, . ~ . + bm_cinc_w15), ds, terms, ff_logit_hc1, C("bC")),
          bD  = fit_keys(update(f0, . ~ . + bm_democ_w15), ds, terms, ff_logit_hc1, C("bD")))
meta <- list(xid = "X052", short_cite = "Daxecker 2007", journal = "EJIR", year = 2007, doi = "10.1177/1354066107083146")
est <- "logit, Stata robust (HC1) SEs; MONADIC state-year"
nt <- paste0("MONADIC (state-year 1950-2000): history adapted -- no dyadic term (b_d NA); b_dm = + own MID-onset history t-1..t-5; b_26 = + t-2..t-6; ",
  "bC/bD are monadic analogues (own lagged CINC; own lagged democracy share). Published Table 2 not accessible (EJIR paywall) -> UNVERIFIED; ",
  "b_repro is from the authors' posted data + do-file (Onset Model 2). DV 'mid' = MID involvement dummy as posted.")
write_p4(meta, "abs_democ", "DEMOCRACY", "KEY", "\"such exacerbating effects are reversed for democratizing states\"", "UNVERIFIED",
  NA, r0["abs_democ", 1], est, n_orig, M, "N (peace-years splines only)", paste(nt, "Key = |regime change| x democratization product term."))
write_p4(meta, "abchg", "DEMOCRACY", "KEY", "\"regime change itself may increase conflict propensities\"", "UNVERIFIED",
  NA, r0["abchg", 1], est, n_orig, M, "N (peace-years splines only)", paste(nt, "abchg = magnitude of Polity change (abschg10 in do-file); constituent at democratization=0 (autocratizing/other)."))
sink()
