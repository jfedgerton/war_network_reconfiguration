# X005 Hensel, Goertz & Diehl 2000 JOP "The Democratic Peace and Rivalries".
# Headline = Table 3 "Overall" panel: dispute propensity in regime-change proto+enduring rivalries by dyadic regime type
#   (jop1.sas: proc freq demdix2*midyr / chisq where demchg=1 and 16>afteryr>-1 and evrivgd>1). Published: joint-dem 47/214 (22.0%),
#   one-dem 104/303, neither 10/31; X2=9.37 p<.01; N=548.
# Authors' estimator is a contingency table. We reparameterize H2 (joint-dem years vs years with >=1 non-democracy) as a bivariate
#   logit, which reproduces the 2x2 log-odds ratio exactly: b_pub (implied) = log((47/167)/(114/220)).
# Key term: jointdem = (demdix2==2). Evidence: "lower dispute propensities while both rivals are democratic".
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p4/P4_common.R"))
WO <- file.path(OUT, "work_p4/X005"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work_p4/X005"); sink(file.path(WO, "log_X005.txt"), split = TRUE)
d <- read.csv(file.path(W, "erdemyr.csv"), header = FALSE,
  col.names = c("gdnumber","dyad","year","evrivgd","demdix2","chgdix","chgdixa","demchg","demchga","midyr","afteryr","rivbegyr","rivendyr","rivdecad","aftdecad"))
d <- d[d$demchg == 1 & d$afteryr > -1 & d$afteryr < 16 & d$evrivgd > 1, ]
tb <- table(d$demdix2, d$midyr); print(tb)
stopifnot(nrow(d) == 548, tb["2", "1"] == 47, tb["1", "1"] == 104, tb["0", "1"] == 10)
d$jointdem <- as.numeric(d$demdix2 == 2)
d$c1 <- d$dyad %/% 1000; d$c2 <- d$dyad %% 1000
stopifnot(all(d$c1 > 0), all(d$c2 > d$c1))
b_pub <- log((47/167) / (114/220)); cat("implied b_pub:", b_pub, "\n")
f0 <- midyr ~ jointdem
r <- fit_keys(f0, d, "jointdem", "iid"); n_orig <- r$n
st <- repro_status(b_pub, TRUE, r$k$jointdem)
cat("years:", range(d$year), "\n")
ds <- common_sample(d, f0, character(0), "c1", "c2", "year")
M <- run_six(f0, ds, "jointdem", "iid", WO)
write_p4("X005", "Hensel, Goertz & Diehl 2000 JOP", "DEMOCRACY", "KEY", "jointdem",
  "abstract: 'significantly lower dispute propensities while both rivals are democratic'", st, b_pub, r$k$jointdem[["b"]],
  "bivariate logit, iid SE (exact reparameterization of Table 3 2x2 contingency table)", n_orig, M,
  "N (no multivariate model; cross-tabs only; rivalry-year sample conditions on being in a Goertz-Diehl rivalry)",
  paste0("Table 3 Overall (regime-change proto+enduring rivalries, rivalry years + 15 yrs after last MID; N=548). b_pub is the log odds ratio implied by the published cell counts (joint-dem vs >=1 non-dem); published test is 3-category chi2=9.37 p<.01. ",
         "Common sample drops ", n_orig - M$m0$n, " obs (pre-1822 years / benchmark missingness). Within-rivalry sample: hist_dyad counts the same rivalry's recent MIDs. bD (lagged joint democracy) NOT a placebo."))
sink()
