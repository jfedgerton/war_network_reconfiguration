# P8_031 Gartzke & Weisiger 2014 ISQ "Under Construction" (10.1111/isqu.12113)
# Data: authors' ISQ Dataverse (doi:10.7910/DVN/26209) gartzke_weisiger_ISQ_2014_table_3_data_012014 (+ do file).
# Paper (author copy, pages.ucsd.edu/~egartzke, draft 17 Jan 2011 -- published table not available here): Table 3, Model 14 ("Difference"),
#   dyadic, fatal MID onset, 1816-2001, King-Zeng rare-events logit, dyad-clustered robust SE, BTSCS peace years + 3 splines. N = 619,104.
#   do file: relogit deadlyl polave pcenerg diff1 demloi engypop dydiff logdist cntgdumy allydumy capratio onemajor deadyrs deadyer*, cluster(dyadid)
# KEY TERMS (a priori from H1-H5 + dyadic democracy): polave (H1 systemic democracy), diff1 (H2 systemic difference), dydiff (H3 dyadic difference),
#   pcenerg (H4 systemic development), engypop (H5 dyadic development), demloi (dyadic democracy low; weak-link democratic peace term).
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/p8_common.R"))
source(file.path(CODE, "work/relogit_kz_ea1.R"))
suppressMessages(library(haven))
WD <- file.path(ROOT, "work_p8/P8_031"); W <- file.path(OUT, "work_p8/P8_031"); dir.create(W, recursive = TRUE, showWarnings = FALSE); sink(file.path(W, "log_P8_031.txt"), split = TRUE)
rhs <- c("polave", "pcenerg", "diff1", "demloi", "engypop", "dydiff", "logdist", "cntgdumy", "allydumy", "capratio", "onemajor",
         "deadyrs", "deadyer1", "deadyer2", "deadyer3")
d <- read_dta(file.path(WD, "t3.dta"), col_select = all_of(c("dyadid", "year", "deadlyl", rhs))) %>% zap_labels() %>% as.data.frame()
# dyad codes: dyadid = ccode1*1000 + ccode2 (checked below against COW state codes)
d$c1 <- d$dyadid %/% 1000; d$c2 <- d$dyadid %% 1000
cs <- unique(peacesciencer::cow_states$ccode)
cat("dyadid range:", range(d$dyadid), " share c1 in COW:", mean(d$c1 %in% cs), " c2 in COW:", mean(d$c2 %in% cs), " c1<c2:", mean(d$c1 < d$c2), "\n")
stopifnot(mean(d$c1 %in% cs) > .99, mean(d$c2 %in% cs) > .99, !anyDuplicated(d[, c("dyadid", "year")]))
# DV timing diagnostic: after an event coded in year t, which year's peace-years counter resets to 0?
d <- d[order(d$dyadid, d$year), ]
nx <- function(x, k) ave(x, d$dyadid, FUN = function(z) c(z[-(1:k)], rep(NA, k)))
ev <- d$deadlyl == 1
r1 <- mean(nx(d$deadyrs, 1)[ev] == 0, na.rm = TRUE); r2 <- mean(nx(d$deadyrs, 2)[ev] == 0, na.rm = TRUE)
cat("P(deadyrs=0 at t+1 | event t) =", round(r1, 3), "; at t+2 =", round(r2, 3), "\n")
lead_dv <- TRUE   # [2026-11-14] data check: deadlyl matches COW MID-years at t+1 (0.97) vs t (0.59); the reset test cannot detect a lead (was r2 > r1)
d$yr_dv <- d$year + as.integer(lead_dv)
cat("DV treated as", if (lead_dv) "LEAD (t+1): history keyed to year+1" else "contemporaneous: history keyed to year", "\n")

terms <- c("polave", "pcenerg", "diff1", "demloi", "engypop", "dydiff")
fitter <- function(r, data) {
  o <- relogit_kz(as.formula(paste("deadlyl ~", paste(r, collapse = " + "))), data, cluster = "dyadid")
  cat("\n---- relogit (KZ-corrected, dyad-clustered) N =", attr(o, "n"), "\n"); print(round(o[c(terms, paste0(terms, "__unc")), ], 5))
  out <- lapply(setNames(terms, terms), function(t) o[t, c("b", "se", "p")]); attr(out, "n") <- attr(o, "n"); out
}
vars <- c("deadlyl", rhs)
d0 <- d[complete.cases(d[, vars]), ]; cat("complete-case N:", nrow(d0), "(published 619,104)\n")
full <- fitter(rhs, d0)
ds <- common_sample(d, vars, "c1", "c2", "yr_dv")
M <- run_six(rhs, ds, fitter)
# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  # [2026-11-10 audit2] hook moved to the new headline term (pcenerg)
  sens_run(id = "P8_031", dat = ds, y = "deadlyl", tie = "pcenerg", ctrl = setdiff(rhs, "pcenerg"), fe = NULL,
           fit = function(dat, add) { r <- fitter(c(rhs, add), dat)$pcenerg; c(b = unname(r[["b"]]), se = unname(r[["se"]]), p = unname(r[["p"]])) })
}
meta <- list(pid = "P8_031", short_cite = "Gartzke & Weisiger 2014", journal = "ISQ", year = 2014, doi = "10.1111/isqu.12113",
  family = "DEMOCRACY", headline_table_model = "Table 3, Model 14 (dyadic, fatal MID onset, 1816-2001)",
  article_estimator = "rare-events logit (King-Zeng relogit), dyad-clustered robust SE, peace years + 3 splines",
  est_used_all_steps = "King-Zeng bias-corrected logit (own implementation, work/relogit_kz_ea1.R), dyad-clustered SE, all six steps",
  author_has_history_control = "N (peace years/splines only)")
pub <- list(polave = c(-0.0530, 0.0338), pcenerg = c(-0.644, 0.130), diff1 = c(0.648, 0.0849), demloi = c(-0.0473, 0.0379),
            engypop = c(-0.0716, 0.0619), dydiff = c(0.105, 0.0233))
ev <- c(polave = "H1: 'Systemic democracy is unrelated to systemic peace'",
        pcenerg = "H4: 'Increases in economic development at the systemic level' reduce conflict; abstract: 'greater systemic development encourages peace'",
        diff1 = "H2: 'Increases in regime type heterogeneity (difference) at the system level' increase conflict",
        demloi = "Dyadic democratic-peace term (weak-link democracy) set against difference: 'Democracy has no significant impact'",
        engypop = "H5: 'Developed dyads are no less prone to experience disputes'",
        dydiff = "H3: 'Regime type difference in dyads should increase disputes'")
for (t in terms) write_p8(meta, t, ev[[t]], pub[[t]][1], pub[[t]][2], 619104, NA, FALSE, M,
  paste0("Authors' data + do file. Published values from the authors' 2011 manuscript (final ISQ table not accessible here). Full-sample KZ b = ",
         signif(full[[t]]["b"], 3), " (SE ", signif(full[[t]]["se"], 3), "). Systemic terms (polave, pcenerg, diff1) vary by year only. DV timing: ",
         if (lead_dv) "lead" else "contemporaneous", " (peace-year reset diagnostic)."))
sink()
