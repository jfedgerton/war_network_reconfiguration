# VERIFY X033 Rasler & Thompson 2001 JPR "Rivalries and the Democratic Peace in the Major Power Subsystem" -- matched to the PDF.
# Headline = Table II (p.674): logit of militarized disputes 1816-1992 on Regime (democratic dyad) + Rivalry + peace years + 3 splines;
#   Regime -.79 (robust SE .31), Rivalry .90 (.15); N = 2,618; Huber/White robust SEs.
# Paper text used for the fixes (no other changes vs rebuild_X033.R):
#  (1) DV = dispute OCCURRENCE: "The general model regressing dispute occurrence on regime, rivalry and peace years" (p.672),
#      "dispute involvement" (p.674); predicted p = .38 at the constant -> occurrence, not onset. -> GML gmlmidongoing;
#      peace years = years since the last dispute-year (onset DV was used before).
#  (2) SEs: "Huber's robust standard errors" (p.672), Table II note "Huber/White" -> heteroskedasticity-robust, not dyad clusters.
#  (3) Major-power years exactly as listed on p.667 (Britain 1816-, France 1816-1940 & 1945-, Russia/USSR 1816-, Austria 1816-1918,
#      Prussia/Germany 1816-1945, Italy 1860-1943, Japan 1895-1945, USA 1899-, China 1950-) instead of peacesciencer COW majors
#      (which have gaps for Germany 1919-24 and Russia 1918-21).
#  (4) Regime dichotomy: democratic dyad = both Polity >= 6 (p.668), else 0; no listwise deletion for missing Polity, as implied by
#      N = 2,618 (>= all major-power dyad-years).
# Not fixable: BTSCS spline knots (Tucker program defaults not reported) -> cubic polynomial kept; Polity III -> polity2;
#   MID 2.1 -> GML dyadic MIDs; polity2 from peacesciencer::ccode_democracy; Thompson rivalry via peacesciencer (Thompson & Dreyer) rather than Thompson 2000a list.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "vb4/vhelp.R"))
source(file.path(CODE, "ea4/est_helpers.R"))
suppressPackageStartupMessages(library(peacesciencer))
WO <- file.path(OUT, "work_p4/X033"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work_p4/X033"); sink(file.path(WO, "log_est_X033.txt"), split = TRUE)
mp <- data.frame(cc = c(200, 220, 220, 365, 300, 255, 325, 740, 2, 710),
                 s  = c(1816, 1816, 1945, 1816, 1816, 1816, 1860, 1895, 1899, 1950),
                 e  = c(9999, 1940, 9999, 9999, 1918, 1945, 1943, 1945, 9999, 9999))
is_mp <- function(cc, y) sapply(seq_along(cc), function(i) any(mp$cc == cc[i] & mp$s <= y[i] & mp$e >= y[i]))
nd <- create_dyadyears(directed = FALSE, subset_years = 1816:1992) %>% as.data.frame()
nd <- nd[nd$ccode1 %in% mp$cc & nd$ccode2 %in% mp$cc, ]
nd <- nd[is_mp(nd$ccode1, nd$year) & is_mp(nd$ccode2, nd$year), ]
cat("major-power dyad-years (paper list):", nrow(nd), " (published N 2,618)\n")
nd <- nd %>% add_gml_mids() %>% add_strategic_rivalries() %>% as.data.frame()
# polity2 joined directly from peacesciencer::ccode_democracy (add_democracy() in v1.2.0 returns NA polity2 for e.g. Austria-Hungary)
pol <- ccode_democracy %>% select(ccode, year, polity2) %>% as.data.frame()
stopifnot(!anyDuplicated(pol[, c("ccode", "year")]))
nd <- nd %>% left_join(pol %>% rename(polity21 = polity2), by = c("ccode1" = "ccode", "year")) %>%
  left_join(pol %>% rename(polity22 = polity2), by = c("ccode2" = "ccode", "year"))
stopifnot("gmlmidongoing" %in% names(nd), "ongoingrivalry" %in% names(nd))
nd <- nd %>% arrange(ccode1, ccode2, year) %>% group_by(ccode1, ccode2) %>%
  mutate(spell = cumsum(lag(gmlmidongoing, default = 0) == 1) + cumsum(c(1, diff(year)) != 1)) %>%
  group_by(ccode1, ccode2, spell) %>% mutate(py = row_number() - 1) %>% ungroup() %>% as.data.frame()
nd$rival <- as.numeric(coalesce(nd$ongoingrivalry, 0) > 0)
nd$jointdem <- as.numeric(coalesce(nd$polity21 >= 6 & nd$polity22 >= 6, FALSE))
# ESTAUDIT CHANGE (only change vs rerun_verify_X033.R): temporal terms = peace years + 3 BTSCS cubic splines, as in the article
# (Table II p.674 reports Peace Years, Spline(1)-(3); p.672 "three cubic spline segments ... derived from Richard Tucker's (1999)
# program"). Basis = Tucker btscs nspline(3): spline_i = (py-k_i)_+^3 - (T-k_i)/T * py^3, knots = quartiles (Stata type 2) of py,
# T = max(py) -- the btscs rule verified against an author log in work_demcap/D1/rerun_D1.R (knots 8,19,36 reproduced).
# Knots computed on the full analysis sample (as btscs is run on the data before estimation). Replaces the cubic polynomial.
kn <- quantile(nd$py, c(.25, .5, .75), type = 2, na.rm = TRUE); Tm <- max(nd$py, na.rm = TRUE)
stopifnot(!anyNA(nd$py), length(unique(kn)) == 3, Tm > max(kn))
for (i in 1:3) nd[[paste0("spline", i)]] <- pmax(nd$py - kn[i], 0)^3 - (Tm - kn[i]) / Tm * nd$py^3
cat("btscs knots:", kn, " T:", Tm, "\n")
print(table(occurrence = nd$gmlmidongoing, jointdem = nd$jointdem)); print(table(rival = nd$rival))
cat("missing polity2 (coded non-democratic):", sum(is.na(nd$polity21) | is.na(nd$polity22)), "\n")
rhs0 <- c("jointdem", "rival", "py", "spline1", "spline2", "spline3")
FF <- function(r, s) fit_bin("gmlmidongoing", r, s, "logit", "hetero", NULL, c("jointdem", "rival"))
rep <- FF(rhs0, nd); cat("full-sample: jointdem b =", rep$jointdem["b"], " N =", attr(rep, "n"), "\n")
s <- common_sample(nd, c("gmlmidongoing", rhs0), "ccode1", "ccode2", "year")
M <- run_six(rhs0, s, FF)

# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "X033", dat = s, y = "gmlmidongoing", tie = "jointdem", ctrl = setdiff(rhs0, "jointdem"), fe = NULL,
           fit = function(dat, add) { r <- FF(c(rhs0, add), dat)$jointdem; c(b = unname(r[["b"]]), se = unname(r[["se"]]), p = unname(r[["p"]])) })
}

R <- lapply(M, function(z) z$jointdem)
f_old <- file.path(OUT, "results_verify/X033__jointdem.csv")  # port: prior result used in notes only; absent in pipeline
old <- if (file.exists(f_old)) read.csv(f_old) else list(b0 = NA_real_, p0 = NA_real_, b_dm = NA_real_, p_dm = NA_real_)
st <- "logit (fixest feglm), HC (Huber/White) SE, peace years + 3 btscs splines"
write_est("X033", "Rasler & Thompson 2001", "jointdem",
  "logit, Huber/White robust SE, peace years + 3 cubic splines (Tucker BTSCS)",
  "pdf p.672 text + Table II p.674",
  rep(st, 6), "N",
  "all six steps used peace years + py^2 + py^3 (cubic polynomial) instead of the article's BTSCS peace years + 3 cubic splines; SE (HC) and logit already matched after verify",
  "Y", R, attr(M$m0, "n"),
  sprintf(paste0("REBUILD (no replication files). btscs knots %s, T=%d (rule from D1 author-log check; R&T's knots unreported). ",
    "Full-sample N=%d; full-sample jointdem b=%.3f, rival b=%.3f (pub -.79, .90). Before fix (cubic): b0 %.4f p0 %.4f, b_dm %.4f p_dm %.4f. Common sample 1822-1992."),
    paste(kn, collapse = "/"), Tm, attr(rep, "n"), rep$jointdem["b"], rep$rival["b"], old$b0, old$p0, old$b_dm, old$p_dm))
sink()
