# P9_088 Xiang 2010, JOP 72(2) "Relevance as a Latent Variable in Dyadic Analysis of Conflict" (10.1017/S0022381609990909). TIE.
#   REBUILD (protocol B; no archive). Estimator: split-population (partial-observability) bivariate probit, ML (work_p9/splitpop_p9.R).
# Article (staged PDF pp. 489-492, Table 1, rasterised): all dyad-years 1870-1992 (Barbieri's trade data), Maoz dyadic MID onset (ongoing
#   MIDs excluded); conflict equation: Trade (Barbieri salience = sqrt of the product of the two trade shares), Major Power, Alliance
#   (defense/neutrality/entente), Contiguity (land or <= 150 miles water), Joint Democracy (product of transformed Polity scores), Relative
#   Capabilities, Distance, peace years + square + cube; relevance equation: Trade, Major Power, Alliance, Contiguity, Distance. EUGene data.
#   Legend * p <= .05, dagger p <= .10.
# HEADLINE (fixed before estimating): Table 1 Model 2 (split population), conflict-equation Trade -0.742 (0.435) dagger (p <= .10), N 137,314.
#   The text calls it "statistically significant at 10% level" -> alpha_pub 0.10. Model 1 (probit, all dyads): Trade 0.114 (0.329).
# Rebuild inputs: as P9_081 (Xiang, Xu & Keteku 2007, same Barbieri sample): DOR/O-R all-dyad file + derived trade file (mzmid1_orig = Maoz MID
#   onset in t+1 against covariates at t; history at outcome_year), COW national trade totals, capital distance ({peacesciencer}).
#   Joint democracy = (polity2_a + 10)(polity2_b + 10); relative capabilities = log CINC ratio; distance = log capital distance; peace years
#   from or_trade (py) with square and cube. History terms enter both equations (heckprob precedent P9_035/P9_045).
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R")); source(file.path(CODE, "work_p9/rebuild_common_p9r.R")); source(file.path(CODE, "work_p9/splitpop_p9.R"))
set.seed(123)
W <- file.path(OUT, "work_p9/P9_088"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_088.txt"), split = TRUE)
splitpop_selftest()
d <- or_dyads(1869:1991)
cd <- suppressMessages(create_dyadyears(system = "cow", directed = FALSE, subset_years = 1869:1991) %>% add_capital_distance()) %>% as.data.frame()
nt_ <- cow_nat_trade(); n0 <- nrow(d)
d <- d %>% left_join(cd[, c("ccode1", "ccode2", "year", "capdist")], by = c("statea" = "ccode1", "stateb" = "ccode2", "year")) %>%
  left_join(nt_ %>% rename(tota = tot), by = c("statea" = "ccode", "year")) %>% left_join(nt_ %>% rename(totb = tot), by = c("stateb" = "ccode", "year"))
stopifnot(nrow(d) == n0)
d <- d %>% mutate(salience = sqrt((trade / tota) * (trade / totb)), jdem = (polity2a + 10) * (polity2b + 10), ldcap = log(capdist + 1),
  py2 = py^2, py3 = py^3, outcome_year = year + 1)
d$salience[!is.finite(d$salience)] <- NA
conf <- c("salience", "majpow", "anyally", "contig150", "jdem", "lcaprat", "ldcap", "py", "py2", "py3")
relv <- c("salience", "majpow", "anyally", "contig150", "ldcap")
vars <- c("mzmid1_orig", unique(c(conf, relv)), "dyadid")
d0 <- d[complete.cases(d[, vars]), ]
cat("complete-case N:", nrow(d0), "(published 137,314)  onsets", sum(d0$mzmid1_orig), "\n")
pm <- glm(as.formula(paste("mzmid1_orig ~", paste(conf, collapse = " + "))), data = d0, family = binomial("probit"))
cat("Model 1 probit check: salience", coef(pm)[["salience"]], "se", sqrt(diag(vcov(pm)))[["salience"]], "(published 0.114 (0.329))\n")
TERM <- "salience"
fitter <- function(r, data) {
  extra <- setdiff(r, conf)
  X <- cbind(`(Intercept)` = 1, as.matrix(data[, c(conf, extra)])); Z <- cbind(`(Intercept)` = 1, as.matrix(data[, c(relv, extra)]))
  f <- fit_splitpop(data$mzmid1_orig, X, Z)
  z <- f$coef / f$se
  cat("\n---- split-population probit | extra:", paste(extra, collapse = " + "), "| N =", f$n, "| LL =", round(f$ll, 2), "| rho =", round(f$rho, 3), "| max|grad| =", signif(f$grad_max, 3), "\n")
  print(round(cbind(b = f$coef, se = f$se, z = z), 5))
  out <- list(c(b = f$coef[["out:salience"]], se = f$se[["out:salience"]], p = 2 * pnorm(-abs(z[["out:salience"]])))); names(out) <- TERM
  attr(out, "n") <- f$n; out
}
full <- fitter(conf, d0)
cat("\nREBUILT vs published: out:salience b", full[[TERM]][["b"]], "(pub -0.742) se", full[[TERM]][["se"]], "(pub 0.435) N", attr(full, "n"), "(pub 137,314)\n")
ds <- common_sample(d, vars, "statea", "stateb", "outcome_year"); rm(d, d0); invisible(gc())
M <- run_six(conf, ds, fitter)
saveRDS(list(M = M, full = full), file.path(W, "M_P9_088.rds"))
b0 <- M$m0[[TERM]][["b"]]; p0 <- M$m0[[TERM]][["p"]]; r_ <- b0 / -0.742; pub_p <- 2 * pnorm(-abs(-0.742 / 0.435))
v <- if (sign(b0) == -1 && (pub_p >= .10 || p0 < .10) && r_ >= 0.5 && r_ <= 2) "REBUILT_MATCH" else "REBUILT_MISMATCH"
meta <- list(pid = "P9_088", cite = "Xiang 2010", journal = "JOP", year = 2010, doi = "10.1017/S0022381609990909", family = "TIE")
nt <- sprintf(paste0("REBUILT (no archive): P9_081 frame (DOR/O-R all dyads, IVs 1869-1991, DV mzmid1_orig onset t+1); split-population bivariate probit by ML ",
  "(work_p9/splitpop_p9.R, Genz BVN, OIM SE); history terms in both equations. Full-sample rebuilt b=%.6f se=%.6f N=%d vs published -0.742 (0.435) N 137,314. ",
  "Common N=%d. Round-2 rule (alpha .10, 0.5-2x)."), full[[TERM]][["b"]], full[[TERM]][["se"]], attr(full, "n"), nrow(ds))
write_p8(meta, TERM, "abstract and Table 1 Model 2: with relevance modelled, trade reduces conflict (significant at 10%)",
  "Table 1 Model 2 (split-population probit, conflict equation; MID onset, all dyads 1870-1992)", "split-population (partial-observability) bivariate probit, ML",
  "ML split-population probit (work_p9/splitpop_p9.R), OIM SE, all six steps", -0.742, 0.435, 137314, FALSE, TRUE, M, "Y (peace years cubic)", nt, verdict = v)
dir.create(file.path(OUT, "results_p9"), recursive = TRUE, showWarnings = FALSE)
src <- file.path(OUT, "results_p8", "P9_088__salience.csv"); stopifnot(file.exists(src), file.rename(src, file.path(OUT, "results_p9", basename(src))))
sink()
