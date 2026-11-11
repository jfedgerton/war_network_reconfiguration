# P9_078 Hewitt 2003, JCR 47(5) "Dyadic Processes and International Crises" (10.1177/0022002703252973). TIE. REBUILD (protocol B).
# Article (staged PDF, Table 4 col c): logistic regression (z statistics; Russett-Oneal Triangulating Peace design), politically relevant
#   dyad-years 1918-1992, DV = involvement in an international crisis (ICB crisis dyads, Hewitt's dyadic data); Kantian variables
#   (dyadic democracy, dyadic economic interdependence, joint IGO memberships) and controls (capability ratio, alliances, geographically
#   noncontiguous, distance, minor-power dyad); one-sided tests (* p<.01, ** p<.001).
# HEADLINE (fixed before estimating): Table 4 col (c), Joint IGO memberships -0.0201 (z -2.849) *, N 35,740 (SE = 0.0201/2.849).
# Rebuild inputs: Russett-Oneal TRIANGLE.tab (PR dyads; demaut, depend, smigoabi, lcaprat2, allies, noncontg, logdstab, minrpwrs; R&O reverse
#   the contiguity and major-power scales: noncontiguous = noncontg + 1, minor-power dyad = minrpwrs + 1) + ICB dyadic data v16
#   (work_p9/P9_078/icb_dyads_v16.csv; Hewitt's dyads, annual observations over each confrontation). Logit, SEs clustered on dyad.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R")); source(file.path(CODE, "work_p9/rebuild_common_p9r.R"))
set.seed(123)
W <- file.path(OUT, "work_p9/P9_078"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_078.txt"), split = TRUE)
t <- read.delim(file.path(ROOT, "recA/or/TRIANGLE.tab"))
stopifnot(nrow(t) == 39996, all(t$statea < t$stateb), all(t$noncontg %in% c(-1, 0)), all(t$minrpwrs %in% c(-1, 0)))
icb <- read.csv(file.path(ROOT, "work_p9/P9_078/icb_dyads_v16.csv"))
icb <- icb %>% transmute(statea = pmin(statea, stateb), stateb = pmax(icb$statea, icb$stateb), year) %>% distinct() %>% mutate(crisis = 1L)
n0 <- nrow(t)
d <- t %>% left_join(icb, by = c("statea", "stateb", "year")) %>% mutate(crisis = coalesce(crisis, 0L))
stopifnot(nrow(d) == n0)
cat("ICB crisis dyad-years matched to TRIANGLE PR rows:", sum(d$crisis), "of", nrow(icb), "\n")
d <- d %>% filter(year >= 1918, year <= 1992) %>% mutate(dyadid = 1000 * statea + stateb, demlo = pmin(demauta, demautb), deplo = pmin(dependa, dependb),
  noncont = noncontg + 1, minor = minrpwrs + 1)
rhs0 <- c("demlo", "deplo", "smigoabi", "lcaprat2", "allies", "noncont", "logdstab", "minor")
vars <- c("crisis", rhs0, "dyadid")
d0 <- d[complete.cases(d[, vars]), ]
cat("complete-case N:", nrow(d0), "(published 35,740)  crisis dyad-years", sum(d0$crisis), "\n")
TERM <- "smigoabi"
# [2026-11-10 audit2] article estimates Table 4 by GEE, AR(1), robust SE (fn 12): geeglm logit, AR(1) spaced by calendar year
#   (as X092), id = dyad, singleton dyads dropped; replaces the pooled logit. Rebuilt full sample: GEE b -0.0352 (0.0076), N 35,740
#   (published N 35,740); 1.75x published -0.0201, inside the round-2 0.5-2x rule.
suppressPackageStartupMessages(library(geepack))
fitter <- function(r, data) {
  g <- data[complete.cases(data[, c("crisis", r, "dyadid", "year")]), ]
  g <- g[order(g$dyadid, g$year), ]
  tb <- table(g$dyadid); g <- g[g$dyadid %in% as.numeric(names(tb)[tb > 1]), ]
  m <- geeglm(as.formula(paste("crisis ~", paste(r, collapse = " + "))), data = g, id = dyadid, waves = year, family = binomial, corstr = "ar1")
  ct <- summary(m)$coefficients; stopifnot(TERM %in% rownames(ct))
  cat("\n---- GEE AR(1) | N =", nrow(g), "\n"); print(round(as.matrix(ct), 6))
  out <- list(); out[[TERM]] <- c(b = ct[TERM, 1], se = ct[TERM, 2], p = ct[TERM, 4]); attr(out, "n") <- nrow(g)
  rm(m); invisible(gc()); out
}
full <- fitter(rhs0, d0)
cat("\nREBUILT vs published: smigoabi b", full[[TERM]][["b"]], "(pub -0.0201) se", full[[TERM]][["se"]], "(pub 0.00706) N", attr(full, "n"), "(pub 35,740)\n")
ds <- common_sample(d, vars, "statea", "stateb", "year")
M <- run_six(rhs0, ds, fitter)
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  # [2026-11-11] dat = the rows the GEE fitter estimates on (complete cases, singleton dyads dropped), so the helper N equals the model N
  dsg <- ds[complete.cases(ds[, c("crisis", rhs0, "dyadid", "year")]), ]; tb <- table(dsg$dyadid); dsg <- dsg[dsg$dyadid %in% as.numeric(names(tb)[tb > 1]), ]
  sens_run(id = "P9_078", dat = dsg, y = "crisis", tie = TERM, ctrl = setdiff(rhs0, TERM), fe = NULL,
           fit = function(dat, add) fitter(c(rhs0, add), dat)[[TERM]])
}
saveRDS(list(M = M, full = full), file.path(W, "M_P9_078.rds"))
meta <- list(pid = "P9_078", cite = "Hewitt 2003", journal = "JCR", year = 2003, doi = "10.1177/0022002703252973", family = "TIE")
nt <- sprintf(paste0("REBUILT (no archive): TRIANGLE PR dyads 1918-92; DV = ICB crisis dyad-year (ICB dyads v16, Hewitt's dyadic data); lower demaut, ",
  "lower dependence, smigoabi, lcaprat2, allies, noncontiguous, logdstab, minor-power dyad; GEE logit AR(1) (calendar year), robust SE by dyad. Full-sample rebuilt b=%.6f se=%.6f N=%d ",
  "vs published -0.0201 (z -2.849; SE 0.00706) N 35,740. Common sample N=%d."), full[[TERM]][["b"]], full[[TERM]][["se"]], attr(full, "n"), nrow(ds))
write_p8(meta, TERM, "abstract: two of three Kantian variables (incl. IGO membership) relate to crisis onset as to MIDs; Table 4 col c",
  "Table 4 col (c) (logit, ICB crisis involvement, PR dyads 1918-1992)", "logistic regression (Triangulating Peace design), one-sided tests",
  "geeglm logit, AR(1) by calendar year, robust SE by dyad, all six steps", -0.0201, 0.0201 / 2.849, 35740, TRUE, TRUE, M, "N (no temporal terms)", nt)
dir.create(file.path(OUT, "results_p9"), recursive = TRUE, showWarnings = FALSE)
src <- file.path(OUT, "results_p8", "P9_078__smigoabi.csv"); stopifnot(file.exists(src), file.rename(src, file.path(OUT, "results_p9", basename(src))))
sink()
