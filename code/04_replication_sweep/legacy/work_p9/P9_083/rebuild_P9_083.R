# P9_083 Gowa 2011, Economics & Politics 23(2) "The Democratic Peace after the Cold War" (10.1111/j.1468-0343.2011.00382.x). DEMOCRACY.
#   REBUILD (protocol B; no archive).
# Article (staged PDF pp. 157-163, Table 2): logit, robust SEs clustered on dyads, duration-dependence controls (peace years + natural cubic
#   spline), all non-directed dyad-years; MID onset; democratic dyad = both Polity2 >= 7; mixed = exactly one >= 7 (non-democratic pairs the
#   base); contiguity = land or <= 12 miles of water; log capital distance; ratio = larger CINC / sum of both; one major power; two major powers;
#   alliance = ATOP defense, neutrality or offense pact.
# HEADLINE (fixed before estimating): Table 2 col (1), 1992-2001: Democratic dyad 0.00 (0.28), p .987 (the article's claim is this null), N 123,769.
#   Verdict rule for a published coefficient of exactly 0.00: REBUILT_MATCH if the rebuilt p0 >= .05 (sign undefined), stated in advance.
# Rebuild inputs: {peacesciencer} COW dyad-years 1980-2001 (peace years need pre-sample history), GML MID onset, COW contiguity (conttype 1-2),
#   capital distance, CINC, COW majors, ATOP 5 dyad-year types; Polity IV 2014 polity2 (work_p8/P8_049). btscs peace years + 3 splines.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R")); source(file.path(CODE, "work_p9/rebuild_common_p9r.R"))
set.seed(123)
W <- file.path(OUT, "work_p9/P9_083"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_083.txt"), split = TRUE)
d <- ps_dyads(1946:2001, trade = FALSE)
d <- d %>% mutate(dyadid = 1000 * ccode1 + ccode2, mid = as.integer(gmlmidonset == 1))
d <- py_calc(d, "dyadid", "mid"); d <- d[d$year >= 1992, ]
p <- polity_cow("polity2"); n0 <- nrow(d)
d <- d %>% left_join(p %>% rename(p2a = polity2), by = c("ccode1" = "ccode", "year")) %>% left_join(p %>% rename(p2b = polity2), by = c("ccode2" = "ccode", "year"))
stopifnot(nrow(d) == n0)
d <- d %>% mutate(demdyad = as.integer(p2a >= 7 & p2b >= 7), mixed = as.integer((p2a >= 7) != (p2b >= 7)), contig = as.integer(conttype %in% 1:2),
  ldist = log(capdist + 1), ratio = pmax(cinc1, cinc2) / (cinc1 + cinc2), onemaj = as.integer(cowmaj1 + cowmaj2 == 1), twomaj = as.integer(cowmaj1 + cowmaj2 == 2),
  ally = as.integer(atop_defense == 1 | atop_neutral == 1 | atop_offense == 1))
d <- btscs_splines(d)
rhs0 <- c("demdyad", "mixed", "contig", "ldist", "ratio", "onemaj", "twomaj", "ally", "py", "spl1", "spl2", "spl3")
vars <- c("mid", rhs0, "dyadid")
cat("dyad-years 1992-2001:", nrow(d), " missing:\n"); print(sapply(d[, vars], function(x) sum(is.na(x))))
d0 <- d[complete.cases(d[, vars]), ]
cat("complete-case N:", nrow(d0), "(published 123,769)  onsets", sum(d0$mid), "\n")
TERM <- "demdyad"
fitter <- function(r, data) fit_bin("mid", r, data, "logit", "~dyadid", terms = TERM)
full <- fitter(rhs0, d0)
cat("\nREBUILT vs published: demdyad b", full[[TERM]][["b"]], "(pub 0.00) se", full[[TERM]][["se"]], "(pub 0.28) N", attr(full, "n"), "(pub 123,769)\n")
ds <- common_sample(d, vars, "ccode1", "ccode2", "year"); rm(d, d0); invisible(gc())
M <- run_six(rhs0, ds, fitter)
saveRDS(list(M = M, full = full), file.path(W, "M_P9_083.rds"))
v <- if (M$m0[[TERM]][["p"]] >= .05) "REBUILT_MATCH" else "REBUILT_MISMATCH"
meta <- list(pid = "P9_083", cite = "Gowa 2011", journal = "EP", year = 2011, doi = "10.1111/j.1468-0343.2011.00382.x", family = "DEMOCRACY")
nt <- sprintf(paste0("REBUILT (no archive): peacesciencer all dyads 1992-2001, GML MID onset; Polity IV 2014 polity2 >= 7; conttype 1-2; log capital distance; ",
  "larger CINC share; one/two COW majors; ATOP defense/neutral/offense; btscs peace years + 3 splines; logit clustered on dyad. Published headline is 0.00 ",
  "(p .987): verdict = match if rebuilt p0 >= .05. Full-sample rebuilt b=%.6f se=%.6f N=%d vs published 0.00 (0.28) N 123,769. Common N=%d."),
  full[[TERM]][["b"]], full[[TERM]][["se"]], attr(full, "n"), nrow(ds))
write_p8(meta, TERM, "abstract: dyadic dispute rates converge after the Cold War (democratic dyad n.s. 1992-2001); Table 2 col 1",
  "Table 2 col (1) (logit, MID onset, all dyads 1992-2001)", "logit, SE clustered on dyad, peace years + spline", "feglm logit, SE clustered on dyad, all six steps",
  0, 0.28, 123769, FALSE, TRUE, M, "Y (peace years + spline)", nt, verdict = v)
dir.create(file.path(OUT, "results_p9"), recursive = TRUE, showWarnings = FALSE)
src <- file.path(OUT, "results_p8", "P9_083__demdyad.csv"); stopifnot(file.exists(src), file.rename(src, file.path(OUT, "results_p9", basename(src))))
sink()
