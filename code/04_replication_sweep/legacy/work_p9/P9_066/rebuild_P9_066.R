# P9_066 Kim & James 2010, Defence and Peace Economics 21(5-6) "Learning from the Great War? The Origin of Democratic Peace"
#   (10.1080/10242694.2010.513493). DEMOCRACY. REBUILD (protocol B; no archive).
# Article (staged PDF, Table II, text pp. 497-498): binary logit with BKT cubic splines (no peace-years term shown), conventional SEs,
#   politically relevant dyads (contiguous or a major power), MID onset (MID 3.0). Democratic dyad = both Polity2 >= 6 (fn. 6).
#   Equal capability = 1 if the dyadic CINC ratio lies in [0.8, 1.2]. Alliance coded 1 when the dyad has NO COW formal alliance
#   (text p. 498). Contiguity = COW contiguity categories 1-4 (land, 12, 24, 150 miles) vs 5-6.
# HEADLINE (fixed before estimating): Table II post-First World War column, Democratic dyad -0.188 (0.049) ***, N 61,217.
# Rebuild inputs: {peacesciencer} COW dyad-years 1919-2000, GML MIDs (stand-in for MID 3.0), Polity (2017), NMC CINC, COW contiguity,
#   COW majors; COW alliance v4.1 dyad-years (work_p8/P8_049). Splines: btscs rule (quartile knots) on peace years within dyad.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R")); source(file.path(CODE, "work_p9/rebuild_common_p9r.R"))
set.seed(123)
W <- file.path(OUT, "work_p9/P9_066"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_066.txt"), split = TRUE)
d <- ps_dyads(1919:2000, trade = FALSE) %>% add_cow_ally()
d <- d[pr_flag(d) == 1, ]
d <- d %>% mutate(dyadid = 1000 * ccode1 + ccode2, mid = as.integer(gmlmidonset == 1),
  demdyad = as.integer(polity21 >= 6 & polity22 >= 6), cr = cinc1 / cinc2, eqcap = as.integer(cr >= 0.8 & cr <= 1.2),
  noally = as.integer(cow_any == 0), contig = as.integer(conttype %in% 1:4))
d <- py_calc(d, "dyadid", "mid") %>% btscs_splines()
rhs0 <- c("demdyad", "eqcap", "noally", "contig", "spl1", "spl2", "spl3")
vars <- c("mid", rhs0, "dyadid")
cat("PR dyad-years 1919-2000:", nrow(d), " missing:\n"); print(sapply(d[, vars], function(x) sum(is.na(x))))
d0 <- d[complete.cases(d[, vars]), ]
cat("complete-case N:", nrow(d0), "(published 61,217)  MID onsets", sum(d0$mid), "\n")
TERM <- "demdyad"
fitter <- function(r, data) fit_bin("mid", r, data, "logit", "iid", terms = TERM)
full <- fitter(rhs0, d0)
cat("\nREBUILT vs published: demdyad b", full[[TERM]][["b"]], "(pub -0.188) se", full[[TERM]][["se"]], "(pub 0.049) N", attr(full, "n"), "(pub 61,217)\n")
ds <- common_sample(d, vars, "ccode1", "ccode2", "year")
M <- run_six(rhs0, ds, fitter)
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "P9_066", dat = ds, y = "mid", tie = TERM, ctrl = setdiff(rhs0, TERM), fe = NULL,
           fit = function(dat, add) fitter(c(rhs0, add), dat)[[TERM]])
}
saveRDS(list(M = M, full = full), file.path(W, "M_P9_066.rds"))
meta <- list(pid = "P9_066", cite = "Kim & James 2010", journal = "DPE", year = 2010, doi = "10.1080/10242694.2010.513493", family = "DEMOCRACY")
nt <- sprintf(paste0("REBUILT (no archive): peacesciencer PR dyads 1919-2000 (contiguity 1-4 or COW major); GML MID onset (paper MID 3.0); ",
  "demdyad = both Polity2 >= 6; eqcap = CINC ratio in [0.8,1.2]; noally = no COW alliance (paper codes alliance=1 for no alliance); ",
  "contig = COW conttype 1-4; btscs splines (quartile knots). Full-sample rebuilt b=%.6f se=%.6f N=%d vs published -0.188 (0.049) N 61,217. ",
  "Common 1822-2011 sample N=%d."), full[[TERM]][["b"]], full[[TERM]][["se"]], attr(full, "n"), nrow(ds))
write_p8(meta, TERM, "abstract and Table II: democratic dyad negative only after the First World War",
  "Table II post-First World War column (logit, MID onset, PR dyads 1919-2000)", "logit, conventional SE, BKT splines",
  "feglm logit, iid vcov, all six steps", -0.188, 0.049, 61217, TRUE, TRUE, M, "N (splines only)", nt)
dir.create(file.path(OUT, "results_p9"), recursive = TRUE, showWarnings = FALSE)
src <- file.path(OUT, "results_p8", "P9_066__demdyad.csv"); stopifnot(file.exists(src), file.rename(src, file.path(OUT, "results_p9", basename(src))))
sink()
