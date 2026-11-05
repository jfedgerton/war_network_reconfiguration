# P9_076 Schramm 2022, Journal of Global Security Studies 7(4) "Tyrants and Terror? Liberal Democracies and Personalist Regimes" /
#   liberal-democratic initiation against personalist regimes (10.1093/jogss/ogac044). DEMOCRACY. REBUILD (protocol B; no archive).
# Article (staged PDF pp. 8-10, Table 1): directed dyad-years (ending 2008), logit with robust SEs and peace years cubic; DV = MID
#   initiation by A (Gibler et al. MID data; col 1 = all MIDs); State B personalist = GWF personalist incl. personalist hybrids;
#   State A liberal democracy = V-Dem liberal democracy index; controls contiguity, trade dependency, CINC A and B, power configuration
#   (maj-maj, maj-min, min-maj; min-min base), defense / neutrality / nonaggression / entente agreements, peace years + square + cube.
#   Legend ***P<0.01 **P<0.05 *P<0.1.
# HEADLINE (fixed before estimating): Table 1 col 1 (MIDs), State B personalist x State A liberal democracy 0.851 (0.487) *, N 275,159.
# Rebuild inputs: {peacesciencer} COW directed dyad-years 1946-2008, GML MIDs (onset with A as initiator, init1), contiguity (conttype 1),
#   CINC, majors, COW trade over PWT GDP of A; COW alliance v4.1 types (work_p8/P8_049); GWF 1.2 (work_p9/gwf); V-Dem v6.2 libdem
#   (work_p9/P9_062 staged file; paper V-Dem v11). Politically relevant directed dyads (contiguous <= 150 miles or a major power) - the
#   article's N (275,159 over 1946-2008) is of that order; the sample rule is not stated (assumption documented).
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R")); source(file.path(CODE, "work_p9/rebuild_common_p9r.R"))
set.seed(123)
W <- file.path(OUT, "work_p9/P9_076"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_076.txt"), split = TRUE)
d <- ps_dyads(1946:2008, directed = TRUE) %>% add_cow_ally()
d <- d %>% mutate(ddyad = 1000 * ccode1 + ccode2, init = as.integer(gmlmidonset == 1 & init1 %in% 1))
d <- py_calc(d, "ddyad", "init")
d <- d[pr_flag(d) == 1, ]
g <- gwf_cy(); v <- vdem_libdem(); n0 <- nrow(d)
d <- d %>% left_join(g %>% transmute(ccode, year, persB = pers), by = c("ccode2" = "ccode", "year")) %>%
  left_join(v %>% rename(libA = libdem), by = c("ccode1" = "ccode", "year"))
stopifnot(nrow(d) == n0)
d <- d %>% mutate(persB = coalesce(persB, 0L), pers_x_lib = persB * libA, contig = as.integer(conttype == 1),
  tdep = (flow1 + flow2) / pwtrgdp1, majmaj = as.integer(cowmaj1 == 1 & cowmaj2 == 1), majmin = as.integer(cowmaj1 == 1 & cowmaj2 == 0),
  minmaj = as.integer(cowmaj1 == 0 & cowmaj2 == 1), py2 = py^2, py3 = py^3)
d$tdep[!is.finite(d$tdep)] <- NA
rhs0 <- c("persB", "libA", "pers_x_lib", "contig", "tdep", "cinc1", "cinc2", "majmaj", "majmin", "minmaj", "cow_def", "cow_neu", "cow_nonagg", "cow_ent", "py", "py2", "py3")
vars <- c("init", rhs0, "ddyad")
cat("PR directed dyad-years 1946-2008:", nrow(d), " missing:\n"); print(sapply(d[, vars], function(x) sum(is.na(x))))
d0 <- d[complete.cases(d[, vars]), ]
cat("complete-case N:", nrow(d0), "(published 275,159)  initiations", sum(d0$init), "\n")
TERM <- "pers_x_lib"
fitter <- function(r, data) fit_bin("init", r, data, "logit", "hetero", terms = TERM)
full <- fitter(rhs0, d0)
cat("\nREBUILT vs published: pers_x_lib b", full[[TERM]][["b"]], "(pub 0.851) se", full[[TERM]][["se"]], "(pub 0.487) N", attr(full, "n"), "(pub 275,159)\n")
ds <- common_sample(d, vars, "ccode1", "ccode2", "year"); rm(d, d0); invisible(gc())
M <- run_six(rhs0, ds, fitter)
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_split(id = "P9_076", dat = ds, fit = function(dat, add) fitter(c(rhs0, add), dat)[[TERM]])
}
saveRDS(list(M = M, full = full), file.path(W, "M_P9_076.rds"))
meta <- list(pid = "P9_076", cite = "Schramm 2022", journal = "JOGSS", year = 2022, doi = "10.1093/jogss/ogac044", family = "DEMOCRACY")
nt <- sprintf(paste0("REBUILT (no archive): peacesciencer PR directed dyads 1946-2008 (sample rule assumed), GML MID onset initiated by A; GWF 1.2 personalist incl. ",
  "hybrids for B; V-Dem v6.2 libdem for A (paper v11); land contiguity, A's COW trade / PWT GDP, CINC A/B, power configuration, COW alliance types, ",
  "peace years cubic; robust (HC) SE. Full-sample rebuilt b=%.6f se=%.6f N=%d vs published 0.851 (0.487) N 275,159. Common sample N=%d. ",
  "Interaction: constituents kept in all steps."), full[[TERM]][["b"]], full[[TERM]][["se"]], attr(full, "n"), nrow(ds))
write_p8(meta, TERM, "abstract: liberal democracies more likely to initiate conflict against personalist regimes; Table 1 col 1",
  "Table 1 col 1 (MIDs; logit, directed dyads 1946-2008)", "logit, robust SE, peace years cubic", "feglm logit, HC vcov, all six steps",
  0.851, 0.487, 275159, FALSE, TRUE, M, "Y (peace years cubic)", nt)
dir.create(file.path(OUT, "results_p9"), recursive = TRUE, showWarnings = FALSE)
src <- file.path(OUT, "results_p8", "P9_076__pers_x_lib.csv"); stopifnot(file.exists(src), file.rename(src, file.path(OUT, "results_p9", basename(src))))
sink()
