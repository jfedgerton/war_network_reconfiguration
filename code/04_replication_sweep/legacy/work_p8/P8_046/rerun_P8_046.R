# P8_046 McDonald 2010 II "Capitalism, Commitment, and Peace" (10.1080/03050621003784911)
# Data: author's Dataverse (doi:10.7910/DVN/IEAASS) CCPIIdirdyadreplicationdata + CCPIIreplicationlog.log.
# Directed-dyad models in the author's log: logit, war targeting of B (wartargetb1970), 1970-1997 (if year<1998), dyad-clustered SE,
#   peace-years + 3 splines. Model A (target attributes): lpublicb = .0459701 (.0111006), N = 360,255.
#   Model B (both sides): lpublicb = .0277735 (.009868), N = 210,314.
# KEY TERM (a priori, abstract/title: public property -> being targeted): lpublicb (target's public property share, lagged).
# Both specifications estimated; headline chosen from the paper (see write step).
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/p8_common.R"))
suppressMessages(library(haven))
WD <- file.path(ROOT, "work_p8/P8_046"); W <- file.path(OUT, "work_p8/P8_046"); dir.create(W, recursive = TRUE, showWarnings = FALSE); sink(file.path(W, "log_P8_046.txt"), split = TRUE)
d <- read_dta(file.path(WD, "dirdyad.dta")) %>% zap_labels() %>% as.data.frame()
d <- d[d$year < 1998, ]
stopifnot(!anyDuplicated(d[, c("ccode1", "ccode2", "year")]))
sp <- c("wartargetbspl", "wartargetb1970sp1", "wartargetb1970sp2", "wartargetb1970sp3")
specs <- list(
  A = c(sp, "lpublicb", "lpolity22", "ldep2b", "cont4", "lndistan", "lncapratio", "majpow2", "s_un_glo", "ally"),
  B = c(sp, "lpublica", "lpublicb", "lpolity21", "lpolity22", "ldep2a", "ldep2b", "cont4", "lndistan", "lncapratio", "majpow1", "majpow2", "s_un_glo", "ally"))
out <- list()
for (s in names(specs)) {
  rhs <- specs[[s]]; vars <- c("wartargetb1970", rhs)
  dd <- d
  if (s == "B") {  # Stata: 'majpow1 != 0 predicts failure perfectly; majpow1 dropped and 8486 obs not used'
    dd <- dd[complete.cases(dd[, vars]), ]; dd <- dd[dd$majpow1 == 0, ]; rhs <- setdiff(rhs, "majpow1"); vars <- setdiff(vars, "majpow1")
  }
  fitter <- function(r, data) fit_bin("wartargetb1970", r, data, "logit", "~dyadid", "lpublicb")
  d0 <- dd[complete.cases(dd[, vars]), ]; cat("\n==== spec", s, "complete-case N:", nrow(d0), "\n")
  full <- fitter(rhs, d0)
  ds <- common_sample(dd, vars, "ccode1", "ccode2", "year")
  out[[s]] <- list(M = run_six(rhs, ds, fitter), full = full)

  # ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
  if (Sys.getenv("SWEEP_SENS") == "1") {
    source(file.path(CODE, "hist/sens_helper.R"))
    if (s == "A")   # headline = spec A (target attributes); rhs, ds, fitter are spec A's in this iteration
      sens_run(id = "P8_046", dat = ds, y = "wartargetb1970", tie = "lpublicb", ctrl = setdiff(rhs, "lpublicb"), fe = NULL,
               fit = function(dat, add) fitter(c(rhs, add), dat)[["lpublicb"]])
  }

}
saveRDS(out, file.path(W, "M_P8_046.rds"))
# ---- port: result-CSV write step for this study, copied verbatim from work_p8/write_pending.R (paths: saved model list read from OUT) ----
{

  o <- readRDS(file.path(W, "M_P8_046.rds"))
  meta <- list(pid = "P8_046", short_cite = "McDonald 2010", journal = "II", year = 2010, doi = "10.1080/03050621003784911",
    family = "CAPITALISM", headline_table_model = "Directed-dyad war-targeting model A (target attributes; author's log, first dyadic model)",
    article_estimator = "logit, dyad-clustered robust SE, peace years + 3 splines (author's replication log)",
    est_used_all_steps = "logit, dyad-clustered SE (fixest), all six steps", author_has_history_control = "N (peace years/splines only)")
  ev <- "Abstract: predominance of public property raises 'the likelihood of being the target of military conflict'"
  write_p8(meta, "lpublicb", ev, 0.0459701, 0.0111006, 360255, NA, FALSE, o$A$M,
    "Author's Dataverse data; published values from the author's own replication log (CCPIIreplicationlog.log); article PDF not accessible, so the table/column label is inferred from the log order. DV = war targeting of B, 1970-97.",
    full = o$A$full)
  meta$headline_table_model <- "Directed-dyad war-targeting model B (both sides' attributes; author's log, second dyadic model)"
  rn <- function(x) { n <- attr(x, "n"); names(x)[names(x) == "lpublicb"] <- "lpublicb_bothsides"; attr(x, "n") <- n; x }
  MB <- lapply(o$B$M, rn); FB <- rn(o$B$full)
  write_p8(meta, "lpublicb_bothsides", ev, 0.0277735, 0.009868, 210314, NA, FALSE, MB,
    "Second dyadic specification (both states' attributes; Stata dropped majpow1 and 8,486 obs as perfect predictors, replicated). Values from the author's replication log.",
    full = FB)
}
sink()
