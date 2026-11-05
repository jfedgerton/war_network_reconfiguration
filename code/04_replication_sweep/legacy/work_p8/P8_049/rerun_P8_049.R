# P8_049 Chatagnier & Castelli 2016 PRQ "A Modern Peace?" (10.1177/1065912916670270)
# Data: authors' replication archive (jtchatagnier.com/Uploads/modern_peace_replication_materials.zipx); dyadic_data_set.dta built by
#   the authors' modern_peace_prepare_data.r (only edit: library(reshape) -> library(reshape2) for melt; current countrycode version).
# Paper (author PDF, air.unipr.it): Table 2, Model 1 (Baseline), politically relevant dyads 1960-2002, fatal MID onset, logit,
#   dyad-clustered robust SE, t, t^2, t^3 (years since last fatal MID). Published IndustryH = -0.05 (0.02), N = 272,682.
# KEY TERM (a priori, H2): ind_hi = industrial share of GDP of the more-industrialized state ("IndustryH").
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/p8_common.R"))
suppressMessages(library(haven))
WD <- file.path(ROOT, "work_p8/P8_049"); W <- file.path(OUT, "work_p8/P8_049"); dir.create(W, recursive = TRUE, showWarnings = FALSE); sink(file.path(W, "log_P8_049.txt"), split = TRUE)
cols <- c("ccode1", "ccode2", "year", "dyadid", "fatal_ongo", "fatal_onset", "ind_hi", "serv_ind_hi", "gdp_lo", "Mrat", "majpow",
          "lodem", "hidem", "contig", "alliance", "s2uni", "lntrade", "t_fatal")
# port: dyadic_data_set.dta is built (into OUT) by the ported authors' prep script from raw inputs under DATA
built <- file.path(W, "dyadic_data_set.dta")
if (!file.exists(built)) {
  rc <- system2(file.path(R.home("bin"), "Rscript"), shQuote(file.path(CODE, "work_p8/P8_049/prep_p8.R")),
                stdout = file.path(W, "prep_p8.out"), stderr = file.path(W, "prep_p8.out"))
  stopifnot(identical(as.integer(rc), 0L), file.exists(built))
}
d <- read_dta(built, col_select = all_of(cols)) %>% zap_labels() %>% as.data.frame()
d <- d[d$fatal_ongo == 0, ]                                    # authors: drop dyad-years with ongoing conflicts
d$lntrade <- ifelse(is.na(d$lntrade), 0, d$lntrade)            # authors: missing trade -> 0
d$t2_fatal <- d$t_fatal^2; d$t3_fatal <- d$t_fatal^3
d$dyad_n <- as.integer(factor(d$dyadid))
ndup <- sum(duplicated(d[, c("ccode1", "ccode2", "year")])); cat("duplicate dyad-years in authors data (kept as in authors models):", ndup, "of", nrow(d), "\n")
rhs <- c("ind_hi", "serv_ind_hi", "gdp_lo", "Mrat", "majpow", "lodem", "hidem", "contig", "alliance", "lntrade", "s2uni",
         "t_fatal", "t2_fatal", "t3_fatal")
vars <- c("fatal_onset", rhs)
fitter <- function(r, data) fit_bin("fatal_onset", r, data, "logit", "~dyad_n", "ind_hi")
d0 <- d[complete.cases(d[, vars]), ]
cat("complete-case N:", nrow(d0), "(published 272,682)\n")
stopifnot(nrow(d0) == 280152)  # port: complete-case N of the original sweep run (log_P8_049.txt)
full <- fitter(rhs, d0)
ds <- common_sample(d, vars, "ccode1", "ccode2", "year")
M <- run_six(rhs, ds, fitter)

# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "P8_049", dat = ds, y = "fatal_onset", tie = "ind_hi", ctrl = setdiff(rhs, "ind_hi"), fe = NULL,
           fit = function(dat, add) fitter(c(rhs, add), dat)[["ind_hi"]])
}

meta <- list(pid = "P8_049", short_cite = "Chatagnier & Castelli 2016", journal = "PRQ", year = 2016, doi = "10.1177/1065912916670270",
  family = "CAPITALISM", headline_table_model = "Table 2, Model 1 (dyadic baseline)",
  article_estimator = "logit, dyad-clustered robust SE, cubic peace-years (t, t2, t3)",
  est_used_all_steps = "logit, dyad-clustered SE (fixest) in all six steps", author_has_history_control = "N")
write_p8(meta, "ind_hi", "H2: 'industrialization of the more-heavily-industrialized state ... negatively related' to violent conflict",
  -0.05, 0.02, 272682, TRUE, FALSE, M,
  paste0("Authors' data rebuilt with authors' prep code (reshape->reshape2; current countrycode may map WB names slightly differently). ",
         "Full-sample b = ", signif(full$ind_hi["b"], 3), " (SE ", signif(full$ind_hi["se"], 3), ", N ", attr(full, "n"), "). ",
         "Published coefficient printed to 2 decimals (-0.05), so the 10% rule is coarse. t_fatal from authors data kept unchanged. Authors data contain ", ndup, " duplicate dyad-years (kept, as in their glm)."))
sink()
