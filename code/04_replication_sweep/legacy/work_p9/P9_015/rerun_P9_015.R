# P9_015 Neumayer & Plumper 2011 JPR 48(1) "Foreign terror on Americans" (10.1177/0022343310390147)
# Data + do-file: Harvard Dataverse doi:10.7910/DVN/DFNTK0 ("Article for JPR (terrorism).dta", "Article for JPR (terrorism).do");
#   staged as np_data.dta / np_main.do. No author log in the archive.
# Directed country dyads (terrorists' home country -> victims' nationality); the article's models use the US-victim slice:
#   countrytarget == "United States" & countryallterrorists != "United States" (one row per origin country-year, 1978-2004 after NA).
# HEADLINE (abstract; text: "US military aid matters most"): Table II model k1 =
#   "nbreg sumuskilled lnpopterr lndistance lngdppcconstterrorists polity2terrorists usmilaid_milexpcomb, cluster(countryallterrorists)".
#   KEY TERM usmilaid_milexpcomb (US military aid relative to the home country's military expenditure).
#   Published (LSE Research Online accepted manuscript, Table II): .055 (SE .019), N = 3360.
# Estimator: NB2 ML, SE clustered on terrorists' home country -> fixest::fenegbin, vcov ~countryallterrorists, all six steps.
# History linkage: dyad (US = 2, origin COW code). Origin names -> COW via countrycode; Germany <= 1990 -> 260,
#   Czech Republic <= 1992 -> 315 (linkage only; the estimation sample and variables are the author's).
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R"))
WD <- file.path(ROOT, "work_p9/P9_015"); W <- file.path(OUT, "work_p9/P9_015"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_015.txt"), split = TRUE)
d <- as.data.frame(haven::zap_labels(haven::read_dta(file.path(WD, "np_data.dta"))))
stopifnot(nrow(d) == 84559)
d <- d[d$countrytarget == "United States" & d$countryallterrorists != "United States", ]
stopifnot(nrow(d) == 8452, !anyDuplicated(d[, c("countryallterrorists", "year")]))
d$ccode_us <- 2
d$ccode_terr <- suppressWarnings(countrycode::countrycode(d$countryallterrorists, "country.name", "cown"))
d$ccode_terr[d$countryallterrorists == "Germany" & d$year <= 1990] <- 260
d$ccode_terr[d$countryallterrorists == "Czech Republic" & d$year <= 1992] <- 315
cat("year range:", range(d$year), " origin names without COW code:", sum(is.na(d$ccode_terr)), "rows\n")

term <- "usmilaid_milexpcomb"
rhs0 <- c("lnpopterr", "lndistance", "lngdppcconstterrorists", "polity2terrorists", term)
fit_nb <- function(rhs, data) {
  f <- as.formula(paste("sumuskilled ~", paste(rhs, collapse = " + ")))
  m <- fenegbin(f, data = data, vcov = ~countryallterrorists, notes = FALSE)
  stopifnot(m$convStatus)
  ct <- coeftable(m); stopifnot(term %in% rownames(ct))
  cat("\n---- nbreg sumuskilled ~", paste(rhs, collapse = " + "), "| N =", nobs(m), "| theta =", round(m$theta, 4), "\n"); print(round(ct, 6))
  out <- list(c(b = ct[term, 1], se = ct[term, 2], p = ct[term, 4])); names(out) <- term
  attr(out, "n") <- nobs(m); out
}
vars <- c("sumuskilled", rhs0, "countryallterrorists")
full <- d[complete.cases(d[, vars]), ]
stopifnot(!anyNA(full$ccode_terr))
r_full <- fit_nb(rhs0, full)
cat("\nREPRO vs article: b", r_full[[term]][["b"]], "(pub .055) se", r_full[[term]][["se"]], "(pub .019) N", attr(r_full, "n"), "(pub 3360)\n")
stopifnot(attr(r_full, "n") == 3360, abs(round(r_full[[term]][["b"]], 3) - 0.055) < 1e-9)
ds <- common_sample(d, vars, "ccode_us", "ccode_terr", "year")
M <- run_six(rhs0, ds, fit_nb)
saveRDS(list(M = M, r_full = r_full, n_used = nrow(ds)), file.path(W, "M_P9_015.rds"))
# no sensitivity hook: negative binomial is outside the hook's estimator families (logit/probit/GEE/OLS)

meta <- list(pid = "P9_015", cite = "Neumayer & Plumper 2011", journal = "JPR", year = 2011, doi = "10.1177/0022343310390147", family = "TIE")
write_p8(meta, term,
  "abstract: more anti-American terrorism from 'countries that receive more US military aid'",
  "Table II model k1 (nbreg, US citizens killed; LSE RO accepted manuscript)",
  "negative binomial (NB2), SE clustered on terrorists' home country",
  "fenegbin (NB2 ML), vcov clustered countryallterrorists, all six steps", 0.055, 0.019, 3360, TRUE, FALSE, M, "N",
  sprintf(paste0("US-victim slice of directed dyads (origin -> US), 1978-2004. Published values from the LSE RO accepted manuscript ",
    "(no author log). Full-sample repro b=%.6f se=%.6f N=%d vs .055 (.019) N=3360. Common 1822-2011 sample N=%d. ",
    "History linked to dyad (2, origin COW); Germany<=1990 -> 260, Czech Republic<=1992 -> 315."),
    r_full[[term]][["b"]], r_full[[term]][["se"]], attr(r_full, "n"), nrow(ds)))
dir.create(file.path(OUT, "results_p9"), recursive = TRUE, showWarnings = FALSE)
src <- file.path(OUT, "results_p8", paste0("P9_015__", term, ".csv")); dst <- file.path(OUT, "results_p9", paste0("P9_015__", term, ".csv"))
stopifnot(file.exists(src), file.rename(src, dst))
sink()
