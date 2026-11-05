# P9_012 Johnson, Leeds & Savun 2024 ISQ "Alliances and Civil War Intervention" (10.1093/isq/sqae128)
# Data + do-file + log: Harvard Dataverse doi:10.7910/DVN/OG2REM (jls_isq.dta, jls_isq.do, jls_isq.smcl).
# HEADLINE: Table 1 (do-file "Table 1"):
#   mlogit typeint consul defense_pure neunon rival rebelties govties capdist colony caprat premilint length, robust cluster(conflict_id)
#   Outcomes 0 None (base), 1 Government_Support, 2 Rebel_Support. Unit = potential intervener (ccode1) x civil-war state (ccode2)
#   x year x conflict, 1975-2017. Published (author's log): N = 194,422, 171 clusters;
#   Government_Support eq: consul 1.066916 (.3691158), defense_pure 1.679701 (.4509977); Rebel_Support eq: neunon -.7833151 (.2176007).
# KEY TERM (a priori, abstract + P8_017 precedent = one term in one equation): the abstract's lead claim is
#   "some kinds of alliances are associated with pro-government intervention"; defense pacts are the alliance type the abstract
#   ties to both troop and non-troop support -> headline = defense_pure in the Government_Support equation (key_term gov_defense_pure).
#   Non-headline rows: consul (Government_Support eq; gov_consul) and neunon (Rebel_Support eq; reb_neunon).
# Estimator (precedents work_p8/P8_017/rerun_P8_017.R and work/52/rebuild_52.R): multinomial logit by Newton-Raphson (exact MLE),
#   base = lowest category, cluster-robust sandwich on conflict_id with Stata's G/(G-1) factor (recB/mlogit_helpers.R mlogit_fit,
#   the helper used by study 52); same estimator in all six steps. History terms enter every equation additively (as in P8_017/52).
#   History attached to the directed (intervener, civil-war state) pair in the observation year.
# No randomness (seed 123 set in common_b3.R).
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R"))
source(file.path(CODE, "recB/mlogit_helpers.R"))   # mlogit_fit(), mlogit_stat()
WD <- file.path(ROOT, "work_p9/P9_012"); W <- file.path(OUT, "work_p9/P9_012"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_012.txt"), split = TRUE)
rhs0 <- c("consul", "defense_pure", "neunon", "rival", "rebelties", "govties", "capdist", "colony", "caprat", "premilint", "length")
d <- haven::read_dta(file.path(WD, "jls_isq.dta"), col_select = all_of(c("ccode1", "ccode2", "year", "conflict_id", "typeint", rhs0)))
d <- as.data.frame(haven::zap_labels(d))
cat("rows", nrow(d), "years", range(d$year), "conflicts", length(unique(d$conflict_id)), "\n"); print(table(d$typeint))
stopifnot(nrow(d) == 194422, all(d$typeint %in% 0:2), length(unique(d$conflict_id)) == 171, min(d$year) == 1975, max(d$year) == 2017)

TERMS <- c(gov_defense_pure = "defense_pure", gov_consul = "consul", reb_neunon = "neunon")
EQ    <- c(gov_defense_pure = 1, gov_consul = 1, reb_neunon = 2)
fit_mnl <- function(rhs, data) {
  stopifnot(all(0:2 %in% data$typeint))
  m <- mlogit_fit(data$typeint, data[, rhs, drop = FALSE], cluster = data$conflict_id)
  stopifnot(identical(as.numeric(m$levels), c(0, 1, 2)))
  cat("\n---- MNL typeint ~", paste(rhs, collapse = " + "), "| N =", m$n, "| clusters =", length(unique(data$conflict_id)),
      "| logLik =", round(m$loglik, 4), "| NR iters =", m$iter, "\n")
  print(signif(m$coeftable, 7))
  out <- lapply(names(TERMS), function(k) mlogit_stat(m, TERMS[[k]], EQ[[k]])); names(out) <- names(TERMS)
  attr(out, "n") <- m$n; out
}
vars <- c("typeint", rhs0, "conflict_id", "ccode1", "ccode2", "year")
full <- d[complete.cases(d[, vars]), ]
PUB <- list(gov_defense_pure = c(b = 1.679701, se = .4509977), gov_consul = c(b = 1.066916, se = .3691158), reb_neunon = c(b = -.7833151, se = .2176007))
# 1. authors' full sample (published N = 194,422)
r_full <- fit_mnl(rhs0, full); n_orig <- attr(r_full, "n")
for (k in names(PUB)) {
  rd <- abs(r_full[[k]][["b"]] - PUB[[k]][["b"]]) / abs(PUB[[k]][["b"]])
  cat("REPRO vs author log:", k, "b", r_full[[k]][["b"]], "(pub", PUB[[k]][["b"]], ") se", r_full[[k]][["se"]], "(pub", PUB[[k]][["se"]], ") rel.diff b", rd, "\n")
  stopifnot(rd <= 0.01)
}
stopifnot(n_orig == 194422)
# 2. common 1822-2011 sample + history (dyadic: intervener x civil-war state)
ds <- common_sample(d, vars, "ccode1", "ccode2", "year")
M <- run_six(rhs0, ds, fit_mnl)
saveRDS(list(M = M, r_full = r_full, n_orig = n_orig, n_used = nrow(ds)), file.path(W, "M_P9_012.rds"))
# ---- result-CSV write step (write_p8 writes to results_p8; the file is moved to results_p9 unchanged) ----
meta <- list(pid = "P9_012", cite = "Johnson, Leeds & Savun 2024", journal = "ISQ", year = 2024, doi = "10.1093/isq/sqae128", family = "TIE")
ev <- "abstract: 'some kinds of alliances are associated with pro-government intervention'"
art <- "multinomial logit (typeint: none [base] / government support / rebel support), robust SE clustered on conflict_id; intervener x civil-war state x year 1975-2017"
est <- "multinomial logit by Newton-Raphson, cluster-robust sandwich on conflict_id x G/(G-1) (recB mlogit_fit), all six steps"
ahc <- "Y (prior military intervention premilint)"
mv <- function(term) {
  src <- file.path(OUT, "results_p8", paste0("P9_012__", term, ".csv")); dst <- file.path(OUT, "results_p9", paste0("P9_012__", term, ".csv"))
  dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE); stopifnot(file.exists(src), file.rename(src, dst))
}
nt <- function(k) sprintf(paste0("Equation-specific MNL coefficient (%s, %s equation). Dyadic intervener x civil-war state; rows are conflict-specific so (ccode1, ccode2, year) repeats. ",
  "Authors' Dataverse data + do-file; published values from the author's log (Table 1). Full-sample repro b=%.6f se=%.6f N=%d. Common 1822-2011 sample N=%d (drops 2012-17)."),
  TERMS[[k]], if (EQ[[k]] == 1) "Government_Support" else "Rebel_Support", r_full[[k]][["b"]], r_full[[k]][["se"]], n_orig, nrow(ds))
write_p8(meta, "gov_defense_pure", ev, "Table 1 (MNL), Government_Support equation: defense_pure", art, est,
  1.679701, .4509977, 194422, TRUE, FALSE, M, ahc, nt("gov_defense_pure")); mv("gov_defense_pure")
write_p8(meta, "gov_consul", ev, "Table 1 (MNL), Government_Support equation: consul", art, est,
  1.066916, .3691158, 194422, TRUE, FALSE, M, ahc, nt("gov_consul")); mv("gov_consul")
write_p8(meta, "reb_neunon", ev, "Table 1 (MNL), Rebel_Support equation: neunon", art, est,
  -.7833151, .2176007, 194422, TRUE, FALSE, M, ahc, nt("reb_neunon")); mv("reb_neunon")
sink()
