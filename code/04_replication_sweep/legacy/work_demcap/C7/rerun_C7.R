# C7: Dafoe & Kelsey 2014 JPR, "Observing the capitalist peace" (qualitative case-selection paper; no coefficient table).
# Key term: capopenl (capital openness, low). Abstract: "Countries with open capital markets tend to have fewer militarized disputes and wars."
# Headline quantitative model = the case-selection model that drives Rule 3 (pr5) in the authors' do-file:
#   logit maoznewl demlo demhi deplo capopenl rgdppclo gdpcontg sun2cati contig logdstab majpdyds allies lncaprt
#         namerica samerica europe africa nafmeast asia, cluster(dyadid)
#   = Gartzke 2007 Model 5 w/o temporal controls + region dummies. b_pub from authors' log: capopenl -.2005327 (p .002), N = 166,152.
#   Stata "allies" abbreviates alliesr in this file.
#   (Figure 1 uses the same spec with a cubic in capopenl on drop==0; not used here because the linear term of a cubic is not interpretable.)
# Data: Gartzke's capitalistpeace_012007 (Dataverse original .dta). Authors merged COW MID + fishing-dispute files first;
#   neither enters this model (their log notes ~12 extra obs from the merge).
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_demcap/common_c3467.R"))
WD <- file.path(ROOT, "work_demcap/C7"); W <- file.path(OUT, "work_demcap/C7"); dir.create(W, recursive = TRUE, showWarnings = FALSE); sink(file.path(W, "log_C7.txt"), split = TRUE)
slim <- file.path(WD, "slim_C7.rds")  # port: staged slim copy of cp012007.dta (read-only input)
V <- c("statea", "stateb", "year", "dyadid", "maoznewl", "demlo", "demhi", "deplo", "capopenl", "rgdppclo", "gdpcontg", "sun2cati", "contig",
       "logdstab", "majpdyds", "alliesr", "lncaprt", "namerica", "samerica", "europe", "africa", "nafmeast", "asia")
if (!file.exists(slim)) {
  d <- haven::read_dta(file.path(WD, "cp012007.dta"), col_select = tidyselect::all_of(V))
  d <- as.data.frame(haven::zap_labels(d)); saveRDS(d, file.path(W, "slim_C7.rds"))
} else d <- readRDS(slim)
stopifnot(all(d$maoznewl %in% c(0, 1, NA)), !anyDuplicated(d[, c("statea", "stateb", "year")]))
f0 <- maoznewl ~ demlo + demhi + deplo + capopenl + rgdppclo + gdpcontg + sun2cati + contig + logdstab + majpdyds + alliesr + lncaprt +
  namerica + samerica + europe + africa + nafmeast + asia
K <- "capopenl"; bp <- structure(-0.2005327, sig = TRUE); n_orig <- 166152
r <- fit_keys(f0, d, K, "dyadid", file.path(W, "cache_repro.rds"))
cat("N repro:", r$n, " vs authors' log:", n_orig, "\n")
ds <- common_sample(d, f0, "dyadid", "statea", "stateb", "year"); rm(d); gc()
M <- run_six(f0, ds, K, "dyadid", W)
# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "C7", dat = ds, y = "maoznewl", tie = "capopenl",
           ctrl = c("demlo", "demhi", "deplo", "rgdppclo", "gdpcontg", "sun2cati", "contig", "logdstab", "majpdyds", "alliesr", "lncaprt", "namerica", "samerica", "europe", "africa", "nafmeast", "asia"), fe = NULL,
           fit = function(dat, add) { ff <- if (length(add)) update(f0, as.formula(paste(". ~ . +", paste(add, collapse = " + ")))) else f0
             fit_keys(ff, dat, "capopenl", "dyadid")$k[["capopenl"]] })
}

write_dc("C7", "Dafoe & Kelsey 2014 (JPR)", "CAPITALISM", K,
  "abstract: 'Countries with open capital markets tend to have fewer militarized disputes and wars.'",
  repro_status(bp, r$k[[K]]), bp, r$k[[K]][["b"]], "logit, dyad-clustered SE (MID onset)", n_orig, M,
  paste0("Case-selection (Rule 3, pr5) MID model: Gartzke 2007 M5 w/o temporal controls + region dummies; NO peace-years/splines in authors' model. ",
         "b_pub from authors' log (paper prints no coefficients). Repro N=", r$n, " vs log ", n_orig, " (authors' COW-MID merge adds ~12 obs). ",
         n_orig - M$m0$n, " obs dropped vs orig in common sample."))
sink()
