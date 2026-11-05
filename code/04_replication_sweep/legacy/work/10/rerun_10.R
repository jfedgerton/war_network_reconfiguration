# Rerun cid 10: Benson 2011 JOP "Unpacking Alliances". Data: author's replication-data.zip (Vanderbilt site).
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work/rerun_helpers_recA.R"))
suppressPackageStartupMessages(library(haven))
WO <- file.path(OUT, "work/10"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)
sink(file.path(WO, "log_10.txt"), split = TRUE)

# 1. Load. Headline = Table 4 Model 5 (logit, cluster ddyad). Key term (abstract, H2): uncondcompel
#    (b_pub 1.285, SE 0.122). Secondary key (H4): condeter (b_pub -0.184, SE 0.078) reported in notes.
X <- "jointdem + cont + caprat2 + sglo + uncondcompel + condcompel + uncondeter + condeter + probabilistic + peaceyears + peace2 + peace3"
keep <- c("ccode1","ccode2","year","ddyad","mzinit", strsplit(X, " \\+ ")[[1]])
d <- as.data.frame(read_dta(file.path(ROOT, "work/10/replication.dta"), col_select = all_of(keep)))
stopifnot(!anyDuplicated(d[, c("ccode1","ccode2","year")]))
cat("rows:", nrow(d), " years:", range(d$year), "\n")
f0 <- as.formula(paste("mzinit ~", X))

# 2. Reproduce on full sample
r_pub <- fit_tie(f0, d, binomial(), "uncondcompel", cluster = "ddyad")
r_pub_cd <- fit_tie(f0, d, binomial(), "condeter", cluster = "ddyad")
b_pub <- 1.285; n_orig <- r_pub$n
stopifnot(n_orig == 172196)
repro <- if (abs(r_pub$b - b_pub)/abs(b_pub) <= .10 && r_pub$p < .05) "REPRODUCED" else "NOT_REPRODUCED"
cat("repro:", repro, r_pub$b, "\n")

# 3. History (DV mzinit is contemporaneous initiation in year t -> merge on year)
d <- add_history(d, "ccode1", "ccode2", "year")
s <- subset(d, year >= 1822 & year <= 2011)
s <- s[complete.cases(s[, c(all.vars(f0), "hist_dyad_w15","hist_mono_w15","hist_dyad_w26","hist_mono_w26","bm_cinc_w15","bm_democ_w15")]), ]
cat("dropped:", nrow(d) - nrow(s), " (pre-1822:", sum(d$year < 1822), ")\n")

# 4. Models on same sample
m0  <- fit_tie(f0, s, binomial(), "uncondcompel", "ddyad")
mD  <- fit_tie(update(f0, . ~ . + hist_dyad_w15), s, binomial(), "uncondcompel", "ddyad")
mDM <- fit_tie(update(f0, . ~ . + hist_dyad_w15 + hist_mono_w15), s, binomial(), "uncondcompel", "ddyad")
m26 <- fit_tie(update(f0, . ~ . + hist_dyad_w26 + hist_mono_w26), s, binomial(), "uncondcompel", "ddyad")
bC  <- fit_tie(update(f0, . ~ . + bm_cinc_w15), s, binomial(), "uncondcompel", "ddyad")
bD  <- fit_tie(update(f0, . ~ . + bm_democ_w15), s, binomial(), "uncondcompel", "ddyad")
# secondary key (conditional deterrent, target side)
c0  <- fit_tie(f0, s, binomial(), "condeter", "ddyad")
cDM <- fit_tie(update(f0, . ~ . + hist_dyad_w15 + hist_mono_w15), s, binomial(), "condeter", "ddyad")

# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "10", dat = s, y = "mzinit", tie = "uncondcompel", ctrl = setdiff(strsplit(X, " \\+ ")[[1]], "uncondcompel"), fe = NULL,
           fit = function(dat, add) { ff <- if (length(add)) update(f0, as.formula(paste(". ~ . +", paste(add, collapse = " + ")))) else f0
             r <- fit_tie(ff, dat, binomial(), "uncondcompel", "ddyad"); c(b = r$b, se = r$se, p = r$p) })
}

notes <- sprintf(paste0("FOUND at author site (cdn.vanderbilt.edu .../2012/08/replication-data.zip). Table 4 Model 5 logit, cluster ddyad (Stata robust). ",
  "Key = uncondcompel (abstract: unconditional compellent +249%%; H2). b_repro %.4f vs pub 1.285 (SE %.3f vs .122). Peace-years cubic (peaceyears, /100, /1000) kept. ",
  "History merged on year (DV contemporaneous). Dropped %d obs (years<1822 or missing benchmarks). ",
  "Secondary key condeter (H4, pub -0.184): repro %.4f; b0 %.4f p0 %.4f -> b_dm %.4f p_dm %.4f (pct %.1f)."),
  r_pub$b, r_pub$se, nrow(d) - nrow(s), r_pub_cd$b, c0$b, c0$p, cDM$b, cDM$p, pct(cDM$b, c0$b))
write_row(10, "Benson 2011", repro, b_pub, r_pub$b, "logit, dyad-clustered SE (Stata robust)", "uncondcompel",
          n_orig, m0$n, m0, mD, mDM, m26, bC, bD, ahc = "Y", notes = notes)
sink()
