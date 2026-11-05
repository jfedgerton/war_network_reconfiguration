# Rerun 122: Chatagnier & Kavakli 2017 JCR -- Table 1 col 1 (All MIDs, main)
# logit mid_init mid_py* contig_dummy lndist l.pwin_cap l.jointdem l.s3un majmaj majmin l.rivalry both_europe dec_1-3
#   l.exports l.lowerdep2015 if ongoing==0, cl(dyad_num)
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work/rerun_common_b8.R"))
library(haven)
WO <- file.path(OUT, "work/122"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
sink(file.path(WO, "log_122.txt"), split = TRUE)

# 1. Load authors' data (nondirected dyad-years; ccode1/ccode2 COW)
v <- c("ccode1","ccode2","year","dyad_num","mid_init","mid_py","mid_py2","mid_py3","contig_dummy","lndist","pwin_cap",
       "jointdem","s3un","majmaj","majmin","rivalry","both_europe","dec_1","dec_2","dec_3","exports","lowerdep2015","ongoing_mid")
d <- as.data.frame(zap_labels(read_dta(file.path(ROOT, "work/122/export_similarity_data.dta"), col_select = all_of(v))))
stopifnot(!any(duplicated(d[, c("dyad_num","year")])))   # tsset dyad_num year must be valid
# Stata l.x under tsset dyad_num year: value at (dyad_num, year-1), missing if that row is absent
lagv <- c("pwin_cap","jointdem","s3un","rivalry","exports","lowerdep2015")
L <- d[, c("dyad_num","year", lagv)]; L$year <- L$year + 1
names(L)[-(1:2)] <- paste0("L_", lagv)
d <- left_join(d, L, by = c("dyad_num","year"))
# `if ongoing==0` abbreviates ongoing_mid (only var starting with "ongoing")
d <- d[!is.na(d$ongoing_mid) & d$ongoing_mid == 0, ]

TIE <- "L_exports"
RHS <- paste("mid_init ~ mid_py + mid_py2 + mid_py3 + contig_dummy + lndist + L_pwin_cap + L_jointdem + L_s3un + majmaj + majmin +",
             "L_rivalry + both_europe + dec_1 + dec_2 + dec_3 + L_exports + L_lowerdep2015")
mv <- all.vars(as.formula(RHS))
d <- d[complete.cases(d[, c(mv, "dyad_num")]), ]
fitter <- function(add, dat) fit_tie_fe(as.formula(paste(RHS, add)), dat, binomial(), TIE, cluster = ~dyad_num)

# 2. Reproduce. b_pub from authors' own saved estimates (beta_maingraph.txt, same covariate list as Table 1 col 1):
#    L.exports = 0.76607997, se = sqrt(vcov) from vcov_maingraph.txt
cat("==== Reproduction ====\n")
r_pub <- fitter("", d)
b_pub <- 0.76607997
stopifnot(abs(r_pub$b - b_pub) / abs(b_pub) <= 0.10, r_pub$p < 0.05)
n_orig <- r_pub$n

# 3. Merge history; restrict to 1822-2011
source(file.path(CODE, "hist/add_history.R"))
d <- add_history(d, "ccode1", "ccode2", "year")
d <- filter(d, year >= 1822, year <= 2011)
cat("\nNA counts in history/benchmark vars:\n"); print(colSums(is.na(d[, HIST_VARS])))
d <- d[complete.cases(d[, c(mv, HIST_VARS)]), ]
cat("Rows after year restriction + complete cases:", nrow(d), " (orig model N", n_orig, ")\n")

# 4. Models on the same sample
R <- run_six(fitter, d)

# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "122", dat = d, y = "mid_init", tie = TIE,
           ctrl = c("mid_py", "mid_py2", "mid_py3", "contig_dummy", "lndist", "L_pwin_cap", "L_jointdem", "L_s3un", "majmaj", "majmin",
                    "L_rivalry", "both_europe", "dec_1", "dec_2", "dec_3", "L_lowerdep2015"), fe = NULL,
           fit = function(dat, add) { r <- fitter(if (length(add)) paste("+", paste(add, collapse = " + ")) else "", dat); c(b = r$b, se = r$se, p = r$p) })
}

# 5. Write row
write_row(122, "Chatagnier & Kavakli 2017", "REPRODUCED", b_pub, r_pub$b, "logit (fixest feglm), SE clustered on dyad_num (Stata-style)",
  "L.exports (export similarity)", n_orig, R$m0$n, R, ahc = "Y",
  notes = sprintf("Table 1 col 1; b_pub from authors' beta_maingraph.txt (identical spec); author controls peace-yr cubic + lagged rivalry (MULTIPLE); ongoing MIDs dropped; dropped %d obs vs orig", n_orig - R$m0$n))
sink()
