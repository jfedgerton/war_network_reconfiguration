# Rerun 12: Kleinberg, Robinson & French 2012 JOP -- Table 1 col 2 (Any MID, Full Model)
# logit mzmid1 lnmindev lower_conc low_share ... numstate pcy pcy2 pcy3, cluster(dyadid)
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work/rerun_common_b8.R"))
library(haven)
WO <- file.path(OUT, "work/12"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)
sink(file.path(WO, "log_12.txt"), split = TRUE)

# 1. Load authors' data (nondirected dyad-years 1948-2001; numa/numb = COW codes). Keep needed columns only (memory).
v <- c("numa","numb","dyadid","year","mzmid1","lnmindev","lower_conc","low_share","low_dem","high_dem",
       "contig_dummy","allies","lncaprat","lndist","majpow","numstate")
d <- as.data.frame(zap_labels(read_dta(file.path(ROOT, "work/12/KleinbergRobinsonFrench2012JOP.dta"), col_select = all_of(v))))
stopifnot(nrow(d) == 1085642)
d$ord <- seq_len(nrow(d))

# btscs mzmid1 year dyadid, gen(pcy): Stata markouts obs with missing mzmid1; counter = years since the
# last event year within dyad (first obs of dyad = 0). Verified empirically: this variant reproduces
# every Table 1 col 2 coefficient to 3 decimals (obs-count variant gives lower_conc -3.035).
d <- d %>% filter(!is.na(mzmid1)) %>% arrange(dyadid, year, ord) %>% group_by(dyadid) %>%
  mutate(.l = coalesce(lag(mzmid1) == 1, FALSE), .sp = cumsum(.l)) %>% group_by(dyadid, .sp) %>%
  mutate(pcy = year - first(year)) %>% ungroup() %>% select(-.l, -.sp) %>% as.data.frame()
d$pcy2 <- d$pcy^2; d$pcy3 <- d$pcy^3
d <- d[complete.cases(d[, v]), ]

TIE <- "lower_conc"
RHS <- "mzmid1 ~ lnmindev + lower_conc + low_share + low_dem + high_dem + contig_dummy + allies + lncaprat + lndist + majpow + numstate + pcy + pcy2 + pcy3"
fitter <- function(add, dat) fit_tie_fe(as.formula(paste(RHS, add)), dat, binomial(), TIE, cluster = ~dyadid)

# 2. Reproduce (published: lower_conc -3.051 (0.665), N = 427,394)
cat("==== Reproduction ====\n")
r_pub <- fitter("", d)
b_pub <- -3.051
stopifnot(r_pub$n == 427394, abs(r_pub$b - b_pub) / abs(b_pub) <= 0.10, r_pub$p < 0.05)
n_orig <- r_pub$n

# 3. Merge history; restrict to 1822-2011 (data are 1948-2001)
source(file.path(CODE, "hist/add_history.R"))
d <- add_history(d, "numa", "numb", "year")
d <- filter(d, year >= 1822, year <= 2011)
cat("\nNA counts in history/benchmark vars:\n"); print(colSums(is.na(d[, HIST_VARS])))
d <- d[complete.cases(d[, c(v, HIST_VARS)]), ]
cat("Rows after year restriction + complete cases:", nrow(d), " (orig model N", n_orig, ")\n")

# 4. Models on the same sample
R <- run_six(fitter, d)

# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "12", dat = d, y = "mzmid1", tie = TIE,
           ctrl = c("lnmindev", "low_share", "low_dem", "high_dem", "contig_dummy", "allies", "lncaprat", "lndist", "majpow", "numstate",
                    "pcy", "pcy2", "pcy3"), fe = NULL,
           fit = function(dat, add) { r <- fitter(if (length(add)) paste("+", paste(add, collapse = " + ")) else "", dat); c(b = r$b, se = r$se, p = r$p) })
}

# 5. Write row
write_row(12, "Kleinberg, Robinson & French 2012", "REPRODUCED", b_pub, r_pub$b, "logit (fixest feglm), SE clustered on dyadid (Stata-style G/(G-1))",
  TIE, n_orig, R$m0$n, R, ahc = "Y",
  notes = sprintf("Table 1 col 2 (Any MID full); author controls BTSCS pcy cubic (PEACE_YEARS); btscs rebuilt in R (year-diff counter, missing-DV rows excluded) matches all coefs; dropped %d obs vs orig (missing benchmark lags)", n_orig - R$m0$n))
sink()
