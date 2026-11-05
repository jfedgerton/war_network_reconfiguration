# Rerun: Gartzke & Li 2003 JPR, Table VI column (2): logit of MID onset, Barbieri (1998a) trade data, all dyads 1950-92,
# robust SE clustered on dyad. Headline tie: dependll (trade dependence, log, t-1): b_pub = -0.4158 (z = 3.317), N = 104,270.
# sharell (trade share) is also in the model (+0.4084); they are r = .98 collinear.
# Data: 405Gartzke.Barbieri.dta (from inner zip file44806). Inner file44804 dta is corrupt (year values garbage) -> unused.
# No do-file in the archive: spec taken from the table (DEMLOW DEMHIGH CONTIG DIST MAJPOW ALLIES CAPRATIO peaceyr + 3 splines).
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work/rerun_helpers_b6.R"))
suppressPackageStartupMessages({library(readstata13); library(fixest)})
WO <- file.path(OUT, "work/50"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work/50"); sink(file.path(WO, "log_50.txt"), split = TRUE)

# 1. load
d <- read.dta13(file.path(W, "file44806_gartzke.barbieri/405Gartzke.Barbieri.dta"), convert.factors = FALSE)
d$year <- as.numeric(format(as.Date(d$year), "%Y"))  # year stored as Stata %t date -> calendar year
stopifnot(all(d$year >= 1816 & d$year <= 1992))
names(d) <- sub("^_", "x_", names(d))
X <- c("sharell","dependll","demlow","demhigh","contig","logdstab","majdyad","alliesr","lncaprt",
       "peaceyr","x_spline1","x_spline2","x_spline3")
d <- d[d$year >= 1950 & complete.cases(d[, c("newmid", X, "dyadid")]), ]
n_orig <- nrow(d); cat("n_orig =", n_orig, " (paper N = 104,270)\n"); print(table(d$sample, useNA = "ifany"))

SSC <- ssc(adj = FALSE, cluster.adj = TRUE)
fit <- function(dat, add) feglm(reformulate(c(X, add), "newmid"), family = binomial, data = dat,
                                cluster = ~dyadid, ssc = SSC)

# 2. reproduce
m_rep <- fit(d, character(0)); print(coeftable(m_rep))
b_pub <- -0.4158; b_repro <- unname(coef(m_rep)["dependll"]); p_repro <- pvalue(m_rep)["dependll"]
repro_ok <- abs(b_repro - b_pub) / abs(b_pub) <= .10 && p_repro < .05
cat("b_repro =", b_repro, " z =", tstat(m_rep)["dependll"], " p =", p_repro, " ok =", repro_ok, "\n")

# 3. history (undirected dyads statea < stateb)
d <- add_history(d, "statea", "stateb", "year")
s <- d[d$year >= 1822 & d$year <= 2011, ]
s <- s[complete.cases(s[, c("hist_dyad_w15","hist_mono_w15","hist_dyad_w26","hist_mono_w26","bm_cinc_w15","bm_democ_w15")]), ]
cat("n_used =", nrow(s), " dropped =", n_orig - nrow(s), "\n")
print(summary(s[, c("hist_dyad_w15","hist_mono_w15","bm_cinc_w15","bm_democ_w15")]))

# 4. model sequence
r <- list()
for (nm in names(ADD_SETS)) {
  m <- fit(s, ADD_SETS[[nm]]); cat("\n====", nm, " N =", nobs(m), "\n"); print(coeftable(m))
  stopifnot(nobs(m) == nrow(s))
  r[[nm]] <- getc(coeftable(m), "dependll")
}

# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "50", dat = s, y = "newmid", tie = "dependll", ctrl = setdiff(X, "dependll"), fe = NULL,
           fit = function(dat, add) getc(coeftable(fit(dat, add)), "dependll"))
}

row <- list(cid = "50", short_cite = "Gartzke & Li 2003", reproduced = ifelse(repro_ok, "REPRODUCED", "NOT_REPRODUCED"),
  b_pub = b_pub, b_repro = b_repro, estimator = "logit, cluster dyadid (fixest)", tie_var = "dependll",
  n_orig = n_orig, n_used = nrow(s), author_has_history_control = "Y",
  notes = paste0("Table VI col 2 (Barbieri data, all dyads 1950-92); no do-file, spec from table. sharell (r=.98 with dependll) kept in model. ",
                 "Dropped ", n_orig - nrow(s), " obs (missing benchmarks). Author ch_control=PEACE_YEARS (BKT splines)."))
row <- fill_seq(row, r)
write_result(row)
sink()
