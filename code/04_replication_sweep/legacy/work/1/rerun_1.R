# Rerun: Leeds 2003 AJPS, Table 1 Model 2 (GEE logit, exchangeable, robust SE clustered on dyad)
# Headline tie: ptargdef (target has defensive ally). b_pub = -.3313126 (authors' log, N=69730)
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work/rerun_helpers_b6.R"))
suppressPackageStartupMessages({library(haven); library(geepack)})
WO <- file.path(OUT, "work/1"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work/1"); sink(file.path(WO, "log_1.txt"), split = TRUE)

# 1. load
d <- read_dta(file.path(W, "LeedsAJPS03.dta")) |> zap_labels() |> as.data.frame()
X <- c("jdem","contig","capprop","atopally","swtglo","ptargdef","pchaloff","pchalneu")
d <- d[complete.cases(d[, c("dispute", X, "dyad", "year")]), ]
d <- d[order(d$dyad, d$year), ]
n_orig <- nrow(d); cat("n_orig =", n_orig, "\n"); stopifnot(n_orig == 69730)

fit <- function(dat, add) {
  f <- reformulate(c(X, add), "dispute")
  m <- geeglm(f, family = binomial, data = dat, id = dyad, corstr = "exchangeable")
  ct <- coef(summary(m)); ct <- cbind(ct[, 1:2], NA, ct[, 4]); m$ct <- ct; m
}

# 2. reproduce
m_rep <- fit(d, character(0)); print(coef(summary(m_rep)))
b_pub <- -.3313126; b_repro <- unname(coef(m_rep)["ptargdef"])
p_repro <- coef(summary(m_rep))["ptargdef", 4]
repro_ok <- abs(b_repro - b_pub) / abs(b_pub) <= .10 && (p_repro < .05) == TRUE
cat("b_repro =", b_repro, " p =", p_repro, " ok =", repro_ok, "\n")

# 3. history: dyad = ccode1*1000 + ccode2
d$ccode1 <- d$dyad %/% 1000; d$ccode2 <- d$dyad %% 1000
stopifnot(all(d$ccode1 >= 2), all(d$ccode2 >= 2))
d <- add_history(d, "ccode1", "ccode2", "year")
s <- d[d$year >= 1822 & d$year <= 2011, ]
s <- s[complete.cases(s[, c("hist_dyad_w15","hist_mono_w15","hist_dyad_w26","hist_mono_w26","bm_cinc_w15","bm_democ_w15")]), ]
s <- s[order(s$dyad, s$year), ]
cat("n_used =", nrow(s), " dropped =", n_orig - nrow(s), "\n")
print(summary(s[, c("hist_dyad_w15","hist_mono_w15","bm_cinc_w15","bm_democ_w15")]))

# 4. model sequence
r <- list()
for (nm in names(ADD_SETS)) {
  m <- fit(s, ADD_SETS[[nm]]); cat("\n====", nm, "\n"); print(coef(summary(m)))
  r[[nm]] <- getc(m$ct, "ptargdef")
}

# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "1", dat = s, y = "dispute", tie = "ptargdef", ctrl = setdiff(X, "ptargdef"), fe = NULL,
           fit = function(dat, add) getc(fit(dat, add)$ct, "ptargdef"))
}

row <- list(cid = "1", short_cite = "Leeds 2003", reproduced = ifelse(repro_ok, "REPRODUCED", "NOT_REPRODUCED"),
  b_pub = b_pub, b_repro = b_repro, estimator = "GEE logit exchangeable, robust (geepack)", tie_var = "ptargdef",
  n_orig = n_orig, n_used = nrow(s), author_has_history_control = "N",
  notes = paste0("Reproduced from authors' log (Table 1 M2). Sample 1816-1944; 1822-2011 restriction + benchmark missingness drop ", n_orig - nrow(s),
                 " obs (", sum(d$year < 1822), " from years 1816-1821, rest missing benchmarks). dyad decoded as ccode1*1000+ccode2. No peace-years in original."))
row <- fill_seq(row, r)
write_result(row)
sink()
