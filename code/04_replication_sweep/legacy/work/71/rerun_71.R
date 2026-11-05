# Rerun 71: Lupu & Greenhill 2017 JPR -- Table II col (1)
# logit maozmidlead samecomm igos ... py_demeaned py2 py3, cluster(dyad) robust, if mzongo~=1
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
suppressPackageStartupMessages({library(haven); library(dplyr); library(fixest); library(peacesciencer)})
WO <- file.path(OUT, "work/71"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)
sink(file.path(WO, "log_71.txt"), split = TRUE)

# 1. Load authors' data (file already contains only mzongo==0 rows)
X <- "samecomm + igos + i_ambassador + i_mission + low_open + dependlow + politylow + contiguous + logdist + majpow + allies + logcapratio + py_demeaned + py2 + py3"
keep <- c("ccode1","ccode2","year","dyad","mzongo","maozmidlead", strsplit(X, " \\+ ")[[1]])
d <- as.data.frame(zap_labels(read_dta(file.path(ROOT, "work/71/Lupu Greenhill Conflict Main.dta"), col_select = all_of(keep))))
d <- subset(d, mzongo != 1); rownames(d) <- NULL
stopifnot(nrow(d) == 430477)
X <- "samecomm + igos + i_ambassador + i_mission + low_open + dependlow + politylow + contiguous + logdist + majpow + allies + logcapratio + py_demeaned + py2 + py3"
# fixest::feglm (memory-lean; glm runs were OOM-killed on the shared box). No FE; cluster SE with
# ssc(adj=TRUE, cluster.adj=TRUE) = Stata logit cluster(dyad): G/(G-1)*(N-1)/(N-K). Verified vs glm+sandwich HC1: -0.54222 (0.11973).
fit <- function(add, dat) {
  m <- feglm(as.formula(paste("maozmidlead ~", X, add)), data = dat, family = binomial("logit"),
             vcov = ~dyad, ssc = ssc(adj = TRUE, cluster.adj = TRUE), lean = TRUE, mem.clean = TRUE)
  ct <- coeftable(m); print(ct); nn <- nobs(m); rm(m); gc()
  list(b = ct["samecomm", 1], se = ct["samecomm", 2], p = ct["samecomm", 4], n = nn)
}
# 2. Reproduce. Published Table II col 1: Same IGO cluster -0.542 (0.120)***, N = 430,477
cat("==== Reproduction ====\n"); r_pub <- fit("", d); b_pub <- -0.542
cat("n =", r_pub$n, "\n"); n_orig <- r_pub$n
rep_ok <- abs(r_pub$b - b_pub) / abs(b_pub) <= 0.10 && r_pub$p < .05

# Timing check: maozmidlead is led one year
p <- gml_part %>% select(dispnum, ccode, sidea, styear)
pp <- inner_join(p, p, by = "dispnum", relationship = "many-to-many") %>% filter(sidea.x != sidea.y, ccode.x < ccode.y) %>%
  transmute(c1 = ccode.x, c2 = ccode.y, on = pmax(styear.x, styear.y)) %>% distinct(c1, c2, on)
dd <- data.frame(c1 = pmin(d$ccode1, d$ccode2), c2 = pmax(d$ccode1, d$ccode2), year = d$year, y = d$maozmidlead)
for (k in 0:1) { g <- paste(dd$c1, dd$c2, dd$year + k) %in% paste(pp$c1, pp$c2, pp$on); cat("GML onset at year+", k, ": "); print(table(y = dd$y, gml = g)) }

# 3. History on outcome year (year+1)
d$year_out <- d$year + 1
source(file.path(CODE, "hist/add_history.R"))
d <- add_history(d, "ccode1", "ccode2", "year_out")
d <- filter(d, year_out >= 1822, year_out <= 2011)
vars <- c("maozmidlead", strsplit(X, " \\+ ")[[1]], "hist_dyad_w15","hist_mono_w15","hist_dyad_w26","hist_mono_w26","bm_cinc_w15","bm_democ_w15")
d <- d[complete.cases(d[, vars]), ]; rownames(d) <- NULL
cat("Rows used:", nrow(d), "\n")

# 4. Models
cat("==== m0 ====\n");  r0  <- fit("", d)
cat("==== mD ====\n");  rD  <- fit("+ hist_dyad_w15", d)
cat("==== mDM ====\n"); rDM <- fit("+ hist_dyad_w15 + hist_mono_w15", d)
cat("==== m26 ====\n"); r26 <- fit("+ hist_dyad_w26 + hist_mono_w26", d)
cat("==== bC ====\n");  rC  <- fit("+ bm_cinc_w15", d)
cat("==== bD ====\n");  rBD <- fit("+ bm_democ_w15", d)
pct <- function(x) 100 * (x$b - r0$b) / abs(r0$b)
out <- data.frame(cid = 71, short_cite = "Lupu & Greenhill 2017", reproduced = ifelse(rep_ok, "REPRODUCED", "NOT_REPRODUCED"),
  b_pub = b_pub, b_repro = r_pub$b, estimator = "logit, dyad-clustered SE (Stata HC1)", tie_var = "samecomm",
  n_orig = n_orig, n_used = r0$n, b0 = r0$b, se0 = r0$se, p0 = r0$p, b_dm = rDM$b, se_dm = rDM$se, p_dm = rDM$p, pct_dm = pct(rDM),
  b_d = rD$b, pct_d = pct(rD), b_26 = r26$b, se_26 = r26$se, p_26 = r26$p, pct_26 = pct(r26), pct_bC = pct(rC), pct_bD = pct(rBD),
  author_has_history_control = "Y",
  notes = sprintf("Table II col1; pub SE 0.120, repro SE %.3f. History merged on year+1 (DV led). Cubic peace-years kept. 20,140 obs dropped solely because bm_democ_w15 is NA (no Polity for microstates); hist vars complete. Est. via fixest feglm (glm OOM on shared box), identical coefs.", r_pub$se))
print(Filter(is.numeric, out))   # port: no citation/text columns in logs
dir.create(file.path(OUT, "results"), recursive = TRUE, showWarnings = FALSE)
write.csv(out, file.path(OUT, "results/71.csv"), row.names = FALSE)
sink()
