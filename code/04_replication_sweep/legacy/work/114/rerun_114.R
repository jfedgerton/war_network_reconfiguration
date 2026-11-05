# Rerun 114: Beardsley, Liu, Mucha, Siegel & Tellez 2020 JOP -- Table 1 col (1) (mod3: spatial-temporal low resolution)
# logit mzmidplusone ~ samecomm*hierlow_close + controls + peaceyears + 3 splines, excl. ongoing MIDs, dyad-clustered SE (no df correction)
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
suppressPackageStartupMessages({library(dplyr); library(fixest); library(peacesciencer)})
WO <- file.path(OUT, "work/114"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)
sink(file.path(WO, "log_114.txt"), split = TRUE)

# 1. Load authors' data (Table 1 col 1 = mod3, per stargazer ordering mod3, mod6, mod2, mod1, mod5, mod4)
d <- readRDS(file.path(ROOT, "work/114/mod3_noongo.rds"))  # built by prep_114.R: mod3 from arms-models.rda, mzongo != 1, model columns
d$pairid <- paste0(d$ccode1, "-", d$ccode2)
stopifnot(nrow(d) == 414158)
X <- "samecomm*hierlow_close + transfer + lngdphigh + lngdplow + demochigh + democlow + igos + contiguous + logdist + majpow + allies + logcapratio + peaceyears + spline1 + spline2 + spline3"
# fixest::feglm (memory-lean; glm run OOM-killed on shared box). Authors: multiwayvcov::cluster.vcov(df_correction=FALSE)
# = no small-sample adjustment -> ssc(adj=FALSE, cluster.adj=FALSE).
fit <- function(add, dat) {
  m <- feglm(as.formula(paste("mzmidplusone ~", X, add)), data = dat, family = binomial("logit"),
             vcov = ~pairid, ssc = ssc(adj = FALSE, cluster.adj = FALSE), lean = TRUE, mem.clean = TRUE)
  ct <- coeftable(m); print(ct); nn <- nobs(m); rm(m); gc()
  list(b = ct["samecomm", 1], se = ct["samecomm", 2], p = ct["samecomm", 4], n = nn,
       bi = ct["samecomm:hierlow_close", 1], sei = ct["samecomm:hierlow_close", 2], pi = ct["samecomm:hierlow_close", 4])
}
# 2. Reproduce. Published Table 1 col 1: samecomm 0.237 (0.246); samecomm x hier -0.004 (0.001)***; N = 390,914
cat("==== Reproduction ====\n")
rownames(d) <- NULL
r_pub <- fit("", d)
b_pub <- 0.237; bi_pub <- -0.004
cat("n =", r_pub$n, "\n")
n_orig <- r_pub$n

# Timing check: mzmidplusone is MID led one year; compare with GML dyadic onsets in year vs year+1
p <- gml_part %>% select(dispnum, ccode, sidea, styear)
pp <- inner_join(p, p, by = "dispnum", relationship = "many-to-many") %>% filter(sidea.x != sidea.y, ccode.x < ccode.y) %>%
  transmute(c1 = ccode.x, c2 = ccode.y, on = pmax(styear.x, styear.y)) %>% distinct(c1, c2, on)
dd <- data.frame(c1 = pmin(d$ccode1, d$ccode2), c2 = pmax(d$ccode1, d$ccode2), year = d$year, y = d$mzmidplusone)
for (k in 0:1) { g <- paste(dd$c1, dd$c2, dd$year + k) %in% paste(pp$c1, pp$c2, pp$on); cat("GML onset at year+", k, ": "); print(table(y = dd$y, gml = g)) }

rm(p, pp, dd); gc()
# 3. History merged on outcome year (year+1) so window t-1..t-5 precedes the led MID
d$year_out <- d$year + 1
source(file.path(CODE, "hist/add_history.R"))
d <- add_history(d, "ccode1", "ccode2", "year_out")
d <- filter(d, year_out >= 1822, year_out <= 2011)
cat("NA counts:\n"); print(colSums(is.na(d[, c("hist_dyad_w15","hist_mono_w15","bm_cinc_w15","bm_democ_w15")])))
# 4b. Sensitivity (not the protocol sample): history vars are never missing, so also run m0/mDM on the full published sample
cat("==== FULL-SAMPLE m0 ====\n");  f0  <- fit("", d)
cat("==== FULL-SAMPLE mDM ====\n"); fDM <- fit("+ hist_dyad_w15 + hist_mono_w15", d)
vars <- c("mzmidplusone","samecomm","hierlow_close","transfer","lngdphigh","lngdplow","demochigh","democlow","igos","contiguous","logdist","majpow","allies","logcapratio","peaceyears","spline1","spline2","spline3",
          "hist_dyad_w15","hist_mono_w15","hist_dyad_w26","hist_mono_w26","bm_cinc_w15","bm_democ_w15")
d <- d[complete.cases(d[, vars]), ]; rownames(d) <- NULL
cat("Rows used:", nrow(d), "\n")

# 4. Models
cat("==== m0 ====\n");  r0  <- fit("", d)
cat("==== mD ====\n");  rD  <- fit("+ hist_dyad_w15", d)
cat("==== mDM ====\n"); rDM <- fit("+ hist_dyad_w15 + hist_mono_w15", d)
cat("==== m26 ====\n"); r26 <- fit("+ hist_dyad_w26 + hist_mono_w26", d)
cat("==== bC ====\n");  rC  <- fit("+ bm_cinc_w15", d)
cat("==== bD ====\n");  rBD <- fit("+ bm_democ_w15", d)
# ---- sensitivity hook (split only): runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  # headline = the interaction samecomm:hierlow_close (bi/sei/pi of the script's fit(), as written to key_114_extract.csv)
  sens_split(id = "114", dat = d,
             fit = function(dat, add) { r <- fit(if (length(add)) paste("+", paste(add, collapse = " + ")) else "", dat); c(b = r$bi, se = r$sei, p = r$pi) })
}
pct <- function(x) 100 * (x$b - r0$b) / abs(r0$b)
pcti <- function(x) 100 * (x$bi - r0$bi) / abs(r0$bi)
rep_ok <- abs(r_pub$b - b_pub) / abs(b_pub) <= 0.10 && (r_pub$p < .05) == FALSE
out <- data.frame(cid = 114, short_cite = "Beardsley et al. 2020", reproduced = ifelse(rep_ok, "REPRODUCED", "NOT_REPRODUCED"),
  b_pub = b_pub, b_repro = r_pub$b, estimator = "logit, dyad-clustered SE (HC0, no df corr.)", tie_var = "samecomm (constituent of samecomm x hierlow_close)",
  n_orig = n_orig, n_used = r0$n, b0 = r0$b, se0 = r0$se, p0 = r0$p, b_dm = rDM$b, se_dm = rDM$se, p_dm = rDM$p, pct_dm = pct(rDM),
  b_d = rD$b, pct_d = pct(rD), b_26 = r26$b, se_26 = r26$se, p_26 = r26$p, pct_26 = pct(r26), pct_bC = pct(rC), pct_bD = pct(rBD),
  author_has_history_control = "Y",
  notes = sprintf("Headline is interaction (H2); per protocol tie=constituent samecomm (n.s. in paper). Interaction: pub -0.004; repro %.5f (p=%.2g); b0 %.5f (p=%.2g) -> b_dm %.5f (p=%.2g), pct_dm %.1f, pct_26 %.1f, pct_bC %.1f. History merged on year+1 (DV led). Table 1 col1 = mod3 (low-res). 19,796 obs dropped solely for NA benchmarks (bm_democ/bm_cinc), which moves samecomm 0.237->0.098. Full-sample sensitivity (N=%d): samecomm %.3f->%.3f; interaction %.5f (p=%.2g) -> %.5f (p=%.2g).",
                  r_pub$bi, r_pub$pi, r0$bi, r0$pi, rDM$bi, rDM$pi, pcti(rDM), pcti(r26), pcti(rC), f0$n, f0$b, fDM$b, f0$bi, f0$pi, fDM$bi, fDM$pi))
print(Filter(is.numeric, out))   # port: no citation/text columns in logs
dir.create(file.path(OUT, "results"), recursive = TRUE, showWarnings = FALSE)
write.csv(out, file.path(OUT, "results/114.csv"), row.names = FALSE)
# port (extraction only, no estimation): the master key term is the interaction samecomm:hierlow_close. Its six-model numbers are
# already estimated above (bi/sei/pi from the same fits; the sweep took them from the identical rerun_key_114.R via write_key_batch.py).
ex <- data.frame(id = "114", key_term = "samecomm:hierlow_close", b0 = r0$bi, se0 = r0$sei, p0 = r0$pi, b_dm = rDM$bi, p_dm = rDM$pi,
  b_26 = r26$bi, p_26 = r26$pi, pct_bC = pcti(rC), pct_bD = pcti(rBD), n_used = r0$n)
stopifnot(all(is.finite(unlist(ex[, -(1:2)]))))
write.csv(ex, file.path(WO, "key_114_extract.csv"), row.names = FALSE)
sink()
