# Rerun 64: Goenner 2010 JPR -- Table III Model 2, logit, dyad-clustered SEs, IVs lagged one year (Stata L. on tsset dyadid year):
#   logit mzonset L.(jntdem2 allies contigkb distance capratio majdyd depend_lev depend_asy lowgrowth) mzpceyrs mzpyrs1 mzpyrs2 mzpyrs3
#         L.(open1_lev..open6_lev open0_lev open1_asy..open6_asy open0_asy), cluster(dyadid)
# KEY-AUDIT COPY: key term = L.open1_lev (energy level). Original headline tie was L.depend_lev ("Interconnectedness level": bilateral trade/GDP, the dyadic trade tie present in both models).
# Disaggregated sector shares (open*_lev/asy) are kept as controls; chemicals level (open3_lev) reported in notes.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work/rerun_common_k6.R"))
WO <- file.path(OUT, "work/64"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work/64")
sink(file.path(WO, "log_key_64.txt"), split = TRUE)
# port: JPRdataToys.dta was unzipped from web-appendix.zip and later deleted; idempotent unzip into OUT
if (!file.exists(file.path(WO, "web-appendix/JPRdataToys.dta")))
  unzip(file.path(W, "web-appendix.zip"), files = "web-appendix/JPRdataToys.dta", exdir = WO)
stopifnot(file.exists(file.path(WO, "web-appendix/JPRdataToys.dta")))

# 1. Authors' analysis data (PRIO JPR 2010/47/5 zip -> web-appendix.zip -> JPRdataToys.dta, 118MB)
LAGV <- c("jntdem2","allies","contigkb","distance","capratio","majdyd","depend_lev","depend_asy","lowgrowth",
          paste0("open", c(1:6, 0), "_lev"), paste0("open", c(1:6, 0), "_asy"))
d <- as.data.frame(zap_labels(read_dta(file.path(WO, "web-appendix/JPRdataToys.dta"),
       col_select = all_of(c("ccode1","ccode2","year","dyadid","mzonset","mzpceyrs","mzpyrs1","mzpyrs2","mzpyrs3", LAGV)))))
stopifnot(nrow(d) == 482240, !any(duplicated(d[, c("dyadid","year")])))

# Build L.x exactly as Stata time-series lag: value of x for the same dyadid at year-1 (NA if that year absent)
lagdf <- d[, c("dyadid","year", LAGV)]; lagdf$year <- lagdf$year + 1
names(lagdf)[-(1:2)] <- paste0("L_", LAGV)
d <- left_join(d, lagdf, by = c("dyadid","year"))
stopifnot(nrow(d) == 482240)
d <- d[, c("ccode1","ccode2","year","dyadid","mzonset","mzpceyrs","mzpyrs1","mzpyrs2","mzpyrs3", paste0("L_", LAGV))]
rm(lagdf); invisible(gc())

LV <- paste0("L_", LAGV)
BASE <- paste("mzonset ~", paste(c(LV[1:9], "mzpceyrs","mzpyrs1","mzpyrs2","mzpyrs3", LV[10:23]), collapse = " + "))
MV <- c("mzonset", LV, "mzpceyrs","mzpyrs1","mzpyrs2","mzpyrs3","dyadid")
# Memory-lean fitter (shared 7GB box; glm+vcovCL on 298k x 28 was OOM-killed): fixest::feglm logit, cluster dyadid,
# ssc(adj=FALSE, cluster.adj=TRUE) = Stata logit cluster() (HC0 sandwich * G/(G-1)).
suppressPackageStartupMessages(library(fixest))
fit_fe <- function(f, dat, tie, quiet = FALSE) {
  m <- feglm(f, data = dat, family = binomial("logit"), cluster = ~dyadid, ssc = ssc(adj = FALSE, cluster.adj = TRUE), notes = FALSE)
  stopifnot(m$convStatus)
  ct <- coeftable(m); stopifnot(tie %in% rownames(ct))
  if (!quiet) { cat("\n---- ", paste(deparse(f, width.cutoff = 500), collapse = ""), " | N =", nobs(m), "\n"); print(signif(ct, 5)) }
  out <- list(b = unname(ct[tie, 1]), se = unname(ct[tie, 2]), p = unname(ct[tie, 4]), n = nobs(m)); rm(m); invisible(gc()); out
}
KEY <- "L_open1_lev"   # key-audit: energy-sector interdependence level (H2/H3; abstract lead commodity)
fitter <- function(add, dat) fit_fe(as.formula(paste0(BASE, add)), dat, tie = KEY)

# 2. Reproduce on the full published sample. Paper Table III Model 2 (PDF text lost minus signs; sign read from text:
#    interconnectedness level "significantly reduce[s] conflict"): Interconnectedness level = -19.4148 (robust se 7.3394, p<.05), N = 298,170
cat("\n==== Reproduction: Table III Model 2 (full sample) ====\n")
d_full <- d[complete.cases(d[, MV]), ]
rp <- fitter("", d_full); rm(d_full); invisible(gc())
b_pub <- 14.5970; cat("\nKEY repro: b_pub", b_pub, "b_repro", rp$b, "rel", abs(rp$b-b_pub)/abs(b_pub), "N", rp$n, "\n")
stopifnot(rp$n == 298170)
n_orig <- rp$n

# 3. Merge history (ccode1<ccode2 undirected; 1962-2000 inside GML range)
d <- d[complete.cases(d[, MV]), ]   # pre-trim to model-complete rows (same rows restrict_sample would keep) to save memory
d <- add_history(d, "ccode1", "ccode2", "year")
d <- restrict_sample(d, MV)

# 4. Six models on one sample
R <- run_six(fitter, d)
# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "64", dat = d, y = "mzonset", tie = KEY,
           ctrl = setdiff(c(LV[1:9], "mzpceyrs", "mzpyrs1", "mzpyrs2", "mzpyrs3", LV[10:23]), KEY), fe = NULL,
           fit = function(dat, add) { r <- fitter(if (length(add)) paste0(" + ", paste(add, collapse = " + ")) else "", dat); c(b = r$b, se = r$se, p = r$p) })
}

out <- do.call(rbind, lapply(names(R), function(k) data.frame(model = k, b = R[[k]]$b, se = R[[k]]$se, p = R[[k]]$p, n = R[[k]]$n)))
out$b_repro <- rp$b; out$n_orig <- n_orig
print(out); write.csv(out, file.path(WO, "key_64_six.csv"), row.names = FALSE)
# port: one-row wide copy of the same six-model numbers (b_hist = mDM) so run_one.R can extract by column; no new estimation
wide <- data.frame(id = "64", key_term = "L.open1_lev", b0 = R$m0$b, p0 = R$m0$p, b_dm = R$mDM$b, p_dm = R$mDM$p, b_26 = R$m26$b, p_26 = R$m26$p,
                   pct_bC = pct(R$bC$b, R$m0$b), pct_bD = pct(R$bD$b, R$m0$b), n_used = R$m0$n)
write.csv(wide, file.path(WO, "key_64_wide.csv"), row.names = FALSE)
sink()
