# Rerun 43: Lupu & Traag 2013 JCR, ERRATUM (2018) Table 2 Model 1:
#   logit mzmidplusone stc_low dependlow maxflow lngdphigh lngdplow demochigh democlow igos contiguous logdist majpow allies
#         logcapratio _spline1 _spline2 _spline3 peaceyears, robust, if mzongo~=1          (published: stc_low 0.117 (0.093), N 394,527)
# Pipeline (authors' Zenodo package 10.5281/zenodo.1343311):
#   (a) lt/src/na_maxflow.py  = authors' network_analysis.py maxflow block  -> lt/results/maxflow.csv
#   (b) lt/src/na_comm.py     = authors' Louvain block (100 runs x 3 resolutions) -> lt/results/dyads_{0.6,1.1,1.7}.csv
#       Deviations forced by current libraries (louvain 0.8.2 / igraph 0.11.9 / pandas 3): pandas operator-precedence fix
#       (non-null & >0, the evident intent), and louvain>=0.7 has no module-level set_rng_seed -> one Optimiser seeded 0,
#       reused for all runs. Partitions therefore cannot be bit-identical to the authors' (louvain 0.6.x) runs.
#   (c) this script = R translation of merge_recode_regress_stata.txt (Stata max()/min() ignore missing; x/0 = missing).
#   btscs: btscs.ado unavailable. peaceyears = years since last mzmidplusone==1 within dyad, 0 at first obs and in the year
#   after an event, counting only rows with non-missing mzmidplusone (as verified by cid 65's rerun). nspline(3) basis
#   s_i = (py-k_i)^3_+ - (T-k_i)/T * py^3, knots = quartiles (type 2) of peaceyears, T = max(peaceyears): the knot rule was
#   checked on btscs output in cid 100's deposited data (knot/T ratios ~ quartiles). Not tuned to the erratum table.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work/rerun_common_g7.R"))
suppressPackageStartupMessages({library(haven); library(data.table)})
WO <- file.path(OUT, "work/43"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work/43"); LT <- file.path(W, "lt")  # port: lt/results/*.csv come from the authors' Python network steps (not ported); read as staged inputs
sink(file.path(WO, "log_43.txt"), split = TRUE)

smax <- function(a, b) { r <- pmax(a, b, na.rm = TRUE); r }   # Stata max(): ignores missing, NA only if both NA
smin <- function(a, b) { r <- pmin(a, b, na.rm = TRUE); r }

# ---- maxflow: mean of both directed maxflows per dyad-year ----
mf <- fread(file.path(LT, "results/maxflow.csv"))
setnames(mf, c("source_ccode","target_ccode"), c("ccode1","ccode2"))
mf[, `:=`(a = pmin(ccode1, ccode2), b = pmax(ccode1, ccode2))]
mf[, maxflow := mean(maxflow), by = .(year, a, b)]
mf <- mf[, .(ccode1, ccode2, year, maxflow)]
stopifnot(!anyDuplicated(mf[, .(ccode1, ccode2, year)]))

# ---- same trading community files ----
rd_stc <- function(res, nm) {
  f <- file.path(LT, sprintf("results/dyads_%s.csv", res))
  if (!file.exists(f)) { cat("MISSING", f, "\n"); return(NULL) }
  x <- fread(f); setnames(x, c("ccode_x","ccode_y","same_comm"), c("ccode1","ccode2", nm))
  x[[nm]] <- as.numeric(as.character(x[[nm]]) %in% c("True","TRUE"))   # insheet True/False -> 1/0
  x
}
stc_low <- rd_stc("0.6", "stc_low"); stc_med <- rd_stc("1.1", "stc_med"); stc_high <- rd_stc("1.7", "stc_high")
stopifnot(!is.null(stc_low))

# ---- merge: igos defines the retained rows (Stata keeps all then 'drop if igos==.'), all files 1:1 on ccode1 ccode2 year ----
K <- c("ccode1","ccode2","year")
igos <- as.data.table(zap_labels(read_dta(file.path(LT, "data/igos.dta"))))
dep  <- as.data.table(zap_labels(read_dta(file.path(LT, "data/dependence.dta"))))[, .(ccode1, ccode2, year, depabgled, depbagled, gdpa, gdpb)]
eu   <- as.data.table(zap_labels(read_dta(file.path(LT, "data/eugene2.dta"))))
d <- igos[!is.na(igos)]
for (x in list(stc_low, stc_med, stc_high, dep, mf, eu)) if (!is.null(x)) d <- merge(d, x, by = K, all.x = TRUE)
stopifnot(!anyDuplicated(d[, ..K]))
cat("merged rows (igos non-missing):", nrow(d), "\n")

# ---- recodes (verbatim logic) ----
d[, contiguous := fifelse(contig < 6, 1, fifelse(contig == 6, 0, NA_real_))]
d[, logdist := fifelse(distance == 0, 0, log(distance))]
d[, allies := fifelse(alliance < 4, 1, fifelse(alliance == 4, 0, NA_real_))]
d[, majpow := fifelse(majpow1 == 1, 1, fifelse(majpow1 == 0, 0, NA_real_))]; d[majpow2 == 1, majpow := 1]
d[, highcap := smax(cap_1, cap_2)][, lowcap := smin(cap_1, cap_2)]
d[, capratio := fifelse(lowcap == 0, NA_real_, highcap / lowcap)][, logcapratio := log(capratio)]
for (v in c("democ1","democ2")) { d[get(v) %in% c(-66, -88), (v) := NA]; d[get(v) == -77, (v) := 0] }
d[, demochigh := smax(democ1, democ2)][, democlow := smin(democ1, democ2)]
d[, dependlow := smin(depabgled, depbagled)]
d[, lngdphigh := log(smax(gdpa, gdpb))][, lngdplow := log(smin(gdpa, gdpb))]
d[is.infinite(lngdplow), lngdplow := NA][is.infinite(lngdphigh), lngdphigh := NA][is.infinite(logcapratio), logcapratio := NA]

# ---- lead DV (tsset id year; F.mzmid needs the row at year+1) ----
setkeyv(d, K)
nx <- d[, .(ccode1, ccode2, year = year - 1, mzmidplusone = mzmid)]
d <- merge(d, nx, by = K, all.x = TRUE)

# ---- btscs mzmidplusone year id, g(peaceyears) nspline(3) ----
setorderv(d, K)
ok <- !is.na(d$mzmidplusone); first <- c(TRUE, (d$ccode1[-1] != d$ccode1[-nrow(d)]) | (d$ccode2[-1] != d$ccode2[-nrow(d)]))
py <- rep(NA_real_, nrow(d)); cur <- -1; lastev <- FALSE; ev <- d$mzmidplusone
for (i in seq_len(nrow(d))) { if (first[i]) { cur <- -1; lastev <- FALSE }
  if (ok[i]) { cur <- if (cur < 0 || lastev) 0 else cur + 1; py[i] <- cur; lastev <- ev[i] == 1 } }
d$peaceyears <- py
kn <- quantile(d$peaceyears, c(.25, .5, .75), type = 2, na.rm = TRUE); Tm <- max(d$peaceyears, na.rm = TRUE)
cat("btscs knots:", kn, " T:", Tm, "\n")
for (i in 1:3) d[[paste0("spline", i)]] <- pmax(d$peaceyears - kn[i], 0)^3 - (Tm - kn[i]) / Tm * d$peaceyears^3

d <- d[ccode1 != 731 & ccode2 != 731]
d <- as.data.frame(d[is.na(mzongo) | mzongo != 1])   # 'if mzongo~=1' (missing counts as ~=1 in Stata)

RHSV <- c("stc_low","dependlow","maxflow","lngdphigh","lngdplow","demochigh","democlow","igos","contiguous","logdist","majpow",
          "allies","logcapratio","spline1","spline2","spline3","peaceyears")
MV <- c("mzmidplusone", RHSV)
BASE <- paste("mzmidplusone ~", paste(RHSV, collapse = " + "))
fitter <- function(add, dat) fit_feglm_tie(as.formula(paste0(BASE, add)), dat, "stc_low", cluster = "ROBUST")

# ---- 2. Reproduce ----
cat("\n==== Reproduction: Erratum Table Model 1 ====\n")
rp <- fitter("", d)
b_pub <- 0.117; n_pub <- 394527
cat("N repro", rp$n, "vs pub", n_pub, "\n")
# published tie is NOT significant (se .093): agreement = both p>=.05
repro <- if (abs(rp$b - b_pub) / abs(b_pub) <= 0.10 && rp$p >= .05) "REPRODUCED" else "NOT_REPRODUCED"
cat("\nb_pub =", b_pub, " b_repro =", rp$b, " p =", rp$p, " status:", repro, "\n")
n_orig <- rp$n
# companion: high-resolution Model 5 (the erratum's only significant STC result, 0.347) if the 1.7 run finished
has_high <- "stc_high" %in% names(d) && any(!is.na(d$stc_high))
if (has_high) {
  B5 <- sub("stc_low", "stc_high", BASE)
  r5 <- fit_feglm_tie(as.formula(B5), d, "stc_high", cluster = "ROBUST")
}

# ---- 3. History merge; sample 1960-2000 so the year restriction drops nothing ----
d <- add_history(d, "ccode1", "ccode2", "year")
d <- d[complete.cases(d[, MV]), c(MV, "year", HIST_VARS, if (has_high) "stc_high")]
SENS <- sens_hist_only(fitter, d, MV)
d <- restrict_sample(d, MV)

# ---- 4. Six models ----
R <- run_six(fitter, d)
C5 <- ""
if (has_high) {
  h0  <- fit_feglm_tie(as.formula(B5), d[!is.na(d$stc_high), ], "stc_high", cluster = "ROBUST", quiet = TRUE)
  hdm <- fit_feglm_tie(as.formula(paste0(B5, " + hist_dyad_w15 + hist_mono_w15")), d[!is.na(d$stc_high), ], "stc_high", cluster = "ROBUST", quiet = TRUE)
  C5 <- sprintf(" Companion M5 stc_high: repro %.3f (pub 0.347), b0=%.3f (p=%.3g) -> b_dm=%.3f (p=%.3g), pct=%.1f.", r5$b, h0$b, h0$p, hdm$b, hdm$p, pct(hdm$b, h0$b))
}

# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "43", dat = d, y = "mzmidplusone", tie = "stc_low", ctrl = setdiff(RHSV, "stc_low"), fe = NULL,
           fit = function(dat, add) { r <- fitter(if (length(add)) paste0(" + ", paste(add, collapse = " + ")) else "", dat); c(b = r$b, se = r$se, p = r$p) })
}

# ---- 5. Results ----
write_row(43, "Lupu & Traag 2013 (erratum)", repro, b_pub, rp$b, "logit, robust SE (fixest::feglm hetero)", "stc_low",
  n_orig, R$m0$n, R, ahc = "Y",
  notes = sprintf("Erratum M1 (low resolution). Full pipeline rebuilt: Louvain (louvain 0.8.2, seed 0 on one optimiser; authors used 0.6.x so partitions differ) + maxflow + R translation of Stata merge; btscs rebuilt (quartile knots %s, T=%g). N repro %d vs pub %d (exact); controls match (peaceyears -0.366 vs -0.368, contiguity 2.71 vs 2.70, const -7.03 vs -7.02) but maxflow -1.34 vs -1.02 and tie 0.075 vs 0.117 (both ns) -> NOT_REPRODUCED attributed to stochastic Louvain partitions under different library versions; not tuned. DV MID at t+1; author controls peace years + splines (Y).%s Dropped %d obs vs orig; %s. Sensitivity %s.",
                  paste(kn, collapse = "/"), Tm, rp$n, n_pub, C5, n_orig - R$m0$n, RS_INFO, SENS))
sink()
