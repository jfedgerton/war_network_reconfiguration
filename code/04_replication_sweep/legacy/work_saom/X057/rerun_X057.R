# X057 Pomeroy, Dasandi & Mikhaylov (2019) PLOS ONE -- Table 1 Model 2 (weak-signal communities), btergm.
# Authors: btergm(mid.nw_full ~ edges + altkstar(2,fixed) + cycle(4) + gwesp(0,fixed) + nodemix("democ",-3)
#   + edgecov(contig, capratio, tradedep, igosec, igoecon, weak_comms) + memory(autoregression, lag 1), R = 2000),
#   seed 1912, fit once per each of 6 embedding settings and coefficients averaged (published: weak comms -0.60 [-1.36,-0.08]).
# Added terms (edgecov, same estimator): mD hist_dyad_w15; mH + hist_mono_w15; m26 hist_dyad_w26 + hist_mono_w26;
#   mP placebo bm_cinc_w15. Model years 1971-1990 (MID onset nets), memory drops 1971.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
WD <- file.path(OUT, "work_saom/X057"); dir.create(WD, recursive = TRUE, showWarnings = FALSE)
setwd(WD)   # res_X057_R<RB>.rds is written here (pipeline output)
DD <- file.path(ROOT, "work_saom/X057")   # author data (read-only)
suppressPackageStartupMessages({library(btergm); library(network); library(countrycode); library(dplyr)})
cat(R.version.string, " btergm", as.character(packageVersion("btergm")), " ergm", as.character(packageVersion("ergm")), "\n")
args <- commandArgs(TRUE); RB <- if (length(args)) as.integer(args[1]) else 2000L
source(file.path(CODE, "hist/add_history.R"))

mid.nw_full <- readRDS(file.path(DD, "data/pauls_cranmer/prepped/system/mid.nw_full.rds"))
directcont.nw_full <- readRDS(file.path(DD, "data/pauls_cranmer/prepped/system/directcont.nw_full.rds"))
capratio_full <- readRDS(file.path(DD, "data/pauls_cranmer/prepped/system/capratio_full.rds"))
tradedep_full <- readRDS(file.path(DD, "data/pauls_cranmer/prepped/system/tradedep_full.rds"))
igosec_full <- readRDS(file.path(DD, "data/pauls_cranmer/prepped/system/igosec_full.rds"))
igoecon_full <- readRDS(file.path(DD, "data/pauls_cranmer/prepped/system/igoecon_full.rds"))
stopifnot(length(mid.nw_full) == 20)
years <- 1971:1990

## ---- history edgecov matrices aligned to each year's node set ----
cm_rev <- c(CSK = 315, DDR = 265, YDYE = 680, YUG = 345, SVM = 816, DEU = 260, YAR = 678, HKG = 999, COK = 998)
hmat <- function(var) lapply(seq_along(years), function(i) {
  nm <- network::get.vertex.attribute(mid.nw_full[[i]], "vertex.names")
  cc <- suppressWarnings(countrycode(nm, "iso3c", "cown", custom_match = cm_rev))
  stopifnot(!any(is.na(cc)))
  g <- expand.grid(a = seq_along(nm), b = seq_along(nm)); g$c1 <- cc[g$a]; g$c2 <- cc[g$b]; g$year <- years[i]
  g <- add_history(g, "c1", "c2", "year")
  v <- g[[var]]; v[is.nan(v)] <- NA
  if (var == "bm_cinc_w15") v[is.na(v)] <- mean(v, na.rm = TRUE)   # few missing CINC (placebo only)
  m <- matrix(v, length(nm), length(nm)); diag(m) <- 0; dimnames(m) <- list(nm, nm); m
})
hd15 <- hmat("hist_dyad_w15"); hm15 <- hmat("hist_mono_w15")
hd26 <- hmat("hist_dyad_w26"); hm26 <- hmat("hist_mono_w26"); bmc <- hmat("bm_cinc_w15")
cat("share dyads with hist_dyad_w15>0 (1990):", mean(hd15[[20]] > 0), "\n")
stopifnot(all(sapply(seq_along(years), function(i) all(dim(hd15[[i]]) == dim(as.matrix(mid.nw_full[[i]]))))))

files <- c("weak_ties_50d_15x.rds","weak_ties_50d_25x.rds","weak_ties_100d_15x.rds",
           "weak_ties_100d_25x.rds","weak_ties_200d_15x.rds","weak_ties_200d_25x.rds")
base_f <- "mid.nw_full ~ edges + altkstar(2, fixed = TRUE) + cycle(4) + gwesp(0, fixed = TRUE) + nodemix('democ', -3) +
  edgecov(directcont.nw_full) + edgecov(capratio_full) + edgecov(tradedep_full) + edgecov(igosec_full) +
  edgecov(igoecon_full) + edgecov(weak_comms) + memory(type = 'autoregression', lag = 1)"
adds <- list(m0 = "", mH = " + edgecov(hd15) + edgecov(hm15)",
             m26 = " + edgecov(hd26) + edgecov(hm26)", mP = " + edgecov(bmc)")
out <- list()
for (f in files) {
  weak_comms <- readRDS(file.path(DD, "data/comm_detection/detected_comms/weak", f))[1970:1989]
  for (m in names(adds)) {
    set.seed(1912)
    t0 <- Sys.time()
    Rm <- if (m %in% c("m0","mH")) RB else min(RB, 100L)   # secondary models: fewer bootstrap reps (point estimate unaffected)
    fit <- btergm(as.formula(paste0(base_f, adds[[m]])), R = Rm, parallel = "multicore", ncpus = 2, verbose = FALSE)
    s <- summary(fit)   # bootstrap percentile CI
    co <- coef(fit)
    bs <- fit@boot$t
    se <- apply(bs, 2, sd)
    pb <- apply(bs, 2, function(x) 2 * min(mean(x <= 0), mean(x >= 0)))   # bootstrap two-sided p
    tn <- sub("\\[\\[i\\]\\]$", "", names(co))
    stopifnot(all(rownames(s) == names(co)), "edgecov.weak_comms" %in% tn)
    out[[paste(f, m)]] <- data.frame(file = f, model = m, term = tn, b = co, se_boot = se, p_boot = pb,
                                     ci_lo = s[, 3], ci_hi = s[, 4], R = Rm, nobs = fit@nobs)
    k <- which(tn == "edgecov.weak_comms")
    cat(sprintf("%s %s: weak_comms b=%.3f [%.3f, %.3f]  (%.1f s)\n", f, m, co[k], s[k, 3], s[k, 4],
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
}
res <- bind_rows(out)
saveRDS(res, sprintf("res_X057_R%d.rds", RB))
avg <- res %>% group_by(model, term) %>% summarise(across(c(b, se_boot, p_boot, ci_lo, ci_hi), mean), .groups = "drop")
print(as.data.frame(avg), digits = 3)
