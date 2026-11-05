# Rerun 15: Kenwick & Vasquez 2017 JOP -- models.do model 1 (main text):
#   logit dispute fbdef nfbdef pchaloff pchalneu ln_distance capprop jdem s_un_glo peaceyrs peaceyrs2 peaceyrs3, robust cluster(ddyad)
# Headline tie: nfbdef (target defense pact formed after the directed dyad's birth; the "caveat emptor" finding).
# fbdef (pact from birth) also tracked in notes. Published (author models.smcl): nfbdef .3461796, fbdef -1.326612, N=1,077,992.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work/rerun_common_g7.R"))
suppressPackageStartupMessages(library(haven))
W <- file.path(ROOT, "work/15")
WO <- file.path(OUT, "work/15"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
sink(file.path(WO, "log_15.txt"), split = TRUE)

# 1. Authors' analysis data (Dataverse 10.7910/DVN/WBWWUO jl_modified.dta)
MV <- c("dispute","fbdef","nfbdef","pchaloff","pchalneu","ln_distance","capprop","jdem","s_un_glo","peaceyrs","peaceyrs2","peaceyrs3","ddyad")
d <- as.data.frame(zap_labels(read_dta(file.path(W, "jl_modified.dta"), col_select = all_of(c(MV, "year")))))
stopifnot(nrow(d) > 1e6)
# ddyad = challenger*1000 + target (Johnson & Leeds coding; verified no a==b, all parts valid COW codes)
d$ccode1 <- d$ddyad %/% 1000; d$ccode2 <- d$ddyad %% 1000
stopifnot(all(d$ccode1 != d$ccode2), all(d$ccode1 >= 2), all(d$ccode2 >= 2))
BASE <- "dispute ~ nfbdef + fbdef + pchaloff + pchalneu + ln_distance + capprop + jdem + s_un_glo + peaceyrs + peaceyrs2 + peaceyrs3"
FB <- list()
# Fits are cached to cache/<tag>_<spec>.rds so an OOM kill (shared box) can resume; tag = sample (N rows) + tie.
dir.create(file.path(WO, "cache"), recursive = TRUE, showWarnings = FALSE)
cached_fit <- function(add, dat, tie) {
  key <- file.path(WO, "cache", paste0(nrow(dat), "_", tie, "_", gsub("[^A-Za-z0-9]", "", add), ".rds"))
  if (file.exists(key)) { r <- readRDS(key); cat("\n[cached]", key, "\n"); print(signif(r$all, 5)); return(r) }
  r <- fit_feglm_tie(as.formula(paste0(BASE, add)), dat, tie, cluster = "ddyad")   # fixest for memory (1M rows)
  saveRDS(r, key); r
}
fitter <- function(add, dat) cached_fit(add, dat, "nfbdef")

# 2. Reproduce
cat("\n==== Reproduction: model 1 (full sample) ====\n")
rp <- fitter("", d)
b_pub <- 0.3461796
stopifnot(rp$n == 1077992)
repro <- if (abs(rp$b - b_pub) / abs(b_pub) <= 0.10 && (rp$p < .05)) "REPRODUCED" else "NOT_REPRODUCED"
cat("\nb_pub =", b_pub, " b_repro =", rp$b, " status:", repro, "\n")
n_orig <- rp$n

# 3. History merge + restriction
HC <- file.path(WO, "cache", "d_hist.rds")                     # cache merged data so an OOM restart skips the merge
if (file.exists(HC)) { rm(d); invisible(gc()); d <- readRDS(HC) } else {
  d <- add_history(d, "ccode1", "ccode2", "year")
  d <- d[complete.cases(d[, MV]), c(MV, "year", HIST_VARS)]  # keep only needed columns (memory; shared 8GB box)
  saveRDS(d, HC)
}
invisible(gc())
SENS <- sens_hist_only(fitter, d, MV)                          # before bm filter (see helper)
d <- restrict_sample(d, MV); invisible(gc())

# 4. Six models (fbdef coefficients visible in the log)
R <- run_six(fitter, d)

# fbdef companion: baseline and +dyad+mono only
fb0  <- R$m0;  fb0$b  <- R$m0$all["fbdef", 1];  fb0$p  <- R$m0$all["fbdef", 4]    # same fits, other coefficient
fbdm <- R$mDM; fbdm$b <- R$mDM$all["fbdef", 1]; fbdm$p <- R$mDM$all["fbdef", 4]
cat("\nfbdef: b0 =", fb0$b, "p0 =", fb0$p, "| b_dm =", fbdm$b, "p_dm =", fbdm$p, "pct =", pct(fbdm$b, fb0$b), "\n")

# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  # calls fit_feglm_tie directly (what cached_fit runs): the cache key (rows + tie + term names) would return a stale fit for the placebo columns
  sens_run(id = "15", dat = d, y = "dispute", tie = "nfbdef",
           ctrl = c("fbdef", "pchaloff", "pchalneu", "ln_distance", "capprop", "jdem", "s_un_glo", "peaceyrs", "peaceyrs2", "peaceyrs3"), fe = NULL,
           fit = function(dat, add) { r <- fit_feglm_tie(as.formula(paste0(BASE, if (length(add)) paste0(" + ", paste(add, collapse = " + ")) else "")), dat, "nfbdef", cluster = "ddyad")
             c(b = r$b, se = r$se, p = r$p) })
}

# 5. Results
write_row(15, "Kenwick & Vasquez 2017", repro, b_pub, rp$b, "logit, SE clustered on directed dyad (fixest::feglm, G/(G-1) adj)", "nfbdef",
  n_orig, R$m0$n, R, ahc = "Y",
  notes = sprintf("models.do M1 (all years); b_pub from author log. Author controls cubic peace years (Y). ccodes decoded from ddyad. Companion fbdef: b0=%.3f (p=%.3g) -> b_dm=%.3f (p=%.3g), pct=%.1f. Dropped %d obs vs orig; %s. Sensitivity %s.",
                  fb0$b, fb0$p, fbdm$b, fbdm$p, pct(fbdm$b, fb0$b), n_orig - R$m0$n, RS_INFO, SENS))
sink()
