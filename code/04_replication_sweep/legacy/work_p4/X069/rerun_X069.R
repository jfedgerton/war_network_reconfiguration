# X069 Tomashevskiy 2025 JPR "Friends and partners: Estimating latent affinity networks with the graphical LASSO".
# Headline = Table 5 (authors' 'gartzke replication.do'), Model 5 = Gartzke (2007) baseline + GLASSO affinity:
#   logit warl demlo demhi deplo contig logdstab majpdyds alliesr lncaprt t t2 t3 friendlassonotime, cluster(dyad)
#   (Model 6 adds all rival affinity measures at once; Model 5 is the first/clean column for the key term.)
# Key term (abstract, fixed before estimating): friendlassonotime. Evidence: "affinity ... was associated with decreased conflict".
# CONTROL: demlo (weak-link democracy).
# finaldata.csv has NO year column -> years recovered by exact match to Gartzke's capitalistpeace_012007 data (work_demcap/C7/slim_C7.rds)
#   on dyad + (demlo, demhi, deplo, lncaprt, logdstab, contig, majpdyds); only unique matches kept.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p4/P4_common.R"))
WO <- file.path(OUT, "work_p4/X069"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work_p4/X069"); sink(file.path(WO, "log_X069.txt"), split = TRUE)
d <- read.csv(file.path(W, "replication/finaldata.csv"))
stopifnot(nrow(d) == 176836)
g <- readRDS(file.path(ROOT, "work_demcap/C7/slim_C7.rds"))
g <- g[g$year >= 1965 & g$year <= 2000, c("dyadid","statea","stateb","year","demlo","demhi","deplo","lncaprt","logdstab","contig","majpdyds","maoznewl","alliesr","sun2cati")]
key <- function(x, id) paste(id, x$demlo, x$demhi, sprintf("%.4f", x$deplo), sprintf("%.4f", x$lncaprt), sprintf("%.4f", x$logdstab), x$contig, x$majpdyds)
d$k <- key(d, d$dyad); g$k <- key(g, g$dyadid)
dupg <- g$k[duplicated(g$k)]
g1 <- g[!(g$k %in% dupg), ]
d$row <- seq_len(nrow(d))
m <- merge(d, g1[, c("k","statea","stateb","year","maoznewl","alliesr","sun2cati")], by = "k", all.x = TRUE, suffixes = c("", ".g"))
m <- m[order(m$row), ]; stopifnot(nrow(m) == nrow(d))
cat("rows:", nrow(d), " matched year:", sum(!is.na(m$year)), " share:", mean(!is.na(m$year)), "\n")
mm <- m[!is.na(m$year), ]
cat("year range:", range(mm$year), "\n"); print(table(mm$year)[1:5])
cat("agreement warl==maoznewl:", mean(mm$warl == mm$maoznewl, na.rm = TRUE), " alliesr:", mean(mm$alliesr == mm$alliesr.g, na.rm = TRUE),
    " sun2cati:", mean(abs(mm$sun2cati - mm$sun2cati.g) < 1e-4, na.rm = TRUE), "\n")
stopifnot(!anyDuplicated(mm[, c("dyad","year")]))
f0 <- warl ~ demlo + demhi + deplo + contig + logdstab + majpdyds + alliesr + lncaprt + t + t2 + t3 + friendlassonotime
K <- c("friendlassonotime", "demlo")
r <- fit_keys(f0, d, K, "dyad", file.path(WO, "cache_repro.rds")); n_orig <- r$n
rm <- fit_keys(f0, mm, K, "dyad", file.path(WO, "cache_repro_matched.rds"))
# DV warl: Gartzke's maoznewl is MID onset in t+1 (lead); history anchored on outcome year (year+1), as in C1.
mm$yr_out <- mm$year + 1
ds <- common_sample(mm, f0, "dyad", "statea", "stateb", "yr_out", min_share = 0.8)
M <- run_six(f0, ds, K, "dyad", WO)
# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "X069", dat = ds, y = "warl", tie = "friendlassonotime",
           ctrl = c("demlo", "demhi", "deplo", "contig", "logdstab", "majpdyds", "alliesr", "lncaprt", "t", "t2", "t3"), fe = NULL,
           fit = function(dat, add) { r <- fit_keys(if (length(add)) update(f0, as.formula(paste(". ~ . +", paste(add, collapse = " + ")))) else f0, dat, "friendlassonotime", "dyad")$k$friendlassonotime
             c(b = unname(r[["b"]]), se = unname(r[["se"]]), p = unname(r[["p"]])) })
}
sc <- "Tomashevskiy 2025 JPR"; est <- "logit, dyad-clustered SE (fixest feglm)"
hc <- "Y (peace years t, t^2, t^3 only)"
nt <- paste0("Table 5 Model 5 (Gartzke 2007 capitalist-peace baseline + GLASSO affinity; 1970-1992). Paper table not accessible (JPR paywall; archive has no log) -> UNVERIFIED. ",
  "finaldata.csv lacks year: years recovered by exact match to Gartzke's data for ", nrow(mm), "/", nrow(d), " rows; b on full file ", round(r$k$friendlassonotime[["b"]], 4),
  " vs matched subset ", round(rm$k$friendlassonotime[["b"]], 4), ". DV = MID onset t+1; history anchored on outcome year. Common sample N ", M$m0$n, ".")
write_p4("X069", sc, "TIE", "KEY", "friendlassonotime", "abstract: 'measure of affinity generated with the three-step GLASSO ... associated with decreased conflict'",
  "UNVERIFIED", NA, r$k$friendlassonotime[["b"]], est, n_orig, M, hc, nt)
write_p4("X069", sc, "DEMOCRACY", "CONTROL", "demlo", "control in Table 5 (weak-link democracy, Gartzke 2007 baseline)",
  "UNVERIFIED", NA, r$k$demlo[["b"]], est, n_orig, M, hc, paste(nt, "bD NOT a placebo for this term."))
sink()
