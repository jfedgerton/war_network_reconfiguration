# D1: Dafoe, Oneal & Russett 2013 ISQ, "The Democratic Peace: Weighing the Evidence and Cautious Inference"
# Headline: Table 1, col. 1 (All MIDs, originators; 1816-2001[data end 2000]) -- the paper's only coefficient table:
#   logit mzmid1_orig jtdem propdem jtdemxpropdem jtaut jtautxpropdem mixed mixedxpropdem contiguity lndstab pwin_lrg
#         lrgcap allies lnNt py _spline1-3, cluster(dyadid)          (Oneal's clarify_3interactions_111412.do)
# b_pub from authors' posted output (clarify_3interactions_111412.docx): jtdem .6292162 (p=.349), jtdemxpropdem -4.036512 (p=.073), N=537603.
# Key terms: jtdem (joint democracy; constituent of JtDem*PropDem, i.e. effect at PropDem=0) and jtdemxpropdem (the interaction
#   that carries the paper's G&W claim). Abstract: "the peacefulness of democratic pairs has actually increased as the proportion of democracies grew".
# Variable construction is an R translation of the do-file (Stata missing semantics handled explicitly).
# btscs: see rule below; knots = quartiles (type 2) of py; checked against the log's spline values at py==3.
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_demcap/D_common.R"))
suppressPackageStartupMessages(library(haven))
WD <- file.path(ROOT, "work_demcap/D1"); W <- file.path(OUT, "work_demcap/D1"); dir.create(W, recursive = TRUE, showWarnings = FALSE); sink(file.path(W, "log_D1.txt"), split = TRUE)

# ---- 1. load ----
d <- read_dta(file.path(WD, "dor.dta"), col_select = c(statea, stateb, year, cap_a, cap_b, alliance, contig, distance, numstate,
                                                       polity2a, polity2b, mzmid1_orig))
d <- as.data.frame(zap_labels(d))
stopifnot(nrow(d) == 656317, !anyDuplicated(d[, c("statea", "stateb", "year")]), all(d$mzmid1_orig %in% c(0, 1, NA)))

# propdem: share of states with polity2>6.9 per year (unique state-years stacked from both sides)
st <- rbind(data.frame(state = d$statea, year = d$year, p = d$polity2a), data.frame(state = d$stateb, year = d$year, p = d$polity2b))
st <- st[order(st$state, st$year, is.na(st$p), st$p), ]; st <- st[!duplicated(st[, c("state", "year")]), ]
st$democ <- ifelse(is.na(st$p), NA, as.numeric(st$p > 6.9))
pd <- aggregate(democ ~ year, data = st, FUN = mean); names(pd)[2] <- "propdem"
d <- merge(d, pd, by = "year", all.x = TRUE)
cat("propdem range:", range(d$propdem), "\n")   # paper: historical range 0.032 - 0.427
stopifnot(abs(min(d$propdem) - 0.032) < 0.002, abs(max(d$propdem) - 0.427) < 0.002)

both <- !is.na(d$polity2a) & !is.na(d$polity2b)
d$smldem <- ifelse(both, pmin(d$polity2a, d$polity2b), NA); d$lrgdem <- ifelse(both, pmax(d$polity2a, d$polity2b), NA)
bc <- !is.na(d$cap_a) & !is.na(d$cap_b)
d$lrgcap <- ifelse(bc, pmax(d$cap_a, d$cap_b), NA)
d$pwin_lrg <- ifelse(bc & (d$cap_a + d$cap_b) > 0, d$lrgcap / (d$cap_a + d$cap_b), NA)   # Stata x/0 = missing
d$dyadid <- d$statea * 1000 + d$stateb
d$lndstab <- ifelse(d$distance > 0, log(d$distance), NA)                               # Stata ln(0) = missing
d$lnNt <- log(d$numstate) - 3.135; d$lnNt[d$contig < 6] <- 0
d$jtdem <- ifelse(is.na(d$smldem), NA, as.numeric(d$smldem > 6)); d$jtdemxpropdem <- d$jtdem * d$propdem
d$jtaut <- ifelse(is.na(d$lrgdem), NA, as.numeric(d$lrgdem < -6)); d$jtautxpropdem <- d$jtaut * d$propdem
d$mixed <- ifelse(is.na(d$lrgdem), NA, as.numeric(d$smldem < -6 & d$lrgdem > 6)); d$mixedxpropdem <- d$mixed * d$propdem
d$contiguity <- abs(d$contig - 6); d$allies <- abs(d$alliance - 4)

# ---- btscs mzmid1_orig year dyadid, gen(py) nspline(3) ----
d <- d[order(d$statea, d$stateb, d$year), ]
# Rule used (diag_btscs*.R): py = years elapsed (time-gap aware, all rows incl. missing-outcome ongoing years) since the dyad's first
# year or since the last event==1 year (reset to 0 the row after an event). This is the only variant tried that reproduces the log's
# knots (8,19,36; T=184); the non-missing-rows-only rule used in earlier reruns gives knots 8,19,35 here.
first <- c(TRUE, d$dyadid[-1] != d$dyadid[-nrow(d)]); ev <- d$mzmid1_orig; y <- d$year
py <- numeric(nrow(d))
for (i in seq_len(nrow(d))) py[i] <- if (first[i]) 0 else if (isTRUE(ev[i - 1] == 1)) 0 else py[i - 1] + (y[i] - y[i - 1])
d$py <- py
kn <- quantile(d$py, c(.25, .5, .75), type = 2, na.rm = TRUE); Tm <- max(d$py, na.rm = TRUE)
for (i in 1:3) d[[paste0("spline", i)]] <- pmax(d$py - kn[i], 0)^3 - (Tm - kn[i]) / Tm * d$py^3
cat("btscs knots:", kn, " T:", Tm, "\n")
s3 <- unlist(d[which(d$py == 3)[1], c("spline1", "spline2", "spline3")]); print(s3)
stopifnot(max(abs(s3 - c(-25.82609, -24.21196, -21.71739))) < 0.01)   # authors' log: summ _spl* if py==3

# ---- 2. reproduce ----
f0 <- mzmid1_orig ~ jtdem + propdem + jtdemxpropdem + jtaut + jtautxpropdem + mixed + mixedxpropdem + contiguity + lndstab +
  pwin_lrg + lrgcap + allies + lnNt + py + spline1 + spline2 + spline3
K <- c("jtdem", "jtdemxpropdem")
r <- fit_keys(f0, d, K, "dyadid", file.path(W, "cache_repro.rds")); n_orig <- r$n
cat("N repro:", n_orig, " (pub 537603)\n")
bp <- c(jtdem = .6292162, jtdemxpropdem = -4.036512); ps <- c(jtdem = .349, jtdemxpropdem = .073)
st_k <- sapply(K, function(k) { b <- bp[[k]]; attr(b, "sig") <- ps[[k]] < .05; repro_status(b, r$k[[k]]) })
print(st_k)

# ---- 3-4. history, common sample (1822-2011), six models ----
d <- d[, c(all.vars(f0), "dyadid", "statea", "stateb", "year")]
ds <- common_sample(d, f0, "dyadid", "statea", "stateb", "year"); rm(d); gc()
M <- run_six(f0, ds, K, "dyadid", W)
# ---- sensitivity hook (split only): runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  # fit_keys with cache = NULL (what run_six runs, uncached): run_six would read the cache_m0/mD/mDM files
  sens_split(id = "D1", dat = ds,
             fit = function(dat, add) { r <- fit_keys(if (length(add)) update(f0, as.formula(paste(". ~ . +", paste(add, collapse = " + ")))) else f0, dat, K, "dyadid")$k$jtdemxpropdem
               c(b = unname(r[["b"]]), se = unname(r[["se"]]), p = unname(r[["p"]])) })
}

sc <- "Dafoe, Oneal & Russett 2013 ISQ"; est <- "logit, dyad-clustered SE (fixest feglm)"
nt <- paste0("Table 1 col.1 (all MID onsets, originators, 1816-2000; Oneal clarify_3interactions do/log). Authors have py + 3 cubic splines. ",
  "Common sample drops ", n_orig - M$m0$n, " obs vs orig (years 1816-1821 + benchmark missingness). N repro ", n_orig, " vs pub 537603 (exact); ",
  "btscs peace-years rebuilt in R (knots match log: 8,19,36) but py coef -.230 vs log -.256, so repro is close not exact. ")
write_dc("D1", sc, "DEMOCRACY", "jtdem",
  "abstract: 'peacefulness of democratic pairs has actually increased as the proportion of democracies grew'",
  st_k[["jtdem"]], bp[["jtdem"]], r$k$jtdem[["b"]], est, n_orig, M,
  paste0(nt, "jtdem is the CONSTITUENT of JtDem*PropDem (effect at PropDem=0; pub .63, n.s.); b0 small/noisy so pct_* are unstable. bD (lagged joint democracy) NOT a placebo."))
write_dc("D1", sc, "DEMOCRACY", "jtdemxpropdem",
  "abstract: 'peacefulness of democratic pairs has actually increased as the proportion of democracies grew'",
  st_k[["jtdemxpropdem"]], bp[["jtdemxpropdem"]], r$k$jtdemxpropdem[["b"]], est, n_orig, M,
  paste0(nt, "Interaction JtDem*PropDem (pub -4.04, p=.073). bD (lagged joint democracy) NOT a placebo."))
sink()
