# P9_062 Hegre, Bernhard & Teorell 2020 JCR 64(1) "Civil Society and the Democratic Peace" (10.1177/0022002719850620). DEMOCRACY.
# DATA BUILD ONLY (not registered; published values pending supplement Table SF-3).
# Inputs (PI-staged 2026-10-10, copied to work_p9/P9_062): V-Dem-DS-CY-v6.2_shortened.dta, V-DemocraticPeace_temp2.dta, the four
#   do-files and the authors' log V-democraticpeace.log (the log stops before the models but records the sample counts).
# Translates AssembleData2_Replication.do + the sample block of Analysis_Replication.do for the headline civil-society model
#   (do-file model 7): logit mzfatal_new cspart_i cspart_j cspart_ij lncapi lncapj lnpopi lnpopj dircont lndstab lnNt py1
#   if include_vdem == 1, cluster(dyadid).
# Only the pieces the headline model and the sample rule need are rebuilt. The factor-analysis indices (multiparty -> schump,
#   mefree) enter only through their MISSINGNESS (include_vdem requires schump_i/j and mefree_i/j non-missing); a factor
#   score from -predict- is missing exactly when one of its inputs is missing, so their values are not computed here.
# Stata's -impute- (cspart filled from v2cscnsult v2csprtcpt v2csgender) is translated as: for each pattern of available
#   predictors, OLS of y on those predictors over rows with y and them non-missing (on the dyad-year data, as in the do-file);
#   rows with no predictor get the constant-only fit (mean of y), which reproduces the article's N of 541,560.
# Validation against the authors' log: include == 1 on 639,194 rows; include_vdem keeps 553,484 rows (264,617 deleted of
#   818,101 = 801,169 dyad-years + 16,932 unmatched V-Dem rows appended by -merge-); dircont fill-forward changes 5,370 rows. Fails loudly if any of these differ.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); stopifnot(nzchar(ROOT), nzchar(OUT))
suppressPackageStartupMessages({library(dplyr); library(fixest)})
WD <- file.path(ROOT, "work_p9/P9_062"); W <- file.path(OUT, "work_p9/P9_062"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_prep_P9_062.txt"), split = TRUE)
inf <- function(x) ifelse(is.na(x), Inf, x)                 # Stata missing = +infinity in comparisons

# ---- 1. V-Dem country-year recodes (AssembleData2 lines 13-30, 52-55) ----
v <- as.data.frame(haven::zap_labels(haven::read_dta(file.path(WD, "V-Dem-DS-CY-v6.2_shortened.dta"))))
stopifnot(nrow(v) == 19858, !anyDuplicated(v[, c("country_id", "year")]))
v <- v[order(v$country_id, v$year), ]
inwin <- v$year >= inf(v$codingstart) & v$year <= inf(v$codingend) & !(v$year > inf(v$gapstart) & v$year < inf(v$gapend))
v$psopp <- pnorm(v$v2psoppaut)
v$psopp[!is.na(v$v2psbars) & is.na(v$v2psoppaut) & inwin] <- 0
cat("v2psoppaut_fill: missing", sum(is.na(pnorm(v$v2psoppaut))), "(log 4,159); set to 0:", sum(!is.na(v$v2psbars) & is.na(v$v2psoppaut) & inwin), "(log 689)\n")
v$elmul <- pnorm(v$v2elmulpar); n0 <- sum(is.na(v$elmul))
ch <- 0L                                                    # replace x = l.x if x == . & window: sequential (cascading) as Stata
for (i in seq_len(nrow(v))[-1]) if (is.na(v$elmul[i]) && inwin[i] && v$country_id[i] == v$country_id[i - 1] && v$year[i] == v$year[i - 1] + 1 &&
                                    !is.na(v$elmul[i - 1])) { v$elmul[i] <- v$elmul[i - 1]; ch <- ch + 1L }
cat("v2elmulpar_fill: missing", n0, "(log 16,513); lag-filled", ch, "(log 8,397)\n")
z <- which(v$v2x_elecreg %in% 0 & (is.na(v$elmul) | v$elmul != 0)); v$elmul[v$v2x_elecreg %in% 0] <- 0
cat("set to 0 where v2x_elecreg == 0:", length(z), "(log 6,743)\n")
v$multiparty_ok <- complete.cases(v[, c("v2psparban", "v2psbars", "psopp", "elmul")])
v$schump_ok <- v$multiparty_ok & !is.na(v$v2x_accex) & !is.na(v$v2xel_frefair)
v$mefree_ok <- complete.cases(v[, c("v2mecenefm", "v2meharjrn", "v2meslfcen", "v2mebias", "v2mecrit", "v2merange")])
vk <- v[, c("COWcode", "year", "v2x_polyarchy", "v2x_egal", "v2x_suffr", "v2x_cspart", "v2cscnsult", "v2csprtcpt", "v2csgender", "schump_ok", "mefree_ok")]
vk <- vk[!is.na(vk$COWcode), ]; stopifnot(!anyDuplicated(vk[, c("COWcode", "year")]))   # merge m:1 requires unique keys

# ---- 2. dyad-year data, merge for state a then b, impute cspart (AssembleData2 lines 96-170) ----
d <- as.data.frame(haven::zap_labels(haven::read_dta(file.path(WD, "V-DemocraticPeace_temp2.dta"))))
cat("temp2 rows", nrow(d), "\n"); stopifnot(nrow(d) == 801169)
imp <- function(y, X, yx = NULL, Xx = NULL) {                 # Stata -impute-: per missing pattern of X, OLS on complete rows
  # yx/Xx: extra rows that are in Stata's data but not in d (V-Dem country-years appended unmatched by -merge-); they enter
  # the regressions only
  out <- y; miss <- which(is.na(y)); if (!length(miss)) return(out)
  if (!is.null(yx)) { y <- c(y, yx); X <- rbind(as.data.frame(X), setNames(as.data.frame(Xx), names(X))) }
  pat <- apply(!is.na(X[miss, , drop = FALSE]), 1, paste, collapse = "")
  for (p in unique(pat)) {
    use <- strsplit(p, "")[[1]] == "TRUE"; rows <- miss[pat == p]
    if (!any(use)) { out[rows] <- mean(y[!is.na(y)]); next }   # no predictor available: constant-only regression (mean of y);
                                                               # this reproduces the authors' N (541,560 in article Table 2)
    Xu <- as.matrix(X[, use, drop = FALSE]); ok <- !is.na(y) & complete.cases(Xu)
    b <- coef(lm.fit(cbind(1, Xu[ok, , drop = FALSE]), y[ok]))
    out[rows] <- drop(cbind(1, Xu[rows, , drop = FALSE]) %*% b)
  }
  out
}
n_um <- integer(0)
for (s in c("a", "b")) {
  m <- vk; names(m)[names(m) != "year"] <- paste0(names(m)[names(m) != "year"], "_", s)
  um <- anti_join(m, d, by = setNames(c(paste0("state", s), "year"), c(paste0("COWcode_", s), "year")))   # _merge == 2 rows
  d <- left_join(d, m, by = setNames(c(paste0("COWcode_", s), "year"), c(paste0("state", s), "year")))
  stopifnot(nrow(d) == 801169); n_um <- c(n_um, nrow(um))
  cs <- paste0("v2x_cspart_", s); sub <- paste0(c("v2cscnsult_", "v2csprtcpt_", "v2csgender_"), s); before <- sum(is.na(d[[cs]]))
  d[[cs]] <- imp(d[[cs]], d[, sub], um[[cs]], um[, sub])
  cat("cspart_", s, ": missing before impute", before, " after", sum(is.na(d[[cs]])), "\n", sep = "")
}
cat("unmatched V-Dem rows appended by the two merges:", sum(n_um), "(log: dyadid missing on 16,932 rows)\n"); stopifnot(sum(n_um) == 16932)
# ---- 3. size ordering (population), i/j terms (AssembleData2 lines 200-300) ----
d$a_largecnt <- ifelse(is.na(d$lnpopa) | is.na(d$lnpopb), NA, ifelse(d$lnpopa >= d$lnpopb, 1, 0))
pick <- function(a, b, first) ifelse(d$a_largecnt == 1, if (first) a else b, if (first) b else a)
d$lnpopi <- pick(d$lnpopa, d$lnpopb, TRUE); d$lnpopj <- pick(d$lnpopa, d$lnpopb, FALSE)
d$lncapi <- pick(d$cap_a, d$cap_b, TRUE);   d$lncapj <- pick(d$cap_a, d$cap_b, FALSE)
d$demi <- pick(d$polity_a, d$polity_b, TRUE); d$demj <- pick(d$polity_a, d$polity_b, FALSE)
for (v0 in c("v2x_cspart", "v2x_polyarchy", "v2x_egal", "v2x_suffr")) {
  d[[paste0(v0, "_i")]] <- pick(d[[paste0(v0, "_a")]], d[[paste0(v0, "_b")]], TRUE)
  d[[paste0(v0, "_j")]] <- pick(d[[paste0(v0, "_a")]], d[[paste0(v0, "_b")]], FALSE)
}
okpick <- function(v0, first) { x <- pick(as.numeric(d[[paste0(v0, "_a")]]), as.numeric(d[[paste0(v0, "_b")]]), first); !is.na(x) & x == 1 }
d$cspart_i <- d$v2x_cspart_i; d$cspart_j <- d$v2x_cspart_j; d$cspart_ij <- d$cspart_i * d$cspart_j
# ---- 4. dircont fill-forward after 2000 (sequential, as Stata) ----
d <- d[order(d$statea, d$stateb, d$year), ]
ch <- 0L; dc <- d$dircont; sa <- d$statea; sb <- d$stateb; yr <- d$year
for (i in seq_along(dc)[-1]) if (yr[i] > 2000 && dc[i] %in% 0 && dc[i - 1] %in% 1 && sa[i] == sa[i - 1] && sb[i] == sb[i - 1]) { dc[i] <- 1; ch <- ch + 1L }
d$dircont <- dc; cat("dircont fill-forward changes:", ch, "(log 5,370)\n"); stopifnot(ch == 5370)
# ---- 5. sample (Analysis_Replication lines 74-91) ----
inc <- !is.na(d$lnpopi) & !is.na(d$lnpopj) & !is.na(d$lncapi) & !is.na(d$lncapj) & !is.na(d$demi) & !is.na(d$demj) &
       !is.na(d$lndstab) & !is.na(d$dircont) & !is.na(d$py1) & !is.na(d$lnNt) & d$year <= 2010
cat("include:", sum(inc), "(log 639,194)\n")
incv <- inc & !is.na(d$v2x_polyarchy_i) & !is.na(d$v2x_polyarchy_j) & !is.na(d$v2x_egal_i) & !is.na(d$v2x_egal_j) &
        okpick("mefree_ok", TRUE) & okpick("mefree_ok", FALSE) & okpick("schump_ok", TRUE) & okpick("schump_ok", FALSE) &
        !is.na(d$v2x_suffr_i) & !is.na(d$v2x_suffr_j)
cat("include_vdem:", sum(incv), "(log keeps 553,484; deletes 264,617)\n")
stopifnot(sum(inc) == 639194)
if (sum(incv) != 553484) cat("WARNING: include_vdem count differs from the authors' log by", sum(incv) - 553484, "\n")
d$dyadid <- d$statea * 1000 + d$stateb
ds <- d[incv, c("statea", "stateb", "year", "dyadid", "mzfatal_new", "cspart_i", "cspart_j", "cspart_ij", "lncapi", "lncapj", "lnpopi", "lnpopj",
                "dircont", "lndstab", "lnNt", "py1")]
saveRDS(ds, file.path(W, "P9_062_model7_data.rds"))
# ---- 6. headline model, full sample (one fit; published values pending) ----
m <- feglm(mzfatal_new ~ cspart_i + cspart_j + cspart_ij + lncapi + lncapj + lnpopi + lnpopj + dircont + lndstab + lnNt + py1,
           data = ds, family = binomial("logit"), vcov = ~dyadid, notes = FALSE, nthreads = 1)
cat("\nmodel 7 (civil society), N =", nobs(m), "(article Table 2: 541,560) AIC =", AIC(m), "(article AICfull 5,854.5)\n"); print(round(coeftable(m), 6))
sink()
