# ESTAUDIT 47 (phase 6): copy of rerun_verify_47.R. ONLY change: duration-dependence terms rebuilt as the Beck-Katz-Tucker
# (btscs) cubic-spline basis the article cites (methods: "peace years since last fatal MID, from which we generated a cubic spline
# with three interior knots (Beck, Katz, and Tucker 1998)"). Before: truncated-power terms (py-k)_+^3 at ad-hoc knots 1,4,7 (not the
# BKT basis). Now: s_i = (py-k_i)_+^3 - (T-k_i)/T * py^3, knots = quartiles (Stata type 2) of py, T = max(py) -- the btscs rule verified
# against an author log in work_demcap/D1/rerun_D1.R. Knots/T computed on py over the full 1949-2011 directed-dyad panel
# (what btscs sees), then evaluated at py1 (t+1), as the paper leads peace years and splines to t+1.
# cid 47 Feldman, Eiran & Rubin 2021 JCR "Naval Power and Effects of Third-Party Trade on Conflict". REBUILD.
# No archive found. Design documented from OA article (SAGE): Table 1 Model 2 ("All states"):
#   directed dyad-years 1950-2010; DV = A takes first militarized action in a new MID with >=1 battle death, coded t+1
#   (peace years also t+1); logit, SE clustered by directed dyad. Key term: 3rd x A's naval (b_pub 0.014, SE 0.006).
# Rebuild choices (deviations flagged in notes): GML MIDs (paper: COW MID presumably), SDP GDP in constant 2011 USD
#   (paper: current USD), allies = COW defense/entente, BKT cubic spline knots (1,4,7) for peace years.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work/rerun_helpers_recA.R"))
source(file.path(CODE, "work_p8/round2_rule.R"))   # rebuild_verdict(): round-2 rule with registry alpha_pub (2026-10-11)
suppressPackageStartupMessages({library(peacesciencer); library(dplyr); library(data.table)})
WO <- file.path(OUT, "work/47"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
sink(file.path(WO, "log_est_47.txt"), split = TRUE)
source(file.path(CODE, "vb5/verify_helpers.R"))
set.seed(123)

# ---- 1. Directed dyad-years 1949-2011 + standard covariates
d <- create_dyadyears(system = "cow", directed = TRUE, subset_years = 1949:2011) %>%
  add_contiguity() %>% add_capital_distance() %>% add_nmc() %>%
  add_democracy() %>% add_sdp_gdp() %>% add_cow_trade()
d <- as.data.table(d)
cat("dyad-years:", nrow(d), "\n")
stopifnot(!anyDuplicated(d[, .(ccode1, ccode2, year)]))

d[, contig_ord := fifelse(conttype == 0, 6, conttype)]          # 1 land ... 5 400mi water, 6 none (paper: ordinal, coef < 0)
d[, ldist := log(capdist + 1)]
# add_cow_alliance() is disabled in this peacesciencer build; merge the bundled cow_alliance object (made symmetric) directly
ca <- as.data.table(cow_alliance)[, .(ccode1, ccode2, year, al = as.integer(cow_defense == 1 | cow_entente == 1))]
ca <- unique(rbind(ca, ca[, .(ccode1 = ccode2, ccode2 = ccode1, year, al)]))[, .(al = max(al)), by = .(ccode1, ccode2, year)]
d <- merge(d, ca, by = c("ccode1","ccode2","year"), all.x = TRUE)
d[, alliance := fifelse(is.na(al), 0L, al)][, al := NULL]
d[, capratio := log(cinc1 / (cinc1 + cinc2))]
d[, democ := pmin(polity21, polity22)]
d[, gdpA := wbgdp2011est1][, gdppcA := wbgdppc2011est1]          # already logged in peacesciencer
d[flow1 < 0, flow1 := NA][flow2 < 0, flow2 := NA]                  # -9 = missing
d[, dyadtrade := flow1 + flow2]                                    # millions current USD
tsy <- as.data.table(cow_trade_sy)[, .(ccode, year, tot = fifelse(imports < 0 | exports < 0, NA_real_, imports + exports))]
d <- merge(d, tsy[, .(ccode1 = ccode, year, totA = tot)], by = c("ccode1", "year"), all.x = TRUE)
d[, third := log1p(pmax(totA - dyadtrade, 0))]                     # ln(A's total trade less dyadic flow)
d[, tdepAB := dyadtrade * 1e6 / exp(gdpA)]                         # A's trade dependence on B (GDP in 2011 USD: deviation)
d[, syssize := .N / 1, by = year]; d[, syssize := (1 + sqrt(1 + 4 * syssize)) / 2]   # n states: N directed = n(n-1)

# ---- 2. Naval power (Crisher & Souva 2014): ln tonnage, 0 for states not in data (paper's "all states" variant)
nv <- fread(file.path(ROOT, "work/47/naval.tab"))[, .(ccode = cowcode, year, totton)]
nv[, lnav := log1p(pmax(totton, 0))]
d <- merge(d, nv[, .(ccode1 = ccode, year, navA = lnav, tonA = totton)], by = c("ccode1", "year"), all.x = TRUE)
d <- merge(d, nv[, .(ccode2 = ccode, year, navB = lnav, tonB = totton)], by = c("ccode2", "year"), all.x = TRUE)
d[is.na(navA), `:=`(navA = 0, tonA = 0)][is.na(navB), `:=`(navB = 0, tonB = 0)]
# allies' naval: ln(1 + summed tonnage of each side's COW defense/entente allies, excluding the other dyad member)
al <- unique(as.data.table(cow_alliance)[(cow_defense == 1 | cow_entente == 1) & year %in% 1949:2011, .(a = ccode1, b = ccode2, year)])
al <- unique(rbind(al, al[, .(a = b, b = a, year)]))                 # both directions
al <- merge(al, nv[, .(b = ccode, year, totton)], by = c("b", "year"), all.x = TRUE)[is.na(totton), totton := 0]
sumal <- al[, .(allsum = sum(totton)), by = .(ccode = a, year)]
d <- merge(d, sumal[, .(ccode1 = ccode, year, sumA = allsum)], by = c("ccode1","year"), all.x = TRUE)
d <- merge(d, sumal[, .(ccode2 = ccode, year, sumB = allsum)], by = c("ccode2","year"), all.x = TRUE)
d <- merge(d, al[, .(ccode1 = a, ccode2 = b, year, isally = 1L)], by = c("ccode1","ccode2","year"), all.x = TRUE)
d[is.na(sumA), sumA := 0][is.na(sumB), sumB := 0][is.na(isally), isally := 0L]
d[, allyA := pmax(sumA - isally * tonB, 0)][, allyB := pmax(sumB - isally * tonA, 0)]   # exclude the dyad partner
d[, `:=`(allynavA = log1p(allyA), allynavB = log1p(allyB))]

# ---- 3. DV: fatal MID initiation. VERIFY FIX: paper p.350 cites Palmer et al. (2015) = COW MID 4 -> use COW MID
#      (peacesciencer cow_mid_dirdisps) instead of GML. A takes first militarized action (orig1 & sidea1) in a new MID
#      (disponset) with >=1 battle death (dispute-level fatality from cow_mid_disps, >=1; -9 treated as not fatal).
cm <- as.data.table(cow_mid_dirdisps)[disponset == 1, .(dispnum, ccode1, ccode2, year, sidea1, orig1)]
cm <- merge(cm, as.data.table(cow_mid_disps)[, .(dispnum, fatality)], by = "dispnum", all.x = TRUE)
stopifnot(!anyNA(cm$fatality))
g <- cm[, .(finit = as.integer(any(sidea1 == 1 & orig1 == 1 & fatality >= 1)),
            fatal_nd = as.integer(any(fatality >= 1))), by = .(ccode1, ccode2, year)]
d <- merge(d, g, by = c("ccode1","ccode2","year"), all.x = TRUE)[is.na(finit), `:=`(finit = 0L, fatal_nd = 0L)]
# nondirected fatal MID onset (either direction) for peace years
d[, c1 := pmin(ccode1, ccode2)][, c2 := pmax(ccode1, ccode2)]
d[, fnd := max(fatal_nd), by = .(c1, c2, year)]
setorder(d, ccode1, ccode2, year)
d[, py := { out <- integer(.N); cnt <- 0L; for (i in seq_len(.N)) { out[i] <- cnt; cnt <- if (fnd[i] == 1L) 0L else cnt + 1L }; out },
  by = .(ccode1, ccode2)]                                          # years since last fatal MID (spell starts 1949)
kn <- quantile(d$py, c(.25, .5, .75), type = 2); Tm <- max(d$py)
stopifnot(length(unique(kn)) == 3, Tm > max(kn)); cat("btscs knots:", kn, " T:", Tm, "\n")
d[, `:=`(dv = shift(finit, -1), py1 = shift(py, -1), yr_next = shift(year, -1)), by = .(ccode1, ccode2)]
d <- d[!is.na(dv) & yr_next == year + 1 & year >= 1950 & year <= 2010]
# BKT/btscs cubic spline basis (ESTAUDIT change)
sp <- function(t, kk) pmax(t - kk, 0)^3 - (Tm - kk) / Tm * t^3
d[, `:=`(s1 = sp(py1, kn[1]), s2 = sp(py1, kn[2]), s3 = sp(py1, kn[3]))]
d[, ddyad := ccode1 * 1000 + ccode2]
d[, dvyear := year + 1]
cat("estimation frame:", nrow(d), " DV events:", sum(d$dv), "\n")

# ---- 4. Model 2 rebuild
X <- "third * navA + third:navB + navB + allynavA + third:allynavA + allynavB + third:allynavB + tdepAB + democ + ldist + contig_ord + alliance + capratio + gdpA + gdppcA + syssize + py1 + s1 + s2 + s3"
f0 <- as.formula(paste("dv ~", X))
keep <- c(all.vars(f0), "ccode1", "ccode2", "year", "dvyear", "ddyad")
d <- d[, ..keep]; gc()
tie <- "third:navA"
r_pub <- fit_tie(f0, as.data.frame(d), binomial(), tie, "ddyad")
b_pub <- 0.014
match <- rebuild_verdict("47", "third:navA", r_pub$b, r_pub$p, b_pub, NA, 0.006)   # round-2 rule (was: sign and p < .05)
cat("rebuild", tie, r_pub$b, "p", r_pub$p, "N", r_pub$n, "| pub 0.014 (0.006) N 995,488 ->", match, "\n")

# ---- 5. History on DV year (t+1)
d <- as.data.table(add_history(as.data.frame(d), "ccode1", "ccode2", "dvyear"))
s <- d[dvyear >= 1822 & dvyear <= 2011]
s <- as.data.frame(s[complete.cases(s[, c(all.vars(f0), "hist_dyad_w15","hist_mono_w15","hist_dyad_w26","hist_mono_w26","bm_cinc_w15","bm_democ_w15"), with = FALSE])])
cat("dropped for history/benchmarks:", r_pub$n - nrow(s), "(relative to rebuild N)\n")
m0  <- fit_tie(f0, s, binomial(), tie, "ddyad")
mD  <- fit_tie(update(f0, . ~ . + hist_dyad_w15), s, binomial(), tie, "ddyad")
mDM <- fit_tie(update(f0, . ~ . + hist_dyad_w15 + hist_mono_w15), s, binomial(), tie, "ddyad")
m26 <- fit_tie(update(f0, . ~ . + hist_dyad_w26 + hist_mono_w26), s, binomial(), tie, "ddyad")
bC  <- fit_tie(update(f0, . ~ . + bm_cinc_w15), s, binomial(), tie, "ddyad")
bD  <- fit_tie(update(f0, . ~ . + bm_democ_w15), s, binomial(), tie, "ddyad")

R <- list(m0 = m0, mD = mD, mDM = mDM, m26 = m26, bC = bC, bD = bD)
cat("rebuild (COW MID) third:navA", r_pub$b, "SE", r_pub$se, "p", r_pub$p, "N", r_pub$n, "\n")
# ---- sensitivity hook (split only): runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_split(id = "47", dat = s,
             fit = function(dat, add) { r <- fit_tie(if (length(add)) update(f0, as.formula(paste(". ~ . +", paste(add, collapse = " + ")))) else f0, dat, binomial(), tie, "ddyad")
               c(b = r$b, se = r$se, p = r$p) })
}
source(file.path(CODE, "ea4/est_helpers.R"))
f_old <- file.path(OUT, "results_verify/47__third_navA.csv")  # port: prior verify-phase output (not ported); feeds notes only
old <- if (file.exists(f_old)) read.csv(f_old) else data.frame(b0 = NA, p0 = NA, b_dm = NA, p_dm = NA)
st <- "logit (fixest feglm), SE clustered on directed dyad, py + 3 BKT/btscs splines (t+1)"
write_est(47, "Feldman, Eiran & Rubin 2021", "third_navA",
  "logit, robust SE clustered on directed dyad, peace years since last fatal MID + BKT cubic spline (3 interior knots), DV/py/splines at t+1",
  "pdf Research Design section + Table 1 note ('robust standard errors clustered on directed dyads')",
  rep(st, 6), "N",
  "all six steps used py + truncated-power terms (py-k)_+^3 at knots 1,4,7, not the BKT (btscs) natural-spline basis the article cites; logit + directed-dyad clustering already matched",
  "Y", lapply(R, function(z) c(b = z$b, se = z$se, p = z$p)), m0$n,
  sprintf(paste0("REBUILD (COW MID DV as in verify). btscs knots %s, T=%d (article's knot positions unstated; btscs default rule). ",
    "Rebuild full-frame b %.4f SE %.4f p %.3f N %d (pub 0.014, SE .006, N 995,488). Before fix: b0 %.5f p0 %.3f, b_dm %.5f p_dm %.3f. Key = product term third:navA."),
    paste(kn, collapse = "/"), Tm, r_pub$b, r_pub$se, r_pub$p, r_pub$n, old$b0, old$p0, old$b_dm, old$p_dm))
sink()
