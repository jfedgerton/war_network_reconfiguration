# ESTAUDIT X059: Maoz, Terris, Kuperman & Talmud 2007 JOP, Table 6 MID column. Copy of rerun_verify_X059.R (same rebuild, spec,
# common sample) with the ESTIMATOR fixed. Article (paper.txt p.~20): "We also use the King and Zeng (2001) procedure for correcting
# for rare-event biases ... We use the prior correction method on the dependent variables, applying the actual proportions of 1's
# ... (MID ... 0.05%)"; Table 6 title "Rare Event Logit Analysis"; robust SE. All previous steps used ordinary logit.
# Now every step: relogit_kz() = ML logit + King-Zeng finite-sample bias correction, SE cluster(dyad) x n/(n+k);
# prior correction with tau = 0.0005 (paper's stated population share) applied to the intercept only (slopes unaffected; the
# sample is the PRD population, not a case-control draw, so it only matters for the constant). Caches -> est_cache/.
# VERIFY X059 Maoz, Terris, Kuperman & Talmud 2007 JOP, Table 6 MID column (rare-event logit, politically relevant dyads,
#   1816-2000, robust SE, N 69,371): Lagged Combined Imbalanced Relations 0.444 (.056)**.
# Copy of rebuild_X059.R with three changes, each justified by the paper text:
#   (1) DV/cumulative MIDs source = Maoz dyadic MID data (p.106-107: MID = "(ongoing) Militarized Interstate Dispute" from the
#       dyadic MID data) -> Maoz Dyadic MID v4.02 (COW site): MID = any dyadic-MID record in the dyad-year; cumulative MIDs = dyadic
#       MID onsets (durindx==1) 1816..t-1. Was GML.
#   (2) "All independent and control variables are lagged one year" (p.108) -> minreg302, distance, caprat301 now at t-1 (were t).
#   (3) "logit analysis, modified by the use of splines and years of peace (Beck, Katz, and Tucker 1998)" (p.108) -> btscs peace
#       years on the DV + 3 splines (quartile knots). Was years since last GML onset + squared + cubed.
# Unchanged: key = |combincons| at t-1; PRD = contig1 | statusa>0 | statusb>0; SE clustered on dyad; ordinary logit
#   (prior correction moves only the intercept). IGO similarity still omitted (no data).
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "vb2/vcommon.R"))
source(file.path(CODE, "work/relogit_kz_ea1.R")); source(file.path(CODE, "work/est_write_ea1.R"))
suppressMessages(library(peacesciencer))
WO <- file.path(OUT, "work_p4/X059"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work_p4/X059"); sink(file.path(WO, "log_est_X059.txt"), split = TRUE)

a <- haven::read_dta(file.path(W, "jopreplicat.dta"))
a <- as.data.frame(lapply(a, as.numeric)) %>% rename(ccode1 = statea, ccode2 = stateb)
stopifnot(all(a$ccode1 < a$ccode2), !anyDuplicated(a[, c("ccode1", "ccode2", "year")]))
L1 <- function(x, yr) ifelse(dplyr::lag(yr) == yr - 1, dplyr::lag(x), NA)
a <- a %>% arrange(ccode1, ccode2, year) %>% group_by(ccode1, ccode2) %>%
  mutate(lag_imbal = abs(L1(combincons, year)), L_minreg = L1(minreg302, year), L_dist = L1(distance, year),
         L_caprat = L1(caprat301, year)) %>% ungroup()

y <- read.csv(file.path(ROOT, "work_p4/X090/dymid/dyadic_mid_4.02.csv"))[, c("statea","stateb","year","durindx")]
y <- y[y$statea < y$stateb, ]
mo <- y %>% distinct(ccode1 = statea, ccode2 = stateb, year) %>% mutate(mid = 1)
on <- y %>% filter(durindx == 1) %>% count(ccode1 = statea, ccode2 = stateb, year, name = "nons")
g <- create_dyadyears(directed = FALSE, subset_years = 1816:2001) %>% select(ccode1, ccode2, year) %>% as.data.frame() %>%
  left_join(mo, by = c("ccode1", "ccode2", "year")) %>% left_join(on, by = c("ccode1", "ccode2", "year")) %>%
  mutate(mid = coalesce(mid, 0), nons = coalesce(nons, 0L)) %>% arrange(ccode1, ccode2, year) %>%
  group_by(ccode1, ccode2) %>% mutate(cummids = cumsum(dplyr::lag(nons, default = 0L))) %>% ungroup() %>% as.data.frame()
g$py <- btscs_py(g$ccode1 * 1000 + g$ccode2, g$year, g$mid)
d <- a %>% inner_join(g %>% select(ccode1, ccode2, year, mid, cummids, py), by = c("ccode1", "ccode2", "year")) %>%
  filter(year <= 2000, contig1 == 1 | statusa > 0 | statusb > 0) %>% mutate(dyad = ccode1 * 1000 + ccode2) %>% as.data.frame()
rm(a, g, y, mo, on); gc()
d <- btscs_splines(d, "py")
cat("PRD dyad-years 1816-2000:", nrow(d), " MID-years:", sum(d$mid), "\n")

f0 <- mid ~ L_minreg + L_dist + L_caprat + cummids + lag_imbal + py + spl1 + spl2 + spl3
terms <- c("lag_imbal", "lag_imbal__unc", "L_minreg")
FF <- function(f, data) { ct <- relogit_kz(f, data, cluster = "dyad", tau = 0.0005)
  cat("prior-corrected intercept:", attr(ct, "prior_intercept"), " KZ intercept:", ct["(Intercept)", "b"], "\n"); ct }
CD <- file.path(WO, "est_cache"); dir.create(CD, recursive = TRUE, showWarnings = FALSE)
rep <- vfit(f0, d[complete.cases(d[, c(all.vars(f0), "dyad")]), ], terms, FF, file.path(CD, "e_repro.rds"))
cat("repro N:", rep$n, "(published 69,371; verify repro 108,396)\n")
ds <- vcommon_sample(d, f0, "dyad", "ccode1", "ccode2", "year"); rm(d); gc()
stopifnot(nrow(ds) == 107459)                      # = results_verify n_used
C <- function(k) file.path(CD, paste0("e_", k, ".rds"))
M <- list(m0  = vfit(f0, ds, terms, FF, C("m0")),
          mD  = vfit(update(f0, . ~ . + hist_dyad_w15), ds, terms, FF, C("mD")),
          mDM = vfit(update(f0, . ~ . + hist_dyad_w15 + hist_mono_w15), ds, terms, FF, C("mDM")),
          m26 = vfit(update(f0, . ~ . + hist_dyad_w26 + hist_mono_w26), ds, terms, FF, C("m26")),
          bC  = vfit(update(f0, . ~ . + bm_cinc_w15), ds, terms, FF, C("bC")),
          bD  = vfit(update(f0, . ~ . + bm_democ_w15), ds, terms, FF, C("bD")))
stopifnot(all(sapply(M, function(m) m$n) == M$m0$n))

# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "X059", dat = ds, y = "mid", tie = "lag_imbal",
           ctrl = c("L_minreg", "L_dist", "L_caprat", "cummids", "py", "spl1", "spl2", "spl3"), fe = NULL,
           fit = function(dat, add) { ct <- FF(if (length(add)) update(f0, as.formula(paste(". ~ . +", paste(add, collapse = " + ")))) else f0, dat)
             c(b = unname(ct["lag_imbal", "b"]), se = unname(ct["lag_imbal", "se"]), p = unname(ct["lag_imbal", "p"])) })
}

RK <- lapply(M, function(m) m$k$lag_imbal); RU <- lapply(M, function(m) m$k$lag_imbal__unc)
print(rbind(KZ = sapply(RK, `[`, "b"), logit = sapply(RU, `[`, "b")))
pc <- function(a, b) 100 * (a - b) / abs(b)
write_est_ea1("X059", "Maoz et al. 2007", "lag_imbal",
  "rare-events logit (King-Zeng; prior correction with population share 0.05%), robust SE, years of peace + splines; PRD 1816-2000",
  "pdf (paper.txt: 'King and Zeng (2001) procedure for correcting for rare-event biases ... prior correction'; Table 6 'Rare Event Logit Analysis', robust SE)",
  "relogit: ML logit + King-Zeng bias correction (work/relogit_kz_ea1.R), prior correction tau=.0005 on intercept, SE cluster dyad x n/(n+k)",
  "N", "all six steps (rerun_verify_X059.R) used ordinary logit; King-Zeng bias/prior correction not applied (verify note: 'prior correction moves only the intercept')",
  "Y", RK, M$m0$n,
  sprintf(paste0("REBUILT_MATCH spec from rerun_verify_X059.R unchanged. KZ b0=%.6g (p=%.3g); uncorrected logit b0=%.6g (p=%.3g) = results_verify; ",
    "uncorrected mDM %.6g (pct_dm %.1f). Full-sample relogit repro %.6g (pub 0.444; logit repro 0.644). Prior correction only shifts the constant. ",
    "Imbalance at t-1 is partly defined by MIDs at t-1 (embeds history)."),
    RK$m0["b"], RK$m0["p"], RU$m0["b"], RU$m0["p"], RU$mDM["b"], pc(RU$mDM["b"], RU$m0["b"]), rep$k$lag_imbal["b"]))
sink()
