# P9_073 Mansfield & Snyder 2002, ISQ 46(4) "Incomplete Democratization and the Outbreak of Military Disputes" (10.1111/1468-2478.00244).
#   DEMOCRACY. REBUILD (protocol B; no archive). Variable construction for the Oneal-Russett controls copied from P9_028 (same base).
# Article (staged PDF pp. 534-537, Table 1): logit, robust SEs, natural spline (3 knots), Oneal-Russett (1997) politically relevant dyads
#   1950-85; DV = a military dispute breaks out in year t. Regime type from the Polity composite REG (democ - autoc): democracy REG > 6,
#   autocracy REG < -6, anocracy otherwise, coded at t-1 and t-6. Incomplete Demtransition = 1 if i or j moves autocracy -> anocracy between
#   t-6 and t-1; Complete Demtransition = autocracy/anocracy -> democracy; Complete Auttransition = democracy/anocracy -> autocracy;
#   Incomplete Auttransition = democracy -> anocracy. Controls (taken from Oneal-Russett 1997): DEM_L, DEM_H, GROWTH_L, ALLIES, CONTIG,
#   CAPRATIO (raw), DEPEND_S (lower dependence), dDEPEND_H. One-tailed tests for the four transition terms.
# HEADLINE (fixed before estimating): Table 1 base model, Incomplete Demtransition_ij 0.449 (0.177) ** (one-tailed), N 18,055.
# Rebuild inputs: Russett-Oneal TRIANGLE.tab (as P9_028: dispute1, demaut, depend lags, allies incl. both-US, contiguity, capratio, SDP
#   growth) + Polity IV 2014 polity2 (stand-in for Polity III REG) at t-1 and t-6 per state (work_p8/P8_049). Splines: btscs on dispute1.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R")); source(file.path(CODE, "work_p9/rebuild_common_p9r.R"))
set.seed(123)
W <- file.path(OUT, "work_p9/P9_073"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_073.txt"), split = TRUE)

# ---- 1. TRIANGLE, all years (lags need pre-1950 rows) ----
t <- read.delim(file.path(ROOT, "recA/or/TRIANGLE.tab"))
stopifnot(nrow(t) == 39996, all(t$statea < t$stateb), anyDuplicated(t[, c("statea", "stateb", "year")]) == 0,
          !anyNA(t[, c("dispute1", "smigoabi", "dependa", "dependb", "demauta", "demautb", "allies", "lcaprat2", "noncontg")]))
t <- t %>% mutate(dyadid = 1000 * statea + stateb)

# ---- 2. lags of IGO and dependence (same dyad, years t-1 .. t-4) ----
lagk <- function(k) t %>% transmute(statea, stateb, year = year + k, igo = smigoabi, da = dependa, db = dependb) %>%
  rename_with(~paste0(.x, "_l", k), c("igo", "da", "db"))
n0 <- nrow(t)
for (k in 1:4) { t <- t %>% left_join(lagk(k), by = c("statea", "stateb", "year")); stopifnot(nrow(t) == n0) }
dtr <- function(l1, l2, l3, l4) ifelse(!is.na(l4), l1 - l4, ifelse(!is.na(l3), l1 - l3, l1 - l2))   # t-1 minus t-4, fallback t-3, t-2
t <- t %>% mutate(igo_l1 = igo_l1, depl_l1 = pmin(da_l1, db_l1),
  dda = dtr(da_l1, da_l2, da_l3, da_l4), ddb = dtr(db_l1, db_l2, db_l3, db_l4),
  ddeph = ifelse(is.na(dda) | is.na(ddb), NA_real_, ifelse(abs(dda) >= abs(ddb), dda, ddb)),
  deml = pmin(demauta, demautb), demh = pmax(demauta, demautb), caprat = exp(lcaprat2), contig = -noncontg)
stopifnot(all(t$contig %in% c(0, 1)))

# ---- 3. ALLIES: direct (TRIANGLE) or both allied with the US (US dyads are all in the PR set) ----
us <- t %>% filter(statea == 2 | stateb == 2, allies == 1) %>% transmute(o = ifelse(statea == 2, stateb, statea), year, usally = 1L) %>% distinct()
t <- t %>% left_join(us %>% rename(usa = usally), by = c("statea" = "o", "year")) %>% left_join(us %>% rename(usb = usally), by = c("stateb" = "o", "year"))
stopifnot(nrow(t) == n0)
t <- t %>% mutate(allies98 = as.integer(allies == 1 | (coalesce(usa, 0L) == 1 & coalesce(usb, 0L) == 1 & statea != 2)))
cat("allies: direct", sum(t$allies), " direct-or-both-US", sum(t$allies98), "\n")

# ---- 4. GROWTH_L: mean annual % change in SDP GDPpc over t-3..t-1 (fallback 1-2 yrs), weak link ----
sy <- tryCatch(create_stateyears(system = "cow", subset_years = 1940:1992) %>% add_sdp_gdp(), error = function(e) { cat("add_sdp_gdp failed:", conditionMessage(e), "\n"); NULL })
stopifnot(!is.null(sy), "wbgdppc2011est" %in% names(sy))
sy <- sy %>% transmute(ccode, year, gpc = exp(wbgdppc2011est))
sy <- sy %>% left_join(sy %>% transmute(ccode, year = year + 1, gpc_p = gpc), by = c("ccode", "year")) %>% mutate(gr = 100 * (gpc / gpc_p - 1))
g3 <- sy %>% select(ccode, year)
for (k in 1:3) g3 <- g3 %>% left_join(sy %>% transmute(ccode, year = year + k, !!paste0("g", k) := gr), by = c("ccode", "year"))
g3 <- g3 %>% mutate(gr3 = ifelse(is.na(g1) & is.na(g2) & is.na(g3), NA_real_, rowMeans(cbind(g1, g2, g3), na.rm = TRUE))) %>% select(ccode, year, gr3)
stopifnot(anyDuplicated(g3[, c("ccode", "year")]) == 0)
t <- t %>% left_join(g3 %>% rename(ga = gr3), by = c("statea" = "ccode", "year")) %>% left_join(g3 %>% rename(gb = gr3), by = c("stateb" = "ccode", "year"))
stopifnot(nrow(t) == n0)
t <- t %>% mutate(growl = pmin(ga, gb)) %>% as.data.frame()

# ---- 5. regime transitions t-6 -> t-1 (Polity IV 2014 polity2 as REG) ----
p <- polity_cow("polity2") %>% mutate(reg = case_when(polity2 > 6 ~ "D", polity2 < -6 ~ "A", !is.na(polity2) ~ "N"))
tr <- p %>% select(ccode, year) %>%
  left_join(p %>% transmute(ccode, year = year + 1, r1 = reg), by = c("ccode", "year")) %>%
  left_join(p %>% transmute(ccode, year = year + 6, r6 = reg), by = c("ccode", "year")) %>%
  mutate(cdem = as.integer(r6 %in% c("A", "N") & r1 == "D"), idem = as.integer(r6 == "A" & r1 == "N"),
         caut = as.integer(r6 %in% c("D", "N") & r1 == "A"), iaut = as.integer(r6 == "D" & r1 == "N")) %>%
  mutate(across(c(cdem, idem, caut, iaut), ~ ifelse(is.na(r1) | is.na(r6), NA_integer_, .x))) %>% select(ccode, year, cdem, idem, caut, iaut)
t <- t %>% left_join(tr %>% rename_with(~paste0(.x, "_a"), c(cdem, idem, caut, iaut)), by = c("statea" = "ccode", "year")) %>%
  left_join(tr %>% rename_with(~paste0(.x, "_b"), c(cdem, idem, caut, iaut)), by = c("stateb" = "ccode", "year"))
stopifnot(nrow(t) == n0)
t <- t %>% mutate(cdem = pmax(cdem_a, cdem_b), idem = pmax(idem_a, idem_b), caut = pmax(caut_a, caut_b), iaut = pmax(iaut_a, iaut_b))
t <- py_calc(t, "dyadid", "dispute1")

# ---- 6. sample 1950-85 ----
d <- t[t$year >= 1950 & t$year <= 1985, ]
d <- btscs_splines(d)
rhs0 <- c("cdem", "idem", "caut", "iaut", "deml", "demh", "growl", "allies98", "contig", "caprat", "depl_l1", "ddeph", "spl1", "spl2", "spl3")
vars <- c("dispute1", rhs0, "dyadid")
cat("1950-85 rows:", nrow(d), "\nmissing by var:\n"); print(sapply(d[, vars], function(x) sum(is.na(x))))
d0 <- d[complete.cases(d[, vars]), ]
cat("complete-case N:", nrow(d0), "(published 18,055); incomplete demtransition dyad-years", sum(d0$idem), "\n")
TERM <- "idem"
fitter <- function(r, data) fit_bin("dispute1", r, data, "logit", "hetero", terms = TERM)
full <- fitter(rhs0, d0)
cat("\nREBUILT vs published: idem b", full[[TERM]][["b"]], "(pub 0.449) se", full[[TERM]][["se"]], "(pub 0.177) N", attr(full, "n"), "(pub 18,055)\n")
ds <- common_sample(d, vars, "statea", "stateb", "year")
M <- run_six(rhs0, ds, fitter)
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "P9_073", dat = ds, y = "dispute1", tie = TERM, ctrl = setdiff(rhs0, TERM), fe = NULL,
           fit = function(dat, add) fitter(c(rhs0, add), dat)[[TERM]])
}
saveRDS(list(M = M, full = full), file.path(W, "M_P9_073.rds"))
meta <- list(pid = "P9_073", cite = "Mansfield & Snyder 2002", journal = "ISQ", year = 2002, doi = "10.1111/1468-2478.00244", family = "DEMOCRACY")
nt <- sprintf(paste0("REBUILT (no archive): TRIANGLE PR dyads 1950-85 (dispute1; O/R 1997 controls as built for P9_028: deml, demh, SDP growth_L, allies incl. both-US, ",
  "contiguity, raw capratio, lower dependence t-1, dDEPEND_H); regime categories from Polity IV 2014 polity2 (paper Polity III REG) at t-1 and t-6; ",
  "btscs splines (no peace-years term); robust (HC) SE. Full-sample rebuilt b=%.6f se=%.6f N=%d vs published 0.449 (0.177) N 18,055. Common sample N=%d."),
  full[[TERM]][["b"]], full[[TERM]][["se"]], attr(full, "n"), nrow(ds))
write_p8(meta, TERM, "abstract: dyads with a state in an incomplete democratic transition are particularly prone to disputes; Table 1 base model",
  "Table 1 base model (composite index; logit, MID outbreak, PR dyads 1950-85)", "logit, robust SE, natural spline (3 knots), one-tailed tests for transitions",
  "feglm logit, HC vcov, all six steps", 0.449, 0.177, 18055, TRUE, TRUE, M, "Y (splines)", nt)
dir.create(file.path(OUT, "results_p9"), recursive = TRUE, showWarnings = FALSE)
src <- file.path(OUT, "results_p8", "P9_073__idem.csv"); stopifnot(file.exists(src), file.rename(src, file.path(OUT, "results_p9", basename(src))))
sink()
