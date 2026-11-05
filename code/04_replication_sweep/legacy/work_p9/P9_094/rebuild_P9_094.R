# P9_094 Krause 2004, International Interactions 30(4) "Hazardous Weapons? Effects of Arms Transfers and Defense Pacts on Militarized
#   Disputes, 1950-1995" (10.1080/03050620490884038). TIE (arms transfers from major powers). REBUILD (protocol B; no archive).
# Article (staged PDF pp. 356-362, Tables 1-2): COW state-years 1950-1995. Table 2 Model 2 = 2SCML second stage (probit) of being a MID
#   initiator (originator on side A of a new MID; ongoing initiator state-years without a new initiation omitted) on arms transfers from
#   major powers (natural log of SIPRI TIV, constant 1990 US$ m), defense pacts with major powers (COW alliance), proportion of major-power
#   transfers from defense-pact allies x defense-pact dummy, national capabilities (CINC), territorial disputes (annual Huth-Allee count),
#   democracy (Polity IV), other MID involvement (any MID other than new initiations), the first-stage residual, and a BTSCS counter + three
#   splines. First stage (Table 1, OLS): arms ~ MID involvement + defense pacts + capabilities + territorial disputes + democracy + lagged
#   arms. SEs robust, clustered by state (both stages; no further correction in the article).
#   Major powers: US, UK, France, USSR/Russia, China 1950-95; Germany and Japan 1991-95.
# HEADLINE (fixed before estimating): Table 2 Model 2, Arms Transfers 0.1241 (0.0153) *, N 5,290.
#   Note: "Significance levels are based on one-tailed tests for variables that rest on specific expectations" (* p<0.05) -> alpha_pub 0.10.
#   Verdict rule (PI, round 2): published sign, significance at alpha_pub, b0 within 0.5x-2x.
# Rebuild inputs and DISCLOSED DEPARTURES:
#   - SIPRI TIV imports by recipient-year from USA, USSR, UK, France and China, delivered 1950-1995 (work_p9/P9_094/
#     SIPRI_TIV_imports_from_5majors_1950_1995.csv; one section per supplier, retrieved 2026-10-11). SIPRI "0" (0-0.5) coded 0.25, blank 0.
#     arms = log(1 + TIV sum). Not in the pull: Russia as a supplier after 1991 (the USSR row is blank from 1992) and Germany/Japan 1991-95.
#     The current SIPRI TIV scale replaces the 2002 constant-1990-dollar figures.
#   - MIDs: GML participant data ({peacesciencer} gml_part) instead of MID 3.02; initiator = side A originator; COW v4.1 defense pacts;
#     NMC CINC and Polity via {peacesciencer}; Huth-Allee territorial dispute-years from work_p8/P8_017/statusquo.dta (datetm = Stata %tm).
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R")); source(file.path(CODE, "work_p9/rebuild_common_p9r.R"))
set.seed(123)
W <- file.path(OUT, "work_p9/P9_094"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_094.txt"), split = TRUE)
YRS <- 1950:1995
# ---- 1. SIPRI file: five supplier sections ----
L <- readLines(file.path(ROOT, "work_p9/P9_094/SIPRI_TIV_imports_from_5majors_1950_1995.csv"), encoding = "UTF-8")
L <- sub("^﻿", "", L)
tot <- grep("^Total imports from ", L); hdr <- grep("^Recipient,", L)
stopifnot(length(hdr) == 5, length(tot) == 5, all(hdr < tot))
sup_names <- sub("^Total imports from ([^,]+),.*$", "\\1", L[tot]); print(sup_names)
stopifnot(identical(sup_names, c("United States", "Soviet Union", "United Kingdom", "France", "China")))
sup_cc <- c("United States" = 2, "Soviet Union" = 365, "United Kingdom" = 200, "France" = 220, "China" = 710)
tiv <- list()
for (s in 1:5) {
  body <- L[(hdr[s] + 1):(tot[s] - 1)]; body <- body[nzchar(body)]
  x <- read.csv(text = c(L[hdr[s]], body), check.names = FALSE, colClasses = "character")
  stopifnot(!anyDuplicated(x$Recipient))                                     # no repeated recipient within a supplier section
  yc <- as.character(YRS); stopifnot(all(yc %in% names(x)))
  v <- as.matrix(x[, yc]); v <- trimws(v)
  num <- ifelse(v == "", 0, ifelse(v == "0", 0.25, suppressWarnings(as.numeric(v)))); stopifnot(!anyNA(num))
  # row sums check against the file's 1950-1995 total column (rounding tolerance)
  tt <- suppressWarnings(as.numeric(x[["1950-1995"]])); ok <- is.finite(tt)
  stopifnot(all(abs(rowSums(num)[ok] - tt[ok]) <= 0.5 * ncol(num) + 1))
  tiv[[s]] <- data.frame(recipient = rep(x$Recipient, times = length(YRS)), year = rep(YRS, each = nrow(x)), supplier = sup_cc[[sup_names[s]]],
                         tiv = as.vector(num))
  cat(sup_names[s], ": recipients", nrow(x), " section total", sum(num), "\n")
}
tiv <- do.call(rbind, tiv); tiv <- tiv[tiv$tiv > 0, ]
# ---- 2. recipient names -> COW codes (year-dependent where needed) ----
manual <- c("Cabo Verde" = 402, "Cote d'Ivoire" = 437, "Congo" = 484, "DR Congo" = 490, "East Germany (GDR)" = 265, "North Korea" = 731,
            "South Korea" = 732, "North Yemen" = 678, "Yemen Arab Republic (North Yemen)" = 678, "South Yemen" = 680, "Yemen" = 679,
            "South Vietnam" = 817, "Viet Nam" = 816, "Turkiye" = 640, "Myanmar" = 775, "Saint Vincent" = 57, "Saint Kitts and Nevis" = 60,
            "Bosnia-Herzegovina" = 346, "Czechoslovakia" = 315, "Yugoslavia" = 345, "Taiwan" = 713, "Soviet Union" = 365)
nm <- unique(tiv$recipient)
cc <- suppressWarnings(countrycode::countrycode(nm, "country.name", "cown"))
cc[nm %in% names(manual)] <- manual[nm[nm %in% names(manual)]]
nonstate <- grepl("\\*|unknown|Biafra|Katanga|Western Sahara", nm)
cc[nonstate] <- NA
map <- data.frame(recipient = nm, ccode = cc)
unm <- map$recipient[is.na(map$ccode)]
cat("UNMAPPED recipient names (non-state or no COW state):", paste(unm, collapse = " | "), "\n")
stopifnot(all(nonstate[match(unm, nm)]))                                      # every unmapped name is a declared non-state recipient
tiv <- tiv %>% left_join(map, by = "recipient") %>% filter(!is.na(ccode))
tiv$ccode[tiv$recipient == "Germany"] <- ifelse(tiv$year[tiv$recipient == "Germany"] <= 1990, 260, 255)
# ---- 3. COW state-years, capabilities, democracy ----
sy <- suppressMessages(create_stateyears(system = "cow", subset_years = YRS) %>% add_nmc()) %>% as.data.frame() %>% dplyr::select(ccode, year, cinc)
p <- polity_cow("polity2"); sy <- sy %>% left_join(p, by = c("ccode", "year"))
miss <- tiv %>% anti_join(sy, by = c("ccode", "year")) %>% distinct(recipient, ccode, year)
cat("TIV rows on mapped names outside the COW state-year frame:", nrow(miss), "\n"); if (nrow(miss)) print(head(miss, 30))
# ---- 4. defense pacts with major powers (COW v4.1), proportion from pact allies ----
mp_cc <- function(y) c(2, 200, 220, 365, 710, if (y >= 1991) c(255, 740))
al <- read.csv(file.path(ROOT, "work_p8/P8_049/alliance_v4_1_by_dyad_yearly.csv")) %>% filter(defense == 1, year %in% YRS) %>%
  transmute(a = ccode1, b = ccode2, year) %>% distinct()
al2 <- rbind(al, transmute(al, a = b, b = a, year)) %>% distinct()
defmp <- al2 %>% rowwise() %>% filter(b %in% mp_cc(year)) %>% ungroup() %>% distinct(ccode = a, supplier = b, year)
def_any <- defmp %>% distinct(ccode, year) %>% mutate(defpact = 1L)
agg <- tiv %>% left_join(mutate(defmp, ally = 1L), by = c("ccode", "supplier", "year")) %>% mutate(ally = coalesce(ally, 0L)) %>%
  group_by(ccode, year) %>% summarise(tiv_mp = sum(tiv), tiv_ally = sum(tiv * ally), .groups = "drop")
sy <- sy %>% left_join(agg, by = c("ccode", "year")) %>% left_join(def_any, by = c("ccode", "year")) %>%
  mutate(tiv_mp = coalesce(tiv_mp, 0), tiv_ally = coalesce(tiv_ally, 0), defpact = coalesce(defpact, 0L),
         arms = log1p(tiv_mp), propdef = ifelse(tiv_mp > 0, tiv_ally / tiv_mp, 0) * defpact)
# ---- 5. Huth-Allee territorial disputes (annual count per state) ----
ha <- as.data.frame(haven::zap_labels(haven::read_dta(file.path(ROOT, "work_p8/P8_017/statusquo.dta"), col_select = c("challenger", "target", "datetm", "dispno"))))
stopifnot(nrow(ha) == 6542); ha$year <- 1960 + floor(ha$datetm / 12); stopifnot(min(ha$year) == 1919, max(ha$year) <= 1995)
ht <- rbind(transmute(ha, ccode = challenger, year, dispno), transmute(ha, ccode = target, year, dispno)) %>% distinct() %>% count(ccode, year, name = "terrdisp")
sy <- sy %>% left_join(ht, by = c("ccode", "year")) %>% mutate(terrdisp = coalesce(terrdisp, 0L))
# ---- 6. MID initiation, involvement (GML participants) ----
gp <- as.data.frame(gml_part) %>% transmute(dispnum, ccode, styear, endyear, sidea, orig)
ini_new <- gp %>% filter(sidea == 1, orig == 1) %>% distinct(ccode, year = styear, dispnum)
inv <- gp %>% rowwise() %>% mutate(year = list(styear:endyear)) %>% tidyr::unnest(year) %>% ungroup() %>% distinct(ccode, year, dispnum, sidea, orig, styear)
I <- sy %>% dplyr::select(ccode, year) %>%
  left_join(ini_new %>% count(ccode, year, name = "n_new_ini"), by = c("ccode", "year")) %>%
  left_join(inv %>% filter(sidea == 1, orig == 1, styear < year) %>% distinct(ccode, year) %>% mutate(ongo_ini = 1L), by = c("ccode", "year")) %>%
  left_join(inv %>% distinct(ccode, year) %>% mutate(mid_inv = 1L), by = c("ccode", "year"))
oth <- inv %>% anti_join(ini_new, by = c("ccode", "year", "dispnum")) %>% distinct(ccode, year) %>% mutate(other_mid = 1L)
sy <- sy %>% left_join(I, by = c("ccode", "year")) %>% left_join(oth, by = c("ccode", "year")) %>%
  mutate(initiator = as.integer(coalesce(n_new_ini, 0L) > 0), ongo_ini = coalesce(ongo_ini, 0L), mid_inv = coalesce(mid_inv, 0L), other_mid = coalesce(other_mid, 0L))
cat("state-years:", nrow(sy), " initiator years:", sum(sy$initiator), " ongoing-initiator years without new initiation (dropped from stage 2):",
    sum(sy$ongo_ini == 1 & sy$initiator == 0), "\n")
sy <- sy %>% arrange(ccode, year) %>% group_by(ccode) %>% mutate(arms_l1 = ifelse(dplyr::lag(year) == year - 1, dplyr::lag(arms), NA)) %>% ungroup() %>% as.data.frame()
# ---- 7. 2SCML ----
s1 <- c("mid_inv", "defpact", "cinc", "terrdisp", "polity2", "arms_l1")
fs <- sy[complete.cases(sy[, c("arms", s1)]), ]
m1 <- feols(as.formula(paste("arms ~", paste(s1, collapse = " + "))), data = fs, vcov = ~ccode)
cat("\n---- first stage (Table 1; published N 5,542: mid_inv 0.2683, defpact 0.2009, terrdisp 0.0372, lagged arms 0.8178)\n"); print(round(coeftable(m1), 4)); cat("N", nobs(m1), "\n")
fs$res <- resid(m1)
d <- sy %>% left_join(fs[, c("ccode", "year", "res")], by = c("ccode", "year"))
d <- d[!(d$ongo_ini == 1 & d$initiator == 0), ]
d <- py_calc(d, "ccode", "initiator"); d <- btscs_splines(d)
rhs0 <- c("arms", "defpact", "propdef", "cinc", "terrdisp", "polity2", "other_mid", "res", "py", "spl1", "spl2", "spl3")
vars <- c("initiator", rhs0, "ccode")
d0 <- d[complete.cases(d[, vars]), ]
cat("complete-case N:", nrow(d0), "(published 5,290)  initiations", sum(d0$initiator), "\n")
TERM <- "arms"
fitter <- function(r, data) fit_bin("initiator", r, data, "probit", "~ccode", terms = TERM)
full <- fitter(rhs0, d0)
cat("\nREBUILT vs published: arms b", full[[TERM]][["b"]], "(pub 0.1241) se", full[[TERM]][["se"]], "(pub 0.0153) N", attr(full, "n"), "(pub 5,290)\n")
ds <- common_sample(d, vars, "ccode", NULL, "year", monadic = TRUE)
M <- run_six(rhs0, ds, fitter, monadic = TRUE)
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "P9_094", dat = ds, y = "initiator", tie = TERM, ctrl = setdiff(rhs0, TERM), fe = NULL,
           fit = function(dat, add) fitter(c(rhs0, add), dat)[[TERM]])
}
saveRDS(list(M = M, full = full), file.path(W, "M_P9_094.rds"))
meta <- list(pid = "P9_094", cite = "Krause 2004", journal = "II", year = 2004, doi = "10.1080/03050620490884038", family = "TIE")
nt <- sprintf(paste0("REBUILT (no archive): COW state-years 1950-95; SIPRI TIV imports from USA/USSR/UK/France/China (PI pull 2026-10-11; '0' = 0.25; ",
  "Russia after 1991 and Germany/Japan 1991-95 not in the pull), arms = log(1+TIV); COW v4.1 defense pacts with major powers; proportion from ",
  "pact allies x pact; CINC; Huth-Allee territorial dispute counts; polity2; GML initiation (side A originators), other MID involvement; 2SCML: OLS ",
  "first stage with lagged arms, probit second stage with residual, BTSCS + 3 splines, SEs clustered by state. Full-sample rebuilt b=%.6f se=%.6f ",
  "N=%d vs published 0.1241 (0.0153) N 5,290. Common sample N=%d."), full[[TERM]][["b"]], full[[TERM]][["se"]], attr(full, "n"), nrow(ds))
b0 <- M$m0[[TERM]][["b"]]; p0 <- M$m0[[TERM]][["p"]]; r_ <- b0 / 0.1241
v <- rebuild_verdict("P9_094", TERM, b0, p0, 0.1241, TRUE)   # round-2 rule, alpha_pub (0.10, one-tailed .05) from registry.csv
cat("round-2 verdict:", v, " ratio", round(r_, 3), " p0", signif(p0, 3), "\n")
write_p8(meta, TERM, "abstract and Table 2 Model 2: increased arms transfers from major powers make states more likely to initiate MIDs",
  "Table 2 Model 2 (2SCML second stage, probit, MID initiator, state-years 1950-95)", "2SCML (OLS first stage, probit second stage), robust SE clustered by state, BTSCS",
  "hand-coded 2SCML, probit stage 2 clustered by state, all six steps", 0.1241, 0.0153, 5290, TRUE, TRUE, M, "Y (BTSCS counter + splines)", nt, verdict = v)
dir.create(file.path(OUT, "results_p9"), recursive = TRUE, showWarnings = FALSE)
src <- file.path(OUT, "results_p8", "P9_094__arms.csv"); stopifnot(file.exists(src), file.rename(src, file.path(OUT, "results_p9", basename(src))))
sink()
