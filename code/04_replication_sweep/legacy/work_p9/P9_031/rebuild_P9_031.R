# P9_031 Mansfield & Pevehouse 2000 IO 54(4) "Trade Blocs, Trade Flows, and International Conflict" (10.1162/002081800551361).
#   REBUILD (protocol B): the authors' 2000 data URLs (Ohio State / Wisconsin personal pages) are not on disk; no author code.
# Article: logit of MID onset, O&R (1997) politically relevant dyads, MIDs 1950-85, Huber SEs clustered on dyad, natural spline
#   (3 knots) of years since a MID last began. Table 1 column 1 (binary logit, MID onset) is the headline equation (1).
# HEADLINE (fixed from the PDF before estimating): Table 1 col 1, key term = TRADE x PTA (product). Published -3.05e-7 (z -3.52,
#   so SE ~ 8.66e-8), p <= .01, N = 20,990. Abstract/intro claim: "hostilities between PTA members are less likely to occur as trade
#   flows rise between them".
# RHS: DEML DEMH GROWTHL ALLIES CONTIG CAPRATIO TRADE GDPL GDPH PTA TRADExPTA GDPLxPTA GDPHxPTA HEGEMONY + spline.
# Rebuild input (user instruction: use data already on hand): Hafner-Burton & Montgomery 2012 PTAMIDTables.dta (sweep study 86,
#   DATA/work/86). That file was built to replicate M&P 2000 Table 1 col 1 (HB&M do-file block "MP2000 Table 1 Column 1"), with
#   M&P-named variables: PTASameGT0 (any shared PTA, GATT-type institutions excluded), traBilMax (bilateral trade), rgdpMin/Max,
#   PTA_* products, polMin/Max, rgdppc_p3Min (3-yr growth, lower), caprat, allies, contigEQ1, gdpHegPerc (hegemon share of world GDP),
#   prdyad, dispon_l1 (MID onset at t+1) and dispon_spl0-3 (years since last onset + 3 spline terms). Rows = covariate year t,
#   DV = onset in t+1 (HB&M "effectively lagged by one year").
#   PTA coding is HB&M's own (WTO / McCall Smith / Schott sources); HB&M fn 13 report r ~ .84 with M&P's PTA on overlapping dyads.
# Choices made a priori: sample prdyad == 1 and DV years 1950-85 (rows 1949-84); spline = dispon_spl* (M&P: years since a MID
#   last BEGAN). For reference only, HB&M's own replication command (disp_spl*, year < 1986, no lower bound) is also fitted and
#   compared with their published Table 1 col 1 (PTA x TRADE -2.14E-04 (1.13E-04), N 25,209).
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R"))   # ROOT, OUT, dplyr, fixest, add_history, common_sample, run_six, write_p8, fit_bin; seed 123
suppressPackageStartupMessages(library(haven))
set.seed(123)
W <- file.path(OUT, "work_p9/P9_031"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_031.txt"), split = TRUE)

# ---- 1. load HB&M PTA/MID file (study 86 data, read-only) ----
keep <- c("year", "dyadid", "dispon_l1", "disp_l1", "prdyad", "PTASame", "PTASameGT0", "traBilMax", "rgdpMin", "rgdpMax",
          "PTA_traBilMax", "PTA_rgdpMin", "PTA_rgdpMax", "polMin", "polMax", "rgdppc_p3Min", "caprat", "allies", "contigEQ1",
          "gdpHegPerc", paste0("dispon_spl", 0:3), paste0("disp_spl", 0:3))
d <- as.data.frame(zap_labels(read_dta(file.path(ROOT, "work/86/PTAMIDTables.dta"), col_select = all_of(keep))))
stopifnot(nrow(d) == 535503, !anyDuplicated(d[, c("dyadid", "year")]))
d$c1 <- d$dyadid %/% 1000; d$c2 <- d$dyadid %% 1000
stopifnot(all(d$c1 < d$c2), all(d$prdyad %in% c(0, 1)), all(d$PTASameGT0 %in% c(0, 1)), all(d$dispon_l1 %in% c(0, 1, NA)))
# product terms in the file equal PTA x component (checked, not recomputed)
for (v in c("traBilMax", "rgdpMin", "rgdpMax"))
  stopifnot(max(abs(d[[paste0("PTA_", v)]] - d$PTASameGT0 * d[[v]]), na.rm = TRUE) < 1e-6 * max(abs(d[[v]]), na.rm = TRUE))
cat("PTASameGT0 == (PTASame > 0) share:", mean(d$PTASameGT0 == as.numeric(d$PTASame > 0)), "\n")

# ---- 2. reference: HB&M's own M&P replication command (Table 1 col 1 of study 86's paper) ----
ref_rhs <- c("PTASameGT0", "PTA_traBilMax", "PTA_rgdpMin", "PTA_rgdpMax", "polMin", "polMax", "rgdppc_p3Min", "traBilMax", "rgdpMin",
             "rgdpMax", "caprat", "allies", "gdpHegPerc", "contigEQ1", paste0("disp_spl", 0:3))
dr <- d[d$prdyad == 1 & d$year < 1986, ]; dr <- dr[complete.cases(dr[, c("dispon_l1", ref_rhs, "dyadid")]), ]
ref <- fit_bin("dispon_l1", ref_rhs, dr, "logit", "~dyadid", terms = "PTA_traBilMax")
cat("REFERENCE (study 86 do-file, P9_031 Table 1 col 1 command): PTA_traBilMax", ref$PTA_traBilMax[["b"]], "se", ref$PTA_traBilMax[["se"]], "N", attr(ref, "n"),
    "| study 86 published -2.14E-04 (1.13E-04) N 25209\n")   # port: ids only, no citation beside estimates

# ---- 3. rebuild sample: PR dyads, DV years 1950-85 ----
d <- d[d$prdyad == 1 & d$year >= 1949 & d$year <= 1984, ]
d$year_dv <- d$year + 1
n_d <- nrow(d); cat("PR dyads, DV years 1950-85:", n_d, " onset mean:", round(mean(d$dispon_l1, na.rm = TRUE), 4),
                    " PTA share:", round(mean(d$PTASameGT0), 4), "\n")

# ---- 4. headline model: Table 1 col 1 ----
rhs0 <- c("polMin", "polMax", "rgdppc_p3Min", "allies", "contigEQ1", "caprat", "traBilMax", "rgdpMin", "rgdpMax", "PTASameGT0",
          "PTA_traBilMax", "PTA_rgdpMin", "PTA_rgdpMax", "gdpHegPerc", paste0("dispon_spl", 0:3))
vars <- c("dispon_l1", rhs0, "dyadid")
fitter <- function(r, data) fit_bin("dispon_l1", r, data, "logit", "~dyadid", terms = c("PTA_traBilMax", "PTASameGT0", "traBilMax"))
d0 <- d[complete.cases(d[, vars]), ]
stopifnot(nrow(d0) > 15000)
full <- fitter(rhs0, d0)
cat("\nREBUILD vs published (Table 1 col 1): PTA_traBilMax b", full$PTA_traBilMax[["b"]], "se", full$PTA_traBilMax[["se"]],
    "p", full$PTA_traBilMax[["p"]], "N", attr(full, "n"), "| published -3.05e-7 (z -3.52) N 20990 (trade in different units)\n")

# history keyed to the OUTCOME year (DV = onset at t+1), as in study 86's port
ds <- common_sample(d, vars, "c1", "c2", "year_dv")
M <- run_six(rhs0, ds, fitter)
saveRDS(list(M = M, full = full, ref = ref, n_used = nrow(ds)), file.path(W, "M_P9_031.rds"))
# Sensitivity hook: not wired. The key term is a product (PTA x TRADE) whose constituents must stay in the model -> "split" class.

meta <- list(pid = "P9_031", cite = "Mansfield & Pevehouse 2000", journal = "IO", year = 2000, doi = "10.1162/002081800551361", family = "TIE")
art <- "logit, Huber SE clustered on dyad, natural spline (3 knots) of years since last MID onset"
used <- "feglm logit, vcov clustered dyadid, HB&M dispon spline terms, all six steps"
ahc <- "N (peace-year spline only)"
nt <- sprintf(paste0("REBUILT from HB&M 2012 PTAMIDTables.dta (sweep study 86 file; not O&R TRIANGLE): PR dyads, DV onset years 1950-85, ",
  "HB&M PTA coding (r ~ .84 with M&P PTA), Gleditsch trade/GDP (units differ from M&P 1985 USD, so b scale is not comparable), ",
  "EUGene allies/contiguity/CINC ratio, Polity IV. Full rebuilt N %d (published 20,990): b=%.4g se=%.4g p=%.3g. ",
  "Reference HB&M replication command (disp spline, year<1986): b=%.4g se=%.4g N %d (HB&M published -2.14e-4 (1.13e-4) N 25,209). ",
  "Common 1822-2011 sample N %d. History keyed to DV year. Interaction constituents kept."),
  attr(full, "n"), full$PTA_traBilMax[["b"]], full$PTA_traBilMax[["se"]], full$PTA_traBilMax[["p"]],
  ref$PTA_traBilMax[["b"]], ref$PTA_traBilMax[["se"]], attr(ref, "n"), nrow(ds))
write_p8(meta, "PTA_traBilMax",
  "abstract: 'hostilities between PTA members are less likely to occur as trade flows rise between them'",
  "Table 1 col 1 (logit MID onset; PR dyads 1950-85; TRADE x PTA)", art, used,
  -3.05e-7, 3.05e-7 / 3.52, 20990, TRUE, TRUE, M, ahc, nt)
dir.create(file.path(OUT, "results_p9"), recursive = TRUE, showWarnings = FALSE)
src <- file.path(OUT, "results_p8", "P9_031__PTA_traBilMax.csv"); dst <- file.path(OUT, "results_p9", "P9_031__PTA_traBilMax.csv")
stopifnot(file.exists(src), file.rename(src, dst))
sink()
