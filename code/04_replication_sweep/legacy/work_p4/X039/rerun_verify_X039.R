# VERIFY X039 Quackenbush & Rudy 2009 CMPS -- copy of rebuild_X039.R matched to the PDF (Table 2 Model 1, p.279).
# Paper: logit of MID occurrence, non-directed politically ACTIVE dyad-years 1816-2001, ongoing dispute-years dropped unless new onset,
#   joiners dropped; both democratic / one democratic (polity2 >= 5), S score, ln distance (contiguous = 0), power parity,
#   peace years + 3 cubic splines; robust SE clustered on dyad. Published: both democratic -0.3527 (0.1279), one democratic
#   0.4586 (0.0776); N = 163,920.
# Only change vs rebuild_X039.R: polity2 joined from peacesciencer::ccode_democracy (COW codes). add_democracy() (v1.2.0) returned NA
#   polity2 for USSR 1923-91, Yugoslavia, Italy etc., silently dropping those dyads (the paper includes them).
# Not fixable: politically ACTIVE dyads (Quackenbush 2006 rule) -> PRD (fn.5: results unchanged with PRD); MID v3 -> GML (joiner
#   exclusion not replicated); BTSCS knots unreported (cubic kept).
suppressMessages({library(peacesciencer)})
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "vb4/vhelp.R"))
WO <- file.path(OUT, "work_p4/X039"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work_p4/X039"); sink(file.path(WO, "log_verify_X039.txt"), split = TRUE)
fp_f <- system.file("extdata", "dyadic_fp_similarity.rds", package = "peacesciencer"); stopifnot(nzchar(fp_f), file.exists(fp_f))  # port
fp <- readRDS(fp_f)
scol <- grep("^srs", names(fp), value = TRUE); scol <- scol[!grepl("w", sub("^srs", "", scol))]; stopifnot(length(scol) >= 1)
scol <- scol[1]; cat("using S column:", scol, "\n"); stopifnot(scol == "srsvas")
d <- create_dyadyears(directed = FALSE, subset_years = 1816:2001) %>% add_gml_mids() %>% add_contiguity() %>% add_cow_majors() %>%
  add_nmc() %>% add_capital_distance()
prdkeys <- d %>% filter_prd() %>% select(ccode1, ccode2, year) %>% as.data.frame()
cat("PRD dyad-years:", nrow(prdkeys), "\n")
pol <- ccode_democracy %>% select(ccode, year, polity2) %>% as.data.frame(); stopifnot(!anyDuplicated(pol[, 1:2]))
d <- d %>% left_join(fp[, c("ccode1","ccode2","year", scol)], by = c("ccode1","ccode2","year")) %>%
  left_join(pol %>% rename(polity21 = polity2), by = c("ccode1" = "ccode", "year")) %>%
  left_join(pol %>% rename(polity22 = polity2), by = c("ccode2" = "ccode", "year"))
d <- d %>% arrange(ccode1, ccode2, year) %>% group_by(ccode1, ccode2) %>%
  mutate(spell = cumsum(lag(gmlmidonset, default = 0) == 1)) %>% group_by(ccode1, ccode2, spell) %>% mutate(py = row_number() - 1) %>% ungroup()
d <- d %>% filter(!(gmlmidongoing == 1 & gmlmidonset == 0))
d <- d %>% semi_join(prdkeys, by = c("ccode1","ccode2","year"))
d <- d %>% mutate(dem1 = as.numeric(polity21 >= 5), dem2 = as.numeric(polity22 >= 5),
                  bothdem = dem1 * dem2, onedem = as.numeric(dem1 + dem2 == 1),
                  sscore = .data[[scol]], dist0 = ifelse(conttype == 1, 0, capdist), lndist = log1p(dist0),
                  parity = pmin(cinc1, cinc2) / pmax(cinc1, cinc2), py2 = py^2 / 10, py3 = py^3 / 100, dyad = ccode1 * 1000 + ccode2) %>%
  as.data.frame()
rhs0 <- c("bothdem","onedem","sscore","lndist","parity","py","py2","py3")
vars <- c("gmlmidonset", rhs0, "dyad")
d <- d[is.finite(d$parity) | is.na(d$parity), ]
FF <- function(r, s) fit_bin("gmlmidonset", r, s, "logit", "~dyad", NULL, c("bothdem", "onedem"))
x <- d[complete.cases(d[, vars]), ]; rep <- FF(rhs0, x); cat("rebuilt N =", attr(rep, "n"), "(pub 163,920; previous rebuild 74,907)\n"); rm(x)
s <- common_sample(d, vars, "ccode1", "ccode2", "year"); rm(d); invisible(gc())
M <- run_six(rhs0, s, FF)
# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "X039", dat = s, y = "gmlmidonset", tie = "bothdem", ctrl = setdiff(rhs0, "bothdem"), fe = NULL,
           fit = function(dat, add) { r <- FF(c(rhs0, add), dat)$bothdem; c(b = unname(r[["b"]]), se = unname(r[["se"]]), p = unname(r[["p"]])) })
}
mism <- "add_democracy() NA polity2 for USSR/Yugoslavia/Italy etc. dropped those dyads; PRD instead of politically active dyads; GML instead of MID v3 (joiners); cubic peace years instead of BTSCS splines"
unf <- "politically active dyad rule (Quackenbush 2006) not rebuilt (PRD used; paper fn.5 says PRD results unchanged); MID v3 joiner exclusion; BTSCS knots unreported"
write_verify("X039", "Quackenbush & Rudy 2009", "DEMOCRACY", "bothdem", "Table 2 Model 1 (Polity, p.279)", -0.3527, 0.1279, 163920,
  "logit, dyad-clustered robust SE, peace years + 3 cubic splines", TRUE, mism, "PARTIAL", unf, TRUE, M,
  sprintf("REBUILD (no author files). Full-sample rebuilt b=%.4g (N=%d). bD (lagged joint democracy) NOT a placebo here.", rep$bothdem["b"], attr(rep, "n")))
write_verify("X039", "Quackenbush & Rudy 2009", "DEMOCRACY", "onedem", "Table 2 Model 1 (Polity, p.279)", 0.4586, 0.0776, 163920,
  "logit, dyad-clustered robust SE, peace years + 3 cubic splines", TRUE, mism, "PARTIAL", unf, TRUE, M,
  sprintf("REBUILD (no author files). Full-sample rebuilt b=%.4g (N=%d).", rep$onedem["b"], attr(rep, "n")))
sink()
