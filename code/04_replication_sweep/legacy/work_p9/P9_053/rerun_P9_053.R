# P9_053 Smith & Vreeland 2023 II 49(4) "UN Security Council membership: Increased security and reduced conflict" (10.1080/03050629.2023.2190589)
# Data + do-file: Harvard Dataverse doi:10.7910/DVN/CEPBYN ("UNSC and MIDS Replication.do", "monadic_working.dta", "working_UNSCdeter.dta");
#   staged as smith.do / monadic_working.dta / working_UNSCdeter.dta. No log in the archive.
# HEADLINE: Table 6 column 1 ("Table 6 presents the main results"; directed dyads, politically relevant = conttype 1-3, 1951-2014,
#   non-P5 pairs). Do-file block "DyadPrElect": logit mids unscA unscB FPrUNSCpooledA PrUNSCpooledA FPrUNSCpooledB PrUNSCpooledB
#   force* demaut* logpopVDEM* e_migdppcln* P5AllyA P5AllyB AnyAlly yr* if $SamCriteria3, cluster(dyad).
#   Stata expansions: force* = forceratio force2 (force2 = forceratio^2 via abbreviation), demaut* = demautA demautB demautAB,
#   logpopVDEM* / e_migdppcln* = A, B; yr* = yr yr2 yr3.
# PrUNSCpooled is generated first, as in the do-file's monadic block: logit t0 on rotate myturn SCratio demaut logpopVDEM e_migdppcln
#   cwar Lcwar (notElection==0 & single==1; 1951-2014, perm==0), predicted for all rows, 0 where notElection==1; FPr = F.Pr (t+1).
# KEY TERM unscA (A on the UNSC; restraint of initiation). SECOND ROW unscB (B on the UNSC; being targeted). Same model, same sample.
#   Published (article Table 6 col 1, author PDF, 3 decimals): unscA -1.240 (0.296), unscB -0.674 (0.240), N = 21,541.
# Estimator: fixest feglm logit, vcov ~dyad, all six steps. DV is contemporaneous (no lead). Rows are directed-dyad-MID records
#   (multiple rows per dyad-year when several MIDs), kept as in the author's model.
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R"))
WD <- file.path(ROOT, "work_p9/P9_053"); W <- file.path(OUT, "work_p9/P9_053"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_053.txt"), split = TRUE)

# ---------------- monadic block: predicted probability of UNSC election ----------------
m <- as.data.frame(haven::zap_labels(haven::read_dta(file.path(WD, "monadic_working.dta"),
  col_select = c("ccode", "year", "perm", "nobs", "t0", "notElection", "rotate", "myturn", "SCratio", "polity2", "logpopVDEM", "e_migdppcln", "cwar"))))
stopifnot(nrow(m) == 24337)
m <- m[!is.na(m$year) & m$year <= 2014 & !is.na(m$perm) & m$perm == 0 & m$year >= 1951, ]
m <- m[order(m$ccode, m$year, m$nobs, na.last = TRUE), ]
m$nn <- ave(seq_len(nrow(m)), m$ccode, m$year, FUN = seq_along)
m$single <- as.numeric(m$nn == 1)
m$demaut <- (10 + m$polity2) / 20                                  # do-file: demaut=(10+polity)/20 (polity abbreviates polity2)
s1 <- m[m$single == 1, c("ccode", "year", "cwar")]
lag1 <- merge(transform(s1, year = year + 1)[, c("ccode", "year", "cwar")], s1[, c("ccode", "year")], by = c("ccode", "year"))
names(lag1)[3] <- "Lcwar"
m <- merge(m, lag1, by = c("ccode", "year"), all.x = TRUE, sort = FALSE)
m <- m[order(m$ccode, m$year, m$nn), ]
ev <- c("rotate", "myturn", "SCratio", "demaut", "logpopVDEM", "e_migdppcln", "cwar", "Lcwar")
es <- m[!is.na(m$notElection) & m$notElection == 0 & m$single == 1 & complete.cases(m[, c("t0", ev)]), ]
ge <- glm(reformulate(ev, "t0"), data = es, family = binomial("logit"))
cat("election logit N =", nobs(ge), "\n"); print(round(summary(ge)$coefficients, 5))
m$PrUNSCpooled <- predict(ge, newdata = m, type = "response")
m$PrUNSCpooled[!is.na(m$notElection) & m$notElection == 1] <- 0
m <- m[m$single == 1, ]
stopifnot(!anyDuplicated(m[, c("ccode", "year")]))
nx <- transform(m[, c("ccode", "year", "PrUNSCpooled")], year = year - 1); names(nx)[3] <- "FPrUNSCpooled"
pr <- merge(m[, c("ccode", "year", "PrUNSCpooled")], nx, by = c("ccode", "year"), all.x = TRUE)
cat("probUNSC rows", nrow(pr), " nonmissing Pr", sum(!is.na(pr$PrUNSCpooled)), " FPr", sum(!is.na(pr$FPrUNSCpooled)), "\n")
rm(m, es, ge); invisible(gc())

# ---------------- dyadic block ----------------
dv <- c("year", "ccodeA", "ccodeB", "unscA", "unscB", "polity2A", "polity2B", "cincA", "cincB", "permA", "permB",
        "logpopVDEMA", "logpopVDEMB", "e_migdppclnA", "e_migdppclnB", "P5AllyA", "P5AllyB", "AnyAlly", "conttype", "mids")
d <- as.data.frame(haven::zap_labels(haven::read_dta(file.path(WD, "working_UNSCdeter.dta"), col_select = all_of(dv))))
cat("dyadic rows", nrow(d), " years", range(d$year, na.rm = TRUE), "\n")
d$mids[is.na(d$mids)] <- 0
d$forceratio <- d$cincA / (d$cincA + d$cincB)
d$demautA <- (d$polity2A + 10) / 20; d$demautB <- (d$polity2B + 10) / 20
d$dyad <- d$ccodeA * 1000 + d$ccodeB
d$yr <- (d$year - 1990) / 10; d$yr2 <- d$yr^2; d$yr3 <- d$yr^3
d$demautAB <- d$demautA * d$demautB
d$force2 <- d$forceratio^2
pA <- pr; names(pA) <- c("ccodeA", "year", "PrUNSCpooledA", "FPrUNSCpooledA")
pB <- pr; names(pB) <- c("ccodeB", "year", "PrUNSCpooledB", "FPrUNSCpooledB")
n0 <- nrow(d); d <- left_join(d, pA, by = c("ccodeA", "year")) %>% left_join(pB, by = c("ccodeB", "year")); stopifnot(nrow(d) == n0)
# $SamCriteria3 (Stata: missing conttype/perm fail the comparisons)
d <- d[!is.na(d$year) & d$year >= 1951 & d$year <= 2014 & !is.na(d$permA) & d$permA == 0 & !is.na(d$permB) & d$permB == 0 &
       !is.na(d$conttype) & d$conttype > 0 & d$conttype <= 3, ]
cat("SamCriteria3 rows", nrow(d), "\n")
terms <- c("unscA", "unscB")
rhs0 <- c("unscA", "unscB", "FPrUNSCpooledA", "PrUNSCpooledA", "FPrUNSCpooledB", "PrUNSCpooledB", "forceratio", "force2",
          "demautA", "demautB", "demautAB", "logpopVDEMA", "logpopVDEMB", "e_migdppclnA", "e_migdppclnB", "P5AllyA", "P5AllyB", "AnyAlly",
          "yr", "yr2", "yr3")
vars <- c("mids", rhs0, "dyad")
pub <- list(unscA = c(b = -1.240, se = 0.296), unscB = c(b = -0.674, se = 0.240)); pub_n <- 21541
fitter <- function(r, data) fit_bin("mids", r, data, "logit", "~dyad", terms = terms)
d0 <- d[complete.cases(d[, vars]), ]
full <- fitter(rhs0, d0)
for (t in terms) cat("\nREPRO vs published:", t, "b", full[[t]][["b"]], "(pub", pub[[t]][["b"]], ") se", full[[t]][["se"]], "(pub", pub[[t]][["se"]], ")")
cat("\nN", attr(full, "n"), "(pub", pub_n, ")\n")
stopifnot(attr(full, "n") == pub_n, all(sapply(terms, function(t) abs(full[[t]][["b"]] - pub[[t]][["b"]]) <= 0.0005)))
ds <- common_sample(d, vars, "ccodeA", "ccodeB", "year")
M <- run_six(rhs0, ds, fitter)
out <- list(M = M, full = full, n_used = nrow(ds))
saveRDS(out, file.path(W, "M_P9_053.rds"))

# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  # headline: unscA enters additively -> full sensitivity (unscB stays among the controls)
  sens_run(id = "P9_053", dat = ds, y = "mids", tie = "unscA", ctrl = setdiff(rhs0, "unscA"), fe = NULL,
           fit = function(dat, add) fitter(c(rhs0, add), dat)[["unscA"]])
}

meta <- list(pid = "P9_053", cite = "Smith 2023", journal = "II", year = 2023, doi = "10.1080/03050629.2023.2190589", family = "TIE")
art <- "logit, SE clustered on directed dyad, cubic time trend; politically relevant non-P5 directed dyads 1951-2014"
used <- "feglm logit, vcov clustered dyad, all six steps"
ahc <- "N (cubic time trend only)"
nt <- function(t, extra) sprintf(paste0("Directed dyad-MID records 1951-2014 (author's data; PrUNSC rebuilt from the do-file's election logit). ",
  "Published values from the author's PDF of the article (Table 6 col 1, 3 decimals). Full-sample repro b=%.6f se=%.6f N=%d vs pub %s (%s) N=%d. ",
  "Common 1822-2011 sample N=%d. %s"), full[[t]][["b"]], full[[t]][["se"]], attr(full, "n"), format(pub[[t]][["b"]]),
  format(pub[[t]][["se"]]), pub_n, nrow(ds), extra)
write_p8(meta, "unscA",
  "abstract: membership 'associated with ... lower rates of initiating conflict'",
  "Table 6 col 1 (dispute, pooled logit)", art, used,
  -1.240, 0.296, pub_n, TRUE, FALSE, M, ahc, nt("unscA", "Headline (initiator on UNSC). unscB kept in all steps."))
write_p8(meta, "unscB",
  "abstract: 'temporary membership on the UNSC is associated with lower rates of being targeted'",
  "Table 6 col 1 (dispute, pooled logit)", art, used,
  -0.674, 0.240, pub_n, TRUE, FALSE, M, ahc, nt("unscB", "Non-headline second row (target on UNSC) from the same model."))
dir.create(file.path(OUT, "results_p9"), recursive = TRUE, showWarnings = FALSE)
for (tm in terms) {
  src <- file.path(OUT, "results_p8", paste0("P9_053__", tm, ".csv")); dst <- file.path(OUT, "results_p9", paste0("P9_053__", tm, ".csv"))
  stopifnot(file.exists(src), file.rename(src, dst))
}
sink()
