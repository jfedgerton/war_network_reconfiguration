# P9_035 Johnson & Wolford 2023 JCR 67(4) "Alliance Reliability and Dispute Escalation" (10.1177/00220027221121140). REBUILD (protocol B).
# Replication materials are on the authors' personal sites only; Harvard Dataverse title search (once, this task) found nothing.
# Spec from the PDF (pp. 627-631, Table 1 rasterized for signs):
#   Directed dyad-years 1816-2000 in which target T has >= 1 ally committed to its defense against challenger C (Johnson & Leeds 2011
#   coding of ATOP). Heckman-style censored probit (Stata heckprob): selection = C initiates a MID against T (Maoz dyadic MID; 997);
#   outcome = T responds with a display or use of force (472). Table 1 Model 3, N = 559,111, uncensored 997, rho = -0.456.
#   Both equations: ally significant power change (any of T's allies' CINC >= 10% above/below its value in the alliance's formation
#   year; Leeds 2003 / Johnson & Joiner 2019), target's probability of winning (capability ratio incl. relevant allies), challenger has a
#   relevant offensive alliance, challenger has a relevant neutrality pact, age of target's alliance. Selection only: ln distance,
#   challenger-target S score (EUGene), peace years + square + cube (Carter & Signorino). "Standard errors" (no robust/cluster wording).
# HEADLINE / KEY (fixed from the PDF before estimating): escalation (outcome) equation, ally significant power change, Model 3:
#   -0.341 (SE 0.138), p < .05. Abstract: less-reliable alliances "decrease the probability of escalation relative to reliable alliances".
# Estimator: heckprob by ML via the phase-8 port work_p8/heckprob_b3.R (precedent P8_003), OIM SEs, history terms in both equations.
# Rebuild inputs (all peacesciencer 1.2.0): COW directed dyad-years, ATOP 5 dyad-year alliances (atop_alliance), COW MID 5.0
#   (add_cow_mids, peace years via add_peace_years), NMC CINC, capital distance, Haege S (valued ATOP alliance, absolute distance,
#   srsvaa; extdata dyadic_fp_similarity.rds).
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R"))      # ROOT, OUT, dplyr, fixest, add_history, common_sample, run_six, write_p8; seed 123
source(file.path(CODE, "work_p8/heckprob_b3.R"))    # fit_heckprob, hp_parts
suppressPackageStartupMessages(library(peacesciencer))
set.seed(123)
W <- file.path(OUT, "work_p9/P9_035"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_035.txt"), split = TRUE)

# ---- 1. directed dyad-years 1816-2000: MIDs, peace years, distance, CINC, S ----
fp <- system.file("extdata", "dyadic_fp_similarity.rds", package = "peacesciencer"); stopifnot(nzchar(fp), file.exists(fp))
d <- create_dyadyears(directed = TRUE, subset_years = 1816:2000) %>% add_cow_mids() %>% add_peace_years() %>%
  add_capital_distance() %>% add_nmc() %>% as.data.frame()
n0 <- nrow(d); stopifnot(anyDuplicated(d[, c("ccode1", "ccode2", "year")]) == 0)
s <- readRDS(fp) %>% as.data.frame() %>% select(ccode1, ccode2, year, srsvaa); stopifnot(anyDuplicated(s[, c("ccode1", "ccode2", "year")]) == 0)
d <- d %>% left_join(s, by = c("ccode1", "ccode2", "year")); stopifnot(nrow(d) == n0)
d <- d %>% select(ccode1, ccode2, year, cowmidonset, sidea1, orig1, hostlev2, cowmidspell, capdist, cinc1, cinc2, srsvaa)
cat("directed dyad-years 1816-2000:", n0, "\n")

# ---- 2. ATOP directed alliance edges (atop_alliance holds both orders) ----
a <- as.data.frame(atop_alliance) %>% filter(year >= 1800, year <= 2000) %>% select(x = ccode1, y = ccode2, year, def = atop_defense, off = atop_offense, neu = atop_neutral)
stopifnot(anyDuplicated(a[, c("x", "y", "year")]) == 0, sum(a$x < a$y) == sum(a$x > a$y))
cinc <- as.data.frame(cow_nmc) %>% select(ccode, year, cinc)
# defense edges T(x) -> ally A(y): formation year = first year of the current consecutive defense spell
de <- a %>% filter(def == 1) %>% arrange(x, y, year) %>% group_by(x, y) %>%
  mutate(sp = cumsum(c(1, diff(year) != 1))) %>% group_by(x, y, sp) %>% mutate(fy = min(year)) %>% ungroup() %>% select(x, y, year, fy)
de <- de %>% left_join(cinc %>% rename(ca = cinc), by = c("y" = "ccode", "year")) %>% left_join(cinc %>% rename(cf = cinc), by = c("y" = "ccode", "fy" = "year")) %>%
  mutate(pc = as.numeric(abs(ca / cf - 1) >= 0.10), age = year - fy) %>% as.data.frame()
cat("defense edges:", nrow(de), " power-change share:", round(mean(de$pc, na.rm = TRUE), 3), " missing pc:", sum(is.na(de$pc)), "\n")
# per target-year aggregates over ALL its defense allies (corrected below when the challenger is one of them)
ta <- de %>% group_by(x, year) %>% arrange(desc(age), .by_group = TRUE) %>%
  summarise(n_al = n(), npc = sum(pc, na.rm = TRUE), scinc = sum(ca, na.rm = TRUE),   # missing CINC: ally contributes 0
             age1 = age[1], arg1 = y[1], age2 = if (n() > 1) age[2] else NA_real_, .groups = "drop")
# challenger offense / neutrality partners and offense partners' CINC
ca_ <- a %>% left_join(cinc, by = c("y" = "ccode", "year")) %>% group_by(x, year) %>%
  summarise(noff = sum(off), nneu = sum(neu), soff = sum(off * coalesce(cinc, 0)), .groups = "drop")
# ---- 3. merge onto directed dyads: T = ccode2, C = ccode1 ----
d <- d %>% left_join(ta, by = c("ccode2" = "x", "year")) %>% left_join(ca_, by = c("ccode1" = "x", "year")) %>%
  left_join(a %>% select(x, y, year, defCT = def, offCT = off, neuCT = neu), by = c("ccode2" = "x", "ccode1" = "y", "year")) %>%   # T's pact with C
  left_join(de %>% select(x, y, year, pcC = pc, ageC = age, caC = ca), by = c("ccode2" = "x", "ccode1" = "y", "year"))
stopifnot(nrow(d) == n0)
d <- d %>% mutate(across(c(n_al, npc, scinc, noff, nneu, soff, defCT, offCT, neuCT), ~coalesce(.x, 0)), isC = coalesce(defCT, 0) == 1) %>%
  mutate(n_rel = n_al - isC,                                          # relevant allies exclude the challenger
         allypc = as.numeric(npc - ifelse(isC, coalesce(pcC, 0), 0) > 0),
         allycinc = scinc - ifelse(isC, coalesce(caC, 0), 0),
         allyage = ifelse(isC & arg1 == ccode1, age2, age1),
         chal_off = as.numeric(noff - offCT > 0), chal_neu = as.numeric(nneu - neuCT > 0),
         offcinc = soff - ifelse(offCT == 1, cinc2, 0),
         pwin_t = (cinc2 + allycinc) / (cinc2 + allycinc + cinc1 + offcinc),
         init = as.numeric(coalesce(cowmidonset, 0) == 1 & coalesce(sidea1, 0) == 1 & coalesce(orig1, 0) == 1),
         esc = ifelse(init == 1, as.numeric(hostlev2 >= 3), NA_real_),
         ldist = log(capdist), sscore = srsvaa, py = cowmidspell, py2 = py^2, py3 = py^3, dyadid = 1000 * ccode1 + ccode2)
d <- d[d$n_rel >= 1, ]
cat("dyad-years with >= 1 relevant defensive ally of the target:", nrow(d), "(paper 559,111 after missing data)\n")
cat("initiations:", sum(d$init), "(paper 997)  escalations:", sum(d$esc, na.rm = TRUE), "(paper 472)  esc missing among inits:", sum(d$init == 1 & is.na(d$esc)), "\n")
cat("means: allypc", round(mean(d$allypc), 3), " pwin_t", round(mean(d$pwin_t, na.rm = TRUE), 3), " chal_off", round(mean(d$chal_off), 3),
    " chal_neu", round(mean(d$chal_neu), 3), " allyage", round(mean(d$allyage, na.rm = TRUE), 2), "\n")
d <- as.data.frame(d)

# ---- 4. heckprob: escalation | initiation ----
outX <- c("allypc", "pwin_t", "chal_off", "chal_neu", "allyage")
selZ <- c(outX, "ldist", "sscore", "py", "py2", "py3")
KEYS <- c("out:allypc", "sel:allypc")
d$esc_chk <- ifelse(d$init == 1, d$esc, 0)
vars <- unique(c("init", "esc_chk", selZ))
# gradient check of the analytic heckprob scores on synthetic data (fail loudly), as P8_003
nn <- 4000; Xs0 <- cbind(a = rnorm(nn), c = 1); Zs0 <- cbind(a = rnorm(nn), b = rnorm(nn), c = 1)
ee <- MASS::mvrnorm(nn, c(0, 0), matrix(c(1, .4, .4, 1), 2)); ss <- as.numeric(Zs0 %*% c(.5, -.3, -.2) + ee[, 2] > 0)
yy <- as.numeric(Xs0 %*% c(.7, .1) + ee[, 1] > 0); yy[ss == 0] <- 0; th <- c(.6, .05, .45, -.25, -.15, .3)
ga <- colSums(hp_parts(th, Xs0, Zs0, yy, ss)$S); gn <- numDeriv::grad(function(t) sum(hp_parts(t, Xs0, Zs0, yy, ss)$ll), th)
cat("synthetic gradient check max abs diff:", max(abs(ga - gn)), "\n"); stopifnot(max(abs(ga - gn)) < 1e-5)
fit_hp <- function(extra, data) {
  X <- cbind(as.matrix(data[, c(outX, extra), drop = FALSE]), `_cons` = 1)
  Z <- cbind(as.matrix(data[, c(selZ, extra), drop = FALSE]), `_cons` = 1)
  f <- fit_heckprob(X, Z, data$esc_chk, data$init)
  z <- f$coef / f$se
  cat("---- heckprob | extra:", paste(extra, collapse = " + "), "| N =", f$n, "selected =", f$n_sel, "| LL =", round(f$ll, 3), "| rho =", round(f$rho, 4), "\n")
  print(round(cbind(b = f$coef, se = f$se, z = z), 5))
  out <- lapply(setNames(KEYS, KEYS), function(k) c(b = f$coef[[k]], se = f$se[[k]], p = 2 * pnorm(-abs(z[[k]]))))
  attr(out, "n") <- f$n; out
}
d0 <- d[complete.cases(d[, vars]), ]
cat("complete-case N:", nrow(d0), " selected:", sum(d0$init), " escalated:", sum(d0$esc_chk), "\n")
full <- fit_hp(character(0), d0)
cat("\nREBUILT vs published (Model 3): out:allypc b", full[["out:allypc"]][["b"]], "se", full[["out:allypc"]][["se"]], "p", full[["out:allypc"]][["p"]],
    "N", attr(full, "n"), "| published -0.341 (0.138), N 559,111, 997 uncensored;  sel:allypc", full[["sel:allypc"]][["b"]], "(pub 0.116 (0.043))\n")
rm(d0); invisible(gc())
ds <- common_sample(d, vars, "ccode1", "ccode2", "year"); rm(d); invisible(gc())
fitter <- function(rhs, data) fit_hp(setdiff(rhs, vars), data)
M <- run_six(vars, ds, fitter)
saveRDS(list(M = M, full = full, n_used = nrow(ds)), file.path(W, "M_P9_035.rds"))
# no sensitivity hook: two-equation selection model (spec: hooks only for single-equation logit/probit)

meta <- list(pid = "P9_035", cite = "Johnson & Wolford 2023", journal = "JCR", year = 2023, doi = "10.1177/00220027221121140", family = "TIE")
art <- "censored (Heckman-style) probit, selection = MID initiation, outcome = target militarized response; standard (OIM) SE"
used <- "heckprob ML (work_p8/heckprob_b3.R port), OIM SE, history terms in both equations, all six steps"
ahc <- "N (peace years + square + cube in the selection equation only)"
nt <- sprintf(paste0("REBUILT (no archive found on Dataverse): directed COW dyad-years 1816-2000 where the target has >= 1 ATOP defense ally other than the challenger ",
  "(paper: Johnson & Leeds 2011 relevance coding vs the specific challenger, not reproducible). Power change = any such ally's CINC >= 10%% from its value in the formation ",
  "year of the current defense spell; age = oldest such spell; challenger offense/neutrality = any ATOP offense/neutrality partner other than the target; ",
  "target win prob = (T + allies CINC)/(T + allies + C + C offense partners). DVs from COW MID 5.0 (add_cow_mids): initiation = onset, C originator on side A; ",
  "escalation = target hostlev >= 3. S = Haege valued-ATOP absolute S (srsvaa; paper EUGene S); ln capital distance; COW peace years. ",
  "Full rebuilt out:allypc b=%.4f se=%.4f p=%.4g N=%d. Common 1822-2011 N=%d."),
  full[["out:allypc"]][["b"]], full[["out:allypc"]][["se"]], full[["out:allypc"]][["p"]], attr(full, "n"), nrow(ds))
write_p8(meta, "out_allypc", "abstract: less-reliable alliances 'decrease the probability of escalation relative to reliable alliances'",
  "Table 1 Model 3, escalation equation (censored probit; ally significant power change)", art, used,
  -0.341, 0.138, 559111, TRUE, TRUE,
  lapply(M, function(m) { o <- list(out_allypc = m[["out:allypc"]]); attr(o, "n") <- attr(m, "n"); o }), ahc, nt)
dir.create(file.path(OUT, "results_p9"), recursive = TRUE, showWarnings = FALSE)
src <- file.path(OUT, "results_p8", "P9_035__out_allypc.csv"); dst <- file.path(OUT, "results_p9", "P9_035__out_allypc.csv")
stopifnot(file.exists(src), file.rename(src, dst))
sink()
