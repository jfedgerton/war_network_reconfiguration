# X040 Bayer & Bernhard 2010 CMPS. No archive found -> REBUILD of Table 1 (1919-1992 politically relevant dyads, replication
# of Russett & Oneal 2001 Table A5.1) on the Russett-Oneal Triangulating Peace PR-dyad data (TRIANGLE.tab, Dataverse
# doi:10.7910/DVN/OMJYE5; already local from cid 20) - the same base data B&B replicate.
# Spec per paper text (author PDF, UMN repository): GEE logit AR(1), robust SE, DV = MID involvement, lagged IVs;
# controls: lower dependence, IGOs, power ratio, alliances, noncontiguity, log distance, only minor powers.
# Key terms (decided from abstract/text BEFORE estimating):
#   lowpol    = lower Polity score (Table 1 col 1, b_pub -0.068 [0.010]); highpol also in model
#   jointdem6 = both Polity >= 6 (Table 1 col 2 'Polity dichotomous', b_pub -1.270 [0.221]); singledem6 also in model
# BNR / Doyle / ACLP codings not available locally -> not rebuilt.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p4/common_p4_x8.R"))
suppressPackageStartupMessages(library(geepack))
WO <- file.path(OUT, "work_p4/X040"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work_p4/X040")
sink(file.path(WO, "log_verify_X040.txt"), split = TRUE)
source(file.path(CODE, "vb5/verify_helpers.R"))
suppressPackageStartupMessages(library(peacesciencer))

d <- read.delim(file.path(ROOT, "recA/or/TRIANGLE.tab"))
cat("TRIANGLE years:", range(d$year), " rows:", nrow(d), "\n")
d <- subset(d, year >= 1919 & year <= 1992)
# VERIFY FIX (B&B 2010 fn 7, p.93): "a version of Polity2 where normal range scores replaced the standardized authority
# scores (-77,-88)" -> use Polity IV polity2 (peacesciencer ccode_democracy) instead of TRIANGLE demaut.
pd <- as.data.frame(ccode_democracy)[, c("ccode", "year", "polity2")]
stopifnot(!anyDuplicated(pd[, c("ccode", "year")]))
d <- merge(d, setNames(pd, c("statea", "year", "p2a")), by = c("statea", "year"), all.x = TRUE)
d <- merge(d, setNames(pd, c("stateb", "year", "p2b")), by = c("stateb", "year"), all.x = TRUE)
cat("polity2 missing a/b:", sum(is.na(d$p2a)), sum(is.na(d$p2b)), " demaut missing a/b:", sum(is.na(d$demauta)), sum(is.na(d$demautb)), "\n")
cat("cor(demauta, p2a):", cor(d$demauta, d$p2a, use = "complete.obs"), "\n")
d$lowpol  <- pmin(d$p2a, d$p2b); d$highpol <- pmax(d$p2a, d$p2b)
d$jointdem6  <- as.numeric(d$p2a >= 6 & d$p2b >= 6)
d$singledem6 <- as.numeric(xor(d$p2a >= 6, d$p2b >= 6))
d$smldep <- pmin(d$dependa, d$dependb)
d$dyadid <- 1000 * d$statea + d$stateb
CTRL <- "smldep + smigoabi + lcaprat2 + allies + noncontg + logdstab + minrpwrs"
SPEC <- list(lowpol    = list(rhs = paste("highpol + lowpol +", CTRL), lab = "Table 1 col 1 (Polity interval)", b_pub = -0.068, se_pub = 0.010),
             jointdem6 = list(rhs = paste("singledem6 + jointdem6 +", CTRL), lab = "Table 1 col 2 (Polity dichotomous)", b_pub = -1.270, se_pub = 0.221))

fit_gee <- function(rhs, extra, s, term, lab) {
  f <- as.formula(paste("dispute1 ~", rhs, extra))
  s <- s[order(s$dyadid, s$year), ]
  m <- geeglm(f, data = s, id = dyadid, waves = year, family = binomial, corstr = "ar1")
  ct <- summary(m)$coefficients; ct <- as.matrix(ct[, c(1, 2, 3, 4)])
  cat("alpha:", round(unlist(m$geese$alpha), 4), "\n")
  grab(ct, term, nrow(s), paste(lab, extra))
}

meta <- list(xid = "X040", short_cite = "Bayer & Bernhard 2010", journal = "CMPS", year = 2010, doi = "10.1177/0738894209352133")
ev <- "abstract: 'dichotomous coding better captures the notion of democracy' (joint democracy vs Polity scale)"
# [2026-11-13 audit] TRIANGLE-based dispute1 is MID at t+1 -> anchor history on the outcome year
d$yr_out <- d$year + 1
d <- add_history(d, "statea", "stateb", "yr_out")
for (k in names(SPEC)) {
  sp <- SPEC[[k]]
  x <- d[complete.cases(d[, all.vars(as.formula(paste("dispute1 ~", sp$rhs)))]), ]
  r1 <- fit_gee(sp$rhs, "", x, k, sp$lab)[[k]]
  rep <- if (sign(r1$b) == sign(sp$b_pub) && r1$p < .05) "REBUILT_MATCH" else "REBUILT_MISMATCH"
  s <- subset(x, year >= 1822 & year <= 2011); s <- s[complete.cases(s[, HISTV]), ]
  cat(k, ": rebuild N", nrow(x), " common N", nrow(s), "\n")
  M <- run_six(function(extra, z) fit_gee(sp$rhs, extra, z, k, sp$lab), s)
  # ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
  if (Sys.getenv("SWEEP_SENS") == "1" && k == "jointdem6") {   # [2026-11-13 PI] headline term is now jointdem6 (authors' preferred dichotomous coding; lowpol is the R&O foil)
    source(file.path(CODE, "hist/sens_helper.R"))
    sens_run(id = "X040", dat = s, y = "dispute1", tie = "jointdem6",
             ctrl = c("singledem6", "smldep", "smigoabi", "lcaprat2", "allies", "noncontg", "logdstab", "minrpwrs"), fe = NULL,
             fit = function(dat, add) {
               r <- fit_gee(SPEC$jointdem6$rhs, if (length(add)) paste("+", paste(add, collapse = " + ")) else "", dat, "jointdem6", SPEC$jointdem6$lab)$jointdem6
               c(b = unname(r$b), se = unname(r$se), p = unname(r$p)) })
  }
  nt <- sprintf(paste0("REBUILT (%s) %s on Russett-Oneal TRIANGLE PR dyads 1919-1991 (authors' base data; B&B restrict to BNR-common sample, pub N=35,052; ours N=%d). ",
    "GEE logit AR(1) robust (geepack, waves=year). b_rebuild %.3f (SE %.3f) vs pub %.3f (%.3f). Polity = demaut (democ-autoc) in TRIANGLE; dichotomous at >=6 per paper. ",
    "BNR/Doyle/ACLP dichotomous codings not rebuilt (data not local). History anchored on outcome year (year+1); common N=%d."),
    rep, sp$lab, nrow(x), r1$b, r1$se, sp$b_pub, sp$se_pub, M$m0[[k]]$n)
  hl <- if (k == "lowpol") "Table 1 col 1 (1919-92 PR, Polity interval)" else "Table 1 col 2 (1919-92 PR, Polity dichotomous)"
  R <- lapply(M, function(z) z[[k]])
  write_verify("X040", "Bayer & Bernhard 2010", "DEMOCRACY", k, hl, sp$b_pub, sp$se_pub, 35052,
    "GEE logit AR(1), robust SE; RO 2001 Table A5.1 design, PR nondirected dyads, DV MID involvement (lagged IVs)",
    "democracy from TRIANGLE demaut (paper: Polity2, fn 7);sample N 35747 vs 35052 (paper restricts to BNR-common dyad-years)",
    "PARTIAL", "BNR-common sample restriction (BNR data not local)",
    verdict_rebuilt(R$m0$b, R$m0$p, sp$b_pub, TRUE, "X040", k, sp$se_pub), R,
    sprintf("REBUILT on TRIANGLE with Polity IV polity2 (peacesciencer) replacing demaut per fn 7. rebuild (full) b %.3f SE %.3f N %d; pub %.3f (%.3f). pct_vs_pub uses m0 on common 1822-2011 sample. bD is NOT a placebo for a democracy term.", r1$b, r1$se, nrow(x), sp$b_pub, sp$se_pub))
}
sink()
