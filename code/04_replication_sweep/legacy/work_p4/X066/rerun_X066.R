# X066 Aydin 2010 JPR "The deterrent effects of economic integration" (10.1177/0022343310370290)
# Data + logs: PRIO JPR archive (Replication - Aydin 47(5).zip): Models1-2.dta (Bennett enduring rivals), Models 3-4.dta (Thompson
#   strategic rivals), Stata logs Models1-2.smcl / Models3-4.smcl.
# Headline: Model 2 (first institutionalized-trade model; enduring rivals), from the author's log:
#   xtgee initiation instrade inscapacity weaklink jointdem development alliance power defensiveally interest, robust
#   (group var = dyad; Stata xtgee defaults: gaussian, identity link, exchangeable; semi-robust SE). N = 1494, 94 groups.
#   Log: instrade -.058131 (se .0230862, p .012); jointdem -.0905828 (p .010); development -.006762 (p .496).
# KEY (fixed from abstract BEFORE estimating): instrade = defenders' trade with the target channelled through regional trade
#   institutions ("extended deterrence success is most likely ... economically integrated through regional trade institutions").
# CONTROL rows: jointdem (DEMOCRACY), development (CAPITALISM family: economic development).
# DATA PROBLEM: the posted files have NO year variable (only directed dyad code = attacker*1000 + target). The history merge needs years.
#   Reconstruction attempted: rows are stored year-sorted (each dyad appears once per block; block sizes grow as rivalries accumulate),
#   so year = base + block index, where blocks split at the first repeated dyad. Base year chosen by maximizing agreement of the
#   author's jointdem with Polity (both polity2 >= 6) and validated; history is run ONLY if agreement >= 0.95.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p4/common_p4_b02.R"))
suppressPackageStartupMessages(library(haven))
WO <- file.path(OUT, "work_p4/X066"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work_p4/X066"); sink(file.path(WO, "log_X066.txt"), split = TRUE)

d <- as.data.frame(zap_labels(read_dta(file.path(W, "Replication/Models1-2.dta"))))
d$.ord <- seq_len(nrow(d))
stopifnot(nrow(d) == 1594, all(d$enduringrivalry == 1))
f0 <- initiation ~ instrade + inscapacity + weaklink + jointdem + development + alliance + power + defensiveally + interest
K <- c("instrade", "jointdem", "development")
FF <- ff_gee_gauss_exch("dyad")
r0 <- FF(f0, d); print(round(r0, 6)); n_orig <- attr(r0, "n")
bp <- c(instrade = -.058131, jointdem = -.0905828, development = -.006762); pp <- c(instrade = .012, jointdem = .010, development = .496)
st <- sapply(K, function(k) if (abs(r0[k, 1] - bp[[k]]) / abs(bp[[k]]) <= .10 && ((r0[k, 3] < .05) == (pp[[k]] < .05))) "REPRODUCED" else "NOT_REPRODUCED")
print(st); cat("repro N:", n_orig, "(log 1494)\n")

# ---- year reconstruction ----
blk <- integer(nrow(d)); seen <- c(); b <- 0L
for (i in seq_len(nrow(d))) { if (d$dyad[i] %in% seen) { b <- b + 1L; seen <- c() }; seen <- c(seen, d$dyad[i]); blk[i] <- b }
d$blk <- blk
d$c1 <- d$dyad %/% 1000; d$c2 <- d$dyad %% 1000
cat("blocks:", max(blk) + 1, "\n")
# consecutive-block check within dyad (rivalry spells can have gaps, so not required to be 100%)
cons <- tapply(d$blk, d$dyad, function(x) all(diff(sort(x)) == 1)); cat("dyads with consecutive blocks:", mean(cons), "\n")
pol <- .H$pol
agree <- sapply(1930:1965, function(b0) {
  x <- data.frame(c1 = d$c1, c2 = d$c2, year = b0 + d$blk, jd = d$jointdem) %>%
    left_join(pol, by = c("c1" = "ccode", "year")) %>% rename(pa = polity2) %>%
    left_join(pol, by = c("c2" = "ccode", "year")) %>% rename(pb = polity2)
  ok <- !is.na(x$pa) & !is.na(x$pb); mean(x$jd[ok] == as.numeric(x$pa[ok] >= 6 & x$pb[ok] >= 6))
})
names(agree) <- 1930:1965; print(round(agree, 3))
base <- as.integer(names(which.max(agree))); best <- max(agree)
cat("best base year:", base, " agreement:", best, " runner-up:", sort(agree, decreasing = TRUE)[2], "\n")
d$year <- base + d$blk

meta <- list(xid = "X066", short_cite = "Aydin 2010", journal = "JPR", year = 2010, doi = "10.1177/0022343310370290")
est <- "GEE gaussian/identity, exchangeable, semi-robust SE (geepack; Stata xtgee ..., robust)"
hc <- "N (no peace years / conflict-history term)"
nt0 <- paste0("Author's PRIO data + log, Model 2 (enduring rivals, Bennett 1998; directed attacker-target dyads, post-1945). b_pub from author's log. ",
  "Posted data lack a year variable: year reconstructed as base + year-block index (file is year-sorted; blocks split at first repeated dyad); ",
  "base year ", base, " chosen by agreement of author's jointdem with Polity polity2>=6 (agreement ", round(best, 3), "). ")
if (best >= 0.95) {
  ds <- common_sample(d, f0, "dyad", "c1", "c2", "year")
  ds <- ds[order(ds$dyad, ds$.ord), ]
  M <- run_six(f0, ds, K, FF, WO)
  # ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
  if (Sys.getenv("SWEEP_SENS") == "1") {
    source(file.path(CODE, "hist/sens_helper.R"))
    sens_run(id = "X066", dat = ds, y = "initiation", tie = "instrade",
             ctrl = c("inscapacity", "weaklink", "jointdem", "development", "alliance", "power", "defensiveally", "interest"), fe = NULL,
             fit = function(dat, add) { ct <- FF(if (length(add)) update(f0, as.formula(paste(". ~ . +", paste(add, collapse = " + ")))) else f0, dat)
               c(b = unname(ct["instrade", 1]), se = unname(ct["instrade", 2]), p = unname(ct["instrade", 3])) })
  }
  nt <- paste0(nt0, "repro N ", n_orig, "; common N ", M$m0$n, ". LPM-GEE coefficients (probability scale).")
  write_p4(meta, "instrade", "TIE", "KEY", "\"extended deterrence success is most likely ... economically integrated through regional trade institutions\"",
    st[["instrade"]], bp[["instrade"]], r0["instrade", 1], est, n_orig, M, hc, paste(nt, "instrade = defenders' RIA-channelled trade with target."))
  write_p4(meta, "jointdem", "DEMOCRACY", "CONTROL", "control in headline Model 2 (attacker-target joint democracy)",
    st[["jointdem"]], bp[["jointdem"]], r0["jointdem", 1], est, n_orig, M, hc, paste(nt, "bD NOT a placebo for this term."))
  write_p4(meta, "development", "CAPITALISM", "CONTROL", "control in headline Model 2 (economic development)",
    st[["development"]], bp[["development"]], r0["development", 1], est, n_orig, M, hc, nt)
} else {
  cat("Year reconstruction NOT validated (agreement < .95): history not merged.\n")
  saveRDS(list(r0 = r0, st = st, agree = agree), file.path(WO, "repro_only.rds"))
}
sink()
