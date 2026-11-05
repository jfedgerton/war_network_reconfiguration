# X009 Gibler & Braithwaite 2013 BJPS "Dangerous Neighbours". Undirected dyad-year 1960-98(ish), DV fatalMID onset.
# Headline: authors' log model 3 (logit, cluster dyadid): lowdem, terr500k10yr, interaction4bjps (+ controls, peace yrs + 3 splines).
# KEY (from abstract, before estimating): interaction4bjps = lowdem x territorial hot spots ("interaction between joint democracy
#   and regional instability confirms..."), and lowdem ("joint democracy is only statistically significant ... once hot spots
#   are jointly estimated"). Both DEMOCRACY family.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p4/xb1/p4c.R"))
suppressPackageStartupMessages(library(readstata13))
WO <- file.path(OUT, "work_p4/X009"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work_p4/X009"); sink(file.path(WO, "log_X009.txt"), split = TRUE)
d <- read.dta13(file.path(W, "Replication Materials/gibler_braithwaite_BJPS_repldata.dta"), convert.factors = FALSE)
names(d)[names(d) == "fatalMID"] <- "fatal"
cat("rows", nrow(d), " years", range(d$year), "\n"); print(summary(d$dyadid))

# decode dyadid -> ccodes (expect c1*1000 + c2, c1 < c2); validate against COW state list
d$c1 <- d$dyadid %/% 1000; d$c2 <- d$dyadid %% 1000
cow <- unique(peacesciencer::cow_ddy$ccode1); nn <- !is.na(d$dyadid)
ok <- d$c1[nn] %in% cow & d$c2[nn] %in% cow; cat("NA dyadid:", sum(!nn), " share decodable to COW codes:", mean(ok), "\n")
print(head(sort(table(d$dyadid[nn][!ok]), decreasing = TRUE)))
stopifnot(mean(ok) > 0.98, all(d$c1[nn] < d$c2[nn]))

rhs <- c("lowdem","terr500k10yr","interaction4bjps","defense","caprat","contig","pnbdem_low","peaceyears","spline1","spline2","spline3")
f0 <- reformulate(rhs, "fatal")
stopifnot(isTRUE(all.equal(d$interaction4bjps, d$lowdem * d$terr500k10yr, check.attributes = FALSE)) ||
          cor(d$interaction4bjps, d$lowdem * d$terr500k10yr, use = "complete.obs") > .999)
terms <- c("interaction4bjps", "lowdem")

d0 <- d[complete.cases(d[, c(all.vars(f0), "dyadid")]), ]; n_orig <- nrow(d0); cat("n_orig", n_orig, "(log: 321041)\n")
rep <- fit_keys(f0, d0, terms, "dyadid", "logit", file.path(WO, "cache_rep.rds"))
PUB <- list(interaction4bjps = c(.0115043, 1), lowdem = c(-.0690097, 1))

ds <- common_sample(d, f0, "dyadid", "c1", "c2", "year"); rm(d, d0); gc()
M <- run_six(f0, ds, terms, "dyadid", WO)
# ---- sensitivity hook (split only): runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  # fit_keys without a cache file (run_six caches each model by file name)
  sens_split(id = "X009", dat = ds,
             fit = function(dat, add) { r <- fit_keys(if (length(add)) update(f0, as.formula(paste(". ~ . +", paste(add, collapse = " + ")))) else f0, dat, terms, "dyadid")$k$lowdem
               c(b = unname(r[["b"]]), se = unname(r[["se"]]), p = unname(r[["p"]])) })
}
EV <- list(interaction4bjps = "abstract: 'interaction between joint democracy and regional instability confirms' regime effects in peaceful regions",
           lowdem = "abstract: 'joint democracy is only statistically significant ... once territorial hot spots are jointly estimated'")
for (t in terms) {
  st <- repro(PUB[[t]][1], PUB[[t]][2] == 1, rep$k[[t]]["b"], rep$k[[t]]["p"])
  write_p4("X009", "DEMOCRACY", "KEY", t, EV[[t]], st, PUB[[t]][1], rep$k[[t]]["b"], "logit, SE clustered on dyadid",
    n_orig, M, "N (peace years + 3 splines; hot-spot count is regional territorial MIDs, not dyadic history)",
    paste0("Authors' log model 3 (paper main interaction model), fatal MID onset, years ", paste(range(ds$year), collapse = "-"),
           ". lowdem is a constituent of the interaction (its b is the effect at zero hot spots). dyadid decoded as c1*1000+c2. ",
           "n_orig=", n_orig, ", dropped ", n_orig - M$m0$n, "."))
}
sink()
