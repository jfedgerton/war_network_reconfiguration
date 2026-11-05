# X021 Crescenzi & Kadera 2015 ISQ "Built to Last" (reanalysis of Gartzke & Weisiger 2014).
# Dyadic headline = Table 3 (authors' do crescenzi-kadera-table3.do; log 'replicating table 3.log'):
#   Model J (their corrected model, "Dropping Dyadic Difference"): logit deadlyl c_polave c_polavesq pcenerg diff1 demloi engypop logdist
#       cntgdumy allydumy capratio onemajor deadyrs deadyer1-3, cluster(dyadid)   -> key: demloi (dyadic democracy, weak link)
#   Model I (Democratic Community): logit deadlyl c_regstrength c_regstrengthsq pcenerg demloi engypop logdist cntgdumy allydumy capratio
#       onemajor deadyrs deadyer1-3, cluster(dyadid)   -> key: c_regstrength, c_regstrengthsq (global democratic community, inverted U)
# Key terms fixed from abstract BEFORE estimation: "dyadic democracy reduces conflict"; "materially strong global democratic community dampens".
# b_pub from authors' log (N=619,104 both). Authors have peace years + 3 splines (deadyrs, deadyer1-3).
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p4/P4_common.R"))
WO <- file.path(OUT, "work_p4/X021"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work_p4/X021"); sink(file.path(WO, "log_X021.txt"), split = TRUE)
V <- c("dyadid","year","deadlyl","c_polave","c_polavesq","pcenerg","diff1","demloi","engypop","logdist","cntgdumy","allydumy","capratio",
       "onemajor","deadyrs","deadyer1","deadyer2","deadyer3","c_regstrength","c_regstrengthsq")
d <- as.data.frame(haven::zap_labels(haven::read_dta(file.path(W, "table3.dta"), col_select = tidyselect::all_of(V))))
stopifnot(!anyDuplicated(d[, c("dyadid","year")]), all(d$deadlyl %in% c(0, 1, NA)))
d$c1 <- d$dyadid %/% 1000; d$c2 <- d$dyadid %% 1000; stopifnot(all(d$c2 > d$c1, na.rm = TRUE))
# Is deadlyl a lead (onset in t+1)? check deadyrs reset in t+1 after deadlyl==1 at t
d <- d[order(d$dyadid, d$year), ]
nx <- c(d$deadyrs[-1], NA); same <- c(d$dyadid[-1] == d$dyadid[-nrow(d)] & diff(d$year) == 1, FALSE)
cat("share deadyrs(t+1)==0 | deadlyl(t)==1:", mean(nx[same & d$deadlyl == 1] == 0, na.rm = TRUE),
    " | deadlyl(t)==0:", mean(nx[same & d$deadlyl == 0] == 0, na.rm = TRUE), "\n")
d$yr_out <- d$year + 1   # DV is fatal MID in t+1 (lead, "deadlyl"); history anchored on the outcome year (as C1 Gartzke convention)
sp <- "deadyrs + deadyer1 + deadyer2 + deadyer3"
fJ <- as.formula(paste("deadlyl ~ c_polave + c_polavesq + pcenerg + diff1 + demloi + engypop + logdist + cntgdumy + allydumy + capratio + onemajor +", sp))
fI <- as.formula(paste("deadlyl ~ c_regstrength + c_regstrengthsq + pcenerg + demloi + engypop + logdist + cntgdumy + allydumy + capratio + onemajor +", sp))
KJ <- c("demloi"); KI <- c("c_regstrength", "c_regstrengthsq", "demloi")
rJ <- fit_keys(fJ, d, KJ, "dyadid", file.path(WO, "cache_J_repro.rds"))
rI <- fit_keys(fI, d, KI, "dyadid", file.path(WO, "cache_I_repro.rds"))
pubJ <- c(demloi = -.0956186); pubI <- c(c_regstrength = .140625, c_regstrengthsq = -.0442074)
dsJ <- common_sample(d, fJ, "dyadid", "c1", "c2", "yr_out")
MJ <- run_six(fJ, dsJ, KJ, "dyadid", WO, "J_"); rm(dsJ); gc()
dsI <- common_sample(d, fI, "dyadid", "c1", "c2", "yr_out")
MI <- run_six(fI, dsI, KI, "dyadid", WO, "I_"); rm(dsI); gc()
# ---- sensitivity hook (split only): runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  # dsI was removed on the estimation line: rebuilt with the same call. fit_keys is called without a cache file (run_six caches by step name)
  sens_split(id = "X021", dat = common_sample(d, fI, "dyadid", "c1", "c2", "yr_out"),
             fit = function(dat, add) { r <- fit_keys(if (length(add)) update(fI, as.formula(paste(". ~ . +", paste(add, collapse = " + ")))) else fI, dat, KI, "dyadid")$k$c_regstrengthsq
               c(b = unname(r[["b"]]), se = unname(r[["se"]]), p = unname(r[["p"]])) })
}
sc <- "Crescenzi & Kadera 2015 ISQ"; est <- "logit, dyad-clustered SE (fixest feglm); authors used logit in place of G&W relogit"
hc <- "Y (fatal-MID peace years + 3 cubic splines only)"
base <- "Table 3 (fatal MID onset, DV = t+1 lead; history anchored on outcome year). "
write_p4("X021", sc, "DEMOCRACY", "KEY", "demloi", "abstract: 'an accurate interpretation ... actually supports the conclusion that dyadic democracy reduces conflict'",
  repro_status(pubJ[["demloi"]], TRUE, rJ$k$demloi), pubJ[["demloi"]], rJ$k$demloi[["b"]], est, rJ$n, MJ, hc,
  paste0(base, "Model J (drops dyadic difference). b_pub from authors' log p=.004. Common sample drops ", rJ$n - MJ$m0$n, " obs. bD NOT a placebo for this term."))
for (k in c("c_regstrength", "c_regstrengthsq"))
  write_p4("X021", sc, "DEMOCRACY", "KEY", k, "abstract: 'a materially strong global democratic community dampens the onset of militarized violence'",
    repro_status(pubI[[k]], TRUE, rI$k[[k]]), pubI[[k]], rI$k[[k]][["b"]], est, rI$n, MI, hc,
    paste0(base, "Model I (centered Democratic Community + square; systemic, constant within year). ", if (k == "c_regstrength") "Linear term (positive; inverted-U with square). " else "Squared term carries the pacifying (inverted-U) claim. ",
           "b_pub from authors' log p<.001. Common sample drops ", rI$n - MI$m0$n, " obs."))
sink()
