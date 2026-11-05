# X062 Leeds & Johnson 2017 JOP, "Theory, Data, and Deterrence: A Response to Kenwick, Vasquez, and Powers"
# Key term (decided from abstract BEFORE estimating): target has a defense pact (ptargdef).
#   Abstract: "claims that defense pacts deter conflict are robust to many of these changes in research design"
# Headline = the JL model the reply defends (Table A1; LeedsJohnsonJOPrep.do first estsimp): probit dispute (MID initiation),
#   directed dyad-years 1816-2000, conventional (non-clustered) SEs as in the authors' command.
# b_pub from authors' log (LeedsJohnsonJOPrep.smcl): ptargdef -.0615538 (p<.001); jdem -.1545372 (p<.001); N = 1,077,992.
# jdem (joint democracy) is recorded as a CONTROL row. Same spec/data as cid 109 (Johnson & Leeds 2011 FPA) in the main sweep.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p4/b6/common_b6.R"))
WO <- file.path(OUT, "work_p4/X062"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work_p4/X062"); sink(file.path(WO, "log_X062.txt"), split = TRUE)
m <- list(xid = "X062", short_cite = "Leeds & Johnson 2017", journal = "JOP", year = 2017, doi = "10.1086/687285")
d <- haven::read_dta(file.path(W, "LJ.dta"), col_select = c(challenger, target, ddyad, year, dispute, ptargdef, pchalally, pchaloff,
  pchalneu, ln_distance, capprop, jdem, s_un_glo, peaceyrs, peaceyrs2, peaceyrs3))
d <- as.data.frame(haven::zap_labels(d)); gc()
cat("rows:", nrow(d), " years:", range(d$year), " dups:", sum(duplicated(d[, c("challenger", "target", "year")])), "\n")
stopifnot(all(d$dispute %in% c(0, 1, NA)))
f0 <- dispute ~ ptargdef + pchalally + pchaloff + pchalneu + ln_distance + capprop + jdem + s_un_glo + peaceyrs + peaceyrs2 + peaceyrs3
K <- c("ptargdef", "jdem")
r <- fit_keys(f0, d, K, "iid", "probit", file.path(WO, "cache_repro.rds"))
bp <- c(ptargdef = -0.0615538, jdem = -0.1545372)
st <- sapply(K, function(k) repro_status(bp[[k]], TRUE, r$k[[k]])); print(st)
ds <- common_sample(d, f0, "ddyad", "challenger", "target", "year"); rm(d); gc()
M <- run_six(f0, ds, K, "iid", "probit", WO)
# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "X062", dat = ds, y = "dispute", tie = "ptargdef",
           ctrl = c("pchalally", "pchaloff", "pchalneu", "ln_distance", "capprop", "jdem", "s_un_glo", "peaceyrs", "peaceyrs2", "peaceyrs3"), fe = NULL,
           fit = function(dat, add) { r <- fit_keys(if (length(add)) update(f0, as.formula(paste(". ~ . +", paste(add, collapse = " + ")))) else f0, dat, "ptargdef", "iid", "probit")$k$ptargdef
             c(b = unname(r[["b"]]), se = unname(r[["se"]]), p = unname(r[["p"]])) })
}
nt <- paste0("Table A1 JL model (LeedsJohnsonJOPrep.do), b_pub from authors' smcl log (N=1,077,992). Probit with conventional SEs as in authors' estsimp call. ",
  "Sample restricted to 1822-2000 for GML history; dropped vs original: ", r$n - M$m0$n, ". Same data/spec as main-sweep cid 109 (Johnson & Leeds 2011).")
write_p4(m, "TIE", "KEY", "ptargdef", "abstract: 'claims that defense pacts deter conflict are robust to many of these changes'",
  st[["ptargdef"]], bp[["ptargdef"]], unname(r$k$ptargdef["b"]), "probit, iid SE (directed dyad-year MID initiation)", r$n, M,
  "Y (peaceyrs, peaceyrs^2, peaceyrs^3)", nt)
write_p4(m, "DEMOCRACY", "CONTROL", "jdem", "control in headline JL model (joint democracy)",
  st[["jdem"]], bp[["jdem"]], unname(r$k$jdem["b"]), "probit, iid SE (directed dyad-year MID initiation)", r$n, M,
  "Y (peaceyrs, peaceyrs^2, peaceyrs^3)", paste(nt, "bD (lagged joint democracy) is not a placebo for this term."))
sink()
