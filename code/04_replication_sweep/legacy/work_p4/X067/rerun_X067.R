# X067 Benson, Bentley & Ray 2013 JPR "Ally provocateur" (10.1177/0022343312454445)
# Data + do-file + log: PRIO JPR archive (Benson, Bentley and Ray replication.zip).
# Headline: Table IV Model 1C (full model), logit, SE clustered on directed dyad (ddyad), directed dyad-years 1816-2000.
# KEY (fixed a priori from abstract H1: "revisionist states with unconditional commitments ... more likely to initiate"):
#   product term idgNcrevA (= idgNcompel x revstateA). CONTROL: jointdem (DEMOCRACY).
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p4/b3/common_b3.R"))
WO <- file.path(OUT, "work_p4/X067"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work_p4/X067"); sink(file.path(WO, "log_X067.txt"), split = TRUE)
m <- list(xid = "X067", short_cite = "Benson, Bentley & Ray 2013", journal = "JPR", year = 2013, doi = "10.1177/0022343312454445")

d <- haven::read_dta(file.path(W, "Ally Provocateur Data/ap_data.dta"))
d <- as.data.frame(lapply(d, function(x) as.numeric(x)))
stopifnot(nrow(d) == 188326)
rhs <- "jointdem + cont + caprat2 + sglo + idgNcompel + idgNcrevA + idiNcompel + idiNcrevA + revstateA + dgNcompel + diNcompel + peaceyears + peace2 + peace3"
f0 <- as.formula(paste("mzinit ~", rhs))
terms <- c("idgNcrevA", "jointdem")

# Reproduction on full author sample
rep <- fit_keys(f0, d, terms, "~ddyad", "logit", file.path(WO, "cache_repro.rds"))
stopifnot(rep$n == 188326)
b_pub <- c(idgNcrevA = 0.8655306, jointdem = -0.3038153)   # author log (ap_log.smcl), Model 1C
sig_pub <- c(idgNcrevA = TRUE, jointdem = TRUE)              # p = .012, .004
st <- sapply(terms, function(t) repro_status(b_pub[t], sig_pub[t], rep$k[[t]]))
print(st)

ds <- common_sample(d, f0, "ddyad", "ccode1", "ccode2", "year")
M <- run_six(f0, ds, terms, "~ddyad", "logit", WO)
# ---- sensitivity hook (split only): runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  # fit_keys is called without a cache file (run_six caches by step name)
  sens_split(id = "X067", dat = ds,
             fit = function(dat, add) { r <- fit_keys(if (length(add)) update(f0, as.formula(paste(". ~ . +", paste(add, collapse = " + ")))) else f0, dat, terms, "~ddyad", "logit")$k$idgNcrevA
               c(b = unname(r[["b"]]), se = unname(r[["se"]]), p = unname(r[["p"]])) })
}
est <- "logit, SE clustered on directed dyad (fixest feglm)"
hc <- "Y (peace years + squared + cubed terms)"
nt <- paste0("Table IV Model 1C; b_pub from authors' log (ap_log.smcl); repro N=", rep$n, "; common sample drops ", attr(ds, "n_cc") - nrow(ds),
             " obs (years 1816-1821 + benchmark missingness). Directed dyad-years 1816-2000, Maoz MID initiation. ")
write_p4(m, "TIE", "KEY", "idgNcrevA", "abstract H1: 'revisionist states with unconditional commitments ... more likely to initiate'",
         st["idgNcrevA"], b_pub["idgNcrevA"], rep$k$idgNcrevA["b"], est, rep$n, M, hc,
         paste0(nt, "Key = product term idgNcompel x revstateA (initiator has 'dg' commitment w/o compellent agreements x revisionist). Mapping of 'dg' to 'unconditional' inferred from variable labels + H1/H2 pattern (idg lincom sig, idi not); idiNcrevA not reported."))
write_p4(m, "DEMOCRACY", "CONTROL", "jointdem", "headline-model control (polity both > 6)",
         st["jointdem"], b_pub["jointdem"], rep$k$jointdem["b"], est, rep$n, M, hc,
         paste0(nt, "Joint-democracy control. bD (lagged joint democracy) is NOT a placebo for this term."))
sink()
