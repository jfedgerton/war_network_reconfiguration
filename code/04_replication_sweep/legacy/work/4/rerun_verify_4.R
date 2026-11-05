# VERIFY 4: Maoz 2009 AJPS -- match Table 4 col 1 (MID Initiation, directed politically relevant dyads, 1870-2001)
# Paper evidence (PDF p.234-236):
#  - Table 4 title: "Directed and Nondirected, Politically Relevant Dyads, 1870-2001"; col 1 "MID Initiation+" (+ = directed dyads)
#  - Col 1 controls: Capability Ratio, Initiator's Regime Score, Contiguity, Size of Initiator's SRG, Strategic Interdep.,
#    Economic Interdep., No. of years without conflict; N = 121,838. Pub strategic interdep b = -0.009 (0.015), n.s.
#  - p.234 "simple logit models with cubic splines and years of peace" -> noinityrs + initspline1-3
#  - p.234 "Contiguity. This is a binary measure" -> contig > 0 (dataset var is 0..7)
# Changes vs rerun_4.R: DV dichmid20 -> dichinit; sample polrel==1 & year>=1870; controls minreg302/distance ->
#  regimea02, contig binary, srga (caprat301 kept); peace years nomidyrs/midsplines -> noinityrs/initsplines.
# Unchanged (paper silent): lagged interdependence vars (lgstrintrdep/lgtrdintrdep), SEs clustered on dyad.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work/rerun_helpers_b118.R"))
source(file.path(CODE, "vb3/vhelp.R"))
library(haven); library(sandwich)
WO <- file.path(OUT, "work/4"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work/4"); logf <- file.path(WO, "log_verify_4.txt"); cat("", file = logf)
sink(logf, split = TRUE, append = TRUE)
dv <- "dichinit"; tie <- "lgstrintrdep"
base_rhs <- c("lgstrintrdep","lgtrdintrdep","caprat301","regimea02","contigb","srga","noinityrs","initspline1","initspline2","initspline3")
v <- c("year","statea","stateb","dyad","polrel","contig","dichinit","lgstrintrdep","lgtrdintrdep","strinterdep","trdinterdep",
       "caprat301","regimea02","srga","noinityrs","initspline1","initspline2","initspline3")
d <- as.data.frame(zap_labels(read_dta(file.path(ROOT, "_parent", "ajpsreplicationdyad.dta"), col_select = all_of(v))))
d <- d[!is.na(d$polrel) & d$polrel == 1 & d$year >= 1870 & d$year <= 2001, ]
d$contigb <- as.numeric(d$contig > 0)
cat("PR dyads 1870-2001 rows:", nrow(d), "\n"); stopifnot(nrow(d) > 1e5)
fitL <- function(rhs, dvv, data) {
  rownames(data) <- NULL
  m <- glm(as.formula(paste(dvv, "~", paste(rhs, collapse = "+"))), data = data, family = binomial("logit"),
           control = glm.control(epsilon = 1e-10, maxit = 100))
  stopifnot(m$converged)
  cl <- data$dyad[as.integer(rownames(model.frame(m)))]
  V <- vcovCL(m, cluster = cl, type = "HC0", cadjust = TRUE)
  b <- coef(m); s <- sqrt(diag(V))
  out <- list(coeftable = cbind(Estimate = b, SE = s, z = b / s, p = 2 * pnorm(-abs(b / s))), n = nobs(m))
  rm(m); invisible(gc()); out }
ts <- function(m, v) c(b = unname(m$coeftable[v, 1]), se = unname(m$coeftable[v, 2]), p = unname(m$coeftable[v, 4]))

cat("\n==== Matched headline, full PR sample ====\n")
m_rep <- fitL(base_rhs, dv, d); print(signif(m_rep$coeftable, 4)); cat("N =", m_rep$n, " (pub 121,838)\n")
cat("\n==== Sensitivity: raw contig (0-7) instead of binary ====\n")
m_s1 <- fitL(sub("contigb", "contig", base_rhs), dv, d); print(signif(m_s1$coeftable[tie, , drop = FALSE], 4))
cat("\n==== Sensitivity: unlagged strinterdep/trdinterdep ====\n")
m_s2 <- fitL(c("strinterdep","trdinterdep", base_rhs[-(1:2)]), dv, d); print(signif(m_s2$coeftable[c("strinterdep","trdinterdep"), ], 4)); cat("N =", m_s2$n, "\n")

d <- add_history(d, "statea", "stateb", "year"); d$year_h <- d$year
ds <- battery_sample(d, c(dv, base_rhs, "dyad")); rm(d); invisible(gc())
ECON <- list()
statf <- function(m) { ECON[[length(ECON) + 1]] <<- ts(m, "lgtrdintrdep"); ts(m, tie) }
R <- run_battery(ds, base_rhs, function(rhs, data) fitL(rhs, dv, data), statf, logf)
names(ECON) <- names(R); cat("lgtrdintrdep by model:\n"); print(do.call(rbind, ECON))

# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "4", dat = ds, y = dv, tie = tie, ctrl = setdiff(base_rhs, tie), fe = NULL,
           fit = function(dat, add) ts(fitL(c(base_rhs, add), dv, dat), tie))
}


mism <- "DV dichmid20 (MID underway) vs pub dichinit (MID initiation);sample all dyads 1816-2001 vs politically relevant 1870-2001;controls minreg302+distance vs initiator regime (regimea02)+contiguity+initiator SRG size (srga);peace years nomidyrs/midsplines vs noinityrs/initsplines"
unf <- "lag of interdependence vars not stated in paper (kept lagged);SE type not stated (kept dyad-clustered);N 1% below pub"
common_notes <- sprintf("Table 4 col1 all coefs incl. controls odd (e.g. contiguity -0.001**); matched full-PR N=%d vs pub 121,838. Sens.: raw contig b=%.3f; unlagged strinterdep b=%.3f. Econ interdep (co-key, pub -0.212**) m0 %.3f (p %.3g) -> mDM %.3f (p %.3g).",
  m_rep$n, m_s1$coeftable[tie, 1], m_s2$coeftable["strinterdep", 1], ECON$m0["b"], ECON$m0["p"], ECON$mDM["b"], ECON$mDM["p"])
write_verify("4", "Maoz 2009", "TIE", tie, "Table 4 col 1 (MID initiation, directed PR dyads)", -0.009, 0.015, 121838,
  "logit, peace years + cubic splines (SE type unstated)", mism, "PARTIAL", unf, FALSE, R, nrow(ds),
  paste("Matched-spec full-sample b =", signif(m_rep$coeftable[tie, 1], 4), "p =", signif(m_rep$coeftable[tie, 4], 3), ".", common_notes))
# co-key economic interdependence row (paper Table 4 col 1: -0.212 (0.063)**)
write_verify("4", "Maoz 2009", "TIE", "lgtrdintrdep", "Table 4 col 1 (MID initiation, directed PR dyads)", -0.212, 0.063, 121838,
  "logit, peace years + cubic splines (SE type unstated)", mism, "PARTIAL", unf, TRUE, ECON, nrow(ds),
  paste("Co-key (economic interdependence). Matched-spec full-sample b =", signif(m_rep$coeftable["lgtrdintrdep", 1], 4), ".", common_notes))
sink()
