# VERIFY X016 Lektzian, Prins & Souva 2010 ISQ. Pre-committed key demternonrival (Table 4 Model 7) reproduces the PDF exactly
# (pub HR 0.825, z 0.23; ours HR 0.825) -> no spec change. PDF check: the abstract's positive claim ("substantially reduce the hazard ...
# outside the context of rivalry") is Table 5 row 4 = Table 4 Model 8, Joint Democracy*Non-Territory*Non-Rivalry HR 0.613 (z 1.30, one-tailed *).
# This script ALSO runs that term (demnonternonriv, Model 8: omitted = mixnonternonriv), same estimator/sample scheme as rerun_X016.R.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p4/common_p4_x8.R"))
source(file.path(CODE, "b1work/v/vhelp.R"))
suppressPackageStartupMessages(library(survival))
WO <- file.path(OUT, "work_p4/X016"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work_p4/X016"); sink(file.path(WO, "log_verify_X016.txt"), split = TRUE)
d <- read.delim(file.path(W, "3185586")); d <- subset(d, X_st == 1 & year > 1900)
CTRL <- "icowsal + recno5 + recmid5 + contig150 + lnrelcap + allies"
ALL8 <- c("demterival","demnonterival","demternonrival","demnonternonriv","mixterival","mixnonterival","mixternonrival","mixnonternonriv")
rhs <- paste(paste(setdiff(ALL8, "mixnonternonriv"), collapse = " + "), "+", CTRL); K <- "demnonternonriv"
fit_cox <- function(extra, s) {
  m <- coxph(as.formula(paste("Surv(X_t0, X_t, X_d) ~", rhs, extra, "+ frailty.gamma(ddyadidclaim2)")), data = s, ties = "breslow")
  sm <- summary(m)$coefficients; ct <- cbind(sm[, "coef"], sm[, "se(coef)"], sm[, "Chisq"], sm[, "p"]); rownames(ct) <- rownames(sm)
  grab(ct[!grepl("frailty", rownames(ct)), , drop = FALSE], K, m$n, paste("Model 8", extra))[[K]]
}
x <- d[complete.cases(d[, all.vars(as.formula(paste("~", rhs)))]), ]; r1 <- fit_cox("", x)
cat("full-sample repro HR:", exp(r1$b), "(pub 0.613) N =", r1$n, "\n")
d <- add_history(d, "ccode1", "ccode2", "year"); s <- subset(d, year >= 1822 & year <= 2011)
s <- s[complete.cases(s[, c(all.vars(as.formula(paste("~", rhs))), HISTV)]), ]
R <- run_six(fit_cox, s); stopifnot(length(unique(sapply(R, `[[`, "n"))) == 1)
vwrite("X016", "Lektzian, Prins & Souva 2010", "DEMOCRACY", K, "Table 4 Model 8 (Table 5 row 4)", log(0.613), log(0.613) / -1.30, 5152,
  "Cox PH, shared frailty on claim dyad, robust z", "none (paper's abstract term, run as extra row)", "N", "", R, pub_sig = FALSE,
  notes = sprintf("EXTRA term flagged from PDF: abstract's 'outside rivalry' pacifying claim = joint dem vs mixed, non-territorial non-rival claims. Published HR 0.613 z 1.30 (one-tailed p<.10 => two-tailed n.s. at .05); pub_se derived from z. Full-sample repro HR %.3f. Not a replacement for pre-committed demternonrival.", exp(r1$b)))
sink()
