# P8_017 Huth & Allee 2002 JCR "Domestic Political Accountability and the Escalation and Settlement of International Disputes"
# (10.1177/002200202237928). Data: Harvard Dataverse doi:10.7910/DVN/XEGO5E, huth_allee_statusquodata (original .dta) + codebook.
# No do-file posted; specification rebuilt from the PDF (Table 1, p. 769) and the codebook variable definitions.
# HEADLINE (conflict-initiation stage; negotiation/settlement stages out of scope per phase-8 screen):
#   Table 1, multinomial logit of challenger decision dvsq1 (0 = no action [base], 1 = call for talks, 2 = threat/use of force),
#   "Political accountability model - comparing differences across regimes", robust (Huber) SEs, N = 6,542 challenger-target
#   decision points (month-specific, 12-month rule for no-action cases), 1919-1995.
#   Column "Force vs. No Action": Challenger democracy level = -.056 (.012), one-tailed p < .01.
# Regressors (Table 1 rows -> codebook names):
#   Challenger democracy level -> netdemoc; Target democracy level -> onetdem; Challenger Democracy x Stalemate -> ntdmrev1
#   (netdem21 x stale1yr); Control for recent stalemate -> stale1yr; Challenger Democracy x Ethnic Ties -> ethvntdm (ethvalue x netdemoc);
#   Control for ethnic ties -> ethvalue; Common security ties -> comopp (Table 4 note: "common TD opponent"); Strategic value -> strvalue;
#   Challenger/Target involved in other dispute -> chgodi/tgtodi; Military balance -> milratio; Months since last challenge -> msafinal.
#   Choices where the PDF is ambiguous (stalemate window 1yr vs 2yr/5yr; comopp vs alliance) are fixed a priori to the first/closest
#   codebook match and NOT tuned; the full-sample fit is compared with the whole Table 1 column as a spec check.
# KEY TERM (a priori, abstract: "the importance of electoral cycles and the strength of opposition parties"; the across-regimes model is
#   the paper's democratic-peace test, "democracies are less likely to use force (hypotheses 1a and 5)"): netdemoc in the Force-vs-No-Action
#   equation. Reported in results as key_term "accountability" (row name fixed in phase 8) = challenger democracy level, force equation.
# Estimator: multinomial logit by Newton-Raphson (exact MLE), Huber sandwich SE with Stata's N/(N-1) factor; same in all six steps.
# History keyed to the calendar year of the decision month (datetm is Stata %tm: months since 1960m1). Directed challenger-target pair.
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/p8_common.R"))
WD <- file.path(ROOT, "work_p8/P8_017"); W <- file.path(OUT, "work_p8/P8_017"); dir.create(W, recursive = TRUE, showWarnings = FALSE); sink(file.path(W, "log_P8_017.txt"), split = TRUE)
rhs0 <- c("netdemoc", "onetdem", "ntdmrev1", "stale1yr", "ethvntdm", "ethvalue", "comopp", "strvalue", "chgodi", "tgtodi", "milratio", "msafinal")
d <- haven::read_dta(file.path(WD, "statusquo.dta"), col_select = all_of(c("challenger", "target", "datetm", "dvsq1", rhs0)))
d <- as.data.frame(haven::zap_labels(d))
d$year <- 1960 + floor(d$datetm / 12)
cat("rows", nrow(d), "years", range(d$year), "\n"); print(table(d$dvsq1))
stopifnot(nrow(d) == 6542, all(d$dvsq1 %in% 0:2), min(d$year) == 1919, max(d$year) == 1995)

# multinomial logit, base = 0; returns list(coef = p x 2 matrix [talks, force], V = robust vcov of stacked (talks, force))
mnl <- function(y, X, tol = 1e-12, maxit = 100) {
  n <- nrow(X); p <- ncol(X); Y <- cbind(y == 1, y == 2) * 1
  b <- rep(0, 2 * p)
  for (it in 1:maxit) {
    B <- matrix(b, p, 2); E <- exp(X %*% B); den <- 1 + rowSums(E); P <- E / den
    g <- c(crossprod(X, Y[, 1] - P[, 1]), crossprod(X, Y[, 2] - P[, 2]))
    H11 <- -crossprod(X * (P[, 1] * (1 - P[, 1])), X); H22 <- -crossprod(X * (P[, 2] * (1 - P[, 2])), X)
    H12 <- crossprod(X * (P[, 1] * P[, 2]), X)
    H <- rbind(cbind(H11, H12), cbind(t(H12), H22))
    step <- solve(H, g); b <- b - step
    if (max(abs(step)) < tol) break
  }
  stopifnot(it < maxit)
  B <- matrix(b, p, 2); E <- exp(X %*% B); P <- E / (1 + rowSums(E))
  S <- cbind(X * (Y[, 1] - P[, 1]), X * (Y[, 2] - P[, 2]))
  Hi <- solve(-H); V <- Hi %*% crossprod(S) %*% Hi * n / (n - 1)
  ll <- sum(log(cbind(1 - rowSums(P), P)[cbind(1:n, y + 1)]))
  list(coef = B, V = V, ll = ll, n = n, p = p, it = it)
}
fit_mnl <- function(rhs, data) {
  X <- model.matrix(as.formula(paste("~", paste(rhs, collapse = " + "))), data); y <- data$dvsq1
  r <- mnl(y, X); p <- ncol(X); se <- sqrt(diag(r$V))
  tab <- data.frame(term = colnames(X), b_talks = r$coef[, 1], se_talks = se[1:p], b_force = r$coef[, 2], se_force = se[(p + 1):(2 * p)])
  # force vs talks (difference of equations), robust SE from the stacked vcov
  dvar <- diag(r$V)[1:p] + diag(r$V)[(p + 1):(2 * p)] - 2 * diag(r$V[1:p, (p + 1):(2 * p)])
  tab$b_fvt <- tab$b_force - tab$b_talks; tab$se_fvt <- sqrt(dvar)
  cat("\n---- MNL dvsq1 ~", paste(rhs, collapse = " + "), "| N =", r$n, "| logLik =", round(r$ll, 2), "| NR iters =", r$it, "\n")
  print(cbind(tab[1], round(tab[-1], 4)), row.names = FALSE)
  k <- which(colnames(X) == "netdemoc"); b <- unname(r$coef[k, 2]); s <- unname(se[p + k])
  out <- list(accountability = c(b = b, se = s, p = 2 * pnorm(-abs(b / s))))
  attr(out, "n") <- r$n; out
}
vars <- c("dvsq1", rhs0)
# 1. authors' full sample (published N = 6,542; published logLik = -4796.27)
r_full <- fit_mnl(rhs0, d); n_orig <- attr(r_full, "n")
# 2. common 1822-2011 sample + history (all decision points are 1919-1995)
ds <- common_sample(d, vars, "challenger", "target", "year")
M <- run_six(rhs0, ds, fit_mnl)
saveRDS(list(M = M, r_full = r_full, n_orig = n_orig, n_used = nrow(ds)), file.path(W, "M_P8_017.rds"))
meta <- list(pid = "P8_017", short_cite = "Huth & Allee 2002", journal = "JCR", year = 2002, doi = "10.1177/002200202237928",
  family = "DEMOCRACY", headline_table_model = "Table 1 (MNL, across-regimes accountability model), Force vs. No Action column",
  article_estimator = "multinomial logit (base = no action), Huber robust SE; challenger-target decision points 1919-95 (12-month rule)",
  est_used_all_steps = "multinomial logit by Newton-Raphson, Huber sandwich SE x N/(N-1), all six steps",
  author_has_history_control = "N (months since last challenge only)")
write_p8(meta, "accountability", "abstract: 'importance of electoral cycles and the strength of opposition parties'; H1a/H5: democracies less likely to use force",
  -0.056, 0.012, 6542, NA, FALSE, M,
  paste0("Key term = challenger democracy level (netdemoc) in the Force-vs-No-Action equation of Table 1 (pub -.056 (.012), one-tailed p<.01). ",
         "Data: Huth & Allee status-quo file (Dataverse XEGO5E); spec rebuilt from Table 1 + codebook (no do-file). Full-sample N=", n_orig,
         " (pub 6,542). Common sample N=", nrow(ds), ". History keyed to the calendar year of the decision month; directed challenger-target pair."),
  full = r_full)
sink()
