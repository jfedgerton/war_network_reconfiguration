# P9_016 Johnson 2016 CMPS "Alliance Treaty Obligations and War Intervention" (10.1177/0738894215577557)
# Data + do-file + log: Harvard Dataverse doi:10.7910/DVN/JUUU8W (johnson_rep.dta, johnson_rep.do, johnson_log.smcl).
# HEADLINE: Table 2 Model 2 (do-file "Table 2", second model, full controls):
#   mlogit intervention offense defense ln_distance icontr tcontr i_win i_demdiff t_demdiff, tech(bfgs)   [default OIM SE]
#   Outcomes 0 Remain_Neutral (base), 1 Join_Initiator, 2 Join_Target. Unit = potential intervener x COW war (95 wars), year = war start.
#   Published (author's log): N = 7,207; Join_Initiator eq offense 5.675436 (.8947592); Join_Target eq defense 3.297405 (.2923497).
# KEY TERM (a priori, abstract + P8_017 precedent = one term in one equation): the abstract states first that "defense pacts are
#   expected to increase the likelihood of members joining the target of a war" and reports the hypothesized effects were found
#   -> headline = defense in the Join_Target equation (key_term tgt_defense). Non-headline row: offense in the Join_Initiator
#   equation (key_term ini_offense).
# Estimator (precedent work_p8/P8_017/rerun_P8_017.R): multinomial logit by Newton-Raphson (exact MLE), base = 0, written locally as in
#   P8_017; SE = OIM (inverse negative Hessian at the optimum), as Stata's mlogit default; same in all six steps.
#   offense is quasi-separated in the Join_Target equation (no offense = 1 row joins the target; Stata stops at -8.46, SE 345):
#   NR stops on Stata's ml rule (relative logLik change < 1e-7 and g'(-H)^-1 g < 1e-5), so separated coefficients stay finite and
#   the OIM is invertible, as in Stata; they are not reported.
#   History terms enter every equation additively (as in P8_017 / study 52).
# History: each row has two dyads (intervener-initiator, intervener-target). History is attached to the dyad of the tie being
#   reported: tgt_defense -> (p_intervener, target); ini_offense -> (p_intervener, initiator). Each term gets its own six-step run.
#   Keyed to the war start year.
# No randomness (seed 123 set in common_b3.R).
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R"))
WD <- file.path(ROOT, "work_p9/P9_016"); W <- file.path(OUT, "work_p9/P9_016"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_016.txt"), split = TRUE)
rhs0 <- c("offense", "defense", "ln_distance", "icontr", "tcontr", "i_win", "i_demdiff", "t_demdiff")
d <- haven::read_dta(file.path(WD, "johnson_rep.dta"), col_select = all_of(c("warnum", "year", "p_intervener", "initiator", "target", "intervention", rhs0)))
d <- as.data.frame(haven::zap_labels(d))
cat("rows", nrow(d), "years", range(d$year), "wars", length(unique(d$warnum)), "\n"); print(table(d$intervention))
stopifnot(nrow(d) == 7829, all(d$intervention %in% 0:2), length(unique(d$warnum)) == 95, all(d$p_intervener != d$target), all(d$p_intervener != d$initiator))

# multinomial logit (J = 2 non-base outcomes), base = 0; Newton-Raphson with step halving; OIM vcov of stacked (eq1, eq2)
mnl_oim <- function(y, X, ltol = 1e-7, nrtol = 1e-5, maxit = 200) {
  n <- nrow(X); p <- ncol(X); Y <- cbind(y == 1, y == 2) * 1
  llf <- function(B) { E <- cbind(0, X %*% B); mx <- apply(E, 1, max); sum(rowSums(cbind(1 - rowSums(Y), Y) * E) - (mx + log(rowSums(exp(E - mx))))) }
  prb <- function(B) { E <- cbind(0, X %*% B); E <- exp(E - apply(E, 1, max)); P <- E / rowSums(E); P[, 2:3] }
  hess <- function(P) {
    H11 <- -crossprod(X * (P[, 1] * (1 - P[, 1])), X); H22 <- -crossprod(X * (P[, 2] * (1 - P[, 2])), X)
    H12 <- crossprod(X * (P[, 1] * P[, 2]), X); rbind(cbind(H11, H12), cbind(t(H12), H22))
  }
  b <- rep(0, 2 * p); ll <- llf(matrix(b, p, 2))
  for (it in 1:maxit) {
    B <- matrix(b, p, 2); P <- prb(B)
    g <- c(crossprod(X, Y[, 1] - P[, 1]), crossprod(X, Y[, 2] - P[, 2]))
    step <- solve(hess(P), g)
    t <- 1; repeat { bn <- b - t * step; lln <- llf(matrix(bn, p, 2)); if ((is.finite(lln) && lln >= ll - 1e-12) || t < 1e-6) break; t <- t / 2 }
    # Stata ml convergence rule: relative logLik change < ltolerance(1e-7) AND scaled gradient g'(-H)^-1 g < nrtolerance(1e-5)
    Pn <- prb(matrix(bn, p, 2)); gn <- c(crossprod(X, Y[, 1] - Pn[, 1]), crossprod(X, Y[, 2] - Pn[, 2]))
    nr <- sum(gn * solve(-hess(Pn), gn))
    conv <- abs(lln - ll) / abs(ll) < ltol && nr < nrtol; b <- bn; ll <- lln
    if (conv) break
  }
  stopifnot(it < maxit)
  B <- matrix(b, p, 2); V <- solve(-hess(prb(B)))
  list(coef = B, V = V, ll = ll, n = n, p = p, it = it)
}
K <- c(tgt_defense = "defense", ini_offense = "offense"); EQ <- c(tgt_defense = 2, ini_offense = 1)
fit_mnl <- function(rhs, data) {
  stopifnot(all(0:2 %in% data$intervention))
  X <- model.matrix(as.formula(paste("~", paste(rhs, collapse = " + "))), data)
  r <- mnl_oim(data$intervention, X); p <- r$p; se <- sqrt(diag(r$V))
  tab <- data.frame(term = colnames(X), b_ini = r$coef[, 1], se_ini = se[1:p], b_tgt = r$coef[, 2], se_tgt = se[(p + 1):(2 * p)])
  cat("\n---- MNL intervention ~", paste(rhs, collapse = " + "), "| N =", r$n, "| logLik =", round(r$ll, 5), "| NR iters =", r$it,
      "| (offense in Join_Target eq is quasi-separated)\n")
  print(cbind(tab[1], signif(tab[-1], 7)), row.names = FALSE)
  out <- lapply(names(K), function(k) { j <- which(colnames(X) == K[[k]]); bb <- r$coef[j, EQ[[k]]]; s <- se[(EQ[[k]] - 1) * p + j]
    stopifnot(is.finite(bb), is.finite(s)); c(b = unname(bb), se = unname(s), p = unname(2 * pnorm(-abs(bb / s)))) })
  names(out) <- names(K); attr(out, "n") <- r$n; out
}
vars <- c("intervention", rhs0, "p_intervener", "initiator", "target", "year")
full <- d[complete.cases(d[, vars]), ]
PUB <- list(tgt_defense = c(b = 3.297405, se = .2923497), ini_offense = c(b = 5.675436, se = .8947592))
# 1. authors' full sample (published N = 7,207, logLik = -598.72466)
r_full <- fit_mnl(rhs0, full); n_orig <- attr(r_full, "n")
for (k in names(PUB)) {
  rd <- abs(r_full[[k]][["b"]] - PUB[[k]][["b"]]) / abs(PUB[[k]][["b"]])
  cat("REPRO vs author log:", k, "b", r_full[[k]][["b"]], "(pub", PUB[[k]][["b"]], ") se", r_full[[k]][["se"]], "(pub", PUB[[k]][["se"]], ") rel.diff b", rd, "\n")
  stopifnot(rd <= 0.01)
}
stopifnot(n_orig == 7207)
# 2. common 1822-2011 sample + history attached to the tie's dyad; one six-step run per term
SIDE <- c(tgt_defense = "target", ini_offense = "initiator")
R <- list()
for (k in names(K)) {
  cat("\n==== term", k, ": history on (p_intervener,", SIDE[[k]], ") ====\n")
  ds <- common_sample(d, vars, "p_intervener", SIDE[[k]], "year")
  M <- run_six(rhs0, ds, fit_mnl)
  R[[k]] <- list(M = M, n_used = nrow(ds))
}
saveRDS(list(R = R, r_full = r_full, n_orig = n_orig), file.path(W, "M_P9_016.rds"))
# ---- result-CSV write step (write_p8 writes to results_p8; the file is moved to results_p9 unchanged) ----
meta <- list(pid = "P9_016", cite = "Johnson 2016", journal = "CMPS", year = 2016, doi = "10.1177/0738894215577557", family = "TIE")
ev <- "abstract: 'defense pacts are expected to increase the likelihood of members joining the target'"
art <- "multinomial logit (remain neutral [base] / join initiator / join target), OIM SE (tech(bfgs)); potential intervener x war, 95 COW wars"
est <- "multinomial logit by Newton-Raphson, OIM SE, all six steps; history on the tie's dyad (intervener-target or intervener-initiator)"
ahc <- "N"
mv <- function(term) {
  src <- file.path(OUT, "results_p8", paste0("P9_016__", term, ".csv")); dst <- file.path(OUT, "results_p9", paste0("P9_016__", term, ".csv"))
  dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE); stopifnot(file.exists(src), file.rename(src, dst))
}
nt <- function(k) sprintf(paste0("Equation-specific MNL coefficient (%s, %s equation). History attached to the (p_intervener, %s) dyad at the war start year. ",
  "Authors' Dataverse data + do-file; published values from the author's log (Table 2 Model 2). Full-sample repro b=%.6f se=%.6f N=%d. Common 1822-2011 sample N=%d. ",
  "offense in Join_Target eq quasi-separated (not reported)."),
  K[[k]], if (EQ[[k]] == 1) "Join_Initiator" else "Join_Target", SIDE[[k]], r_full[[k]][["b"]], r_full[[k]][["se"]], n_orig, R[[k]]$n_used)
write_p8(meta, "tgt_defense", ev, "Table 2 Model 2 (MNL), Join_Target equation: defense", art, est,
  3.297405, .2923497, 7207, TRUE, FALSE, R$tgt_defense$M, ahc, nt("tgt_defense")); mv("tgt_defense")
write_p8(meta, "ini_offense", ev, "Table 2 Model 2 (MNL), Join_Initiator equation: offense", art, est,
  5.675436, .8947592, 7207, TRUE, FALSE, R$ini_offense$M, ahc, nt("ini_offense")); mv("ini_offense")
sink()
