# ADD_D11 Joo, Chan, ... et al. 2025 JCR "Right-Wing Populist Leaders, Nationalist Rhetoric, and Dispute Initiation in International
#   Politics" (10.1177/00220027241247041). PORT of the observational component from the authors' Dataverse deposit
#   (doi:10.7910/DVN/MW0PSJ: FP_data.dta, LargeNTables.do, Tab1.txt). The survey experiments are out of scope.
# Design (LargeNTables.do): leader-year panel (xtset ltid year), DV pinit = MID initiation, Stata xtlogit (random intercept by
#   leader-term ltid; default 12-point adaptive quadrature), vce(cluster ccode).
# HEADLINE (fixed before estimating): Table 1 Model 1, rightparti2 = right-wing populist leader (right2) x participatory democracy
#   (partipdem; product term, constituents in the model): 2.101 (robust SE 1.049), N = 11,656, 2,302 leader-terms (Tab1.txt).
#   Model 1 is the paper's first and most parsimonious specification (Models 2-5 add leader controls; Model 5 2.899 (1.402) N 9,711).
# Estimator here: random-intercept logit by ML with 40-node Gauss-Hermite quadrature (log-sum-exp), analytic scores, numerical
#   Hessian of the scores; Stata-style cluster-robust V on ccode (G/(G-1)). Cross-checked against lme4::glmer nAGQ = 12.
# Unit is the leader-year (monadic): history = the state's own MID onsets (add_history_monadic on ccode, year); mD is NA.
# Author history control: pinit_lag (lagged MID initiation) -> A_monadic.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R"))   # dplyr, fixest, add_history_monadic, common_sample, run_six, write_p8; set.seed(123)
set.seed(123)
suppressMessages(library(haven))
W <- file.path(OUT, "work_add20261110/ADD_D11"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_ADD_D11.txt"), split = TRUE)

# ---- 1. authors' leader-year data ----
d <- read_dta(file.path(ROOT, "_incoming/Joo2025/FP_data.dta")) %>% zap_labels() %>% as.data.frame()
stopifnot(nrow(d) == 16330, all(d$pinit %in% c(0, 1, NA)))
rhs <- c("pinit_lag", "rightparti2", "partipdem", "right2", "cinc", "majpower", "atopally")
vars <- c("pinit", rhs, "ltid", "ccode", "year")
KEY <- "rightparti2"
x0 <- d[complete.cases(d[, vars]), ]
stopifnot(max(abs(x0$rightparti2 - x0$right2 * x0$partipdem)) < 1e-6,                 # product term of its two constituents
          all(tapply(x0$ccode, x0$ltid, function(z) length(unique(z))) == 1))         # leader-terms nested in countries

# ---- 2. random-intercept logit, GH quadrature, cluster-robust on ccode ----
GH <- statmod::gauss.quad(40, kind = "hermite")
re_parts <- function(th, X, y, gid, G) {            # th = (beta incl. constant, lnsig2u)
  k <- ncol(X); b <- th[1:k]; sig <- exp(th[k + 1] / 2)
  u <- sqrt(2) * sig * GH$nodes; K <- length(u)
  eta <- drop(X %*% b) + matrix(u, length(y), K, byrow = TRUE)
  lo <- y * eta - log1p(exp(eta)); lo[!is.finite(lo)] <- (y * eta - eta)[!is.finite(lo)]
  lg <- rowsum(lo, gid, reorder = TRUE) + matrix(log(GH$weights / sqrt(pi)), G, K, byrow = TRUE)
  mx <- apply(lg, 1, max); lL <- mx + log(rowSums(exp(lg - mx)))
  om <- exp(lg - lL)                                  # posterior node weights per group
  res <- y - plogis(eta); omo <- om[gid, , drop = FALSE]
  sb <- rowSums(omo * res); st <- rowSums(omo * res * matrix(u / 2, length(y), K, byrow = TRUE))
  list(ll = lL, S = rowsum(cbind(X * sb, st), gid, reorder = TRUE))
}
fit_re <- function(rhs_v, data) {
  X <- cbind(as.matrix(data[, rhs_v, drop = FALSE]), `_cons` = 1); y <- data$pinit
  gid <- as.integer(factor(data$ltid)); G <- max(gid)
  cl <- tapply(data$ccode, gid, function(z) z[1])     # country of each leader-term (nested)
  sx <- apply(X, 2, sd); sx[sx == 0] <- 1; Xs <- sweep(X, 2, sx, "/")
  st <- c(coef(glm(y ~ Xs - 1, family = binomial)), 0)
  o <- nlminb(st, function(t) -sum(re_parts(t, Xs, y, gid, G)$ll), function(t) -colSums(re_parts(t, Xs, y, gid, G)$S),
              control = list(eval.max = 5000, iter.max = 3000, rel.tol = 1e-13))
  th <- o$par
  hfun <- function(t) { H <- numDeriv::jacobian(function(v) colSums(re_parts(v, Xs, y, gid, G)$S), t, method.args = list(r = 2)); (H + t(H)) / 2 }
  for (it in 1:10) { g <- colSums(re_parts(th, Xs, y, gid, G)$S); if (max(abs(g)) < 1e-7) break; th <- th + solve(-hfun(th), g) }
  P <- re_parts(th, Xs, y, gid, G); gmax <- max(abs(colSums(P$S))); stopifnot(gmax < 1e-4)
  Ainv <- solve(-hfun(th)); Sc <- rowsum(P$S, cl); NC <- nrow(Sc)
  V <- Ainv %*% crossprod(Sc) %*% Ainv * NC / (NC - 1)
  D <- c(1 / sx, 1); coef <- th * D; V <- V * outer(D, D)
  names(coef) <- c(colnames(X), "lnsig2u"); se <- sqrt(diag(V)); z <- coef / se
  cat("\n---- RE logit (40-node GH) | rhs:", paste(rhs_v, collapse = " + "), "| N =", length(y), "| groups =", G, "| clusters =", NC,
      "| LL =", round(sum(P$ll), 4), "| max|score|", signif(gmax, 3), "\n")
  print(round(cbind(b = coef, se = se, z = z), 6))
  out <- list(c(b = coef[[KEY]], se = se[[KEY]], p = 2 * pnorm(-abs(z[[KEY]])))); names(out) <- KEY
  attr(out, "n") <- length(y); attr(out, "coef") <- coef; out
}

# ---- 3. reproduction on the published sample ----
cat("complete-case N:", nrow(x0), "(published 11,656; 2,302 leader-terms)\n"); stopifnot(nrow(x0) == 11656, length(unique(x0$ltid)) == 2302)
full <- fit_re(rhs, x0)
b_pub <- 2.101; se_pub <- 1.049
cat("REPRO", KEY, "b", full[[KEY]][["b"]], "(pub", b_pub, ") se", full[[KEY]][["se"]], "(pub", se_pub, ") lnsig2u",
    attr(full, "coef")[["lnsig2u"]], "(pub 0.173)\n")
stopifnot(abs(full[[KEY]][["b"]] - b_pub) / abs(b_pub) <= 0.10, abs(full[[KEY]][["se"]] - se_pub) / se_pub <= 0.10)
# cross-check the quadrature against lme4 adaptive quadrature (12 points, the Stata default)
g12 <- lme4::glmer(as.formula(paste("pinit ~", paste(rhs, collapse = " + "), "+ (1 | ltid)")), data = x0, family = binomial, nAGQ = 12,
                   control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))
cat("glmer nAGQ=12", KEY, lme4::fixef(g12)[[KEY]], " vs GH40", full[[KEY]][["b"]], "\n")
stopifnot(abs(lme4::fixef(g12)[[KEY]] - full[[KEY]][["b"]]) < 1e-3 * max(1, abs(full[[KEY]][["b"]])))

# ---- 4. history window + six models on one common sample (monadic) ----
fitter <- function(r, data) fit_re(r, data)
ds <- common_sample(d, vars, "ccode", NULL, "year", monadic = TRUE)
M <- run_six(rhs, ds, fitter, monadic = TRUE)
saveRDS(list(M = M, full = full), file.path(W, "M_ADD_D11.rds"))

# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
# headline is a product term (split class), but the design is monadic: sens_helper.R needs a dyadic-history column, which is not
# defined for leader-years. Like the other monadic studies (P9_075, P9_094) this id is NOT proposed for sens_studies.csv; the hook
# stops loudly instead of running a split with a substituted column.
if (Sys.getenv("SWEEP_SENS") == "1") {
  stop("ADD_D11 is monadic (leader-years): sens_helper.R requires hist_dyad_w15; not in sens_studies.csv")
}

meta <- list(pid = "ADD_D11", cite = "ADD_D11", journal = "JCR", year = 2025, doi = "10.1177/00220027241247041", family = "DEMOCRACY")
nt <- sprintf(paste0("PORT of Dataverse doi:10.7910/DVN/MW0PSJ (observational component). Table 1 Model 1 xtlogit RE by leader-term, SE clustered on ccode. ",
  "Full-sample b=%.6f se=%.6f N=%d vs published 2.101 (1.049) N 11,656. Monadic leader-year: mD not defined; mDM = + own w15 history; m26 = + own w26. ",
  "Common 1822-2011 sample N=%d."), full[[KEY]][["b"]], full[[KEY]][["se"]], attr(full, "n"), nrow(ds))
write_p8(meta, KEY, "Title/Table 1: right-wing populist leaders in participatory democracies initiate more MIDs (right2 x partipdem)",
  "Table 1 Model 1 (xtlogit RE, leader-years 1887-2012)", "random-effects logit (xtlogit, re; 12-pt adaptive quadrature), SE clustered on ccode",
  "RE logit ML, 40-node Gauss-Hermite (R port), cluster-robust on ccode, all six steps (monadic)", b_pub, se_pub, 11656, TRUE, FALSE, M,
  "Y (pinit_lag, lagged MID initiation)", nt)
dir.create(file.path(OUT, "results_add"), recursive = TRUE, showWarnings = FALSE)
stopifnot(file.rename(file.path(OUT, "results_p8", "ADD_D11__rightparti2.csv"), file.path(OUT, "results_add", "ADD_D11__rightparti2.csv")))
sink()
