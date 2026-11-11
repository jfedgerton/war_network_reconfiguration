# ADD_D10 Zellman, Wajner, Lutmar & Ozdamar 2026 ISQ "Assessing Populist Belligerence: Evidence of Populism and Democratic
#   Governance's Interaction on Interstate Violence" (10.1093/isq/sqag085). PORT of the authors' Dataverse deposit
#   (doi:10.7910/DVN/IXWNFI: Assessing-Populist-Belligerence-1970-2015.dta, .do, .log).
# Design (do-file + log): directed politically relevant dyad-years 1970-2014, DV hihosta_3 = highest hostility of StateA toward
#   StateB (0 none, 1 threat/display, 2 use of force, 3 war; Dyadic MID 4.02), Stata ologit, year dummies, vce(cluster dyad).
# HEADLINE (fixed before estimating; title = the populism x democratic-governance interaction): Table 1 Model 2,
#   statea_populist_strict # c.statea_v2x_polyarchy = -1.534721 (robust SE 0.6288662, p .015), N = 101,468 (log lines 159-182).
#   Model 3's interaction (edem01 x continuous populism) is -0.072 (0.459), n.s., so Model 2 carries the interaction claim.
# Estimator here: ordered logit by ML (analytic scores, numerical Hessian of the scores), Stata-style cluster-robust V on dyad
#   (G/(G-1)), year dummies, all six steps. History/benchmarks: directed dyad -> add_history(statea, stateb, year).
# Author history control: none (no peace years, no prior-MID terms) -> C_none.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R"))   # dplyr, fixest, add_history, common_sample, run_six, write_p8; set.seed(123)
set.seed(123)
suppressMessages(library(haven))
W <- file.path(OUT, "work_add20261110/ADD_D10"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_ADD_D10.txt"), split = TRUE)

# ---- 1. authors' analysis data ----
d <- read_dta(file.path(ROOT, "_incoming/Zellman2026/Assessing-Populist-Belligerence-1970-2015.dta")) %>% zap_labels() %>% as.data.frame()
stopifnot(nrow(d) == 138046, anyDuplicated(d[, c("statea", "stateb", "year")]) == 0, all(d$hihosta_3 %in% 0:3))
d$pop_x_poly <- d$statea_populist_strict * d$statea_v2x_polyarchy
rhs <- c("statea_populist_strict", "statea_v2x_polyarchy", "pop_x_poly", "statea_gdp_pcap_2015USD", "statea_inflation_wb",
         "stateb_v2x_polyarchy", "statea_cinc_share", "s_un_atop", "statea_tradebalance", "conttype")
vars <- c("hihosta_3", rhs, "year", "dyad")
KEY <- "pop_x_poly"

# ---- 2. ordered logit with year dummies, cluster-robust on dyad ----
ologit_parts <- function(th, X, y, J) {           # th = (beta, cut_1..cut_{J-1}); y in 0..J-1
  k <- ncol(X); b <- th[1:k]; cut <- th[k + 1:(J - 1)]
  xb <- drop(X %*% b); cc <- c(-Inf, cut, Inf)
  au <- cc[y + 2] - xb; al <- cc[y + 1] - xb
  Fu <- plogis(au); Fl <- plogis(al); fu <- dlogis(au); fl <- dlogis(al)
  p <- pmax(Fu - Fl, 1e-300)
  Sb <- X * (-(fu - fl) / p)
  Sc <- matrix(0, length(y), J - 1)
  for (j in 1:(J - 1)) Sc[, j] <- (y + 1 == j) * fu / p - (y == j) * fl / p   # cut j is the upper bound for y = j-1, lower for y = j
  list(ll = log(p), S = cbind(Sb, Sc))
}
fit_ologit <- function(rhs_v, data) {
  yr <- factor(data$year); YD <- model.matrix(~ yr)[, -1, drop = FALSE]             # year dummies (first year = base, as Stata i.year)
  X <- cbind(as.matrix(data[, rhs_v, drop = FALSE]), YD); y <- data$hihosta_3; J <- 4
  stopifnot(!anyNA(X), all(y %in% 0:(J - 1)))
  sx <- apply(X, 2, sd); sx[sx == 0] <- 1; Xs <- sweep(X, 2, sx, "/")
  q <- cumsum(table(factor(y, levels = 0:(J - 1))))[1:(J - 1)] / length(y)
  st <- c(rep(0, ncol(X)), qlogis(q))
  o <- nlminb(st, function(t) -sum(ologit_parts(t, Xs, y, J)$ll), function(t) -colSums(ologit_parts(t, Xs, y, J)$S),
              control = list(eval.max = 5000, iter.max = 3000, rel.tol = 1e-13))
  th <- o$par
  hfun <- function(t) { H <- numDeriv::jacobian(function(u) colSums(ologit_parts(u, Xs, y, J)$S), t, method.args = list(r = 2)); (H + t(H)) / 2 }
  for (it in 1:10) { g <- colSums(ologit_parts(th, Xs, y, J)$S); if (max(abs(g)) < 1e-7) break; th <- th + solve(-hfun(th), g) }
  P <- ologit_parts(th, Xs, y, J); gmax <- max(abs(colSums(P$S)))
  stopifnot(gmax < 1e-4)
  Ainv <- solve(-hfun(th)); Sc <- rowsum(P$S, data$dyad); G <- nrow(Sc)
  V <- Ainv %*% crossprod(Sc) %*% Ainv * G / (G - 1)
  D <- c(1 / sx, rep(1, J - 1)); coef <- th * D; V <- V * outer(D, D)
  names(coef) <- c(colnames(X), paste0("cut", 1:(J - 1))); se <- sqrt(diag(V)); z <- coef / se
  cat("\n---- ologit | rhs:", paste(rhs_v, collapse = " + "), "+ year dummies | N =", length(y), "| clusters =", G,
      "| LL =", round(sum(P$ll), 4), "| max|score|", signif(gmax, 3), "\n")
  print(round(cbind(b = coef, se = se, z = z)[c(rhs_v, paste0("cut", 1:(J - 1))), ], 6))
  out <- list(c(b = coef[[KEY]], se = se[[KEY]], p = 2 * pnorm(-abs(z[[KEY]])))); names(out) <- KEY
  attr(out, "n") <- length(y); out
}

# ---- 3. reproduction on the published sample ----
d0 <- d[complete.cases(d[, vars]), ]
cat("complete-case N:", nrow(d0), "(published 101,468)\n"); stopifnot(nrow(d0) == 101468)
full <- fit_ologit(rhs, d0)
b_pub <- -1.534721; se_pub <- 0.6288662
cat("REPRO", KEY, "b", full[[KEY]][["b"]], "(pub", b_pub, ") se", full[[KEY]][["se"]], "(pub", se_pub, ")\n")
stopifnot(abs(full[[KEY]][["b"]] - b_pub) / abs(b_pub) <= 0.10, abs(full[[KEY]][["se"]] - se_pub) / se_pub <= 0.10)

# ---- 4. history window + six models on one common sample ----
fitter <- function(r, data) fit_ologit(r, data)
ds <- common_sample(d, vars, "statea", "stateb", "year")
M <- run_six(rhs, ds, fitter)
saveRDS(list(M = M, full = full), file.path(W, "M_ADD_D10.rds"))

# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
# headline is a product term -> split only (sens_split)
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_split(id = "ADD_D10", dat = ds, fit = function(dat, add) fitter(c(rhs, add), dat)[[KEY]])
}

meta <- list(pid = "ADD_D10", cite = "ADD_D10", journal = "ISQ", year = 2026, doi = "10.1093/isq/sqag085", family = "DEMOCRACY")
nt <- sprintf(paste0("PORT of Dataverse doi:10.7910/DVN/IXWNFI. Table 1 Model 2 ologit, year dummies, SE clustered on dyad. Full-sample b=%.6f se=%.6f N=%d ",
  "vs published -1.534721 (0.6288662) N 101,468 (authors' log). Common 1822-2011 sample N=%d (2012-2014 outside the history window)."),
  full[[KEY]][["b"]], full[[KEY]][["se"]], attr(full, "n"), nrow(ds))
write_p8(meta, KEY, "Title/Model 2: populism x democratic governance interaction (statea_populist_strict # c.statea_v2x_polyarchy)",
  "Table 1 Model 2 (ologit, directed PR dyad-years 1970-2014)", "ordered logit, year dummies, SE clustered on dyad (Stata ologit vce(cluster dyad))",
  "ordered logit ML (R port, analytic scores), cluster-robust on dyad, all six steps", b_pub, se_pub, 101468, TRUE, FALSE, M, "N (no prior-conflict controls)", nt)
dir.create(file.path(OUT, "results_add"), recursive = TRUE, showWarnings = FALSE)
stopifnot(file.rename(file.path(OUT, "results_p8", "ADD_D10__pop_x_poly.csv"), file.path(OUT, "results_add", "ADD_D10__pop_x_poly.csv")))
sink()
