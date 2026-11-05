# Rerun 58: Goldsmith 2007 JPR, "A liberal peace in Asia?"
# Headline: Model 3, cdsimeq (Keshk 2SPLS: OLS lntradeab <-> probit cwfataldmy), 1950-2000.
# Tie: lngleddeplow (lower trade dependence) in the probit equation; enters with x Asia interaction ->
# constituent term = non-Asia effect. b_pub from authors' smcl log: lngleddeplow = -90.52128 (corrected SE 15.16587), N = 403,092.
# cdsimeq is ported below line-for-line from cdsimeq.ado v2.0 (downloaded, in this folder):
#   stage 1: OLS cont_dep on all exog; probit dich_dep on all exog (V_o = its vcov)
#   stage 2: probit dich_dep ~ I_cont + dich_iv ; Amemiya/Maddala corrected VCE
#   alpha2_vce = dd + MA_d * dd G'V0^-1 (X'X)^-1 V0^-1 G dd,  dd = (G'V0^-1 G)^-1, G = [pie_1, J2]
#   MA_d = gamma2^2 sigma_11 - 2 gamma2 sigma_12; sigma_12 = mean(y2 * v1 / dnorm(xb_probit1))
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work/rerun_common_b7.R"))
WO <- file.path(OUT, "work/58"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work/58"); sink(file.path(WO, "log_58.txt"), split = TRUE)

# Stata probit default VCE = inverse OBSERVED information (OIM); glm.fit gives expected information,
# which differs for probit. oim_vcov() computes -H^-1 with H = -sum lambda(lambda + xb) x x'.
oim_vcov <- function(X, y, b) {
  xb <- drop(X %*% b); q <- 2 * y - 1
  lam <- q * dnorm(q * xb) / pnorm(q * xb)
  w <- lam * (lam + xb)
  V <- solve(crossprod(X * sqrt(w)))
  dimnames(V) <- list(colnames(X), colnames(X)); V
}

# Lean Newton-Raphson probit (glm.fit copies X several times and got OOM-killed on the shared box).
# Converges to the same MLE; start from OLS-scaled values; stop when max|step| < 1e-10.
probit_nr <- function(X, y) {
  b <- rep(0, ncol(X)); b[ncol(X)] <- qnorm(mean(y))
  for (it in 1:100) {
    xb <- drop(X %*% b); q <- 2 * y - 1
    lam <- q * dnorm(q * xb) / pnorm(q * xb)
    g <- crossprod(X, lam)
    H <- crossprod(X * sqrt(lam * (lam + xb)))
    step <- drop(solve(H, g)); b <- b + step
    if (max(abs(step)) < 1e-10) break
  }
  stopifnot(max(abs(step)) < 1e-8)
  names(b) <- colnames(X); list(coefficients = b, iter = it)
}

C <- function(k) file.path(WO, paste0("cache_", k, ".rds"))
cdsimeq_probit <- function(d, cont_dep, cont_iv, dich_dep, dich_iv, tie, cache = NULL) {
  if (!is.null(cache) && file.exists(cache)) { o <- readRDS(cache); stopifnot(!is.null(o$tab)); cat(o$txt, sep = "\n"); return(o) }  # port: cache must carry full table (extraction)
  exog <- unique(c(cont_iv, dich_iv))
  d <- d[complete.cases(d[, c(cont_dep, dich_dep, exog)]), ]           # marksample over both equations
  X <- cbind(as.matrix(d[, exog]), `_cons` = 1)                          # Stata puts _cons last
  y1 <- d[[cont_dep]]; y2 <- d[[dich_dep]]
  # stage 1 OLS
  XtX <- crossprod(X); stopifnot(rcond(XtX) > 1e-14)
  pie_1 <- drop(solve(XtX, crossprod(X, y1))); names(pie_1) <- colnames(X); v1 <- y1 - drop(X %*% pie_1)
  sigma_11 <- sum(v1^2) / (nrow(X) - ncol(X))
  # stage 1 probit
  p1 <- probit_nr(X, y2)
  pie_2 <- p1$coefficients
  xb1 <- drop(X %*% pie_2)
  V_o <- oim_vcov(X, y2, pie_2)
  sigma_12 <- mean(y2 * v1 / dnorm(xb1))
  # stage 2 probit
  I_cont <- drop(X %*% pie_1)
  Z <- cbind(I_cont = I_cont, as.matrix(d[, dich_iv]), `_cons` = 1)
  p2 <- probit_nr(Z, y2)
  alpha_2 <- p2$coefficients; names(alpha_2) <- colnames(Z)
  gamma_2 <- alpha_2["I_cont"]
  se_uncorr <- sqrt(diag(oim_vcov(Z, y2, alpha_2)))
  # corrected VCE
  MA_d <- unname(gamma_2^2 * sigma_11 - 2 * gamma_2 * sigma_12)
  J2 <- matrix(0, ncol(X), ncol(Z) - 1, dimnames = list(colnames(X), colnames(Z)[-1]))
  for (v in colnames(Z)[-1]) J2[v, v] <- 1
  G <- cbind(pie_1, J2)
  XTXinv <- solve(XtX); aa <- solve(V_o)
  dd <- solve(t(G) %*% aa %*% G)
  vce <- dd + MA_d * dd %*% (t(G) %*% aa %*% XTXinv %*% aa %*% G) %*% dd
  se <- sqrt(diag(vce)); z <- alpha_2 / se
  tab <- cbind(b = alpha_2, se_corr = se, z = z, p = 2 * pnorm(-abs(z)), se_uncorr = se_uncorr)
  txt <- c("", paste("---- cdsimeq probit eq:", dich_dep, "~ I_", cont_dep, "+", paste(dich_iv, collapse = " + "), "| N =", nrow(X)),
           paste("sigma_11", sigma_11, "sigma_12", sigma_12, "gamma_2", gamma_2, "MA_d", MA_d), capture.output(print(round(tab, 5))))
  cat(txt, sep = "\n")
  out <- list(b = unname(alpha_2[tie]), se = unname(se[tie]), p = unname(2 * pnorm(-abs(z[tie]))), n = nrow(X), txt = txt,
              tab = tab)  # port: keep full-precision coefficient table for extraction of the product term
  if (!is.null(cache)) saveRDS(out, cache)
  rm(X, Z, p1, p2); gc(); out
}

# 1. Load authors' data
keep <- c("ccode1", "ccode2", "year", "lntradeab", "lnpolity21", "lnpolity22", "lncgdpa", "lncgdpb", "lndistance",
          "lntradeablag", "cwfataldmy", "pol2low", "pol2lowAsia", "lngleddeplow", "lngleddeplowAsia", "lnnumigofill",
          "lnnumigofillAsia", "Asia", "minpwrdyad", "paritysqrt", "colcont24", "contig24", "lcwfataldmy")
# slim_58.csv = the keep-columns of LPAgoldsmith.dta for 1950-2000, extracted with pandas read_stata (chunked)
# because reading the 90MB dta in R was OOM-killed on the shared box; values are exact (float32 -> %.17g).
# port: slim_58 is rebuilt into OUT from the authors' LPAgoldsmith.dta (staged under DATA): keep-columns, years
# 1950-2000, whole-number columns stored as integer exactly as the old CSV round trip (read.csv) typed them.
slim_out <- file.path(WO, "slim_58.rds")
if (!file.exists(slim_out)) {
  raw <- haven::read_dta(file.path(W, "LPAgoldsmith.dta"), col_select = dplyr::all_of(keep))
  raw <- as.data.frame(haven::zap_formats(haven::zap_labels(haven::zap_label(raw))))[, keep]
  raw <- raw[!is.na(raw$year) & raw$year >= 1950 & raw$year <= 2000, ]
  for (v in keep) { x <- raw[[v]]; if (all(is.na(x) | (x == round(x) & abs(x) < .Machine$integer.max))) raw[[v]] <- as.integer(x) }
  rownames(raw) <- NULL
  stopifnot(nrow(raw) == 581464)  # slim_58 row count in the original run
  saveRDS(raw, slim_out); rm(raw); gc()
}
d <- readRDS(slim_out)
stopifnot(identical(sort(names(d)), sort(keep)))
d <- d[!is.na(d$year) & d$year > 1949 & d$year < 2001, ]
stopifnot(all(d$cwfataldmy %in% c(0, 1, NA)))

cont_dep <- "lntradeab"; cont_iv <- c("lnpolity21", "lnpolity22", "lncgdpa", "lncgdpb", "lndistance", "lntradeablag")
dich_dep <- "cwfataldmy"
dich_iv <- c("pol2low", "pol2lowAsia", "lngleddeplow", "lngleddeplowAsia", "lnnumigofill", "lnnumigofillAsia", "Asia",
             "minpwrdyad", "paritysqrt", "colcont24", "contig24", "lndistance", "lcwfataldmy")
T <- "lngleddeplow"; b_pub <- -90.52128; n_orig <- 403092

# 2. Reproduce
r <- cdsimeq_probit(d, cont_dep, cont_iv, dich_dep, dich_iv, T, C("repro"))
stopifnot(r$n == n_orig)
repro_ok <- abs(r$b - b_pub) / abs(b_pub) <= 0.10 && r$p < .05
cat("b_pub", b_pub, "b_repro", r$b, "rel diff", abs(r$b - b_pub) / abs(b_pub), "\n")

# 3. History (undirected ccode1<ccode2?) and restriction 1822-2011 (sample 1950-2000 already inside)
cat("ccode1<ccode2 share:", mean(d$ccode1 < d$ccode2, na.rm = TRUE), "\n")
d <- d[!is.na(d$ccode1) & !is.na(d$ccode2), ]
d <- add_history(d, "ccode1", "ccode2", "year")
d <- d[complete.cases(d[, c(cont_dep, dich_dep, cont_iv, dich_iv, "hist_dyad_w15", "hist_mono_w15", "hist_dyad_w26", "hist_mono_w26")]), ]
cat("common sample N:", nrow(d), "\n")

# 4. Added controls enter the probit (conflict) equation, hence the exogenous list of stage 1 as well
m0  <- cdsimeq_probit(d, cont_dep, cont_iv, dich_dep, dich_iv, T, C("m0"))
mD  <- cdsimeq_probit(d, cont_dep, cont_iv, dich_dep, c(dich_iv, "hist_dyad_w15"), T, C("mD"))
mDM <- cdsimeq_probit(d, cont_dep, cont_iv, dich_dep, c(dich_iv, "hist_dyad_w15", "hist_mono_w15"), T, C("mDM"))
m26 <- cdsimeq_probit(d, cont_dep, cont_iv, dich_dep, c(dich_iv, "hist_dyad_w26", "hist_mono_w26"), T, C("m26"))
bC  <- cdsimeq_probit(d, cont_dep, cont_iv, dich_dep, c(dich_iv, "bm_cinc_w15"), T, C("bC"))
bD  <- cdsimeq_probit(d, cont_dep, cont_iv, dich_dep, c(dich_iv, "bm_democ_w15"), T, C("bD"))

# 5. Results
write_row(58, "Goldsmith 2007", if (repro_ok) "REPRODUCED" else "NOT_REPRODUCED", b_pub, r$b,
  "cdsimeq 2SPLS probit eq (Keshk corrected SE, ported to R)", "lngleddeplow (constituent; x Asia)", n_orig, m0$n,
  m0, mD, mDM, m26, bC, bD, ahc = "Y",
  notes = paste0("Model 3 (fatal MID 1950-2000); b and Keshk-corrected SE (15.16587) match authors log exactly. Tie is constituent term of lngleddeplow x Asia (non-Asia dyads). ",
                 "Author already controls lagged DV (lcwfataldmy). History added to probit eq (thus also to stage-1 exog). ",
                 n_orig - m0$n, " obs dropped for history merge. bC N=", bC$n, ", bD N=", bD$n, "."))

# port: extraction-only writer. The master row is the product term lngleddeplowAsia (probit eq), read from the
# same fitted models' full-precision coefficient tables (no re-estimation).
K <- "lngleddeplowAsia"
kt <- function(m, col) { stopifnot(K %in% rownames(m$tab)); unname(m$tab[K, col]) }
key58 <- data.frame(id = "58", key_term = K,
  b0 = kt(m0, "b"), se0 = kt(m0, "se_corr"), p0 = kt(m0, "p"),
  b_dm = kt(mDM, "b"), se_dm = kt(mDM, "se_corr"), p_dm = kt(mDM, "p"),
  b_26 = kt(m26, "b"), p_26 = kt(m26, "p"),
  pct_bC = pct(kt(bC, "b"), kt(m0, "b")), pct_bD = pct(kt(bD, "b"), kt(m0, "b")),
  n_used = m0$n, stringsAsFactors = FALSE)
stopifnot(nrow(key58) == 1, all(is.finite(unlist(key58[, -(1:2)]))))
write.csv(key58, file.path(WO, "key_58_extract.csv"), row.names = FALSE)
print(key58, digits = 10)
sink()
