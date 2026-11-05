# Verify-phase (batch 2) helpers. Sequential, fail-loud. Caches go to <W>/verify_cache/ (never the old phase caches).
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
suppressPackageStartupMessages({library(dplyr); library(fixest)})
source(file.path(CODE, "hist/add_history.R"))
set.seed(123)
HV <- c("hist_dyad_w15", "hist_mono_w15", "hist_dyad_w26", "hist_mono_w26", "bm_cinc_w15", "bm_democ_w15")

# fitfun(f, data) -> matrix with columns (b, se, p), rownames = terms, attr "n"
vfit <- function(f, data, terms, fitfun, cache) {
  if (file.exists(cache)) { o <- readRDS(cache); cat("\n---- [cached]", o$txt, sep = "\n"); return(o) }
  ct <- fitfun(f, data)
  stopifnot(all(terms %in% rownames(ct)), all(is.finite(ct[terms, 1])), all(is.finite(ct[terms, 2])))
  txt <- c(paste("----", paste(deparse(f, width.cutoff = 500), collapse = " "), "| N =", attr(ct, "n")),
           capture.output(print(round(ct, 6))))
  cat("\n", txt, sep = "\n")
  out <- list(k = lapply(setNames(terms, terms), function(t) c(b = ct[t, 1], se = ct[t, 2], p = ct[t, 3])), n = attr(ct, "n"), txt = txt)
  saveRDS(out, cache); gc()
  out
}

vcommon_sample <- function(d, f0, extra, c1, c2, yr) {
  d <- d[complete.cases(d[, unique(c(all.vars(f0), extra))]), ]
  n_cc <- nrow(d)
  d <- add_history(d, c1, c2, yr)
  d <- d[d[[yr]] >= 1822 & d[[yr]] <= 2011, ]
  d <- d[complete.cases(d[, HV]) & is.finite(d$bm_cinc_w15) & is.finite(d$bm_democ_w15), ]
  cat("complete-case N:", n_cc, " common sample N:", nrow(d), " dropped:", n_cc - nrow(d), "\n")
  stopifnot(nrow(d) > 0.8 * n_cc, sd(d$hist_dyad_w15) > 0)
  attr(d, "n_cc") <- n_cc
  d
}

vrun_six <- function(f0, ds, terms, fitfun, W) {
  CD <- file.path(W, "verify_cache"); dir.create(CD, showWarnings = FALSE)
  C <- function(k) file.path(CD, paste0("v_", k, ".rds"))
  M <- list(m0  = vfit(f0, ds, terms, fitfun, C("m0")),
            mD  = vfit(update(f0, . ~ . + hist_dyad_w15), ds, terms, fitfun, C("mD")),
            mDM = vfit(update(f0, . ~ . + hist_dyad_w15 + hist_mono_w15), ds, terms, fitfun, C("mDM")),
            m26 = vfit(update(f0, . ~ . + hist_dyad_w26 + hist_mono_w26), ds, terms, fitfun, C("m26")),
            bC  = vfit(update(f0, . ~ . + bm_cinc_w15), ds, terms, fitfun, C("bC")),
            bD  = vfit(update(f0, . ~ . + bm_democ_w15), ds, terms, fitfun, C("bD")))
  stopifnot(all(sapply(M, function(m) m$n) == M$m0$n))   # one common sample
  M
}

# long table of the six models (+ optional repro) for the key terms -> <W>/verify_six.csv
vsave <- function(M, terms, W, rep = NULL) {
  L <- if (is.null(rep)) M else c(list(repro = rep), M)
  out <- do.call(rbind, lapply(names(L), function(m) do.call(rbind, lapply(terms, function(t)
    data.frame(model = m, term = t, b = L[[m]]$k[[t]][["b"]], se = L[[m]]$k[[t]][["se"]], p = L[[m]]$k[[t]][["p"]], n = L[[m]]$n)))))
  write.csv(out, file.path(W, "verify_six.csv"), row.names = FALSE); print(out)
  invisible(out)
}

ff_logit_cl <- function(cl) function(f, data) {
  m <- feglm(f, data = data, family = binomial("logit"), vcov = as.formula(paste0("~", cl)), notes = FALSE, mem.clean = TRUE, lean = TRUE)
  stopifnot(m$convStatus); ct <- coeftable(m)[, c(1, 2, 4), drop = FALSE]; attr(ct, "n") <- nobs(m); ct
}
ff_logit_cl_fe <- function(cl, fe) function(f, data) {
  f2 <- as.formula(paste(paste(deparse(f, width.cutoff = 500), collapse = " "), "|", fe))
  m <- feglm(f2, data = data, family = binomial("logit"), vcov = as.formula(paste0("~", cl)), notes = FALSE, mem.clean = TRUE, lean = TRUE)
  stopifnot(m$convStatus); ct <- coeftable(m)[, c(1, 2, 4), drop = FALSE]; attr(ct, "n") <- nobs(m); ct
}

# Tucker btscs rule used in earlier reruns (cid 43/65/100): py = years since last event within unit (0 at first obs and in
# the year after an event); nspline(3) basis s_i = (py-k_i)^3_+ - (T-k_i)/T * py^3, knots = quartiles (type 2) of py, T = max(py).
btscs_py <- function(unit, year, ev) {       # data must be sorted by unit, year
  stopifnot(!is.unsorted(order(unit, year)) || TRUE)
  o <- order(unit, year); stopifnot(identical(o, seq_along(o)))
  first <- c(TRUE, unit[-1] != unit[-length(unit)])
  py <- numeric(length(ev)); cur <- -1; lastev <- FALSE
  for (i in seq_along(ev)) { if (first[i]) { cur <- -1; lastev <- FALSE }
    cur <- if (cur < 0 || lastev) 0 else cur + 1; py[i] <- cur; lastev <- !is.na(ev[i]) && ev[i] == 1 }
  py
}
btscs_splines <- function(d, py = "py") {
  kn <- quantile(d[[py]], c(.25, .5, .75), type = 2, na.rm = TRUE); Tm <- max(d[[py]], na.rm = TRUE)
  cat("btscs knots:", kn, " T:", Tm, "\n")
  for (i in 1:3) d[[paste0("spl", i)]] <- (pmax(d[[py]] - kn[i], 0)^3 - (Tm - kn[i]) / Tm * d[[py]]^3) / 1000
  attr(d, "knots") <- c(kn, T = Tm); d
}
