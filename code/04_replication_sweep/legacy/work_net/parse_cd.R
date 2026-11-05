# R port of work_net/parse_cd.py (CD2011): parse the btergm log written by CD2011/est_cd2011.R (via run_all.sh) and write
# results_net/CD2011__{SEq,MinRegime,SEq_LDV}.csv. Numeric logic is the Python original's (values taken from the printed log,
# pct = round(100 * (x - b0) / |b0|, 1), pct_placebo against the re-estimated placebo-sample baseline 'pbase').
# The long free-text columns (study, convergence, notes) of the original are not reproduced.
# Addition for the pipeline: the parsed per-model coefficient table is also saved to work_net/CD2011/log_cd2011_parsed.rds,
# because two master rows (kstar2, triangle) were read directly off this log and have no results_net CSV.
OUT <- Sys.getenv("SWEEP_OUT"); stopifnot(nzchar(OUT))
W <- file.path(OUT, "work_net/CD2011"); O <- file.path(OUT, "results_net"); dir.create(O, recursive = TRUE, showWarnings = FALSE)
txt <- paste(readLines(file.path(W, "log_cd2011.txt")), collapse = "\n")
M <- list()
for (blk in strsplit(txt, "====", fixed = TRUE)[[1]][-1]) {
  lines <- strsplit(trimws(blk), "\n", fixed = TRUE)[[1]]
  name <- strsplit(trimws(lines[1]), "\\s+")[[1]][1]
  M[[name]] <- list()
  for (l in lines[-(1:2)]) {
    p <- strsplit(trimws(l), "\\s+")[[1]]
    if (length(p) == 7 && grepl("^[0-9]+$", p[1])) {
      v <- as.numeric(p[3:7]); stopifnot(!anyNA(v))
      M[[name]][[p[2]]] <- list(b = v[1], sd = v[2], lo = v[3], hi = v[4], p = v[5])
    }
  }
}
stopifnot(all(c("m0", "hist", "m26", "ldv0", "ldvhist", "pbase", "placebo") %in% names(M)))
saveRDS(M, file.path(W, "log_cd2011_parsed.rds"))
pct <- function(x, b0) round(100 * (x - b0) / abs(b0), 1)
row <- function(key, term, pub, pubci, base, hist, m26, pb, pl, verdict) {
  a <- M[[base]][[term]]; h <- M[[hist]][[term]]; stopifnot(!is.null(a), !is.null(h))
  c26 <- if (!is.null(m26)) M[[m26]][[term]] else NULL; q <- if (!is.null(pl)) M[[pl]][[term]] else NULL
  qb <- if (!is.null(pb)) M[[pb]][[term]] else NULL
  data.frame(tag = "CD2011", key_effect = key, model_term = term, pub_b = pub, pub_ci = pubci,
    b0 = a$b, se0 = a$sd, p0 = a$p, ci0_excl0 = !(a$lo <= 0 & 0 <= a$hi), b_hist = h$b, se_hist = h$sd, p_hist = h$p,
    ci_hist_excl0 = !(h$lo <= 0 & 0 <= h$hi), pct_hist = pct(h$b, a$b),
    b_26 = if (!is.null(c26)) c26$b else "", pct_26 = if (!is.null(c26)) pct(c26$b, a$b) else "",
    b_placebo = if (!is.null(q)) q$b else "", pct_placebo = if (!is.null(q)) pct(q$b, qb$b) else "",
    placebo_base_b0 = if (!is.null(qb)) qb$b else "",
    convergence = "MPLE via btergm (fastglm), 1000 network-year bootstrap reps", verdict = verdict, check.names = FALSE)
}
rows <- list(
  SEq       = row("Integrated SEq (ERGM col.)", "edgecov.seqv[[i]]", -0.511, "[-0.682, -0.352]", "m0", "hist", "m26", "pbase", "placebo", "CLOSE"),
  MinRegime = row("Minimum Regime Score (ERGM col.)", "edgecov.minreg[[i]]", -0.001, "[-0.0035, 0.0001]", "m0", "hist", "m26", "pbase", "placebo", "CLOSE"),
  SEq_LDV   = row("Integrated SEq (ERGM LDV col.)", "edgecov.seqv[[i]]", -0.344, "[-0.544, -0.171]", "ldv0", "ldvhist", NULL, NULL, NULL, "NOT_REPRODUCED"))
for (fn in names(rows)) {
  r <- rows[[fn]]; write.csv(r, file.path(O, paste0("CD2011__", fn, ".csv")), row.names = FALSE)
  cat(fn, r$b0, r$b_hist, r$pct_hist, r$b_26, r$pct_26, r$b_placebo, r$pct_placebo, "\n")
}
