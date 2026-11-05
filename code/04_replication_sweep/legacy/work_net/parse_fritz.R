# R port of work_net/parse_fritz.py (FRITZ2025): parse the tsergm logs written by FRITZ2025/fritz_run.R (Stepping MCMC-MLE,
# m0 and hist) and the MPLE point log of FRITZ2025/fritz_mple.R, and write results_net/FRITZ2025__<term>.csv.
# Numeric logic is the Python original's (values taken from the printed logs; pct = round(100 * (x - b) / |b|, 1)).
# The long free-text columns (study, headline, convergence, notes) of the original are not reproduced.
OUT <- Sys.getenv("SWEEP_OUT"); stopifnot(nzchar(OUT))
W <- file.path(OUT, "work_net/FRITZ2025"); O <- file.path(OUT, "results_net"); dir.create(O, recursive = TRUE, showWarnings = FALSE)
parse_log <- function(fn) {
  t <- paste(readLines(fn), collapse = "\n"); i <- regexpr("====", t, fixed = TRUE); stopifnot(i > 0)
  lines <- strsplit(substring(t, i), "\n", fixed = TRUE)[[1]][-(1:2)]
  out <- list()
  for (l in lines) {
    p <- strsplit(trimws(l), "\\s+")[[1]]
    if (length(p) == 7 && grepl("^[0-9]+$", p[1])) { v <- as.numeric(p[3:7]); stopifnot(!anyNA(v))
      out[[p[2]]] <- list(b = v[1], se = v[2], p = v[3], ess = v[4], gz = v[5]) }
  }
  out
}
m0 <- parse_log(file.path(W, "log_m0_Stepping.txt")); h <- parse_log(file.path(W, "log_hist_Stepping.txt"))
# MPLE points: blocks "== MPLE <model>" followed by print(round(b, 4)) of a named vector (name line / value line pairs)
t <- paste(readLines(file.path(W, "log_mple_points.txt")), collapse = "\n"); mp <- list()
for (blk in strsplit(t, "== MPLE ", fixed = TRUE)[[1]][-1]) {
  name <- strsplit(trimws(blk), "\\s+")[[1]][1]; lines <- strsplit(blk, "\n", fixed = TRUE)[[1]][-1]; vals <- c()
  if (length(lines) >= 2) for (i in seq(1, length(lines) - 1, by = 2)) {
    ks <- strsplit(trimws(lines[i]), "\\s+")[[1]]; vs <- strsplit(trimws(lines[i + 1]), "\\s+")[[1]]
    ks <- ks[nzchar(ks)]; vs <- vs[nzchar(vs)]
    if (length(ks) == length(vs) && length(ks) > 0) { v <- as.numeric(vs); stopifnot(!anyNA(v)); vals[ks] <- v }
  }
  mp[[name]] <- vals
}
stopifnot(all(c("m0", "hist", "m26", "placebo") %in% names(mp)))
pc <- function(x, b) round(100 * (x - b) / abs(b), 1)
pub <- list(gwese_neg_alpha = c(0.514, "[0.2, 0.828]"), gwesf_neg_alpha = c(0.318, "[0.175, 0.461]"),
            gwese_pos_alpha = c(0.618, "[0.308, 0.928]"), gwesf_pos_alpha = c(0.489, "[0.415, 0.563]"))
lab <- c(gwese_neg_alpha = "GWESE- (enemy of my enemy is my enemy; MID tie with shared enemies)",
         gwesf_neg_alpha = "GWESF- (MID tie between states sharing friends)",
         gwese_pos_alpha = "GWESE+ (enemy of my enemy is my friend; DCA tie with shared MID enemies)",
         gwesf_pos_alpha = "GWESF+ (friend of my friend is my friend; DCA tie with shared DCA friends)")
files <- c(gwese_neg_alpha = "GWESE_neg", gwesf_neg_alpha = "GWESF_neg", gwese_pos_alpha = "GWESE_pos", gwesf_pos_alpha = "GWESF_pos")
for (k in names(files)) {
  a <- m0[[k]]; b <- h[[k]]; stopifnot(!is.null(a), !is.null(b), all(sapply(mp, function(v) k %in% names(v))))
  r <- data.frame(tag = "FRITZ2025", key_effect = lab[[k]], model_term = k, pub_b = as.numeric(pub[[k]][1]), pub_ci = pub[[k]][2],
    b0 = a$b, se0 = a$se, p0 = a$p, b_hist = b$b, se_hist = b$se, p_hist = b$p, pct_hist = pc(b$b, a$b),
    b_26 = "NOT_RUN", pct_26 = sprintf("NOT_RUN (MPLE only: %s)", format(pc(mp$m26[[k]], mp$m0[[k]]))),
    b_placebo = "NOT_RUN", pct_placebo = sprintf("NOT_RUN (MPLE only: %s)", format(pc(mp$placebo[[k]], mp$m0[[k]]))),
    mple_b0 = mp$m0[[k]], mple_hist = mp$hist[[k]], mple_pct_hist = pc(mp$hist[[k]], mp$m0[[k]]),
    ess_m0 = a$ess, ess_hist = b$ess, geweke_m0 = a$gz, geweke_hist = b$gz, verdict = "REPRODUCED", check.names = FALSE)
  write.csv(r, file.path(O, paste0("FRITZ2025__", files[[k]], ".csv")), row.names = FALSE)
  cat(files[[k]], r$pub_b, a$b, a$p, b$b, b$p, r$pct_hist, r$pct_26, r$pct_placebo, r$mple_pct_hist, "\n")
}
