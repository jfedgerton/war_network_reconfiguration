# Round-2 rebuild verdict (PI rule, 2026-10-11). Sourced by every legacy helper that writes a rebuild verdict
# (work_p8/common_b3.R, common_b2.R, p8_common.R, vb3/vhelp.R, vb4/vhelp.R, vb5/verify_helpers.R) and by the registered
# rebuild scripts that compute their verdict inline.
# A protocol-B rebuild MATCHES when all three hold:
#   (1) the rebuilt baseline b0 has the published sign;
#   (2) b0 is significant at the paper's alpha_pub whenever the published estimate is significant at alpha_pub;
#   (3) b0 lies within 0.5x-2x of the published coefficient.
# alpha_pub is read from registry.csv (row matching id and key_term; else the id's single alpha). Published significance is
# judged at alpha_pub from pub_b/pub_se when pub_se is given, otherwise the caller's pub_sig is used.
# Ports of author data (REPRODUCED / CLOSE / NOT_REPRODUCED) keep their own rule in each helper; no magnitude gate there.
.R2_REG <- NULL
alpha_pub_of <- function(id, term = NULL) {
  if (is.null(.R2_REG)) {
    f <- file.path(dirname(normalizePath(Sys.getenv("SWEEP_CODE"))), "registry.csv"); stopifnot(file.exists(f))
    .R2_REG <<- read.csv(f, colClasses = "character")
  }
  r <- .R2_REG[.R2_REG$id == as.character(id), ]
  if (!is.null(term) && any(r$key_term == term)) r <- r[r$key_term == term, ]
  a <- unique(r$alpha_pub)
  if (length(a) == 0) return(NA_real_)
  stopifnot(length(a) == 1, a %in% c("0.05", "0.10"))
  as.numeric(a)
}
rebuild_verdict <- function(id, term, b0, p0, pub_b, pub_sig = NA, pub_se = NA) {
  if (is.na(pub_b)) return("NO_PUBLISHED")
  a <- alpha_pub_of(id, term)
  if (is.na(a)) { a <- 0.05; cat("NOTE: rebuild", id, "is not registered; alpha_pub 0.05 assumed for its verdict\n") }
  ps <- if (!is.na(pub_se) && pub_se > 0) 2 * pnorm(-abs(pub_b / pub_se)) < a else isTRUE(pub_sig)
  ratio <- b0 / pub_b
  ok <- isTRUE(sign(b0) == sign(pub_b)) && (!ps || isTRUE(p0 < a)) && isTRUE(ratio >= 0.5 && ratio <= 2)
  cat(sprintf("round-2 rebuild verdict %s/%s: alpha_pub %.2f, published sig %s, b0/pub %.3f, p0 %.3g -> %s\n", id, term, a, ps, ratio, p0,
              if (ok) "REBUILT_MATCH" else "REBUILT_MISMATCH"))
  if (ok) "REBUILT_MATCH" else "REBUILT_MISMATCH"
}
