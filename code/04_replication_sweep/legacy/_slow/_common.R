# Shared by legacy/_slow/<id>.R (slow network-model studies): environment checks, one step runner, the final writer.
# Each <id>.R runs the recorded chain (estimation -> parse/aggregate) and writes SWEEP_OUT/_slow/<id>.csv with columns
#   id,key_term,b0,p0,b_dm,p_dm,b_26,p_26,pct_bC,pct_bD,n_used   (one row per master row of that id).
# Env: SWEEP_CODE, SWEEP_DATA, SWEEP_OUT, SWEEP_HIST (as for run_one.R); SWEEP_LIB_AMEN / SWEEP_LIB_CHEN where noted.
#      SLOW_SKIP_EST=1 skips every estimation step and only re-runs parse/aggregate on fits already under SWEEP_OUT
#      (used to test the extraction against the cached estimates without spending hours).
# The steps call Rscript directly (no flock): whoever launches <id>.R is expected to hold the R lock for the whole chain.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE), dir.exists(ROOT), dir.exists(CODE))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
SKIP_EST <- Sys.getenv("SLOW_SKIP_EST") == "1"
RSCRIPT <- file.path(R.home("bin"), "Rscript")
REG <- read.csv(file.path(dirname(normalizePath(CODE)), "registry.csv"), colClasses = "character")
stopifnot(file.exists(Sys.getenv("SWEEP_HIST")), file.exists(file.path(CODE, "hist/add_history.R")))

# key_term of a registry row (asserts the row belongs to `id`)
kt <- function(row_id, id) { k <- which(REG$row_id == row_id); stopifnot(length(k) == 1, REG$id[k] == id); REG$key_term[k] }

# Run one ported R script: cd wd && <env> Rscript CODE/<script> <args> > log 2>&1  (">>" when append = TRUE).
# Stops on a non-zero exit. ok_marker = TRUE touches <log>.ok on success (the queue scripts' convention; agg_mine.R needs it).
run_r <- function(script, args = character(0), log, wd, env = character(0), append = FALSE, ok_marker = FALSE) {
  f <- file.path(CODE, script); stopifnot(file.exists(f))
  dir.create(wd, recursive = TRUE, showWarnings = FALSE); dir.create(dirname(log), recursive = TRUE, showWarnings = FALSE)
  e <- c(SWEEP_CODE = CODE, SWEEP_DATA = ROOT, SWEEP_OUT = OUT, SWEEP_HIST = Sys.getenv("SWEEP_HIST"), env)
  cmd <- paste("cd", shQuote(wd), "&&", paste0(names(e), "=", shQuote(e), collapse = " "), shQuote(RSCRIPT), shQuote(f),
               paste(shQuote(args), collapse = " "), if (append) ">>" else ">", shQuote(log), "2>&1")
  cat(sprintf("[%s] start %s %s -> %s\n", format(Sys.time()), script, paste(args, collapse = " "), log))
  t0 <- Sys.time(); rc <- system(cmd)
  cat(sprintf("[%s] end   %s rc=%d (%.1f min)\n", format(Sys.time()), script, rc, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  if (rc != 0) stop("step failed (rc ", rc, "): ", script, " ", paste(args, collapse = " "), "; see ", log)
  if (ok_marker) file.create(paste0(log, ".ok"))
  invisible(rc)
}

# read one result CSV written by a legacy parser/aggregator (all columns character, as in run_one.R)
read_res <- function(rel) { f <- file.path(OUT, rel); if (!file.exists(f)) stop("missing result file: ", f)
  x <- read.csv(f, check.names = FALSE, colClasses = "character"); stopifnot(nrow(x) == 1); x }
num <- function(x) suppressWarnings(as.numeric(x))   # "NOT_RUN ..." / "" -> NA

# final per-id CSV; asserts exactly the registry's master rows of `id`, finite b0 and b_dm
write_slow <- function(id, d) {
  cols <- c("id", "key_term", "b0", "p0", "b_dm", "p_dm", "b_26", "p_26", "pct_bC", "pct_bD", "n_used")
  stopifnot(identical(names(d), cols), all(d$id == id), !anyDuplicated(d$key_term),
            setequal(d$key_term, REG$key_term[REG$id == id]), all(is.finite(d$b0)), all(is.finite(d$b_dm)))
  f <- file.path(OUT, "_slow", paste0(id, ".csv")); dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
  write.csv(d, f, row.names = FALSE)
  cat("wrote", f, "(", nrow(d), "rows )\n")
  invisible(f)
}
