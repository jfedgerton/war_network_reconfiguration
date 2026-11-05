# 01_fetch_data.R -- recreate data/raw/replication_sweep/ from the public sources in data_sources.csv.
#
# One row of code/04_replication_sweep/data_sources.csv per file (hist/ is rebuilt by 00_build_history.R and is not listed).
#   source_type dataverse   source = "doi:... fileid=<id> format=<original|tab|file>"
#                             original = ingested file, fetch the original upload (?format=original)
#                             tab      = ingested file, fetch the .tab export (no query)
#                             file     = non-ingested file, plain download
#   source_type url         source = direct URL
#   archive_member non-empty: the source is an archive; it is downloaded once to data/raw/_fetch_cache/ and that member
#                             is extracted to <rel_path>
#   prio_login / unknown    not fetched (listed at the end for manual download)
#   derived / package       not fetched (built by a pipeline or prep step; listed at the end only if missing)
# Every fetched file is checked against data_sources.csv md5 before it is moved into place. Existing files with the
# right md5 are skipped, so the script can be re-run after a failure. Any md5 mismatch stops the script.
#
# Needs: R >= 3.6 (download.file headers=), base R + utils (+ tools::md5sum, which ships with base R), internet.
# RAR and .zipx archives (X031, P8_049) and zips R's internal unzip cannot read are handed to the first external tool
# found among unzip, bsdtar, 7z, unar; without one those files are reported as failed.
# The Dataverse 2.4 GB MINHAS2022 archive and the 988 MB P8_033 member make a full run take a while.
# Run from the project root:  Rscript code/04_replication_sweep/01_fetch_data.R   (or set PIPE_ROOT)

options(timeout = 3600)
PIPE <- Sys.getenv("PIPE_ROOT", unset = "")
if (!nzchar(PIPE)) PIPE <- here::here()
DATA <- file.path(PIPE, "data/raw/replication_sweep")
SRC <- file.path(PIPE, "code/04_replication_sweep/data_sources.csv")
CACHE <- file.path(PIPE, "data/raw/_fetch_cache")   # downloaded archives; delete after a successful run
stopifnot(dir.exists(PIPE), file.exists(SRC))
dir.create(DATA, recursive = TRUE, showWarnings = FALSE)
dir.create(CACHE, recursive = TRUE, showWarnings = FALSE)

# Hard-coded special cases
GUNZIP <- c("work_net/CHEN2021/data")        # Dataverse original is a gzip-compressed .RData; the pipeline file is gunzipped
BROWSER_UA_HOSTS <- c("saramitchell.org")    # serves an HTML page unless the request looks like a browser
BROWSER_UA <- "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120 Safari/537.36"
DV_ACCESS <- "https://dataverse.harvard.edu/api/access/datafile/"

# ---- read and validate the manifest --------------------------------------------------------------
src <- read.csv(SRC, colClasses = "character", check.names = FALSE, na.strings = character(0))
cols <- c("rel_path", "bytes", "md5", "id", "source_type", "source", "archive_member", "retrieval_note", "verified")
stopifnot(identical(names(src), cols), nrow(src) > 0, !anyDuplicated(src$rel_path))
stopifnot(all(src$source_type %in% c("dataverse", "url", "prio_login", "derived", "package", "unknown")))
stopifnot(all(grepl("^[0-9a-f]{32}$", src$md5)), all(grepl("^[0-9]+$", src$bytes)))
stopifnot(!any(grepl("^/|(^|/)\\.\\.(/|$)", src$rel_path)), !any(grepl("^hist/", src$rel_path)))
is_dv <- src$source_type == "dataverse"
is_url <- src$source_type == "url"
stopifnot(all(grepl("^doi:10\\.7910/DVN/[A-Z0-9]+ fileid=[0-9]+ format=(original|tab|file)$", src$source[is_dv])))
stopifnot(all(grepl("^https?://[^ ]+$", src$source[is_url])))
stopifnot(!any(grepl("::", src$archive_member[is_dv | is_url], fixed = TRUE)))   # nested members only occur in prio_login rows
stopifnot(all(src$rel_path[src$rel_path %in% GUNZIP] %in% src$rel_path[is_dv & !nzchar(src$archive_member)]))
cat("data_sources.csv:", nrow(src), "rows\n"); print(table(src$source_type))

# ---- helpers ---------------------------------------------------------------------------------------
md5_of <- function(f) unname(tools::md5sum(f))

fetch_url <- function(row) {
  if (row$source_type == "url") return(row$source)
  fid <- sub("^.* fileid=([0-9]+) .*$", "\\1", row$source)
  fmt <- sub("^.* format=", "", row$source)
  paste0(DV_ACCESS, fid, if (fmt == "original") "?format=original" else "")
}

# download url to dest via dest.part; returns "" on success, otherwise the error text
download_to <- function(url, dest) {
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  part <- paste0(dest, ".part")
  unlink(part)
  hdr <- if (any(vapply(BROWSER_UA_HOSTS, grepl, logical(1), x = url, fixed = TRUE))) c("User-Agent" = BROWSER_UA) else NULL
  cat("    downloading", url, "\n")
  rc <- tryCatch(download.file(url, part, mode = "wb", quiet = TRUE, method = "libcurl", headers = hdr),
                 error = function(e) conditionMessage(e), warning = function(w) conditionMessage(w))
  if (!identical(rc, 0L) && !identical(rc, 0)) { unlink(part); return(paste("download failed:", rc)) }
  if (!file.exists(part) || file.size(part) == 0) { unlink(part); return("download failed: empty file") }
  if (!file.rename(part, dest)) return(paste("could not rename", part))
  ""
}

archive_kind <- function(path) {
  con <- file(path, "rb"); b <- readBin(con, "raw", 300); close(con)
  if (length(b) >= 4 && identical(b[1:4], as.raw(c(0x50, 0x4b, 0x03, 0x04)))) return("zip")
  if (length(b) >= 4 && identical(b[1:4], charToRaw("Rar!"))) return("rar")
  if (length(b) >= 262 && identical(b[258:262], charToRaw("ustar"))) return("tar")
  "other"
}

run_tool <- function(tool, archive, member, exdir) {
  args <- switch(tool,
    unzip  = c("-o", "-qq", shQuote(archive), shQuote(member), "-d", shQuote(exdir)),
    bsdtar = c("-x", "-f", shQuote(archive), "-C", shQuote(exdir), shQuote(member)),
    "7z"   = c("x", "-y", shQuote(paste0("-o", exdir)), shQuote(archive), shQuote(member)),
    unar   = c("-q", "-f", "-o", shQuote(exdir), shQuote(archive), shQuote(member)))
  rc <- suppressWarnings(system2(tool, args, stdout = FALSE, stderr = FALSE))
  identical(as.integer(rc), 0L) && file.exists(file.path(exdir, member))
}

# extract one member to out_file; returns "" on success, otherwise the error text
extract_member <- function(archive, member, out_file) {
  exdir <- tempfile("x_", tmpdir = CACHE)
  dir.create(exdir)
  on.exit(unlink(exdir, recursive = TRUE), add = TRUE)
  kind <- archive_kind(archive)
  got <- file.path(exdir, member)
  ok <- FALSE
  if (kind == "zip") ok <- tryCatch({ utils::unzip(archive, files = member, exdir = exdir); TRUE }, error = function(e) FALSE, warning = function(w) FALSE)
  if (kind == "tar") ok <- tryCatch({ utils::untar(archive, files = member, exdir = exdir, tar = "internal"); TRUE }, error = function(e) FALSE, warning = function(w) FALSE)
  if (!ok) unlink(got)   # drop any partial output before trying external tools
  if (!file.exists(got)) {
    for (tool in c("unzip", "bsdtar", "7z", "unar")) {
      if (tool == "unzip" && kind != "zip") next
      if (!nzchar(Sys.which(tool))) next
      cat("    extracting with", tool, "\n")
      if (run_tool(tool, archive, member, exdir)) break
      unlink(got)
    }
  }
  if (!file.exists(got)) return(sprintf("could not extract '%s' from %s (%s archive; tried utils and unzip/bsdtar/7z/unar)", member, basename(archive), kind))
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  if (!file.copy(got, out_file, overwrite = TRUE)) return(paste("could not copy extracted member to", out_file))
  ""
}

gunzip_file <- function(gz, out_file) {
  con_in <- gzfile(gz, "rb"); con_out <- file(out_file, "wb")
  repeat { b <- readBin(con_in, "raw", n = 1e7); if (!length(b)) break; writeBin(b, con_out) }
  close(con_in); close(con_out)
}

cache_path <- function(url) {
  nm <- gsub("[^A-Za-z0-9._-]+", "_", sub("^https?://", "", url))
  file.path(CACHE, substr(nm, max(1, nchar(nm) - 149), nchar(nm)))
}

# ---- main loop -------------------------------------------------------------------------------------
# whole-file rows first, so an archive that is itself a pipeline file (X003 weeks zip, P8_033 sup001.zip) is reused
src <- src[order(nzchar(src$archive_member), src$rel_path), ]
status <- setNames(rep(NA_character_, nrow(src)), src$rel_path)
problem <- setNames(rep("", nrow(src)), src$rel_path)

for (i in seq_len(nrow(src))) {
  row <- src[i, ]
  dest <- file.path(DATA, row$rel_path)
  cat(sprintf("[%d/%d] %s (%s)\n", i, nrow(src), row$rel_path, row$source_type))

  if (file.exists(dest)) {
    got <- md5_of(dest)
    if (identical(got, row$md5)) { status[row$rel_path] <- "present"; next }
    stop(sprintf("existing file has the wrong md5: %s\n  got %s, expected %s (%s bytes on disk, %s expected). Inspect or delete it, then re-run.",
                 dest, got, row$md5, file.size(dest), row$bytes))
  }
  if (row$source_type %in% c("prio_login", "unknown")) { status[row$rel_path] <- "manual"; next }
  if (row$source_type %in% c("derived", "package")) { status[row$rel_path] <- "not_fetched_missing"; next }

  url <- fetch_url(row)
  tmp <- paste0(dest, ".fetch")
  unlink(tmp)
  err <- ""
  if (!nzchar(row$archive_member)) {
    if (row$rel_path %in% GUNZIP) {
      gz <- paste0(dest, ".gz")
      err <- download_to(url, gz)
      if (!nzchar(err)) { gunzip_file(gz, tmp); unlink(gz) }
    } else {
      err <- download_to(url, tmp)
    }
  } else {
    # reuse a pipeline file that is the same archive, else the cache, else download into the cache
    same <- src$rel_path[src$source == row$source & !nzchar(src$archive_member)]
    arch <- if (length(same) && file.exists(file.path(DATA, same[1]))) file.path(DATA, same[1]) else cache_path(url)
    if (!file.exists(arch)) err <- download_to(url, arch)
    if (!nzchar(err)) err <- extract_member(arch, row$archive_member, tmp)
  }
  if (nzchar(err)) {
    unlink(tmp)
    status[row$rel_path] <- "failed"; problem[row$rel_path] <- err
    cat("    FAILED:", err, "\n")
    next
  }

  got <- md5_of(tmp)
  if (!identical(got, row$md5)) {
    bad <- paste0(dest, ".bad")
    file.rename(tmp, bad)
    stop(sprintf(paste0("md5 mismatch for %s\n  got %s (%s bytes), expected %s (%s bytes)\n  source: %s %s\n",
                        "  downloaded copy kept at %s. If it came from a cached archive, delete the archive in %s and re-run."),
                 row$rel_path, got, file.size(bad), row$md5, row$bytes, row$source, row$archive_member, bad, CACHE))
  }
  stopifnot(file.rename(tmp, dest))
  status[row$rel_path] <- "fetched"
}

# ---- summary ---------------------------------------------------------------------------------------
stopifnot(!anyNA(status))
cat("\nSummary:\n"); print(table(status))
todo <- src[status[src$rel_path] %in% c("manual", "failed", "not_fetched_missing"), ]
if (nrow(todo)) {
  cat("\n", nrow(todo), " file(s) could not be fetched automatically. Get each one by hand and place it at\n",
      "  ", DATA, "/<rel_path>\n(for an archive, extract archive_member; a::b means member b inside the nested zip a).\n\n", sep = "")
  for (j in seq_len(nrow(todo))) {
    r <- todo[j, ]
    why <- switch(status[[r$rel_path]], manual = r$source_type, failed = paste("failed:", problem[[r$rel_path]]),
                  not_fetched_missing = paste(r$source_type, "file missing; build it with the step named in source"))
    cat(sprintf("- %s  [id %s]\n    %s\n    source: %s\n", r$rel_path, r$id, why, r$source))
    if (nzchar(r$archive_member)) cat("    archive_member:", r$archive_member, "\n")
  }
  cat("\nThen re-run this script to verify the md5 of every file.\n")
  stop(nrow(todo), " file(s) still missing; see the list above. Stage 5 cannot run until every file is present.")
} else {
  cat("All", nrow(src), "files present with the expected md5.\n")
}
