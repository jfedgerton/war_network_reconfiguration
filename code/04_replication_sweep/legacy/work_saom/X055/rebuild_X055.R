# X055 Li, Bradshaw, Clary & Cranmer (2017) Sci Adv -- Table 1 Model 1 TERGM (btergm, bootstrapped MPLE, R = 1000).
# No replication archive (data 'available upon request'). REBUILD from:
#   * Pauls & Cranmer (2017) network data as posted in Pomeroy et al.'s Dataverse (DVN/45GXK5, data/pauls_cranmer/original):
#     36 yearly networks 1965-2000, fixed 125-state node set, MID onset network with 'democ' vertex attribute,
#     direct contiguity, capability ratio, trade dependence -- the same covariates Li et al. list (Cranmer co-authors both).
#   * k-degree alliance matrices computed here from ATOP defense pacts (peacesciencer::atop_alliance, ATOP v5) as
#     shortest-path distance in each year's undirected defense-pact network (1,2,3,4 degrees), plus contiguity interactions.
# Model 1: edges + altkstar(2) + cycle(4) + gwesp(0) + joint democracy + contiguity + capratio + tradedep
#          + deg1..deg4 + contig x deg1..deg4. Seed 123 (authors' seed unknown).
# History (edgecov): mH = hist_dyad_w15 + hist_mono_w15; m26 = w26 versions; mP = bm_cinc_w15 placebo.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
WD <- file.path(OUT, "work_saom/X055"); dir.create(WD, recursive = TRUE, showWarnings = FALSE)
setwd(WD)   # res_X055_R<RB>.rds is written here (pipeline output)
suppressPackageStartupMessages({library(btergm); library(network); library(igraph); library(dplyr); library(peacesciencer)})
cat(R.version.string, " btergm", as.character(packageVersion("btergm")), " peacesciencer", as.character(packageVersion("peacesciencer")), "\n")
args <- commandArgs(TRUE); RB <- if (length(args)) as.integer(args[1]) else 1000L
source(file.path(CODE, "hist/add_history.R"))

PC <- file.path(ROOT, "work_saom/X055/pc")
mid <- readRDS(file.path(PC, "mid_nw.rds")); contig_nw <- readRDS(file.path(PC, "directcont_nw.rds"))
capratio <- readRDS(file.path(PC, "capratio.rds")); tradedep <- readRDS(file.path(PC, "tradedep.rds"))
years <- 1965:2000; stopifnot(length(mid) == 36)
NM <- lapply(mid, function(x) network::get.vertex.attribute(x, "vertex.names"))   # node set varies by year
for (i in seq_along(years)) stopifnot(all(rownames(capratio[[i]]) == NM[[i]]), all(rownames(tradedep[[i]]) == NM[[i]]),
                                      all(network::get.vertex.attribute(contig_nw[[i]], "vertex.names") == NM[[i]]))
cat("nodes per year:", lengths(NM), "\n")
contig <- lapply(contig_nw, function(x) { m <- as.matrix(x); m[is.na(m)] <- 0; m })

## ---- k-degree alliance matrices from ATOP defense pacts ----
data(atop_alliance, package = "peacesciencer")
atop <- atop_alliance %>% filter(atop_defense == 1, year %in% years)
deg <- lapply(1:4, function(k) vector("list", length(years)))
for (i in seq_along(years)) {
  e <- atop %>% filter(year == years[i]) %>% transmute(a = as.character(ccode1), b = as.character(ccode2))
  nm <- NM[[i]]; e <- e %>% filter(a %in% nm, b %in% nm)
  g <- graph_from_data_frame(e, directed = FALSE, vertices = data.frame(name = nm))
  g <- simplify(g)
  D <- distances(g, v = nm, to = nm)
  for (k in 1:4) { m <- (D == k) * 1; diag(m) <- 0; dimnames(m) <- list(nm, nm); deg[[k]][[i]] <- m }
}
a1 <- deg[[1]]; a2 <- deg[[2]]; a3 <- deg[[3]]; a4 <- deg[[4]]
ca1 <- Map(`*`, contig, a1); ca2 <- Map(`*`, contig, a2); ca3 <- Map(`*`, contig, a3); ca4 <- Map(`*`, contig, a4)
cat("mean share of dyads at degree 1..4 (1980):", sapply(1:4, function(k) mean(deg[[k]][[16]][upper.tri(deg[[k]][[16]])])), "\n")
cat("MID onset ties per year:", sapply(mid, network.edgecount), "\n")

## ---- history ----
hmat <- function(var) lapply(seq_along(years), function(i) {
  nm <- NM[[i]]; cc <- as.numeric(nm); N <- length(cc)
  g <- expand.grid(a = seq_len(N), b = seq_len(N)); g$c1 <- cc[g$a]; g$c2 <- cc[g$b]; g$year <- years[i]
  g <- add_history(g, "c1", "c2", "year"); v <- g[[var]]; v[is.nan(v)] <- NA
  if (var == "bm_cinc_w15") v[is.na(v)] <- mean(v, na.rm = TRUE)
  m <- matrix(v, N, N); diag(m) <- 0; dimnames(m) <- list(nm, nm); m
})
hd15 <- hmat("hist_dyad_w15"); hm15 <- hmat("hist_mono_w15"); hd26 <- hmat("hist_dyad_w26"); hm26 <- hmat("hist_mono_w26"); bmc <- hmat("bm_cinc_w15")

base_f <- "mid ~ edges + altkstar(2, fixed = TRUE) + cycle(4) + gwesp(0, fixed = TRUE) + nodemix('democ', -3) +
  edgecov(contig) + edgecov(capratio) + edgecov(tradedep) + edgecov(a1) + edgecov(a2) + edgecov(a3) + edgecov(a4) +
  edgecov(ca1) + edgecov(ca2) + edgecov(ca3) + edgecov(ca4)"
adds <- list(m0 = "", mH = " + edgecov(hd15) + edgecov(hm15)", m26 = " + edgecov(hd26) + edgecov(hm26)", mP = " + edgecov(bmc)")
out <- list()
for (m in names(adds)) {
  set.seed(123); t0 <- Sys.time()
  Rm <- if (m %in% c("m0","mH")) RB else min(RB, 100L)
  fit <- btergm(as.formula(paste0(base_f, adds[[m]])), R = Rm, parallel = "multicore", ncpus = 2, verbose = FALSE)
  s <- summary(fit); co <- coef(fit); bs <- fit@boot$t
  tn <- sub("\\[\\[i\\]\\]$", "", names(co))
  out[[m]] <- data.frame(model = m, term = tn, b = co, se_boot = apply(bs, 2, sd),
                         p_boot = apply(bs, 2, function(x) 2 * min(mean(x <= 0), mean(x >= 0))),
                         ci_lo = s[, 3], ci_hi = s[, 4], R = Rm, nobs = fit@nobs)
  cat(sprintf("== %s (R=%d, %.0f s)\n", m, Rm, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  print(out[[m]][, c("term","b","ci_lo","ci_hi")], digits = 3, row.names = FALSE)
}
saveRDS(bind_rows(out), sprintf("res_X055_R%d.rds", RB))
