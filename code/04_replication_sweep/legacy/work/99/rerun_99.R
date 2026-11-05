# Rerun 99: Gartzke & Hewitt 2010 II, "International Crises and the Capitalist Peace"
# Headline: crisis-onset Model A5 (full controls; model used for Figure 2 simulations):
#   logit crisis1 ~ demlo2 demhi2 deplo capopnl2 gdppcln gdplncon sun2cati contig logdstab majpdyds alliesr
#                   lncaprt nafmeast crisyrs _spline1-3, cluster(dyadid)
# Tie: sun2cati (UN affinity / interest similarity). b_pub from authors' log_01.txt: -.7368035 (se .2432785), N = 254,018.
# Uses authors' temp_data_01.dta (saved by their do-file after btscs, so crisyrs/_spline* are theirs).
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work/rerun_common_b7.R"))
WO <- file.path(OUT, "work/99"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work/99"); sink(file.path(WO, "log_99.txt"), split = TRUE, append = TRUE)
C <- function(k) file.path(WO, paste0("cache_", k, ".rds"))

# 1. Load (cache a slim copy; the dta is 146MB)
# port: slim_99.rds is a staged input under DATA; if absent it is rebuilt from the dta into OUT
slim <- file.path(W, "slim_99.rds"); if (!file.exists(slim)) slim <- file.path(WO, "slim_99.rds")
if (!file.exists(slim)) {
  d <- haven::read_dta(file.path(W, "temp_data_01.dta"),
    col_select = c(statea, stateb, year, dyadid, crisis1, demlo2, demhi2, deplo, capopnl2, gdppcln, gdplncon, sun2cati,
                   contig, logdstab, majpdyds, alliesr, lncaprt, nafmeast, crisyrs, `_spline1`, `_spline2`, `_spline3`))
  d <- as.data.frame(haven::zap_labels(d)); names(d) <- sub("^_spline", "spline", names(d))
  saveRDS(d, slim)
} else d <- readRDS(slim)
stopifnot(all(d$crisis1 %in% c(0, 1, NA)))

f0 <- crisis1 ~ demlo2 + demhi2 + deplo + capopnl2 + gdppcln + gdplncon + sun2cati + contig + logdstab + majpdyds +
  alliesr + lncaprt + nafmeast + crisyrs + spline1 + spline2 + spline3
LG <- binomial(link = "logit"); T <- "sun2cati"; b_pub <- -0.7368035; n_orig <- 254018

# 2. Reproduce
r <- fit_tie(f0, d, LG, T, "dyadid", C("repro"))
stopifnot(r$n == n_orig)
repro_ok <- abs(r$b - b_pub) / abs(b_pub) <= 0.10 && r$p < .05
cat("b_pub", b_pub, "b_repro", r$b, "rel diff", abs(r$b - b_pub) / abs(b_pub), "\n")

# 3. History + restriction (statea<stateb undirected)
cat("statea<stateb share:", mean(d$statea < d$stateb), " year range:", range(d$year), "\n")
d <- d[complete.cases(d[, c(all.vars(f0), "dyadid")]), ]
d <- add_history(d, "statea", "stateb", "year")
ds <- d[d$year >= 1822 & d$year <= 2011, ]
ds <- ds[complete.cases(ds[, c("hist_dyad_w15", "hist_mono_w15", "hist_dyad_w26", "hist_mono_w26")]), ]
cat("common sample N:", nrow(ds), " dropped vs orig:", n_orig - nrow(ds), "\n"); rm(d); gc()

# 4. Models
m0  <- fit_tie(f0, ds, LG, T, "dyadid", C("m0"))
mD  <- fit_tie(update(f0, . ~ . + hist_dyad_w15), ds, LG, T, "dyadid", C("mD"))
mDM <- fit_tie(update(f0, . ~ . + hist_dyad_w15 + hist_mono_w15), ds, LG, T, "dyadid", C("mDM"))
m26 <- fit_tie(update(f0, . ~ . + hist_dyad_w26 + hist_mono_w26), ds, LG, T, "dyadid", C("m26"))
bC  <- fit_tie(update(f0, . ~ . + bm_cinc_w15), ds, LG, T, "dyadid", C("bC"))
bD  <- fit_tie(update(f0, . ~ . + bm_democ_w15), ds, LG, T, "dyadid", C("bD"))

# 5. Results
write_row(99, "Gartzke & Hewitt 2010", if (repro_ok) "REPRODUCED" else "NOT_REPRODUCED", b_pub, r$b,
  "logit, dyad-clustered SE", T, n_orig, m0$n, m0, mD, mDM, m26, bC, bD, ahc = "Y",
  notes = paste0("Crisis-onset Model A5 (ICB onset, not MID; history vars are MID-based). Author has crisis peace-years + splines. ",
                 n_orig - m0$n, " obs dropped vs orig. bC N=", bC$n, ", bD N=", bD$n, "."))
sink()
