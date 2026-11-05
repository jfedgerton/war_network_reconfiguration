# Rerun 100: Peterson & Quackenbush 2010 II -- Table 1 Model 1:
#   stcox imposed lowerdep lowerdepxdurmo higherdep interdep bothdem2 bothdem2xdurmo relpow relpowdif lndist lndistxdurmo strival
#     if unclearsettle==0, cluster(dyadid) nohr        (data already stset: _t0 _t _d, time in months since prior MID)
# Headline tie: lowerdep (lower of the two trade-dependence values); enters with lowerdep x duration (PH correction) -> constituent term.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work/rerun_common_g7.R"))
suppressPackageStartupMessages({library(haven); library(survival)})
W <- file.path(ROOT, "work/100")
WO <- file.path(OUT, "work/100"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)
sink(file.path(WO, "log_100.txt"), split = TRUE)

# 1. Authors' data (Dataverse 10.7910/DVN/KUUULZ rep_data_all_mids.dta, original format, stset stored)
d <- as.data.frame(zap_labels(read_dta(file.path(W, "allmids.dta"))))
stopifnot(nrow(d) == 30765, all(d$`_st` == 1))
names(d)[names(d) == "_t0"] <- "t0"; names(d)[names(d) == "_t"] <- "t1"; names(d)[names(d) == "_d"] <- "dd"
d <- d[d$unclearsettle == 0, ]
RHS <- "imposed + lowerdep + lowerdepxdurmo + higherdep + interdep + bothdem2 + bothdem2xdurmo + relpow + relpowdif + lndist + lndistxdurmo + strival"
MV <- c("t0","t1","dd","imposed","lowerdep","lowerdepxdurmo","higherdep","interdep","bothdem2","bothdem2xdurmo","relpow","relpowdif","lndist","lndistxdurmo","strival","dyadid")
TIE <- "lowerdep"

# Cox PH, Breslow ties (Stata default), cluster-robust SE with Stata's G/(G-1) small-sample factor
fitter <- function(add, dat) {
  dat <- dat[complete.cases(dat[, MV]), ]
  if (add != "") dat <- dat[complete.cases(dat[, trimws(strsplit(sub("^ \\+ ", "", add), "\\+")[[1]])]), ]
  f <- as.formula(paste0("Surv(t0, t1, dd) ~ ", RHS, add))
  m <- coxph(f, data = dat, ties = "breslow", cluster = dyadid)
  G <- length(unique(dat$dyadid))
  V <- m$var * G / (G - 1)
  b <- coef(m); se <- sqrt(diag(V)); names(se) <- names(b); z <- b / se; p <- 2 * pnorm(-abs(z))
  cat("\n---- ", deparse(f, width.cutoff = 500), "| N =", m$n, " events =", m$nevent, " clusters =", G, "\n")
  print(signif(cbind(b = b, se = se, z = z, p = p), 5))
  list(b = unname(b[TIE]), se = unname(se[TIE]), p = unname(p[TIE]), n = m$n)
}

# 2. Reproduce
cat("\n==== Reproduction: Table 1 Model 1 ====\n")
rp <- fitter("", d)
b_pub <- NA   # article closed-access, no author log in deposit
n_orig <- rp$n

# 3. History merge; sample years 1885-2000 so 1822-2011 restriction drops nothing
d <- add_history(d, "ccode1", "ccode2", "year")
d <- restrict_sample(d, MV)

# 4. Six models
R <- run_six(fitter, d)

# 5. Results
write_row(100, "Peterson & Quackenbush 2010", "UNVERIFIED", b_pub, rp$b, "Cox PH (Breslow), SE clustered on dyadid x G/(G-1)", "lowerdep (constituent; model has lowerdep x durmo)",
  n_orig, R$m0$n, R, ahc = "Y",
  notes = sprintf("Table 1 M1 recurrence design (post-MID peace spells; author has strival rivalry = Y, and sample conditioned on prior MID). hist_dyad_w15 mechanically includes the spell-starting MID when within 5 yrs. b_pub not obtainable (closed article, no log). Tie insignificant in reproduction (p~.56); coxph warns interdep scale (tiny values) but estimates finite. Dropped %d obs vs orig; %s.", n_orig - R$m0$n, RS_INFO))
sink()
