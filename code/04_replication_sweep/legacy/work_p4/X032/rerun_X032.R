# X032 Boehmer & Sobek 2005 JPR "Violent Adolescence". Authors' data + do-file from Boehmer's UTEP page
# (http://utminers.utep.edu/crboehmer/ -> "Violent Adolescence JPR_logit and GEE.zip").
# MONADIC state-year design. Headline = first model in do-file: logit MID initiation, robust SE, BTSCS splines.
# Key terms (decided from abstract before estimating): logdev + logdev2 (curvilinear development).
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p4/common_p4_x8.R"))
WO <- file.path(OUT, "work_p4/X032"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work_p4/X032")
sink(file.path(WO, "log_X032.txt"), split = TRUE)

d <- haven::zap_labels(haven::read_dta(file.path(W, "Violent Adolescence Public Data.dta")))
d <- as.data.frame(d)
X <- "logdev + logdev2 + dem + popdens + popma10 + cap + openness + inityrs + initspl1 + initspl2 + initspl3"
f0 <- paste("initiate ~", X)
KEY <- c("logdev", "logdev2")
fit <- function(extra, s) fit_fx(as.formula(paste(f0, extra)), s, KEY, "logit", "hetero")

# 1. reproduce on authors' full sample (Stata: logit ..., robust)
full <- d[complete.cases(d[, all.vars(as.formula(f0))]), ]
r_full <- fit("", full)
n_orig <- r_full$logdev$n

# 2. monadic history (own-state MID onsets w15/w26) + own benchmarks; restrict 1822-2011
d <- add_history_monadic(d, "ccode", "year")
s <- subset(d, year >= 1822 & year <= 2011)
s <- s[complete.cases(s[, c(all.vars(as.formula(f0)), "hist_mono_w15", "hist_mono_w26", "bm_cinc_w15", "bm_democ_w15")]), ]
cat("orig N", n_orig, " common sample N", nrow(s), " dropped", n_orig - nrow(s), "\n")
stopifnot(nrow(s) > 0.8 * n_orig)

m0  <- fit("", s)
mD  <- fit("+ hist_mono_w15", s)       # monadic: no dyadic history; own-state history only
mDM <- mD
m26 <- fit("+ hist_mono_w26", s)
bC  <- fit("+ bm_cinc_w15", s)
bD  <- fit("+ bm_democ_w15", s)
M <- list(m0 = m0, mD = mD, mDM = mDM, m26 = m26, bC = bC, bD = bD)

meta <- list(xid = "X032", short_cite = "Boehmer & Sobek 2005", journal = "JPR", year = 2005, doi = "10.1177/0022343305049664")
ev <- "abstract: 'non-linear relationship where both extreme poverty and prosperity reduce this risk'"
nt <- sprintf(paste0("MONADIC state-year (history = own-state MID onsets t-1..t-5 / t-2..t-6; bC = own CINC mean, bD = own polity>=6 mean; mD=mDM). ",
  "Headline = do-file model 1 (logit initiate, robust SE, inityrs + 3 splines); key claim is inverted-U so logdev and logdev2 both reported. ",
  "Paper PDF not accessible -> b_pub unavailable, no author log -> UNVERIFIED. Full-sample repro N=%d (1816-1994); common 1822-1994 sample N=%d. ",
  "Authors control peace-years splines only (no MID-count history). dem is monadic democracy (not joint) and openness is trade openness -> not recorded as CONTROL rows."),
  n_orig, m0$logdev$n)
for (k in KEY) write_p4(meta, k, "CAPITALISM", "KEY", ev, "UNVERIFIED", NA, r_full[[k]]$b,
                        "logit, robust SE (fixest hetero), MONADIC", n_orig, M, "N (peace-year splines only)", nt)
sink()
