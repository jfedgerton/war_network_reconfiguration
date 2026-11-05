# X093 Martin, Mayer & Thoenig 2008 RESTUD "Make Trade Not War?"
# Headline = Table 3 col. 4 (authors' mid_regs_tab3.do; log mid_regs_tab3.txt): full-sample logit with all controls
#   logit mid ltrl ztl lopenl muldis bildis peace votcol pol totwar diswar larea allied ldis comlng contig colony comcol fta_f gatt
#         yr1-yr51 con1-con20, robust cluster(dyad)      N = 223,788
# Key terms fixed from abstract before estimating: ltrl (bilateral openness, t-4) and lopenl (multilateral openness, t-4).
#   NOTE both are constituents of distance interactions (bildis = ldis x ltrl, muldis = ldis x lopenl) -> effect at ln distance = 0.
# CONTROL: pol ('sum of democracy indexes' - additive dyadic democracy, not a joint-democracy dummy).
# Year dummies -> year fixed effects in feglm (identical ML for logit with dummies). Authors already control lagged conflict:
#   con1-con20 = 'war [MID] in t-k' dummies, peace years, totwar.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p4/P4_common.R"))
WO <- file.path(OUT, "work_p4/X093"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work_p4/X093"); sink(file.path(WO, "log_X093.txt"), split = TRUE)
cons <- paste0("con", 1:20)
V <- c("year","ccodea","ccodeb","dyad","mid","ltrl","ztl","lopenl","muldis","bildis","peace","votcol","pol","totwar","diswar","larea","allied",
       "ldis","comlng","contig","colony","comcol","fta_f","gatt", cons)
d <- as.data.frame(haven::zap_labels(haven::read_dta(file.path(W, "mid_base_tab3.dta"), col_select = tidyselect::all_of(V))))
stopifnot(all(d$mid %in% c(0, 1, NA)), !anyDuplicated(d[, c("ccodea","ccodeb","year")]))
cat("years:", range(d$year), " N:", nrow(d), "\n")
rhs <- paste(c("ltrl","ztl","lopenl","muldis","bildis","peace","votcol","pol","totwar","diswar","larea","allied","ldis","comlng","contig",
               "colony","comcol","fta_f","gatt", cons), collapse = " + ")
f0 <- as.formula(paste("mid ~", rhs, "+ factor(year)"))   # year dummies (update()-safe; same ML as yr1-yr51)
K <- c("ltrl", "lopenl", "pol")
r <- fit_keys(f0, d, K, "dyad", file.path(WO, "cache_repro.rds")); cat("N repro", r$n, "vs pub 223788\n")
pub <- c(ltrl = -.2361011, lopenl = 1.520434, pol = -.283361); sig <- c(ltrl = FALSE, lopenl = TRUE, pol = TRUE)
ds <- common_sample(d, f0, "dyad", "ccodea", "ccodeb", "year"); rm(d); gc()
M <- run_six(f0, ds, K, "dyad", WO)
# ---- sensitivity hook (split only): runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  # fit_keys without a cache file (run_six caches each model by file name)
  sens_split(id = "X093", dat = ds,
             fit = function(dat, add) { r <- fit_keys(if (length(add)) update(f0, as.formula(paste(". ~ . +", paste(add, collapse = " + ")))) else f0, dat, K, "dyad")$k$lopenl
               c(b = unname(r[["b"]]), se = unname(r[["se"]]), p = unname(r[["p"]])) })
}
sc <- "Martin, Mayer & Thoenig 2008 RESTUD"; est <- "logit + year dummies, dyad-clustered SE (fixest feglm)"
hc <- "Y (con1-con20 lagged MID dummies t-1..t-20, peace years, # other wars)"
nt <- paste0("Table 3 col 4 (dyad-years 1950-2000). b_pub from authors' log. Common sample drops ", r$n - M$m0$n, " obs. ",
  "Authors already include 20 annual lagged-conflict dummies, so hist_* terms are largely redundant here. ")
write_p4("X093", sc, "TIE", "KEY", "ltrl", "abstract: 'probability of escalation is lower for countries that trade more bilaterally'",
  repro_status(pub[["ltrl"]], sig[["ltrl"]], r$k$ltrl), pub[["ltrl"]], r$k$ltrl[["b"]], est, r$n, M, hc,
  paste0(nt, "ltrl is the constituent of bildis (ltrl x ln distance): effect at ln distance = 0; published p=.074 (n.s. at .05)."))
write_p4("X093", sc, "TIE", "KEY", "lopenl", "abstract: 'countries more open to global trade have a higher probability of war'",
  repro_status(pub[["lopenl"]], sig[["lopenl"]], r$k$lopenl), pub[["lopenl"]], r$k$lopenl[["b"]], est, r$n, M, hc,
  paste0(nt, "lopenl is the constituent of muldis (lopenl x ln distance): effect at ln distance = 0."))
write_p4("X093", sc, "DEMOCRACY", "CONTROL", "pol", "control in all Table 3 models ('sum of democracy indexes')",
  repro_status(pub[["pol"]], sig[["pol"]], r$k$pol), pub[["pol"]], r$k$pol[["b"]], est, r$n, M, hc,
  paste0(nt, "pol = sum of the two states' democracy indexes (additive, not a joint-democracy dummy). bD NOT a placebo."))
sink()
