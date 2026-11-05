# Builder (run once, 2026-10-11): column subset of the authors' SV_rep.dta (Peterson & Drury 2011 JCR replication zip, Dropbox;
#   zip md5 2734f15b602ac71af5ae106d79e5eca6; member "JCR - sanct viol - replication/SV_rep.dta", 89,290,439 bytes) restricted to the
#   do-file's estimation filter (fcwongo == 0 & active == 1) and the Table 3 Model 1 variables. The 89 MB .dta is not kept (disk).
# Usage: Rscript build_subset_P9_079.R <path to SV_rep.dta> <output .rds>
a <- commandArgs(TRUE); stopifnot(length(a) == 2, file.exists(a[1]))
x <- haven::read_dta(a[1], col_select = c("dyad", "dyadid", "year", "fcwinit", "anysanctB", "sanctionAtoB", "allyAsBdpsio", "caprat", "jdem",
  "lndist", "s_un_glo", "coldwar", "fcwpceyrs", "fpysq", "fpycu", "fcwongo", "active"))
x <- as.data.frame(haven::zap_labels(x)); stopifnot(nrow(x) == 1223026)
x <- x[x$fcwongo %in% 0 & x$active %in% 1, ]; x$fcwongo <- NULL; x$active <- NULL
cat("rows kept:", nrow(x), "\n"); saveRDS(x, a[2], compress = "xz")
