# X016 Lektzian, Prins & Souva 2010 ISQ. Authors' archive: Dataverse doi:10.7910/DVN/ULAVUA (LPSDataArchive.tab + Table do-files).
# Design: ICOW claim-dyad-years (Western Hemisphere), Cox PH with shared (gamma) frailty on claim-dyad (stcox ..., shared()).
# Key terms (decided from abstract BEFORE estimating): joint democracy conditional on territory x rivalry:
#   demterival  (Table 4 Model 5, omitted = mixterival)     -> joint dem vs mixed, territorial claim + rivalry (Table 5 row 1: HR 1.05, n.s.)
#   demternonrival (Table 4 Model 7, omitted = mixternonrival) -> joint dem vs mixed, territorial claim, no rivalry
# CONTROL row: jntd6 in Table 2 Model 2 (territory x rivalry headline model).
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p4/common_p4_x8.R"))
suppressPackageStartupMessages(library(survival))
WO <- file.path(OUT, "work_p4/X016"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work_p4/X016")
sink(file.path(WO, "log_X016.txt"), split = TRUE)

d <- read.delim(file.path(W, "3185586"))
stopifnot(all(d$time == d$X_t, na.rm = TRUE))
d <- subset(d, X_st == 1 & year > 1900)
CTRL <- "icowsal + recno5 + recmid5 + contig150 + lnrelcap + allies"
ALL8 <- c("demterival","demnonterival","demternonrival","demnonternonriv","mixterival","mixnonterival","mixternonrival","mixnonternonriv")
rhs_t4 <- function(omit) paste(setdiff(ALL8, omit), collapse = " + ")
SPEC <- list(
  demterival     = list(rhs = paste(rhs_t4("mixterival"), "+", CTRL), lab = "Table 4 Model 5"),
  demternonrival = list(rhs = paste(rhs_t4("mixternonrival"), "+", CTRL), lab = "Table 4 Model 7"),
  jntd6          = list(rhs = paste("territory + thomriva + terival +", CTRL, "+ jntd6"), lab = "Table 2 Model 2"))

fit_cox <- function(rhs, extra, s, term, lab) {
  f <- as.formula(paste("Surv(X_t0, X_t, X_d) ~", rhs, extra, "+ frailty.gamma(ddyadidclaim2)"))
  m <- coxph(f, data = s, ties = "breslow")
  sm <- summary(m)$coefficients
  ct <- cbind(sm[, "coef"], sm[, "se(coef)"], sm[, "Chisq"], sm[, "p"]); rownames(ct) <- rownames(sm)
  ct <- ct[!grepl("frailty", rownames(ct)), , drop = FALSE]
  cat("frailty theta:", m$history[[1]]$theta, "\n")
  grab(ct, term, m$n, paste(lab, extra))
}

# 1. reproduce on authors' full sample (year > 1900)
R1 <- lapply(names(SPEC), function(k) { sp <- SPEC[[k]]; x <- d[complete.cases(d[, all.vars(as.formula(paste("~", sp$rhs)))]), ]
  fit_cox(sp$rhs, "", x, k, sp$lab)[[k]] }); names(R1) <- names(SPEC)
cat("demterival HR (pub 1.05, n.s.):", exp(R1$demterival$b), " p =", R1$demterival$p, "\n")

# 2. history + common sample (1901-2001 -> <=2011 fine)
d <- add_history(d, "ccode1", "ccode2", "year")
s_all <- subset(d, year >= 1822 & year <= 2011)

meta <- list(xid = "X016", short_cite = "Lektzian, Prins & Souva 2010", journal = "ISQ", year = 2010, doi = "10.1111/j.1468-2478.2010.00627.x")
ev <- "abstract: 'jointly democratic ... do not appear to reduce the hazard ... within the context of rivalry'"
for (k in names(SPEC)) {
  sp <- SPEC[[k]]
  s <- s_all[complete.cases(s_all[, c(all.vars(as.formula(paste("~", sp$rhs))), HISTV)]), ]
  M <- run_six(function(extra, x) fit_cox(sp$rhs, extra, x, k, sp$lab), s)
  b_pub <- if (k == "demterival") log(1.05) else NA
  rep <- if (k == "demterival") {
    if (abs(R1[[k]]$b - b_pub) / abs(b_pub) <= 0.10 && R1[[k]]$p >= .05) "REPRODUCED" else "NOT_REPRODUCED"
  } else "UNVERIFIED"
  nt <- sprintf(paste0("%s; Cox PH, Breslow ties, gamma shared frailty on claim-dyad (R survival::frailty.gamma ~ Stata stcox shared()). ",
    "b_pub only available for demterival (do-file Table5 note: HR 1.05 n.s.); repro HR %.3f p=%.3f. Paper PDF not accessible (OUP bot wall). ",
    "Sample year>1900 claim-dyad-years; full N=%d, common N=%d. Authors control recmid5 (MIDs over the claim in prior 5 yrs) and recno5 (recent failed negotiations) = claim-level history. ",
    "History merged on claim-year (dyad ccode1/ccode2)."), sp$lab, exp(R1[[k]]$b), R1[[k]]$p, R1[[k]]$n, M$m0[[k]]$n)
  write_p4(meta, k, "DEMOCRACY", if (k == "jntd6") "CONTROL" else "KEY", ev, rep, b_pub, R1[[k]]$b,
           "Cox PH shared gamma frailty (coxph)", R1[[k]]$n, M, "Y (recmid5: claim MIDs prior 5 yrs)", nt)
}
sink()
