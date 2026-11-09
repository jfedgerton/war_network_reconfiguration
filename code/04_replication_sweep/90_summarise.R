# 90_summarise.R -- pooled statistics from est_<id>.rds.
# One headline term per study (registry study_headline_term == "Y"); rows with in_headline == "N" are not pooled
# (methods papers' reanalyses, non-comparable ports, headline terms not significant at baseline, out-of-scope outcomes;
# reason in methods_note). Writes results/replication_sweep/summary_stats.rds (aggregates, incl. $topline) and
# results/replication_sweep/study_level_by_id.csv (one row per pooled study, keyed by study id only; no names or citations).
# Significance uses each paper's own threshold for its headline claim (registry alpha_pub: 0.10 for one-tailed .05 tests or
# headline claims made at .10, else 0.05; coding and evidence in alpha_pub_coding.csv).
# Topline ("changed") = significant at baseline and, with history, either not significant or significant with the opposite
# sign. It is reported for all pooled studies and for studies whose authors did not already control for prior conflict
# (registry author_history in B_duration / C_none), with the already-controlled group (A_dyadic / A_monadic) alongside.
# Benchmarks are reported on the same scale as the history change: 100 * (|b_bench| / |b0| - 1) for same-sign estimates,
# i.e. 100 * (b_bench / b0 - 1) = sign(b0) * pct, since pct = 100 * (b_bench - b0) / |b0|.
suppressPackageStartupMessages(library(dplyr))
set.seed(123)
PIPE <- Sys.getenv("PIPE_ROOT", unset = here::here())
SW <- file.path(PIPE, "code/04_replication_sweep"); EST <- file.path(PIPE, "results/replication_sweep")
RUN_APPENDIX_FAMILIES <- as.logical(Sys.getenv("RUN_APPENDIX_FAMILIES", "TRUE"))
reg <- read.csv(file.path(SW, "registry.csv"), colClasses = "character")
ids <- unique(reg$id[reg$family == "TIE" | RUN_APPENDIX_FAMILIES])
f <- file.path(EST, paste0("est_", ids, ".rds"))
if (!all(file.exists(f))) stop("missing estimates for: ", paste(ids[!file.exists(f)], collapse = ", "))
extra <- setdiff(sub("^est_(.*)\\.rds$", "\\1", list.files(EST, pattern = "^est_.*\\.rds$")), ids)
if (length(extra)) stop("est files for studies not in this run: ", paste(extra, collapse = ", "))
est <- bind_rows(lapply(f, function(x) { e <- readRDS(x); e$row_id <- as.character(e$row_id); e }))
stopifnot(!anyDuplicated(est$row_id), setequal(est$row_id, reg$row_id[reg$id %in% ids]))

stopifnot(all(reg$alpha_pub %in% c("0.05", "0.10")))   # per-paper threshold for the headline claim (alpha_pub_coding.csv)
# [2026-11-08] pool_status (headline rows): main = pooled (in_headline Y); the others are reported only in the all-attempts
# counts below and the appendix all-attempts figure (91_figures.R). Values set by apply_pool_status_20261108.py.
POOL_STATUS <- c("main", "author_null", "rebuild_outside_rule", "not_reproduced", "null_common_sample",
                 "not_comparable", "methods_reanalysis", "out_of_scope", "duplicate")
hr <- reg$study_headline_term == "Y"
stopifnot("pool_status" %in% names(reg), all(reg$pool_status[hr] %in% POOL_STATUS), all(reg$pool_status[!hr] == ""),
          all((reg$pool_status[hr] == "main") == (reg$in_headline[hr] == "Y")))
d <- reg %>% filter(id %in% ids, study_headline_term == "Y", in_headline == "Y") %>%
  select(row_id, id, family, tie_type, citations_openalex, alpha_pub) %>%
  inner_join(est %>% select(-id, -key_term), by = "row_id") %>%
  mutate(alpha_pub = as.numeric(alpha_pub),
         cites = suppressWarnings(as.numeric(citations_openalex)),
         chg = 100 * (abs(b_hist) - abs(b0)) / abs(b0),
         placebo_chg = sign(b0) * placebo_pct, jointdem_chg = sign(b0) * jointdem_pct,
         # significance judged at the paper's own threshold, at baseline and with history
         sig0 = p0 < alpha_pub, sigh = p_hist < alpha_pub,
         # standard error implied by the estimate and its two-sided p-value (|b| / |z|); exact for Wald z tests, close for
         # t tests with many degrees of freedom. The estimators' own standard errors for the hooked studies are in sens_<id>.rds.
         se_ratio = (abs(b_hist) / -qnorm(pmax(p_hist, 1e-300) / 2)) / (abs(b0) / -qnorm(pmax(p0, 1e-300) / 2)),   # lower tail: finite for tiny p
         lose = sig0 & !sigh, flip = sig0 & sigh & sign(b_hist) != sign(b0))
stopifnot(!anyDuplicated(d$id), all(is.finite(d$b0)), all(is.finite(d$b_hist)), !anyNA(d$p0), !anyNA(d$p_hist))

boot_med <- function(x, R = 2000) { stopifnot(all(is.finite(x))); q <- quantile(replicate(R, median(sample(x, replace = TRUE))), c(.025, .975)); unname(q) }
summ <- function(x) {
  ci <- boot_med(x$chg)
  data.frame(studies = nrow(x), median_chg = median(x$chg), ci_lo = ci[1], ci_hi = ci[2],
             share_shrink = mean(abs(x$b_hist) < abs(x$b0)),
             se_n = sum(is.finite(x$se_ratio)), median_se_ratio = median(x$se_ratio[is.finite(x$se_ratio)]),
             share_se_up5 = mean(x$se_ratio[is.finite(x$se_ratio)] > 1.05), share_se_down5 = mean(x$se_ratio[is.finite(x$se_ratio)] < 0.95),
             sig_baseline = sum(x$sig0), lose_sig = sum(x$lose),
             share_lose = sum(x$lose) / sum(x$sig0), sign_flips = sum(x$flip),
             placebo_n = sum(!is.na(x$placebo_chg)), placebo_median_chg = median(x$placebo_chg, na.rm = TRUE),
             jointdem_n = sum(!is.na(x$jointdem_chg)), jointdem_median_chg = median(x$jointdem_chg, na.rm = TRUE),
             cites_n = sum(!is.na(x$cites[x$sig0])),
             cite_share_lose = sum(x$cites[x$lose], na.rm = TRUE) / sum(x$cites[x$sig0], na.rm = TRUE))
}
S <- bind_rows(cbind(group = "All", summ(d)),
               bind_rows(lapply(split(d, d$family), function(x) cbind(group = x$family[1], summ(x)))))
saveRDS(list(by_family = S, by_tie_type = d %>% filter(family == "TIE") %>% group_by(tie_type) %>%
               summarise(studies = n(), median_chg = median(chg), lose_sig = sum(lose), sig_baseline = sum(sig0), .groups = "drop")),
        file.path(EST, "summary_stats.rds"))
print(S, digits = 3, row.names = FALSE)

# ---- topline: share of significant headline terms that change (null or opposite-sign) once history is added ----
ah <- reg %>% filter(study_headline_term == "Y") %>% select(row_id, author_history, interaction_hand)
d <- d %>% left_join(ah, by = "row_id") %>%
  mutate(changed = sig0 & (!sigh | sign(b_hist) != sign(b0)),
         author_ctrl = author_history %in% c("A_dyadic", "A_monadic"))
stopifnot(all(d$author_history %in% c("A_dyadic", "A_monadic", "B_duration", "C_none")),
          all(d$interaction_hand %in% c("Y", "N")), all(d$sig0))   # pooled set = significant at baseline (registry rule)
top1 <- function(x, set) data.frame(set = set, sig_baseline = sum(x$sig0), changed = sum(x$changed),
                                    to_null = sum(x$lose), sig_flip = sum(x$flip), share_changed = mean(x$changed),
                                    median_chg = median(x$chg))
topf <- function(x, g) rbind(cbind(group = g, top1(x, "all pooled")),
                             cbind(group = g, top1(x[!x$author_ctrl, ], "authors without history control")),
                             cbind(group = g, top1(x[x$author_ctrl, ], "authors with history control")))
TOP <- bind_rows(topf(d, "All"), bind_rows(lapply(split(d, d$family), function(x) topf(x, x$family[1])))) %>%
  filter(sig_baseline > 0)   # e.g. no capitalism study had an author history control
stopifnot(all(TOP$changed == TOP$to_null + TOP$sig_flip))
S_list <- readRDS(file.path(EST, "summary_stats.rds")); S_list$topline <- TOP
saveRDS(S_list, file.path(EST, "summary_stats.rds"))
cat("\nTOPLINE (changed = significant at baseline -> null, or significant with opposite sign)\n")
print(TOP, digits = 3, row.names = FALSE)

# ---- [2026-11-08] all attempts: pooled studies plus the attempts kept out of the pool ----
# Included: main, rebuild_outside_rule (rebuild fails the round-2 rule: sign, significance or 0.5-2x magnitude), not_reproduced, null_common_sample (headline term not
# significant on the common sample before history is added), author_null (author argues a null; asked only whether the term is
# significant with history). Not included: not_comparable, methods_reanalysis, out_of_scope, duplicate.
# "changed" uses the same rule as the topline; "gain" = not significant at baseline, significant with history.
ATT <- c("main", "rebuild_outside_rule", "not_reproduced", "null_common_sample", "author_null")
da <- reg %>% filter(id %in% ids, study_headline_term == "Y", pool_status %in% ATT) %>%
  select(row_id, id, family, pool_status, alpha_pub) %>%
  inner_join(est %>% select(row_id, b0, p0, b_hist, p_hist), by = "row_id") %>%
  mutate(alpha_pub = as.numeric(alpha_pub), sig0 = p0 < alpha_pub, sigh = p_hist < alpha_pub,
         lose = sig0 & !sigh, flip = sig0 & sigh & sign(b_hist) != sign(b0), changed = lose | flip, gain = !sig0 & sigh)
stopifnot(!anyDuplicated(da$id), all(is.finite(da$b0)), all(is.finite(da$b_hist)), !anyNA(da$p0), !anyNA(da$p_hist),
          all(!da$sig0[da$pool_status == "null_common_sample"]),
          setequal(da$row_id[da$pool_status == "main"], d$row_id))
att1 <- function(x) data.frame(attempts = nrow(x), sig_baseline = sum(x$sig0), changed = sum(x$changed), to_null = sum(x$lose),
                               sig_flip = sum(x$flip), gain_sig = sum(x$gain))
AT <- bind_rows(
  da %>% group_by(family, pool_status) %>% group_modify(~ att1(.x)) %>% ungroup(),
  # every attempt whose author claims an effect (author_null left out): main + failed rebuilds + not reproduced + null at baseline
  da %>% filter(pool_status != "author_null") %>% group_by(family) %>% group_modify(~ att1(.x)) %>% ungroup() %>%
    mutate(pool_status = "all attempts (excl. author null)")) %>%
  mutate(share_changed = ifelse(sig_baseline > 0, changed / sig_baseline, NA_real_)) %>% arrange(family, pool_status)
S_list <- readRDS(file.path(EST, "summary_stats.rds")); S_list$all_attempts <- AT
saveRDS(S_list, file.path(EST, "summary_stats.rds"))
cat("\nALL ATTEMPTS (share_changed = changed / significant at baseline, each paper's own threshold)\n")
print(as.data.frame(AT), digits = 3, row.names = FALSE)

# ---- study-level file for checking (keyed by id only) ----
L <- d %>% transmute(row_id, id, family, tie_type, author_history, interaction = interaction_hand, alpha_pub,
                     b0, p0, b_hist, p_hist, b_26, p_26, chg_pct = chg, placebo_chg, jointdem_chg, se_ratio,
                     sig_baseline = sig0, sig_hist = sigh, sign_flip = sign(b_hist) != sign(b0), changed, n) %>%
  arrange(family, id)
stopifnot(!any(c("short_cite", "doi", "methods_note") %in% names(L)))
write.csv(L, file.path(EST, "study_level_by_id.csv"), row.names = FALSE)
cat("wrote study_level_by_id.csv:", nrow(L), "rows\n")
cat("se_ratio = implied standard error with history / implied standard error at baseline (from b and p);",
    "share_se_up5 / share_se_down5 = share of studies whose standard error rises / falls by more than 5%\n")
cat("placebo/joint-democracy columns: *_n = studies with that benchmark (the rest could not estimate it);",
    "cites_n = significant-at-baseline studies with a citation count (cite_share_lose excludes the others)\n")
