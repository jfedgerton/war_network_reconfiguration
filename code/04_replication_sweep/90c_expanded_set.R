# 90c_expanded_set.R -- [2026-11-14] appendix footnote: results for the expanded set (pooled studies plus the studies whose
# outcome codes every year of a dispute, pool_status == "outcome_overlap"). Same rule as 90_summarise.R / 90b_broader_rule.R:
# significant at baseline at the paper's alpha_pub; changed = not significant with history, or significant with the opposite sign.
# Writes results/replication_sweep/expanded_set_count.csv (counts only).
suppressPackageStartupMessages(library(dplyr))
PIPE <- Sys.getenv("PIPE_ROOT", unset = here::here())
SW <- file.path(PIPE, "code/04_replication_sweep"); EST <- file.path(PIPE, "results/replication_sweep")
reg <- read.csv(file.path(SW, "registry.csv"), colClasses = "character")
hd  <- reg %>% filter(study_headline_term == "Y", pool_status %in% c("main", "outcome_overlap"))
stopifnot(!anyDuplicated(hd$id), sum(hd$pool_status == "outcome_overlap") > 0)
e <- bind_rows(lapply(hd$id, function(i) { x <- readRDS(file.path(EST, paste0("est_", i, ".rds"))); x$row_id <- as.character(x$row_id); x }))
d <- hd %>% select(row_id, id, family, pool_status, alpha_pub) %>% inner_join(e %>% select(row_id, b0, p0, b_hist, p_hist), by = "row_id") %>%
  mutate(alpha = as.numeric(alpha_pub), sig0 = p0 < alpha, changed = sig0 & (p_hist >= alpha | sign(b_hist) != sign(b0)),
         chg = 100 * (abs(b_hist) - abs(b0)) / abs(b0))
stopifnot(nrow(d) == nrow(hd), all(is.finite(d$b0)), all(is.finite(d$b_hist)))
one <- function(x) data.frame(studies = nrow(x), sig_baseline = sum(x$sig0), changed = sum(x$changed),
                              share_changed = round(sum(x$changed) / sum(x$sig0), 3), median_chg = round(median(x$chg[x$sig0]), 2))
out <- bind_rows(
  d %>% filter(pool_status == "main") %>% group_by(family) %>% group_modify(~ one(.x)) %>% mutate(set = "pooled"),
  d %>% filter(pool_status == "outcome_overlap") %>% group_by(family) %>% group_modify(~ one(.x)) %>% mutate(set = "overlap only"),
  d %>% group_by(family) %>% group_modify(~ one(.x)) %>% mutate(set = "expanded"),
  d %>% mutate(family = "ALL") %>% group_by(family) %>% group_modify(~ one(.x)) %>% mutate(set = "expanded")) %>%
  ungroup() %>% select(set, family, everything())
print(as.data.frame(out), row.names = FALSE)
write.csv(out, file.path(EST, "expanded_set_count.csv"), row.names = FALSE)
cat("wrote expanded_set_count.csv\n")
