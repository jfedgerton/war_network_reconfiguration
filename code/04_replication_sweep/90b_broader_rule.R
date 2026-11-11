# 90b_broader_rule.R -- appendix count under a broader inclusion rule (PI decision 2026-11-11).
# Adds to the pooled tie studies the 10 tie studies excluded as out of scope (refugee flows, P9_008, stays out):
#   N01, N02 (terms built from the MID network), P9_012, P9_036 (intervention in civil wars), P9_016 (joining ongoing wars),
#   P9_023, P9_065 (duration / termination), P9_041 (event-data outcome), 104 (claim-level unit), X088 (sequence model).
# Same rule as 90_summarise.R: significant at baseline at the paper's alpha_pub; changed = null with history or significant
# with the opposite sign. Writes results/replication_sweep/broader_rule_count.csv (counts only; no names).
suppressPackageStartupMessages(library(dplyr))
PIPE <- Sys.getenv("PIPE_ROOT", unset = here::here())
SW <- file.path(PIPE, "code/04_replication_sweep"); EST <- file.path(PIPE, "results/replication_sweep")
BROADER <- c("N01", "N02", "P9_012", "P9_036", "P9_016", "P9_023", "P9_065", "P9_041", "104", "X088")
reg <- read.csv(file.path(SW, "registry.csv"), colClasses = "character")
hd  <- reg %>% filter(study_headline_term == "Y", family == "TIE")
stopifnot(all(BROADER %in% hd$id), all(hd$pool_status[hd$id %in% BROADER] == "out_of_scope"))
sets <- list(pooled = hd$id[hd$pool_status == "main"], broader = c(hd$id[hd$pool_status == "main"], BROADER))
one <- function(ids) {
  h <- hd %>% filter(id %in% ids)
  e <- bind_rows(lapply(ids, function(i) { x <- readRDS(file.path(EST, paste0("est_", i, ".rds"))); x$row_id <- as.character(x$row_id); x }))
  d <- h %>% select(row_id, id, alpha_pub) %>% inner_join(e %>% select(row_id, b0, p0, b_hist, p_hist), by = "row_id") %>%
    mutate(alpha = as.numeric(alpha_pub), sig0 = p0 < alpha,
           changed = sig0 & (p_hist >= alpha | sign(b_hist) != sign(b0)),
           chg = 100 * (abs(b_hist) - abs(b0)) / abs(b0))   # % change in magnitude, as in 90_summarise.R
  stopifnot(nrow(d) == length(ids), all(is.finite(d$b0)), all(is.finite(d$b_hist)))
  data.frame(studies = nrow(d), sig_baseline = sum(d$sig0), changed = sum(d$changed),
             share_changed = round(sum(d$changed) / sum(d$sig0), 3), median_chg = round(median(d$chg[d$sig0]), 2))
}
out <- bind_rows(lapply(names(sets), function(k) cbind(set = k, one(sets[[k]]))))
print(out, row.names = FALSE)
write.csv(out, file.path(EST, "broader_rule_count.csv"), row.names = FALSE)
cat("wrote broader_rule_count.csv\n")
