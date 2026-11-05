# 94_study_flow_table.R -- appendix study-flow (attrition) table for the replication sweep.
# Input:  code/04_replication_sweep/study_flow.csv (one row per article that reached "screened eligible"; final_status and reason
#         per article) and registry.csv (pooled set = study_headline_term == "Y" & in_headline == "Y").
# Output: results/appendix/tab_study_flow.tex and results/appendix/tab_study_flow.csv (aggregate counts only; no study names).
# Checks: every pooled registry id maps to exactly one pooled article, and every pooled article has at least one pooled id.
suppressPackageStartupMessages(library(dplyr))
set.seed(123)
PIPE <- Sys.getenv("PIPE_ROOT", unset = here::here())
SW <- file.path(PIPE, "code/04_replication_sweep"); OUT <- file.path(PIPE, "results/appendix")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
fl <- read.csv(file.path(SW, "study_flow.csv"), colClasses = "character", na.strings = character(0))
reg <- read.csv(file.path(SW, "registry.csv"), colClasses = "character")
stopifnot(nrow(fl) > 0, !anyDuplicated(fl$key), all(fl$family %in% c("TIE", "DEMOCRACY", "CAPITALISM")))

# ---- reconcile pooled articles with the registry ----
pooled_ids <- reg$id[reg$study_headline_term == "Y" & reg$in_headline == "Y"]
stopifnot(!anyDuplicated(pooled_ids))
fl_ids <- strsplit(fl$sweep_id, ";", fixed = TRUE)
n_pooled_ids <- vapply(fl_ids, function(x) sum(x %in% pooled_ids), integer(1))
stopifnot(all(n_pooled_ids[fl$final_status == "pooled"] >= 1),     # pooled article -> a pooled registry term
          all(n_pooled_ids[fl$final_status != "pooled"] == 0))     # no pooled term hiding under another status
hits <- unlist(fl_ids[fl$final_status == "pooled"])
stopifnot(setequal(intersect(hits, pooled_ids), pooled_ids), !anyDuplicated(hits[hits %in% pooled_ids]))
terms_by_family <- table(reg$family[reg$study_headline_term == "Y" & reg$in_headline == "Y"])

# ---- status groups (order = table rows) ----
grp <- c(pooled = "Pooled in the sweep",
         excluded_not_significant_baseline = "Headline not significant at baseline (paper's own threshold)",
         excluded_out_of_scope = "Out of scope (outcome, unit or key term)",
         rebuilt_no_match = "Rebuilt from published sources; baseline did not match the article",
         no_data = "No replication data located",
         infeasible_design = "Design cannot carry a conflict-history control",
         no_ids_for_history = "Design cannot carry a conflict-history control",
         custom_measure_unavailable = "Design cannot carry a conflict-history control",
         excluded_estimand_not_portable = "Design cannot carry a conflict-history control",
         duplicate = "Duplicate or repeat specification",
         excluded_pi_decision = "Excluded by researcher decision (documented)",
         pending_pdf = "Pending (article or data in hand, port not finished)",
         pending_data = "Pending (article or data in hand, port not finished)",
         other = "Not yet adjudicated")
stopifnot(all(fl$final_status %in% names(grp)))
fl$group <- factor(grp[fl$final_status], levels = unique(grp))

tab <- fl %>% count(group, family) %>%
  tidyr::pivot_wider(names_from = family, values_from = n, values_fill = 0) %>%
  select(group, TIE, DEMOCRACY, CAPITALISM) %>% mutate(Total = TIE + DEMOCRACY + CAPITALISM) %>% arrange(group)
tot <- tab %>% summarise(across(c(TIE, DEMOCRACY, CAPITALISM, Total), sum)) %>% mutate(group = "Articles screened eligible", .before = 1)
tab <- bind_rows(mutate(tab, group = as.character(group)), tot)
stopifnot(tot$Total == nrow(fl))
print(tab, n = Inf)
write.csv(tab, file.path(OUT, "tab_study_flow.csv"), row.names = FALSE)

# ---- LaTeX ----
row_tex <- function(r) paste0(r$group, " & ", r$TIE, " & ", r$DEMOCRACY, " & ", r$CAPITALISM, " & ", r$Total, " \\\\")
body <- vapply(seq_len(nrow(tab) - 1), function(i) row_tex(tab[i, ]), character(1))
foot <- sprintf(paste0("Unit is the article. Pooled articles contribute one headline term per family; the pooled analyses use %d tie, ",
                       "%d democracy and %d capitalism terms. Significance at baseline is judged at the threshold each article ",
                       "uses for its headline claim."),
                terms_by_family[["TIE"]], terms_by_family[["DEMOCRACY"]], terms_by_family[["CAPITALISM"]])
tex <- c("\\begin{table}[htbp]\\centering", "\\caption{Replication sweep: study flow}", "\\label{tab:study_flow}",
         "\\begin{tabular}{lrrrr}", "\\hline",
         "Final status & Ties & Democracy & Capitalism & Total \\\\", "\\hline",
         body, "\\hline", row_tex(tab[nrow(tab), ]), "\\hline", "\\end{tabular}",
         paste0("\\par\\vspace{0.5em}\\begin{minipage}{0.9\\textwidth}\\footnotesize ", foot, "\\end{minipage}"),
         "\\end{table}")
writeLines(tex, file.path(OUT, "tab_study_flow.tex"))
cat("wrote", file.path(OUT, "tab_study_flow.tex"), "\n")
