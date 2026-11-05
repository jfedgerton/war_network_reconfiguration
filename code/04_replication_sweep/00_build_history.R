# 00_build_history.R -- conflict-history lookup tables used by every study (GML MIDs via peacesciencer).
# Writes data/processed/replication_sweep/hist/history_tables.rds (generated data; never written under data/raw):
# dyadic and monadic MID-onset counts in t-1..t-5 and t-2..t-6, the opponent-exclusion counts used for the monadic term,
# plus lagged CINC and polity2 for the placebo/benchmark controls. Deterministic; no seed needed.
PIPE <- Sys.getenv("PIPE_ROOT", unset = here::here())
OUTF <- Sys.getenv("SWEEP_HIST", unset = file.path(PIPE, "data/processed/replication_sweep/hist/history_tables.rds"))
dir.create(dirname(OUTF), recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages({library(dplyr); library(peacesciencer)})
stopifnot(packageVersion("peacesciencer") >= "1.2.0")

p <- gml_part %>% select(dispnum, ccode, sidea, styear)
# dyadic onsets: opposite-side pairs within a dispute; onset year = later participant start
pp <- inner_join(p, p, by = "dispnum", suffix = c("1", "2"), relationship = "many-to-many") %>%
  filter(sidea1 != sidea2, ccode1 < ccode2) %>%
  transmute(dispnum, c1 = ccode1, c2 = ccode2, year = pmax(styear1, styear2)) %>% distinct()
dy <- pp %>% count(c1, c2, year, name = "d_on")
st <- p %>% distinct(dispnum, ccode, styear) %>% count(ccode, year = styear, name = "s_on")
stopifnot(nrow(dy) > 3000, nrow(st) > 3000)
win <- function(tab, keys, var, a, b) {
  out <- NULL
  for (k in a:b) out <- bind_rows(out, tab %>% mutate(year = year + k))
  out %>% group_by(across(all_of(keys)), year) %>% summarise(v = sum(.data[[var]]), .groups = "drop")
}
dyh <- full_join(win(dy, c("c1", "c2"), "d_on", 1, 5) %>% rename(d15 = v),
                 win(dy, c("c1", "c2"), "d_on", 2, 6) %>% rename(d26 = v), by = c("c1", "c2", "year")) %>%
  mutate(across(c(d15, d26), ~ coalesce(.x, 0)))
# opponent exclusion for the monadic term of dyad (a, b): each state's own onsets (counted at its own start year, as in
# `st`) in disputes where the other state of the dyad is an opponent. Monadic history of (a, b) = s(a) + s(b) - e(a, b).
# (Before 2026-10-01 the code subtracted 2 * d, the dyadic onset count at the later start year, which misallocates
# disputes whose participants joined in different years; Codex audit, area 4.)
own <- p %>% distinct(dispnum, ccode, sidea, styear)
ex <- inner_join(own, own %>% select(dispnum, opp = ccode, sidea_o = sidea), by = "dispnum", relationship = "many-to-many") %>%
  filter(sidea != sidea_o) %>% distinct(dispnum, ccode, styear, opp) %>%
  transmute(c1 = pmin(ccode, opp), c2 = pmax(ccode, opp), year = styear) %>% count(c1, c2, year, name = "e_on")
stopifnot(sum(ex$e_on) >= 2 * nrow(pp))   # one own onset per side of each opposing dispute-dyad (more with re-entries)
exh <- full_join(win(ex, c("c1", "c2"), "e_on", 1, 5) %>% rename(e15 = v),
                 win(ex, c("c1", "c2"), "e_on", 2, 6) %>% rename(e26 = v), by = c("c1", "c2", "year")) %>%
  mutate(across(c(e15, e26), ~ coalesce(.x, 0)))
sth <- full_join(win(st, "ccode", "s_on", 1, 5) %>% rename(s15 = v),
                 win(st, "ccode", "s_on", 2, 6) %>% rename(s26 = v), by = c("ccode", "year")) %>%
  mutate(across(c(s15, s26), ~ coalesce(.x, 0)))
nmc <- cow_nmc %>% select(ccode, year, cinc)
pol <- ccode_democracy %>% select(ccode, year, polity2)
H <- list(dyh = dyh, sth = sth, exh = exh, nmc = nmc, pol = pol, gml_range = range(p$styear))
stopifnot(!anyDuplicated(dyh[, c("c1", "c2", "year")]), !anyDuplicated(sth[, c("ccode", "year")]),
          !anyDuplicated(exh[, c("c1", "c2", "year")]))
# the exclusion can never exceed the two states' own onsets
chk <- exh %>% left_join(sth %>% rename(a = ccode, sa15 = s15, sa26 = s26), by = c("c1" = "a", "year")) %>%
  left_join(sth %>% rename(b = ccode, sb15 = s15, sb26 = s26), by = c("c2" = "b", "year"))
stopifnot(all(chk$e15 <= coalesce(chk$sa15, 0) + coalesce(chk$sb15, 0)), all(chk$e26 <= coalesce(chk$sa26, 0) + coalesce(chk$sb26, 0)))

# If a reference copy exists (the tables the published sweep used), the rebuild must match it exactly.
REF <- Sys.getenv("SWEEP_HIST_REF", unset = "")
if (nzchar(REF) && file.exists(REF)) {
  R0 <- readRDS(REF)
  for (k in setdiff(intersect(names(H), names(R0)), "gml_range")) {   # exh is new (2026-10-01); older references lack it
    a <- as.data.frame(H[[k]]); b <- as.data.frame(R0[[k]])
    a <- a[do.call(order, a), ]; b <- b[do.call(order, b[names(a)]), names(a)]; rownames(a) <- rownames(b) <- NULL
    if (!isTRUE(all.equal(a, b, check.attributes = FALSE))) stop("history table '", k, "' differs from reference ", REF)
  }
  cat("history tables match reference\n")
}
saveRDS(H, OUTF)
cat("wrote", OUTF, "\n")
