# source() this. add_history(df, c1, c2, year) -> df + hist/benchmark columns
suppressPackageStartupMessages(library(dplyr))
.H <- readRDS(Sys.getenv("SWEEP_HIST"))   # data/processed/replication_sweep/hist/history_tables.rds (00_build_history.R)
stopifnot(!is.null(.H$exh))
add_history <- function(df, c1, c2, year) {
  stopifnot(all(c(c1, c2, year) %in% names(df)))
  key <- data.frame(.a = as.numeric(df[[c1]]), .b = as.numeric(df[[c2]]), year = as.numeric(df[[year]]))
  key$c1 <- pmin(key$.a, key$.b); key$c2 <- pmax(key$.a, key$.b)
  k <- key %>% left_join(.H$dyh, by = c("c1","c2","year")) %>%
    left_join(.H$sth %>% rename(sa15 = s15, sa26 = s26), by = c(".a" = "ccode", "year")) %>%
    left_join(.H$sth %>% rename(sb15 = s15, sb26 = s26), by = c(".b" = "ccode", "year")) %>%
    left_join(.H$exh, by = c("c1","c2","year")) %>%
    mutate(across(c(d15,d26,sa15,sa26,sb15,sb26,e15,e26), ~coalesce(.x, 0)))
  stopifnot(nrow(k) == nrow(df))
  df$hist_dyad_w15 <- log1p(k$d15); df$hist_dyad_w26 <- log1p(k$d26)
  # monadic: both states' own onsets minus their onsets in disputes against each other (participant level, own start year)
  stopifnot(all(k$sa15 + k$sb15 >= k$e15), all(k$sa26 + k$sb26 >= k$e26))
  df$hist_mono_w15 <- log1p(k$sa15 + k$sb15 - k$e15)
  df$hist_mono_w26 <- log1p(k$sa26 + k$sb26 - k$e26)
  # benchmarks: mean over t-1..t-5 of CINC ratio (stronger share) and joint democracy
  bm <- NULL
  for (L in 1:5) {
    x <- key %>% mutate(yl = year - L) %>%
      left_join(.H$nmc, by = c(".a" = "ccode", "yl" = "year")) %>% rename(ca = cinc) %>%
      left_join(.H$nmc, by = c(".b" = "ccode", "yl" = "year")) %>% rename(cb = cinc) %>%
      left_join(.H$pol, by = c(".a" = "ccode", "yl" = "year")) %>% rename(pa = polity2) %>%
      left_join(.H$pol, by = c(".b" = "ccode", "yl" = "year")) %>% rename(pb = polity2)
    bm[[L]] <- data.frame(cr = pmax(x$ca, x$cb) / (x$ca + x$cb), jd = as.numeric(x$pa >= 6 & x$pb >= 6))
  }
  df$bm_cinc_w15  <- rowMeans(sapply(bm, `[[`, "cr"), na.rm = TRUE)
  df$bm_democ_w15 <- rowMeans(sapply(bm, `[[`, "jd"), na.rm = TRUE)
  df$hist_in_gml_range <- key$year >= 1817 & key$year <= 2016
  df
}
