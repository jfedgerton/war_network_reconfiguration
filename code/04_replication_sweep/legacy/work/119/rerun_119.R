# Rerun 119: Maoz, Kuperman, Terris & Talmud 2006 JCR "Structural Equivalence and International Conflict"
# Headline: Table 2, MIDs Model 3: logit dichmid20 (MID under way) on lagged integrated SEq (lgintegse)
# + minreg302 + caprat301 + distance + peace years + 3 splines. Published (paper text via academia.edu):
# lgintegse -0.330 (0.067), constant -2.057 (0.113), n = 449,969. No do-file/log in the archive.
# Authors' spline vars midspline1-3 are kept exactly as supplied (they are built on the no-war counter, not peaceyrs).
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work/rerun_helpers_b118.R"))
library(haven)
WO <- file.path(OUT, "work/119"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work/119"); logf <- file.path(WO, "log_119.txt"); cat("", file = logf)
sink(logf, split = TRUE, append = TRUE)

dv <- "dichmid20"; tie <- "lgintegse"
base_rhs <- c("lgintegse","minreg302","caprat301","distance","peaceyrs","midspline1","midspline2","midspline3")
d <- as.data.frame(read_dta(file.path(ROOT, "_parent", "jcr2006replication.dta"),
                            col_select = all_of(c("year","statea","stateb","dyad", dv, base_rhs))))
stopifnot(nrow(d) == 758756)

# ---- 2. Reproduce (converged ML; SE type chosen by rule: closest match to published SEs of tie & constant) ----
f <- as.formula(paste(dv, "~", paste(base_rhs, collapse = "+")))
# fixest's IRLS does not settle here (fitted probs hit 0/1; peaceyrs=0 whenever a MID is under way ->
# quasi-separation). R's glm converges cleanly (strict epsilon), so estimation uses glm; SEs via sandwich,
# scaled like Stata (robust: n/(n-1); cluster: G/(G-1)).
library(sandwich)
glm_fit <- function(f, data, se = "hetero") {
  m <- glm(f, data = data, family = binomial("logit"), control = glm.control(epsilon = 1e-10, maxit = 100))
  stopifnot(m$converged)
  n <- nobs(m)
  V <- switch(se, iid = vcov(m), hetero = sandwich(m) * n / (n - 1),
              cl_dyad = vcovCL(m, cluster = data$dyad[as.integer(rownames(model.frame(m)))], type = "HC0", cadjust = TRUE))
  b <- coef(m); s <- sqrt(diag(V)); ct <- cbind(Estimate = b, SE = s, z = b / s, p = 2 * pnorm(-abs(b / s)))
  list(coeftable = ct, n = n, b = b) }
m_rep <- glm_fit(f, d, "iid")
stopifnot(m_rep$n == 449969)
rownames(d) <- NULL
ses <- sapply(c("iid", "hetero", "cl_dyad"), function(v) glm_fit(f, d, v)$coeftable[c(tie, "(Intercept)"), 2])
cat("\n==== REPRODUCTION ====\n"); print(m_rep$b); cat("SEs by type (pub: tie .067, const .113):\n"); print(ses)
dist <- abs(ses[1, ] - 0.067) / 0.067 + abs(ses[2, ] - 0.113) / 0.113
se_type <- names(which.min(dist)); cat("SE type chosen:", se_type, "\n")
vc <- se_type
b_pub <- -0.330
m_rep <- glm_fit(f, d, vc); print(m_rep$coeftable)
b_repro <- unname(m_rep$b[tie]); p_repro <- m_rep$coeftable[tie, 4]
rel <- abs(b_repro - b_pub) / abs(b_pub)
reproduced <- ifelse(rel <= 0.10 & p_repro < .05, "REPRODUCED", "NOT_REPRODUCED")
cat(sprintf("b_pub=%g b_repro=%g rel.diff=%.4f p=%.3g -> %s\n", b_pub, b_repro, rel, p_repro, reproduced))

# ---- 3-4. History + battery ----
d <- add_history(d, "statea", "stateb", "year"); d$year_h <- d$year
ds <- battery_sample(d, c(dv, base_rhs, "dyad"))
n_pre1822 <- sum(complete.cases(d[, c(dv, base_rhs)]) & d$year < 1822)
fitter <- function(rhs, data) glm_fit(as.formula(paste(dv, "~", paste(rhs, collapse = "+"))), data, vc)
rm(m_rep); invisible(gc())
rownames(ds) <- NULL
R <- run_battery(ds, base_rhs, fitter, function(m) c(b = unname(m$coeftable[tie, 1]), se = unname(m$coeftable[tie, 2]), p = unname(m$coeftable[tie, 4])), logf)
# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "119", dat = ds, y = dv, tie = tie, ctrl = setdiff(base_rhs, tie), fe = NULL,
           fit = function(dat, add) { m <- fitter(c(base_rhs, add), dat)
             c(b = unname(m$coeftable[tie, 1]), se = unname(m$coeftable[tie, 2]), p = unname(m$coeftable[tie, 4])) })
}

notes <- sprintf("Headline lgintegse (lagged integrated SEq) Tbl2 MID Model3; b & N=449,969 match paper exactly; SE type=%s by rule (pub se .067). No archive code; spec from paper text. %d obs pre-1822 dropped; rest from benchmark missingness. Authors' DV is MID under way (dichmid20).",
                 se_type, n_pre1822)
write_result(119, "Maoz et al. 2006", reproduced, b_pub, b_repro, paste0("logit (glm), SE ", se_type), tie,
             449969, nrow(ds), R, "N", notes)
sink()
