# 101 full chain (Bayesian dyadic model, custom MCMC sampler xiang_mcmc.cpp): 16 segments x 25,000 iterations per model
# (m0 and mDM interleaved), then the write step (burn-in 100,000). Logged runtime ~7.4 h in the sweep.
# Runs the ported shell driver work/101/est/run_chains.sh, which writes results_estaudit/101__smigoabi.csv under SWEEP_OUT.
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE), nzchar(Sys.getenv("SWEEP_DATA")), nzchar(Sys.getenv("SWEEP_OUT")), nzchar(Sys.getenv("SWEEP_HIST")))
Sys.setenv(RSCRIPT = file.path(R.home("bin"), "Rscript"))   # RLOCK optional (see run_chains.sh)
rc <- system2("bash", shQuote(file.path(CODE, "work/101/est/run_chains.sh")))
if (!identical(as.integer(rc), 0L)) stop("101 chain failed; see SWEEP_OUT/work/101/est/runner.log")
