## =============================================================================
## 06_section5_numbers.R -- the numbers quoted in the text of Section 5, computed
## from the analysis object and the eight Figure 6 fits. Run from mec/ after
## 03_fit_mec.R has produced all eight fits.
## Output: <results_dir>/section5_numbers.txt (also printed)
## =============================================================================
source("config.R")

dl <- readRDS(file.path(data_dir, "MEC_genus_mediation.rds"))
pri_dir <- file.path(results_dir, "primary")
sen_dir <- file.path(results_dir, "sensitivity")
pri_file <- function(a, e)
  sprintf("%s/pfat_ethall_metfRetained_adj%s_eta1over%d_gibbs_4chains_iter5000.rds", pri_dir, a, e)
fits <- list(
  "primary (eta = 1/20, diabetes)" = pri_file("diabet", 20),
  "eta = 1/3, diabetes"            = pri_file("diabet", 3),
  "eta = 1/5, diabetes"            = pri_file("diabet", 5),
  "eta = 1/20, metformin"          = pri_file("metformin", 20),
  "eta = 1/20, neither"            = pri_file("none", 20),
  "no adiposity"                   = sprintf("%s/pfat_ethall_metfRetained_adjdiabet_covnoadip_eta1over20_gibbs_4chains_iter5000.rds", sen_dir),
  "metformin excluded"             = sprintf("%s/pfat_ethall_noMetf_adjdiabet_eta1over20_gibbs_4chains_iter5000.rds", sen_dir),
  "diabetes excluded"              = sprintf("%s/pfat_ethall_noT2D_adjnone_eta1over20_gibbs_4chains_iter5000.rds", sen_dir))
miss <- unlist(fits)[!file.exists(unlist(fits))]
if (length(miss)) stop("missing fit(s):\n  ", paste(miss, collapse = "\n  "))
obj <- lapply(fits, readRDS)

out <- character(0)
say <- function(...) { line <- sprintf(...); out <<- c(out, line); cat(line, "\n") }
ci  <- function(v) sprintf("%.4f (95%% CrI %.4f to %.4f)", mean(v), quantile(v, 0.025), quantile(v, 0.975))
draws <- function(o, what) unlist(lapply(o$fits, function(f) switch(what,
  gamma = as.numeric(f$gamma), alpha = as.numeric(f$alpha[, 1]), beta = as.numeric(f$beta[, 1]),
  ie = as.numeric(f$alpha[, 1] * f$beta[, 1]))))

## sample and data
known <- !is.na(dl$metformin)
say("Analytic sample: n = %d (participants with known metformin use)", sum(known))
say("Metformin users: %d, of whom %d have diabetes", sum(dl$metformin[known] == 1),
    sum(dl$metformin[known] == 1 & dl$covariates$OQ3_diabet[known] == 1))
say("Metformin users excluded: n = %d; participants with diabetes excluded: n = %d",
    sum(dl$metformin[known] == 0), sum(dl$covariates$OQ3_diabet[known] == 0))
say("Genera retained: %d; zeros in the filtered abundance matrix: %.1f%%",
    ncol(dl$exposures), 100 * mean(dl$exposures == 0.5))

## primary model
p <- obj[[1]]
de <- draws(p, "gamma"); ie <- draws(p, "ie")
say("Primary: gamma %s", ci(de))
say("Primary: alpha %s", ci(draws(p, "alpha")))
say("Primary: beta %s", ci(draws(p, "beta")))
say("Primary: alpha*beta %s", ci(ie))
say("Primary: mediated proportion IE / (DE + IE) %s", ci(ie / (de + ie)))

## stability across the eight fits
gm <- sapply(obj, function(o) mean(draws(o, "gamma")))
am <- sapply(obj, function(o) mean(draws(o, "alpha")))
say("gamma posterior means across the 8 fits: %.3f to %.3f; all CrIs exclude zero: %s",
    min(gm), max(gm), all(sapply(obj, function(o) quantile(draws(o, "gamma"), 0.025) > 0)))
say("alpha posterior means across the 8 fits: %.3f to %.3f; all CrIs exclude zero: %s",
    min(am), max(am), all(sapply(obj, function(o) quantile(draws(o, "alpha"), 0.025) > 0)))
adj <- obj[names(obj) != "no adiposity"]
iem <- sapply(adj, function(o) mean(draws(o, "ie")))
ppos <- sapply(adj, function(o) mean(draws(o, "ie") > 0))
say("Adiposity-adjusted fits: IE posterior means %.4f to %.4f; minimum P(IE > 0) = %.3f", min(iem), max(iem), min(ppos))
for (nm in names(obj)) {
  v <- draws(obj[[nm]], "ie")
  say("  %-32s IE %s; P(IE > 0) = %.3f; CrI includes zero: %s", nm, ci(v), mean(v > 0),
      quantile(v, 0.025) <= 0)
}
na <- obj[["no adiposity"]]
say("No adiposity covariate: beta %s", ci(draws(na, "beta")))
say("No adiposity covariate: alpha*beta %s", ci(draws(na, "ie")))

## convergence
rh <- sapply(obj, function(o) max(o$diagnostics$rhat))
say("Largest split-Rhat over the 8 fits and 4 reported quantities: %.3f", max(rh))
say("Smallest cross-chain correlation of inclusion probabilities: %.3f", min(sapply(obj, `[[`, "pip_cross_cor")))

## selected taxa
t0 <- p$taxa
sel <- t0$taxon[t0$P_incl > 0.5]
n_fits <- sapply(sel, function(tx) sum(sapply(obj, function(o) o$taxa$P_incl[o$taxa$taxon == tx] > 0.5)))
say("Genera selected in the primary model (PIP > 0.5): %d (%d numerator, %d denominator); %d selected in all 8 fits",
    length(sel), sum(t0$P_incl > 0.5 & t0$P_plus > t0$P_minus), sum(t0$P_incl > 0.5 & t0$P_plus < t0$P_minus),
    sum(n_fits == 8))

writeLines(out, file.path(results_dir, "section5_numbers.txt"))
