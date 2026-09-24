## =============================================================================
## 03_fit_mec.R -- fit the balance mediation model to the MEC-APS data (Section 5).
##
##   balance B(z, X) --alpha--> LBP (mediator) --beta--> log percent liver fat
##          \----------------------- gamma -----------------------/
##
## Usage (from mec/):
##   Rscript 03_fit_mec.R [eta_denom] [eth] [adjust] [n_chains] [n_iter] [burn_in] \
##                        [sampler] [mode] [chain] [covars] [sample]
## All arguments are optional; pass "" to skip one. Defaults give the primary model.
##   [1] eta_denom  prior eta = 1/eta_denom for each side of the balance: 20 (primary), 5, 3
##   [2] eth        "all" (primary); "J" or "L" restricts to one ethnic group (not used in the paper)
##   [3] adjust     "diabet" (primary): adjust for diabetes; "metformin": adjust for
##                  metformin use; "none": adjust for neither
##   [4] n_chains   number of chains (4)
##   [5] n_iter     iterations per chain, including burn-in (5000)
##   [6] burn_in    burn-in iterations per chain (1000)
##   [7] sampler    "gibbs" (primary) or "mh"
##   [8] mode       "all": run all chains in parallel in this process;
##                  "chain": run a single chain and save it (for a SLURM array);
##                  "aggregate": combine saved chains into diagnostics and tables
##   [9] chain      chain index for mode = "chain"
##  [10] covars     "new" (primary): total body fat percentage is the adiposity covariate;
##                  "noadip": no adiposity covariate; "bmi": body mass index instead (not in the paper)
##  [11] sample     "retained" (primary): metformin users kept; "nometf": metformin users
##                  excluded; "not2d": participants with diabetes excluded
## The eight fits shown in Figure 6 are listed in submit_mec.sh.
##
## Covariates (all fits): age, sex, race and ethnicity, place of birth, education,
## AHEI-2010 score, smoking pack-years, sedentary hours, season of stool collection,
## total body fat percentage (unless covars = "noadip"), plus diabetes or metformin
## according to `adjust`. Continuous covariates are standardized. LBP is scaled to
## unit root mean square. Participants with unknown metformin use are excluded from
## every fit so that the adjustment settings share one sample.
##
## Chain 1 starts from the leading principal component of the clr-transformed
## genera (loadings above 0.1 in absolute value); chains 2 to 4 start at random.
## Chain c is seeded with 1000 + c.
##
## Input:  <data_dir>/MEC_genus_mediation.rds (02_build_mediation_inputs.R)
## Output: <results_dir>/primary/ when covars = "new", sample = "retained" and
##         eth = "all"; <results_dir>/sensitivity/ otherwise. Each fit writes
##         <tag>.rds (draws, diagnostics, taxa table, settings), taxa_<tag>.csv and
##         diag_<tag>.pdf; mode = "chain" writes parts/<tag>/chain_<c>.rds.
## =============================================================================

args <- commandArgs(trailingOnly = TRUE)
print(args)

get_arg <- function(i, default) {
  v <- if (length(args) >= i) args[i] else NA
  if (is.na(v) || v == "") default else v
}
eta_denom  <- as.numeric(get_arg(1, 20))
eth_subset <- as.character(get_arg(2, "all"))
adjust     <- tolower(as.character(get_arg(3, "diabet")))
n_chains   <- as.numeric(get_arg(4, 4))
n_iter     <- as.numeric(get_arg(5, 5000))
burn_in    <- as.numeric(get_arg(6, 1000))
sampler    <- tolower(as.character(get_arg(7, "gibbs")))
mode       <- tolower(as.character(get_arg(8, "all")))
chain_arg  <- get_arg(9, NA)
covars     <- tolower(as.character(get_arg(10, "new")))
sample_arg <- tolower(as.character(get_arg(11, "retained")))
eta <- c(1, 1) / eta_denom
if (!adjust %in% c("diabet", "metformin", "none"))
  stop("adjust (arg 3) must be diabet|metformin|none (got '", adjust, "')")
if (!covars %in% c("new", "bmi", "noadip"))
  stop("covars (arg 10) must be new|bmi|noadip (got '", covars, "')")
if (!sample_arg %in% c("retained", "nometf", "not2d"))
  stop("sample (arg 11) must be retained|nometf|not2d (got '", sample_arg, "')")
if (sample_arg == "nometf" && adjust == "metformin")
  stop("sample=nometf leaves no metformin users -- use adjust=diabet or none")
if (sample_arg == "not2d" && adjust == "diabet")
  stop("sample=not2d leaves no participants with diabetes -- use adjust=metformin or none")
if (!sampler %in% c("gibbs", "mh"))
  stop("sampler (arg 7) must be 'gibbs' or 'mh' (got '", sampler, "')")
if (!mode %in% c("all", "chain", "aggregate"))
  stop("mode (arg 8) must be all|chain|aggregate (got '", mode, "')")

## ---- sampler and packages -----------------------------------------------------
source("config.R")
source("../R/balance_mediation_sampler.R")
sampler_fun <- if (sampler == "mh") run_mcmc_MH else run_mcmc_Gibbs
cat("sampler:", sampler, "\n")

suppressPackageStartupMessages({
  library(dplyr); library(parallel); library(coda)
})

## One BLAS thread per chain avoids oversubscribing cores when chains run in parallel.
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
try(suppressWarnings(RhpcBLASctl::blas_set_num_threads(1)), silent = TRUE)

## ---- data -----------------------------------------------------------------------
DATA_TAG  <- "genus2024"
data_path <- file.path(data_dir, "MEC_genus_mediation.rds")
if (!file.exists(data_path)) stop(data_path, " not found; run 02_build_mediation_inputs.R first")
data.list <- readRDS(data_path)
cat("data:", data_path, "|", nrow(data.list$response), "participants,",
    ncol(data.list$exposures), "features\n")

outcome_var <- "OMRI_pct_liver_fat_corr"

## ---- sample selection -----------------------------------------------------------
apply_select <- function(dl, keep) {
  dl$response   <- dl$response[keep, , drop = FALSE]
  dl$covariates <- dl$covariates[keep, , drop = FALSE]
  dl$exposures  <- dl$exposures[keep, , drop = FALSE]
  dl$mediator   <- dl$mediator[keep]
  if (!is.null(dl$metformin)) dl$metformin <- dl$metformin[keep]
  dl
}
if (eth_subset %in% c("J", "L")) {
  data.list <- apply_select(data.list, as.character(data.list$covariates$Q1_eth) == eth_subset)
  data.list$covariates <- droplevels(data.list$covariates)
}
data.list <- apply_select(data.list, !is.na(data.list$metformin))
if (sample_arg == "nometf")
  data.list <- apply_select(data.list, data.list$metformin == 0)
if (sample_arg == "not2d")
  data.list <- apply_select(data.list, data.list$covariates$OQ3_diabet == 0)

n <- nrow(data.list$response)
cat(sprintf("n = %d (%d metformin users, %d non-users) | eth=%s | adjust=%s | covars=%s | sample=%s | eta=1/%g | sampler=%s\n",
            n, sum(data.list$metformin == 1), sum(data.list$metformin == 0),
            eth_subset, adjust, covars, sample_arg, eta_denom, sampler))

## ---- covariates -----------------------------------------------------------------
vars <- c("OST_sample_age", "OQ3_DP_AHEI2010_TOTSCORE",
          "ODXA_pfat_tot_corr", "OQ3_packyrs", "OQ3_ac_hrsit")
if (covars == "noadip") vars <- setdiff(vars, "ODXA_pfat_tot_corr")
if (covars == "bmi")    vars <- c(setdiff(vars, "ODXA_pfat_tot_corr"), "OCL_anthro_BMI")
vars_all <- c(vars, "Q1_CORR_SEX", "Q1_eth", "Q1_POB", "edu_cat", "season")
if (adjust == "diabet")          vars_all <- c(vars_all, "OQ3_diabet")
if (adjust == "metformin") {
  data.list$covariates$metformin <- data.list$metformin
  vars_all <- c(vars_all, "metformin")
}
if (eth_subset %in% c("J", "L")) vars_all <- setdiff(vars_all, "Q1_eth")

data.list$covariates <- data.list$covariates %>%
  dplyr::mutate(across(all_of(vars), ~ as.numeric(scale(.)))) %>%
  dplyr::select(all_of(vars_all))
W <- model.matrix(~ ., data = data.list$covariates)

X     <- as.matrix(data.list$exposures)
X.rel <- sweep(X, 1, rowSums(X), FUN = "/")     # zeros were replaced by 0.5 in step 2
p     <- ncol(X)
m     <- as.matrix(scale(data.list$mediator, center = FALSE, scale = TRUE))

y <- unlist(data.list$response[, outcome_var]); names(y) <- NULL

hyper <- list(
  h_alpha = 1e-6, h_beta = 1e-6, h_gamma = 1e-6, hm = 1e-6, hy = 1e-6,
  nu = n, lambda_e = as.numeric(var(m)), lambda_eps = as.numeric(var(y)),
  tau_a = 1/3, tau_d = 1/3, tau_s = 1/3,
  eta1 = eta[1], eta2 = eta[2]
)

## ---- starting values: chain 1 from the leading clr principal component ----------
X.clr <- t(apply(X.rel, 1, function(a) log(a) - mean(log(a))))
v     <- svd(X.clr)$v[, 1]
z_pca <- sign(v * (abs(v) > 0.1))
if (!any(z_pca == 1) || !any(z_pca == -1)) {
  warning("PCA start lacks both signs at threshold 0.1; chain 1 uses a random start")
  z_pca <- NULL
} else {
  cat(sprintf("PCA start for chain 1: %d in numerator, %d in denominator\n",
              sum(z_pca == 1), sum(z_pca == -1)))
}

## ---- output paths ---------------------------------------------------------------
smp_lab <- c(retained = "metfRetained", nometf = "noMetf", not2d = "noT2D")[sample_arg]
cov_lab <- if (covars == "new") "" else paste0("_cov", covars)
tag <- sprintf("pfat_eth%s_%s_adj%s%s_eta1over%g_%s_%dchains_iter%d",
               eth_subset, smp_lab, adjust, cov_lab, eta_denom, sampler, n_chains, n_iter)
out_dir <- file.path(results_dir,
                     if (covars == "new" && sample_arg == "retained" && eth_subset == "all")
                       "primary" else "sensitivity")
parts_dir <- file.path(out_dir, "parts", tag)
dir.create(parts_dir, showWarnings = FALSE, recursive = TRUE)

run_one <- function(chain) {
  set.seed(1000 + chain)
  z0 <- if (chain == 1) z_pca else NULL
  sampler_fun(X = X.rel, C = W, M = m, Y = y, hyper = hyper,
              n_iter = n_iter, burn_in = burn_in, z_init = z0)
}

## ---- mode = chain: run one chain, save it, exit -----------------------------------
if (mode == "chain") {
  k <- if (!is.na(chain_arg)) as.integer(chain_arg) else
       as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "0")) + 1L
  if (k < 1L || k > n_chains) stop("chain index ", k, " out of [1, ", n_chains, "]")
  fp <- file.path(parts_dir, sprintf("chain_%d.rds", k))
  if (file.exists(fp)) { cat("chain", k, "already done:", fp, "\n"); quit(save = "no") }
  cat(sprintf("chain %d/%d | %s | n_iter=%d burn_in=%d\n", k, n_chains, tag, n_iter, burn_in))
  t0 <- Sys.time(); fit <- run_one(k)
  cat("elapsed:", round(difftime(Sys.time(), t0, units = "secs"), 1), "s\n")
  tmp <- paste0(fp, ".tmp"); saveRDS(fit, tmp); file.rename(tmp, fp)
  cat("saved:", fp, "\n"); quit(save = "no")
}

## ---- assemble chains: aggregate reads saved chains, all runs them now ----------
if (mode == "aggregate") {
  cf <- file.path(parts_dir, sprintf("chain_%d.rds", seq_len(n_chains)))
  if (!all(file.exists(cf)))
    stop("missing chain file(s): ", paste(which(!file.exists(cf)), collapse = ","),
         " -- run those chain tasks first (mode=chain)")
  fits <- lapply(cf, readRDS)
  cat("aggregating", n_chains, "chains from", parts_dir, "\n")
} else {
  n_cores <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "")))
  if (is.na(n_cores) || n_cores < 1L) n_cores <- detectCores()
  n_cores <- min(n_chains, n_cores)
  cat("running", n_chains, "chains on", n_cores, "cores | n_iter =", n_iter, "burn_in =", burn_in, "\n")
  t0 <- Sys.time()
  fits <- mclapply(seq_len(n_chains), run_one, mc.cores = n_cores)
  cat("elapsed:", round(difftime(Sys.time(), t0, units = "secs"), 1), "s\n")
}
ok <- vapply(fits, function(f) is.list(f) && !is.null(f$gamma), logical(1))
if (!all(ok)) stop("chain(s) failed/missing: ", paste(which(!ok), collapse = ","))

## ---- convergence diagnostics and effects ------------------------------------------
DE    <- lapply(fits, function(f) as.numeric(f$gamma))
IE    <- lapply(fits, function(f) as.numeric(f$alpha[, 1] * f$beta[, 1]))
Araw  <- lapply(fits, function(f) as.numeric(f$alpha[, 1]))
Braw  <- lapply(fits, function(f) as.numeric(f$beta[, 1]))

## split-Rhat: each chain split in half, giving 2 * n_chains sub-chains
split_rhat <- function(chain_list) {
  n <- min(lengths(chain_list)); half <- floor(n / 2)
  subs <- unlist(lapply(chain_list, function(x) {
    x <- x[seq_len(2 * half)]
    list(x[1:half], x[(half + 1):(2 * half)])
  }), recursive = FALSE)
  M <- length(subs); nn <- half
  means <- sapply(subs, mean); vars <- sapply(subs, var)
  B <- nn / (M - 1) * sum((means - mean(means))^2)
  Wv <- mean(vars)
  Vhat <- (nn - 1) / nn * Wv + B / nn
  sqrt(Vhat / Wv)
}
ess_of <- function(chain_list) {
  n <- min(lengths(chain_list))
  ml <- as.mcmc.list(lapply(chain_list, function(x) as.mcmc(x[seq_len(n)])))
  as.numeric(effectiveSize(ml))
}
qsum <- function(chain_list) {
  v <- unlist(chain_list); c(mean = mean(v), q2.5 = quantile(v, .025), q97.5 = quantile(v, .975))
}

diag_tbl <- data.frame(
  param = c("direct (gamma)", "indirect (alpha*beta)", "alpha (B->LBP)", "beta (LBP->y)"),
  rhat  = round(c(split_rhat(DE), split_rhat(IE), split_rhat(Araw), split_rhat(Braw)), 3),
  ess   = round(c(ess_of(DE), ess_of(IE), ess_of(Araw), ess_of(Braw))),
  t(sapply(list(DE, IE, Araw, Braw), qsum)), row.names = NULL, check.names = FALSE)
cat("\n================ convergence & effects (pooled over chains) ================\n")
print(diag_tbl, digits = 4)

## ---- posterior inclusion probabilities ------------------------------------------
Pplus  <- sapply(fits, function(f) colMeans(f$z ==  1))   # p x n_chains
Pminus <- sapply(fits, function(f) colMeans(f$z == -1))
Pincl  <- Pplus + Pminus
taxa_tbl <- data.frame(
  taxon   = colnames(X),
  P_plus  = rowMeans(Pplus),
  P_minus = rowMeans(Pminus),
  P_incl  = rowMeans(Pincl),
  P_incl_sd_across_chains = apply(Pincl, 1, sd),
  direction = ifelse(rowMeans(Pplus) >= rowMeans(Pminus), "+", "-")
)
if (!is.null(data.list$taxonomy)) {
  tx <- data.list$taxonomy
  taxa_tbl$family <- tx$family[match(taxa_tbl$taxon, tx$feature)]
  taxa_tbl$genus  <- tx$genus[match(taxa_tbl$taxon, tx$feature)]
  taxa_tbl$n_asv  <- tx$n_asv[match(taxa_tbl$taxon, tx$feature)]
}
taxa_tbl <- taxa_tbl[order(-taxa_tbl$P_incl), ]
pip_cross_cor <- if (n_chains > 1) round(min(cor(Pincl)[upper.tri(cor(Pincl))]), 3) else NA
cat(sprintf("\nTaxon-inclusion cross-chain agreement: min pairwise cor(P_incl) = %s;  max across-chain sd = %.3f\n",
            pip_cross_cor, max(taxa_tbl$P_incl_sd_across_chains)))
cat("Top 15 taxa by mean P_incl:\n"); print(head(taxa_tbl, 15), row.names = FALSE, digits = 3)

## ---- outputs ----------------------------------------------------------------------
saveRDS(list(fits = fits, diagnostics = diag_tbl, taxa = taxa_tbl,
             config = list(data = basename(data_path), data_tag = DATA_TAG,
                           outcome = outcome_var, n = n, eta = eta, n_chains = n_chains,
                           n_iter = n_iter, burn_in = burn_in,
                           eth_subset = eth_subset, sample = sample_arg,
                           covars = covars, adjust = adjust,
                           covariates = vars_all, sampler = sampler),
             pip_cross_cor = pip_cross_cor),
        file = file.path(out_dir, paste0(tag, ".rds")))
write.csv(taxa_tbl, file.path(out_dir, paste0("taxa_", tag, ".csv")), row.names = FALSE)

## trace and running-mean plots of the direct and indirect effects (base PDF device,
## which also works on headless cluster nodes)
draw_diag <- function(fits, file) {
  DE <- lapply(fits, function(f) as.numeric(f$gamma))
  IE <- lapply(fits, function(f) as.numeric(f$alpha[, 1] * f$beta[, 1]))
  cols <- seq_len(length(fits)) + 1
  pdf(file, width = 12, height = 8.5)
  op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  on.exit({ par(op); dev.off() }, add = TRUE)
  matplot(do.call(cbind, DE), type = "l", lty = 1, col = cols,
          xlab = "iteration (post burn-in)", ylab = expression(gamma), main = "Direct effect trace")
  abline(h = 0, lty = 3)
  matplot(do.call(cbind, IE), type = "l", lty = 1, col = cols,
          xlab = "iteration (post burn-in)", ylab = expression(alpha * beta), main = "Indirect effect trace")
  abline(h = 0, lty = 3)
  rmean <- function(x) cumsum(x) / seq_along(x)
  matplot(sapply(DE, rmean), type = "l", lty = 1, col = cols,
          xlab = "iteration", ylab = expression(bar(gamma)), main = "Direct effect running mean")
  matplot(sapply(IE, rmean), type = "l", lty = 1, col = cols,
          xlab = "iteration", ylab = expression(bar(alpha * beta)), main = "Indirect effect running mean")
}
fig <- file.path(out_dir, paste0("diag_", tag, ".pdf"))
figok <- tryCatch({ draw_diag(fits, fig); TRUE },
                  error = function(e) { while (length(dev.list())) dev.off()
                    message("diagnostic figure failed (", conditionMessage(e), ")"); FALSE })

cat("\nsaved to", normalizePath(out_dir), ":\n  ", tag, ".rds / taxa_", tag, ".csv",
    if (figok) paste0(" / diag_", tag, ".pdf") else "", "\n", sep = "")
cat("Rhat should be < 1.01 for DE and IE; if not, increase n_iter / burn_in.\n")
