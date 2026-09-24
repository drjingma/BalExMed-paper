# =============================================================================
# sim_lib.R -- simulation library for Section 4 of the paper.
#
# Provides: data generation (Scenarios I/II/III with PVE equalization),
# fits for the proposed method + three competitors (oracle, principal balance
# via constrained PCs [Martin-Fernandez et al. 2018], regDOC), and evaluation
# metrics (IE estimation, taxon selection, balance recovery). A driver sources
# this and loops over the grid.
#
# The FITTED proposed model uses compute_balance() (equal weights); the DGP
# uses weighted_balance() -- equal weights reproduce compute_balance exactly
# (Scenario I, correct), non-equal weights do not (Scenario II, misspecified).
# =============================================================================

suppressPackageStartupMessages({
  library(MASS); library(PRROC); library(coda)
})
source("../R/balance_mediation_sampler.R")  # run_mcmc_Gibbs/MH, compute_balance
source("../R/regDOC.R")                     # regDOC_sampler, compute_pip

`%||%` <- function(a, b) if (is.null(a)) b else a

## ---- weighted balance for the DGP (standalone copy of the arXiv helper) -----
weighted_balance <- function(X, z, weights = NULL) {
  if (is.null(weights)) weights <- rep(1, ncol(X))
  z[is.na(z)] <- 0
  bplus <- z == 1; bminus <- z == -1
  mplus <- sum(weights[bplus]); mminus <- sum(weights[bminus])
  if (mplus == 0 || mminus == 0)
    stop("weighted_balance needs both numerator and denominator taxa")
  apply(X, 1, function(x) {
    lx <- weights * log(x)
    (sum(lx[bplus]) / mplus - sum(lx[bminus]) / mminus) *
      sqrt(1 / (1 / mplus + 1 / mminus))
  })
}

## ---- 1. compositional exposure (log-normal -> round -> 0.5 -> TSS) ----------
gen_composition <- function(n, d, rho = 0.2, var_log = 9) {
  S <- rho^abs(outer(seq_len(d), seq_len(d), "-"))
  sdl <- sqrt(var_log)
  S <- diag(sdl, d) %*% S %*% diag(sdl, d)            # Sigma_jj' = var_log * rho^|j-j'|
  logX <- MASS::mvrnorm(n = n, mu = rep(0, d), Sigma = S)
  Xc <- round(exp(logX)); Xc[Xc == 0] <- 0.5
  Xr <- sweep(Xc, 1, rowSums(Xc), FUN = "/")
  colnames(Xr) <- paste0("v", seq_len(d))
  Xr
}

## ---- 2. mediation data (Scenarios I / II / III) ----------------------------
# scenario: "I"  equal-weight balance (correctly specified)
#           "II" weighted log-contrast, PVE-matched to I (balance-shape misspec.)
#           "III" equal-weight balance + hidden confounder u (delta_m, delta_y)
# Returns m, y, cc (covariate), z_true, B (the balance actually used).
gen_mediation <- function(Xr, scenario = "I", alpha = 1/3, beta = 2, gamma = 1/3,
                          delta = c(0, 0), a_plus = 3, a_minus = 3) {
  n <- nrow(Xr); d <- ncol(Xr)
  if (scenario == "II") { a_plus <- 3; a_minus <- 3 }    # weighted contrast is fixed on 6 taxa
  z_true <- c(rep(1, a_plus), rep(-1, a_minus), rep(0, d - a_plus - a_minus))
  B_equal <- compute_balance(Xr, z_true)                 # = the sampler's balance
  if (scenario == "II") {
    w <- c(1, 0.4, 1.2, -1.5, -0.8, -0.3, rep(0, d - 6)) # signs match z_true
    B <- weighted_balance(Xr, sign(w), abs(w))
    B <- B * sd(B_equal) / sd(B)                         # PVE equalization -> same Var as I
  } else {
    B <- B_equal
  }
  cc <- rnorm(n)                                         # one covariate, psi = 1
  u  <- if (scenario == "III") rnorm(n) else rep(0, n)   # hidden confounder
  e  <- rnorm(n); eps <- rnorm(n)
  m <- alpha * B + cc + delta[1] * u + e
  y <- gamma * B + beta * m + cc + delta[2] * u + eps
  list(m = m, y = y, cc = cc, z_true = z_true, B = B,
       IE_true = alpha * beta, DE_true = gamma)
}

## ---- 3a. proposed method (BalMed-Gibbs / MH) -------------------------------
fit_proposed <- function(Xr, m, y, cc, eta = 1/10, n_iter = 2e4, burn_in = 1e4,
                         sampler = c("gibbs", "mh"), z_init = NULL) {
  sampler <- match.arg(sampler)
  C <- cbind(1, cc)                                      # intercept + covariate
  hyper <- list(h_alpha = 1e-6, h_beta = 1e-6, h_gamma = 1e-6, hm = 1e-6, hy = 1e-6,
                nu = length(y), lambda_e = stats::var(m), lambda_eps = stats::var(y),
                tau_a = 1/3, tau_d = 1/3, tau_s = 1/3, eta1 = eta, eta2 = eta)
  fun <- if (sampler == "mh") run_mcmc_MH else run_mcmc_Gibbs
  tt <- system.time(
    fit <- fun(X = Xr, C = C, M = matrix(m, ncol = 1), Y = y,
               hyper = hyper, n_iter = n_iter, burn_in = burn_in, z_init = z_init))
  ie <- fit$alpha[, 1] * fit$beta[, 1]
  p_plus  <- colMeans(fit$z ==  1)
  p_minus <- colMeans(fit$z == -1)
  # efficiency: effective sample size of the IE (alpha*beta) and gamma traces,
  # and the sampler CPU time (user+sys, robust to core contention under mclapply).
  ess_ie    <- as.numeric(coda::effectiveSize(coda::as.mcmc(ie)))
  ess_gamma <- as.numeric(coda::effectiveSize(coda::as.mcmc(fit$gamma)))
  fit_cpu   <- unname(tt["user.self"] + tt["sys.self"])
  mix_rate  <- if (!is.null(fit$accept_rate)) fit$accept_rate else fit$move_rate  # MH accept / Gibbs move
  list(IE = mean(ie), IE_lo = quantile(ie, .025), IE_hi = quantile(ie, .975),
       DE = mean(fit$gamma), DE_lo = quantile(fit$gamma, .025), DE_hi = quantile(fit$gamma, .975),
       pip = p_plus + p_minus, p_plus = p_plus, p_minus = p_minus,
       ess_ie = ess_ie, ess_gamma = ess_gamma, fit_cpu = fit_cpu, mix_rate = mix_rate,
       fit = fit)
}

## ---- 3b. oracle: 2SLS product-of-coefficients at the TRUE z -----------------
fit_oracle <- function(Xr, z_true, m, y, cc) {
  B <- compute_balance(Xr, z_true)
  am <- lm(m ~ B + cc); ay <- lm(y ~ B + m + cc)
  a <- coef(am)["B"]; b <- coef(ay)["m"]
  sea <- summary(am)$coefficients["B", "Std. Error"]
  seb <- summary(ay)$coefficients["m", "Std. Error"]
  ie <- a * b
  se <- sqrt((b * sea)^2 + (a * seb)^2)                  # delta method
  de  <- coef(ay)["B"]; sede <- summary(ay)$coefficients["B", "Std. Error"]
  c(IE = unname(ie), IE_lo = unname(ie - 1.96 * se), IE_hi = unname(ie + 1.96 * se),
    DE = unname(de), DE_lo = unname(de - 1.96 * sede), DE_hi = unname(de + 1.96 * sede))
}

## ---- 3c. principal balance via CONSTRAINED PCs (threshold-free) -------------
# Martin-Fernandez, Pawlowsky-Glahn, Egozcue & Tolosana-Delgado (2018),
# "Advances in principal balances for compositional data", Math Geosci 50:273-298,
# Sec. 4.2 (constrained PC approach of Chipman & Gu 2005). The first principal
# balance is the balance whose clr-coefficient vector alpha (Eq. 4: values in
# {-c1, 0, c2}) is closest to the leading clr-PC gamma_1, i.e. maximizes
##  |cos(gamma_1, alpha)|. The number of parts is chosen by the search (candidates
# use the top-k parts by |loading|, grouped by loading sign), removing the
# subjective magnitude threshold of the old ad-hoc version.
pb_constrained_pc <- function(Xr) {
  clr <- t(apply(Xr, 1, function(a) log(a) - mean(log(a))))
  g <- svd(scale(clr, center = TRUE, scale = FALSE))$v[, 1]   # leading clr-PC loading (unit norm)
  D <- length(g)
  pos <- order(g, decreasing = TRUE);  pos <- pos[g[pos] > 0]
  neg <- order(g, decreasing = FALSE); neg <- neg[g[neg] < 0]
  if (!length(pos) || !length(neg)) return(list(z = rep(0L, D), loading = g))
  seed <- c(pos[1], neg[1])                                    # largest +, largest - (2-part start)
  add  <- setdiff(order(abs(g), decreasing = TRUE), seed)      # remaining parts by |loading|
  best_ip <- -Inf; best_z <- rep(0L, D); incl <- seed
  for (k in 2:D) {                                             # candidate balances: top-k parts
    if (k > 2) incl <- c(incl, add[k - 2])
    z <- integer(D); z[incl[g[incl] > 0]] <- 1L; z[incl[g[incl] < 0]] <- -1L
    r <- sum(z == 1); s <- sum(z == -1)
    if (r == 0 || s == 0) next
    nf <- sqrt(r * s / (r + s)); a <- numeric(D)               # balance clr-coeffs (Eq. 4), unit norm
    a[z == 1] <- nf / r; a[z == -1] <- -nf / s
    ip <- abs(sum(g * a))                                      # |cos angle| to the PC (g, a both unit)
    if (ip > best_ip) { best_ip <- ip; best_z <- z }
  }
  list(z = best_z, loading = g)
}

fit_pb <- function(Xr, m, y, cc) {
  pb <- pb_constrained_pc(Xr); z_pb <- pb$z
  ie <- NA; de <- NA
  if (any(z_pb == 1) && any(z_pb == -1)) {
    B <- compute_balance(Xr, z_pb)
    a <- coef(lm(m ~ B + cc))["B"]; fy <- lm(y ~ B + m + cc)
    ie <- unname(a * coef(fy)["m"]); de <- unname(coef(fy)["B"])
  }
  list(IE = ie, DE = de, score = abs(pb$loading), z = z_pb)    # |PC loading| = selection score
}

## ---- 3d. regDOC (selection only) -------------------------------------------
fit_regdoc <- function(Xr, m, y, n_iter = 1e4, burn_in = 5e3, epsilon = 0.01) {
  clr <- t(apply(Xr, 1, function(a) log(a) - mean(log(a))))
  s <- regDOC_sampler(X = clr, M = as.numeric(m), Y = as.numeric(y),
                      n_iter = n_iter, burn_in = burn_in)
  list(pip = compute_pip(s$alpha, epsilon = epsilon))
}

## ---- 4. metrics -------------------------------------------------------------
auprc <- function(score, label) {                        # label: 1 = truly active
  if (length(unique(label)) < 2) return(NA_real_)
  PRROC::pr.curve(scores.class0 = score[label == 1],
                  scores.class1 = score[label == 0])$auc.integral
}
f1_at <- function(pip, label, thr = 0.5) {
  sel <- as.integer(pip > thr)
  tp <- sum(sel == 1 & label == 1); fp <- sum(sel == 1 & label == 0)
  fn <- sum(sel == 0 & label == 1)
  if (tp == 0) return(0)
  prec <- tp / (tp + fp); rec <- tp / (tp + fn)
  2 * prec * rec / (prec + rec)
}
# sign-aligned balance recovery for the proposed fit
balance_recovery <- function(Xr, p_plus, p_minus, B_true) {
  z_hat <- ifelse(p_plus > 0.5 & p_plus >= p_minus, 1,
                  ifelse(p_minus > 0.5, -1, 0))
  if (!any(z_hat == 1) || !any(z_hat == -1)) return(NA_real_)
  abs(cor(compute_balance(Xr, z_hat), B_true))
}

## ---- 5. one replicate over all methods for a single config -----------------
# cfg: list(n, d, scenario, beta, eta, rho, delta, sampler, n_iter, burn_in,
#           competitors = TRUE/FALSE)
run_replicate <- function(cfg, seed) {
  set.seed(seed)
  Xr  <- gen_composition(cfg$n, cfg$d, rho = cfg$rho %||% 0.2)
  dat <- gen_mediation(Xr, scenario = cfg$scenario, alpha = cfg$alpha %||% (1/3),
                       beta = cfg$beta, gamma = cfg$gamma %||% (1/3),
                       delta = cfg$delta %||% c(0, 0),
                       a_plus = cfg$a_plus %||% 3, a_minus = cfg$a_minus %||% 3)
  lab <- as.integer(dat$z_true != 0)

  pr <- fit_proposed(Xr, dat$m, dat$y, dat$cc, eta = cfg$eta %||% (1/10),
                     n_iter = cfg$n_iter %||% 2e4, burn_in = cfg$burn_in %||% 1e4,
                     sampler = cfg$sampler %||% "gibbs")
  out <- data.frame(
    seed = seed, method = "BalMed", IE = pr$IE, IE_lo = pr$IE_lo, IE_hi = pr$IE_hi,
    DE = pr$DE, DE_lo = pr$DE_lo, DE_hi = pr$DE_hi, IE_true = dat$IE_true, DE_true = dat$DE_true,
    cover = pr$IE_lo <= dat$IE_true & dat$IE_true <= pr$IE_hi,
    width = pr$IE_hi - pr$IE_lo, reject0 = !(pr$IE_lo <= 0 & 0 <= pr$IE_hi),
    auprc = auprc(pr$pip, lab), f1 = f1_at(pr$pip, lab),
    brecov = balance_recovery(Xr, pr$p_plus, pr$p_minus, dat$B),
    row.names = NULL)

  if (isTRUE(cfg$competitors)) {
    orc <- fit_oracle(Xr, dat$z_true, dat$m, dat$y, dat$cc)
    out <- rbind(out, data.frame(
      seed = seed, method = "Oracle", IE = orc["IE"], IE_lo = orc["IE_lo"], IE_hi = orc["IE_hi"],
      DE = orc["DE"], DE_lo = orc["DE_lo"], DE_hi = orc["DE_hi"],
      IE_true = dat$IE_true, DE_true = dat$DE_true,
      cover = orc["IE_lo"] <= dat$IE_true & dat$IE_true <= orc["IE_hi"],
      width = orc["IE_hi"] - orc["IE_lo"], reject0 = !(orc["IE_lo"] <= 0 & 0 <= orc["IE_hi"]),
      auprc = NA, f1 = NA, brecov = NA, row.names = NULL))
    pb <- fit_pb(Xr, dat$m, dat$y, dat$cc)
    out <- rbind(out, data.frame(
      seed = seed, method = "PrinBal", IE = pb$IE, IE_lo = NA, IE_hi = NA,
      DE = pb$DE, DE_lo = NA, DE_hi = NA,
      IE_true = dat$IE_true, DE_true = dat$DE_true, cover = NA, width = NA, reject0 = NA,
      auprc = auprc(pb$score, lab), f1 = NA, brecov = NA, row.names = NULL))
    rd <- fit_regdoc(Xr, dat$m, dat$y)
    out <- rbind(out, data.frame(
      seed = seed, method = "regDOC", IE = NA, IE_lo = NA, IE_hi = NA,
      DE = NA, DE_lo = NA, DE_hi = NA,
      IE_true = dat$IE_true, DE_true = dat$DE_true, cover = NA, width = NA, reject0 = NA,
      auprc = auprc(rd$pip, lab), f1 = f1_at(rd$pip, lab), brecov = NA, row.names = NULL))
  }
  # direct-effect interval calibration (methods with a DE interval; NA otherwise)
  out$DE_cover <- out$DE_lo <= out$DE_true & out$DE_true <= out$DE_hi
  out$DE_width <- out$DE_hi - out$DE_lo
  # efficiency columns (BalMed only; NA for the competitors)
  out$ess_ie <- NA_real_; out$ess_gamma <- NA_real_; out$fit_cpu <- NA_real_
  out$ess_ie_per_sec <- NA_real_; out$mix_rate <- NA_real_
  bm <- out$method == "BalMed"
  out$ess_ie[bm] <- pr$ess_ie; out$ess_gamma[bm] <- pr$ess_gamma
  out$fit_cpu[bm] <- pr$fit_cpu
  out$ess_ie_per_sec[bm] <- pr$ess_ie / pr$fit_cpu
  out$mix_rate[bm] <- pr$mix_rate
  # per-taxon selection scores (list-columns; averaged over replicates at cell level).
  # score_vec = PIP/score used for AUPRC; signed_vec = P(+)-P(-) orientation (BalMed).
  out$score_vec  <- vector("list", nrow(out))
  out$signed_vec <- vector("list", nrow(out))
  out$score_vec[[which(out$method == "BalMed")]]  <- pr$pip
  out$signed_vec[[which(out$method == "BalMed")]] <- pr$p_plus - pr$p_minus
  if (isTRUE(cfg$competitors)) {
    out$score_vec[[which(out$method == "PrinBal")]] <- pb$score
    out$score_vec[[which(out$method == "regDOC")]]  <- rd$pip
  }
  out
}

## replicate-averaged per-taxon scores per method (for a component-wise PIP figure);
## returns list: method -> list(score = mean PIP/score, signed = mean P(+)-P(-) for BalMed).
avg_pip <- function(res) {
  stack_mean <- function(vs) {
    vs <- vs[!vapply(vs, is.null, logical(1))]
    if (!length(vs)) return(NULL)
    L <- max(lengths(vs))
    M <- do.call(rbind, lapply(vs, function(v) { length(v) <- L; v }))
    colMeans(M, na.rm = TRUE)
  }
  o <- list()
  for (mth in unique(res$method)) {
    s <- stack_mean(res$score_vec[res$method == mth])
    if (!is.null(s)) o[[mth]] <- list(score = s)
  }
  sg <- stack_mean(res$signed_vec[res$method == "BalMed"])
  if (!is.null(sg)) o[["BalMed"]]$signed <- sg
  o
}
