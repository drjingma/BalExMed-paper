## =============================================================================
## Samplers for the balance mediation model, as used for every result in the
## paper (simulations in Section 4 and the MEC-APS analysis in Section 5).
##
##   m_i = alpha * B(z, x_i) + psi_m' c_i + e_i
##   y_i = gamma * B(z, x_i) + beta * m_i + psi_y' c_i + eps_i
##
## run_mcmc_Gibbs()  collapsed Gibbs sampler (Algorithm 1)
## run_mcmc_MH()     Metropolis-Hastings sampler (Section 3.3)
##
## Arguments of both samplers
##   X        n x d composition (rows sum to one, strictly positive)
##   M        n x 1 mediator matrix
##   C        n x k covariate matrix; include an intercept column yourself
##   Y        length-n outcome
##   hyper    list(h_alpha, h_beta, h_gamma, hm, hy, nu, lambda_e, lambda_eps,
##                 eta1, eta2, tau_a, tau_d)
##   n_iter, burn_in, z_init (NULL = random valid start)
## Both return the post-burn-in draws of z, alpha, beta, gamma and the error
## variances, plus the move rate (Gibbs) or acceptance rate (MH) and the rate at
## which draws were relabeled to enforce gamma > 0.
##
## This file is the single-mediator code path of the research library used for
## the paper, with an unpublished extension to multiple correlated mediators removed.
##
## Version note. The log posterior of z weights log|Lambda_y| by 1/2, as in the
## paper. Earlier versions used (k + 2) / 2, where k is the number of covariate
## columns. The MEC-APS results in the paper were computed after this correction;
## the simulation results in simulations/results were computed before it.
## =============================================================================

library(MASS)   # mvrnorm

gibbs_update_z_site <- function(z_curr, j, logX, M, C, Y, hyper) {
  candidate_states <- c(-1, 0, 1)
  log_probs <- rep(-Inf, 3)
  for (k in seq_along(candidate_states)) {
    z_prop <- z_curr
    z_prop[j] <- candidate_states[k]
    if (!any(z_prop == 1) || !any(z_prop == -1)) next   # keep the balance defined
    log_probs[k] <- log_posterior_z(z_prop, logX, M, C, Y, hyper)
  }
  if (all(!is.finite(log_probs))) return(list(z = z_curr, moved = 0))
  maxlp <- max(log_probs[is.finite(log_probs)])
  probs <- exp(log_probs - maxlp)
  probs <- probs / sum(probs)
  new_state <- sample(candidate_states, size = 1, prob = probs)
  z_new <- z_curr
  z_new[j] <- new_state
  list(z = z_new, moved = as.integer(z_new[j] != z_curr[j]))
}

gibbs_update_z_sweep <- function(z_curr, logX, M, C, Y, hyper) {
  d <- length(z_curr)
  visit_order <- sample.int(d, size = d, replace = FALSE)
  moved_total <- 0
  z_work <- z_curr
  for (j in visit_order) {
    upd <- gibbs_update_z_site(z_work, j, logX, M, C, Y, hyper)
    z_work <- upd$z
    moved_total <- moved_total + upd$moved
  }
  list(z = z_work, n_moved = moved_total)
}

## Balance B(z, X): normalized log-ratio of the geometric means of the parts with
## z = 1 and z = -1 (equation 1 of the paper).
compute_balance <- function(X, z) balance_from_log(log(X), z)

## The same balance from log(X), which the samplers compute once per run.
balance_from_log <- function(logX, z) {
  pos <- which(z == 1)
  neg <- which(z == -1)
  a_p <- length(pos); a_m <- length(neg)
  if (a_p == 0 || a_m == 0)
    stop("Balance undefined unless both numerator (z=1) and denominator (z=-1) nonempty")
  num_mean <- rowMeans(logX[, pos, drop = FALSE])
  den_mean <- rowMeans(logX[, neg, drop = FALSE])
  norm <- sqrt((a_p * a_m) / (a_p + a_m))
  as.numeric(norm * (num_mean - den_mean))
}

## Log posterior of z up to an additive constant, with the coefficients and
## error variances integrated out (equation 5 of the paper). logX = log(X).
log_posterior_z <- function(z, logX, M, C, Y, hyper) {
  n <- nrow(logX)
  q <- ncol(M)
  k <- ncol(C)
  Bz <- balance_from_log(logX, z)
  if (k == 0) {
    Dm <- matrix(Bz, ncol = 1)
    H2 <- matrix(hyper$h_alpha, 1, 1)
  } else {
    Dm <- cbind(Bz, C)
    H2 <- diag(c(hyper$h_alpha, rep(hyper$hm, k)), nrow = 1 + k)
  }
  ## mediator model
  Lambda_m <- crossprod(Dm) + H2
  SS_m <- 0
  for (j in 1:q) {
    mj <- M[, j]
    mu_j <- solve(Lambda_m, crossprod(Dm, mj))
    SS_m <- SS_m + (crossprod(mj) - t(mu_j) %*% Lambda_m %*% mu_j)
  }
  logdet_Lm <- as.numeric(determinant(Lambda_m, logarithm = TRUE)$modulus)
  log_fM <- -(q / 2) * logdet_Lm - ((hyper$nu + n * q) / 2) * log(hyper$lambda_e + SS_m)
  ## outcome model
  Bsc <- seq_len(q)
  qB <- length(Bsc)
  Dy <- cbind(M[, Bsc, drop = FALSE], Bz, C)
  H1_diag <- c(rep(hyper$h_beta, qB), hyper$h_gamma, rep(hyper$hy, k))
  Lambda_y <- crossprod(Dy) + diag(H1_diag, nrow = qB + 1 + k)
  mu_y <- solve(Lambda_y, crossprod(Dy, Y))
  SS_y <- as.numeric(crossprod(Y) - t(mu_y) %*% Lambda_y %*% mu_y)
  logdet_Ly <- as.numeric(determinant(Lambda_y, logarithm = TRUE)$modulus)
  log_fY <- -0.5 * logdet_Ly - ((hyper$nu + n) / 2) * log(hyper$lambda_eps + SS_y)
  ## prior on z: weights (eta2, 1 - eta1 - eta2, eta1) for (-1, 0, 1)
  w <- c(`-1` = hyper$eta2, `0` = 1 - hyper$eta1 - hyper$eta2, `1` = hyper$eta1)
  log_prior_z <- sum(log(w[as.character(z)]))
  log_fM + log_fY + log_prior_z
}

sample_inv_gamma <- function(n, shape, rate) 1 / rgamma(n, shape = shape, rate = rate)

make_valid_start <- function(p, prob_active = 0.1) {
  z <- rep(0, p)
  k <- max(2, ceiling(prob_active * p))
  active <- sample.int(p, size = k, replace = FALSE)
  z[active] <- sample(c(-1, 1), size = k, replace = TRUE)
  if (!any(z == 1)) z[sample(active, 1)] <- 1
  if (!any(z == -1)) z[sample(active, 1)] <- -1
  z
}

resolve_z_start <- function(z_init, d) {
  if (is.null(z_init)) return(make_valid_start(d))
  z_init <- as.numeric(z_init)
  if (length(z_init) != d)
    stop("z_init must have length d = ", d, " (got ", length(z_init), ")")
  if (!all(z_init %in% c(-1, 0, 1)))
    stop("z_init entries must all be in {-1, 0, 1}")
  if (!any(z_init == 1) || !any(z_init == -1))
    stop("z_init must contain at least one +1 and one -1 so the balance is defined")
  z_init
}

## Per-iteration update shared by both samplers: sign relabeling (gamma > 0),
## then draws of the mediator and outcome regression parameters given z.
mcmc_common_step <- function(z_curr, logX, M, C, Y, hyper, k, q, n) {
  Bz <- balance_from_log(logX, z_curr)

  ## sign check from a draw of the outcome coefficients
  Dy <- cbind(M, Bz, C)
  H1_diag <- c(rep(hyper$h_beta, q), hyper$h_gamma, rep(hyper$hy, k))
  Lambda_y <- crossprod(Dy) + diag(H1_diag, nrow = q + 1 + k)
  Lambda_y_inv <- solve(Lambda_y)
  mu_y <- as.numeric(Lambda_y_inv %*% crossprod(Dy, Y))
  SS_y <- as.numeric(crossprod(Y) - t(mu_y) %*% Lambda_y %*% mu_y)
  sig2_y_i <- sample_inv_gamma(1, shape = (hyper$nu + n) / 2, rate = (hyper$lambda_eps + SS_y) / 2)
  coef_y <- as.numeric(MASS::mvrnorm(1, mu_y, sig2_y_i * Lambda_y_inv))
  flip <- 0
  if (coef_y[q + 1] < 0) {
    z_curr <- -z_curr
    Bz <- -Bz
    flip <- 1
  }

  ## mediator model: M_j = alpha_j * Bz + C theta_j + e_j
  Dm <- cbind(Bz, C)
  DmtDm <- crossprod(Dm)
  mu_list <- vector("list", q)
  Linv_list <- vector("list", q)
  SS_m <- 0
  for (j in seq_len(q)) {
    Lambda_m_j <- DmtDm + diag(c(hyper$h_alpha, rep(hyper$hm, k)), nrow = 1 + k)
    Linv_j <- solve(Lambda_m_j)
    mu_j <- as.numeric(Linv_j %*% crossprod(Dm, M[, j]))
    SS_m <- SS_m + as.numeric(crossprod(M[, j]) - t(mu_j) %*% Lambda_m_j %*% mu_j)
    mu_list[[j]] <- mu_j
    Linv_list[[j]] <- Linv_j
  }
  sig2_m_i <- sample_inv_gamma(1, shape = (hyper$nu + n * q) / 2, rate = (hyper$lambda_e + SS_m) / 2)
  alpha_i <- numeric(q)
  for (j in seq_len(q)) {
    coef_j <- as.numeric(MASS::mvrnorm(1, mu_list[[j]], sig2_m_i * Linv_list[[j]]))
    alpha_i[j] <- coef_j[1]
  }

  ## outcome model: Y = M beta + gamma * Bz + C theta_y + eps
  Dy <- cbind(M, Bz, C)
  H1_diag <- c(rep(hyper$h_beta, q), hyper$h_gamma, rep(hyper$hy, k))
  Lambda_y <- crossprod(Dy) + diag(H1_diag, nrow = q + 1 + k)
  Lambda_y_inv <- solve(Lambda_y)
  mu_y <- as.numeric(Lambda_y_inv %*% crossprod(Dy, Y))
  SS_y <- as.numeric(crossprod(Y) - t(mu_y) %*% Lambda_y %*% mu_y)
  sig2_y_i <- sample_inv_gamma(1, shape = (hyper$nu + n) / 2, rate = (hyper$lambda_eps + SS_y) / 2)
  coef_y <- as.numeric(MASS::mvrnorm(1, mu_y, sig2_y_i * Lambda_y_inv))

  list(z = z_curr, flip = flip, alpha = alpha_i, beta = coef_y[seq_len(q)],
       gamma = coef_y[q + 1], sig2_m = sig2_m_i, sig2_y = sig2_y_i)
}

run_mcmc_Gibbs <- function(X, M, C, Y, hyper, n_iter = 20000, burn_in = 5000, z_init = NULL) {
  n <- nrow(X); d <- ncol(X); q <- ncol(M)
  C <- if (is.null(C)) matrix(nrow = n, ncol = 0) else as.matrix(C)
  k <- ncol(C)
  n_keep <- n_iter - burn_in
  z_samps <- matrix(0, nrow = n_keep, ncol = d)
  alpha_samps <- matrix(0, nrow = n_keep, ncol = q)
  beta_samps <- matrix(0, nrow = n_keep, ncol = q)
  gamma_samps <- sigma2_m <- sigma2_y <- numeric(n_keep)
  n_moved_trace <- numeric(n_iter)
  logX <- log(X)
  z_curr <- resolve_z_start(z_init, d)
  flip <- 0
  for (iter in seq_len(n_iter)) {
    gibbs_out <- gibbs_update_z_sweep(z_curr, logX, M, C, Y, hyper)
    z_curr <- gibbs_out$z
    n_moved_trace[iter] <- gibbs_out$n_moved
    step <- mcmc_common_step(z_curr, logX, M, C, Y, hyper, k, q, n)
    z_curr <- step$z
    flip <- flip + step$flip
    if (iter > burn_in) {
      idx <- iter - burn_in
      z_samps[idx, ] <- z_curr
      alpha_samps[idx, ] <- step$alpha
      sigma2_m[idx] <- step$sig2_m
      gamma_samps[idx] <- step$gamma
      beta_samps[idx, ] <- step$beta
      sigma2_y[idx] <- step$sig2_y
    }
  }
  cat("Mean sites changed per sweep:", mean(n_moved_trace), "\n")
  cat("Fraction of sweeps with any change:", mean(n_moved_trace > 0), "\n")
  cat("Flip rate:", flip / n_iter, "\n")
  list(z = z_samps, alpha = alpha_samps, sigma2_m = sigma2_m, gamma = gamma_samps,
       beta = beta_samps, sigma2_y = sigma2_y, n_moved_trace = n_moved_trace,
       move_rate = mean(n_moved_trace > 0), flip_rate = flip / n_iter)
}

run_mcmc_MH <- function(X, M, C, Y, hyper, n_iter = 20000, burn_in = 5000, z_init = NULL) {
  n <- nrow(X); d <- ncol(X); q <- ncol(M)
  C <- if (is.null(C)) matrix(nrow = n, ncol = 0) else as.matrix(C)
  k <- ncol(C)
  tau_a <- hyper$tau_a
  tau_d <- hyper$tau_d
  tau_s <- 1 - tau_a - tau_d
  value_pairs <- list(c(0, 1), c(0, -1), c(1, -1))
  pair_prob <- rep(1 / 3, 3)
  n_keep <- n_iter - burn_in
  z_samps <- matrix(0, nrow = n_keep, ncol = d)
  alpha_samps <- matrix(0, nrow = n_keep, ncol = q)
  beta_samps <- matrix(0, nrow = n_keep, ncol = q)
  gamma_samps <- sigma2_m <- sigma2_y <- numeric(n_keep)
  logX <- log(X)
  z_curr <- resolve_z_start(z_init, d)
  lp_curr <- log_posterior_z(z_curr, logX, M, C, Y, hyper)
  accept <- 0
  flip <- 0
  for (iter in seq_len(n_iter)) {
    pair_id <- sample(seq_along(value_pairs), 1, prob = pair_prob)
    l <- value_pairs[[pair_id]][1]
    lp <- value_pairs[[pair_id]][2]
    elig <- which(z_curr == l | z_curr == lp)
    if (length(elig) > 0) {
      j <- sample(elig, 1)
      curr_val <- z_curr[j]
      prop_val <- if (curr_val == l) lp else l
      z_prop <- z_curr
      z_prop[j] <- prop_val
      if (any(z_prop == 1) && any(z_prop == -1)) {
        if (curr_val == 0 && prop_val != 0) {
          tau_f <- tau_a; tau_b <- tau_d
        } else if (curr_val != 0 && prop_val == 0) {
          tau_f <- tau_d; tau_b <- tau_a
        } else {
          tau_f <- tau_s; tau_b <- tau_s
        }
        q_f <- (1 / 3) * tau_f * (1 / length(elig))
        elig_rev <- which(z_prop == l | z_prop == lp)
        q_b <- (1 / 3) * tau_b * (1 / length(elig_rev))
        lp_prop <- log_posterior_z(z_prop, logX, M, C, Y, hyper)
        log_r <- lp_prop - lp_curr + log(q_b) - log(q_f)
        if (log(runif(1)) < log_r) {
          z_curr <- z_prop
          lp_curr <- lp_prop
          accept <- accept + 1
        }
      }
    }
    step <- mcmc_common_step(z_curr, logX, M, C, Y, hyper, k, q, n)
    z_curr <- step$z
    flip <- flip + step$flip
    lp_curr <- log_posterior_z(z_curr, logX, M, C, Y, hyper)
    if (iter > burn_in) {
      idx <- iter - burn_in
      z_samps[idx, ] <- z_curr
      alpha_samps[idx, ] <- step$alpha
      sigma2_m[idx] <- step$sig2_m
      gamma_samps[idx] <- step$gamma
      beta_samps[idx, ] <- step$beta
      sigma2_y[idx] <- step$sig2_y
    }
  }
  cat("Acceptance rate:", accept / n_iter, "\n")
  cat("Flip rate:", flip / n_iter, "\n")
  list(z = z_samps, alpha = alpha_samps, sigma2_m = sigma2_m, gamma = gamma_samps,
       beta = beta_samps, sigma2_y = sigma2_y, accept_rate = accept / n_iter,
       flip_rate = flip / n_iter)
}
