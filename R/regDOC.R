## regDOC (Wang et al. 2019): Bayesian regularized difference-of-coefficients mediation
## analysis with multiple exposures. Used in Section 4 as a taxon-selection competitor.

#' Function for performing Bayesian mediation analysis with high-dimensional exposures.
#' @references Wang, Y. B., Chen, Z., Goldstein, J. M., Buck Louis, G. M., & Gilman, S. E. (2019). A Bayesian regularized mediation analysis with multiple exposures. Statistics in medicine, 38(5), 828-843.
#' @param X n by p exposures
#' @param M n by 1 mediator
#' @param Y n x 1 outcome
#' @param lambda1_init the scale parameter in the Laplace prior for alpha (direct effects)
#' @param lambda2_init the scale parameter in the Laplace prior for delta (indirect effects)
#' @param sigma2_0 variance of the regression coefficient of the mediator in the outcome regression model
#' @param theta1 shape parameter in the Gamma prior for lambda1
#' @param nu1 scale parameter in the Gamma prior for lambda1
#' @param theta2 shape parameter in the Gamma prior for lambda2
#' @param nu2 scale parameter in the Gamma prior for lambda2
#' @return posterior means of alpha (direct) and delta (indirect) effects
#' @details
#' This method uses the difference of coefficients approach to estimate indirect effects. Note Laplace priors
#' can be viewed as scale mixtures of normal distributions (a scale mixture of normal distributions with exponential mixing
#' distribution).  
#' 
regDOC_sampler <- function(X, M, Y, n_iter = 10000, burn_in = 5000,
                           lambda1_init = 0.1, lambda2_init = 0.1, sigma2_0 = 100,
                           theta1 = 1, nu1 = 0.1, theta2 = 1, nu2 = 0.1, verbose = TRUE) {
  n <- nrow(X)
  p <- ncol(X)
  
  # Initialize
  alpha <- rep(0, p)
  alpha_star <- rep(0, p)
  delta <- alpha_star - alpha
  beta <- 0
  mu1 <- mean(Y)
  mu3 <- mean(Y)
  sigma2_1 <- var(Y)
  eta <- 2
  tauA <- rep(1, p)
  tauD <- rep(1, p)
  lambda1 <- lambda1_init
  lambda2 <- lambda2_init
  
  samples <- list(alpha = matrix(0, n_iter - burn_in, p),
                  delta = matrix(0, n_iter - burn_in, p),
                  beta = numeric(n_iter - burn_in),
                  sigma2_1 = numeric(n_iter - burn_in),
                  eta = numeric(n_iter - burn_in))
  
  for (iter in 1:n_iter) {
    # Sample alpha
    VA <- solve(t(X) %*% X + diag(1 / tauA, p))
    mA <- VA %*% t(X) %*% (Y - mu1 - beta * M)
    alpha <- as.numeric(rmvnorm(1, mA, sigma2_1 * VA))
    
    # Sample alpha_star
    VD <- solve(t(X) %*% X + diag(1 / tauD, p))
    mD <- VD %*% t(X) %*% (Y - mu3)
    alpha_star <- as.numeric(rmvnorm(1, mD, eta * sigma2_1 * VD))
    
    # Update delta
    delta <- alpha_star - alpha
    
    # Sample beta
    V_beta <- 1 / (sum(M^2) + 1 / sigma2_0)
    m_beta <- V_beta * sum(M * (Y - mu1 - X %*% alpha))
    beta <- rnorm(1, m_beta, sqrt(sigma2_1 * V_beta))
    
    # Update mu1 and mu3
    mu1 <- rnorm(1, mean(Y - X %*% alpha - beta * M), sqrt(sigma2_1 / n))
    mu3 <- rnorm(1, mean(Y - X %*% alpha_star), sqrt(eta * sigma2_1 / n))
    
    # Update tauA and tauD
    tauA <- 1 / rinvgauss(p, mean = sqrt(lambda1^2 * sigma2_1 / (alpha^2 + 1e-8)), shape = lambda1^2)
    tauD <- 1 / rinvgauss(p, mean = sqrt(lambda2^2 * eta * sigma2_1 / (delta^2 + 1e-8)), shape = lambda2^2)
    
    # Update lambda1 and lambda2
    lambda1 <- sqrt(rgamma(1, shape = theta1 + p, rate = nu1 + sum(tauA) / 2))
    lambda2 <- sqrt(rgamma(1, shape = theta2 + p, rate = nu2 + sum(tauD) / 2))
    
    # Update sigma2_1
    residual1 <- Y - mu1 - X %*% alpha - beta * M
    residual3 <- Y - mu3 - X %*% alpha_star
    shape_sigma <- (2 * n + 2 * p + 1) / 2
    rate_sigma <- (sum(residual1^2) + sum(residual3^2) / eta + sum(alpha^2 / tauA) + sum(delta^2 / tauD) / eta + beta^2 / sigma2_0) / 2
    sigma2_1 <- 1 / rgamma(1, shape = shape_sigma, rate = rate_sigma)
    
    # Update eta (truncated InvGamma with eta >= 1)
    shape_eta <- (n + p) / 2 - 1
    rate_eta <- (sum(residual3^2) / sigma2_1 + sum(delta^2 / tauD) / sigma2_1) / 2
    repeat {
      eta <- 1 / rgamma(1, shape = shape_eta, rate = rate_eta)
      if (eta >= 1) break
    }
    
    # Store samples
    if (iter > burn_in) {
      idx <- iter - burn_in
      samples$alpha[idx, ] <- alpha
      samples$delta[idx, ] <- delta
      samples$beta[idx] <- beta
      samples$sigma2_1[idx] <- sigma2_1
      samples$eta[idx] <- eta
    }
    
    if (verbose && iter %% 10000 == 0) {
      cat("Iteration:", iter, "\n")
    }
  }
  
  return(samples)
}

library(mvtnorm)   # rmvnorm
library(statmod)   # rinvgauss

## functions to derive the precision recall curve

#' @param beta_samples: S x p matrix, S samples, p coefficients
# Score can also be defined by binary inclusion of the credible interval
credible_interval_score <- function(beta_samples) {
  apply(beta_samples, 2, function(samples) {
    ci <- quantile(samples, probs = c(0.025, 0.975))
    # score is defined as minimum absolute distance from 0 to the nearest bound
    if (ci[1] > 0 | ci[2] < 0) {
      # If 0 is excluded
      score = min(abs(ci[1]), abs(ci[2]))
    } else {
      score = 0
    }
    return(score)
  })
}

#' Function to compute posterior inclusion probability from the posterior samples
#' and a pre-specified cutoff
#' @param beta_samples S x p matrix, S samples, p coefficients
#' @param epsilon a prespecified cutoff for determining whether a coefficient is zero. 
compute_pip <- function(beta_samples, epsilon = 0.01) {
  apply(beta_samples, 2, function(samples) {
    mean(abs(samples) > epsilon)
  })
}



