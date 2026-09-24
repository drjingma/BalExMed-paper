## =============================================================================
## 07_convergence_figures.R -- convergence diagnostics for the primary MEC-APS
## fit: traces and running means of the direct and indirect effects, and mixing
## of the balance configuration across chains. Run from mec/ after 03_fit_mec.R.
## Output (in <results_dir>):
##   convergence_trace.pdf    traces and running means, one colour per chain
##   convergence_mixing.pdf   cross-chain inclusion probabilities and sweep activity
##   convergence_numbers.txt  numbers quoted in the supplement, including split-Rhat
##                            and effective sample size per quantity
## =============================================================================
source("config.R")
suppressPackageStartupMessages({ library(ggplot2); library(cowplot) })
theme_set(theme_cowplot(font_size = 11))

pri <- file.path(results_dir, "primary",
                 "pfat_ethall_metfRetained_adjdiabet_eta1over20_gibbs_4chains_iter5000.rds")
if (!file.exists(pri)) stop("primary fit not found: ", pri)
o <- readRDS(pri)
fits <- o$fits
nch <- length(fits)
burn_in <- o$config$burn_in
cols <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7")[seq_len(nch)]

draws <- do.call(rbind, lapply(seq_len(nch), function(i) {
  g  <- as.numeric(fits[[i]]$gamma)
  ie <- as.numeric(fits[[i]]$alpha[, 1] * fits[[i]]$beta[, 1])
  data.frame(chain = factor(i), draw = seq_along(g), gamma = g, ie = ie)
}))
run_mean <- function(x) cumsum(x) / seq_along(x)
draws$gamma_rm <- ave(draws$gamma, draws$chain, FUN = run_mean)
draws$ie_rm    <- ave(draws$ie,    draws$chain, FUN = run_mean)

## traces are thinned for plotting only; running means use every draw
thin_plot <- 4
panel <- function(y, ylab, ttl, thin = thin_plot) {
  d <- draws[draws$draw %% thin == 0, ]
  ggplot(d, aes(draw, .data[[y]], colour = chain)) +
    geom_line(linewidth = 0.2, alpha = 0.65) +
    scale_colour_manual(values = cols, name = "chain") +
    labs(x = "retained draw", y = ylab, title = ttl)
}
p1 <- panel("gamma",    expression(gamma),          "Direct effect")
p2 <- panel("ie",       expression(alpha * beta),   "Indirect effect")
p3 <- panel("gamma_rm", expression(bar(gamma)),     "Direct effect, running mean")
p4 <- panel("ie_rm",    expression(bar(alpha * beta)), "Indirect effect, running mean")
leg <- get_legend(p1 + theme(legend.position = "bottom"))
body <- plot_grid(plotlist = lapply(list(p1, p2, p3, p4), function(p) p + theme(legend.position = "none")),
                  ncol = 2, align = "hv")
ggsave(file.path(results_dir, "convergence_trace.pdf"),
       plot_grid(body, leg, ncol = 1, rel_heights = c(1, 0.06)), width = 9, height = 6)

## mixing over the balance configuration
pip <- vapply(fits, function(f) colMeans(f$z != 0), numeric(ncol(fits[[1]]$z)))
pip_df <- do.call(rbind, lapply(2:nch, function(i)
  data.frame(chain = factor(i), x = pip[, 1], y = pip[, i])))
cc <- cor(pip)
pA <- ggplot(pip_df, aes(x, y, colour = chain)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60") +
  geom_point(size = 1, alpha = 0.7) +
  scale_colour_manual(values = cols[-1], name = "chain") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "inclusion probability, chain 1", y = "chains 2 to 4",
       title = "Taxon inclusion probabilities")
mv <- do.call(rbind, lapply(seq_len(nch), function(i)
  data.frame(chain = factor(i), iteration = seq_along(fits[[i]]$n_moved_trace),
             moved = fits[[i]]$n_moved_trace)))
pB <- ggplot(mv[mv$iteration %% 2 == 0, ], aes(iteration, moved, colour = chain)) +
  geom_vline(xintercept = burn_in, linetype = 3, colour = "grey40") +
  geom_line(linewidth = 0.2, alpha = 0.7) +
  scale_colour_manual(values = cols, name = "chain") +
  labs(x = "iteration", y = "coordinates changed", title = "Gibbs sweep activity")
ggsave(file.path(results_dir, "convergence_mixing.pdf"),
       plot_grid(pA + theme(legend.position = "none"), pB + theme(legend.position = "none"),
                 ncol = 2, rel_widths = c(1, 1.25), labels = c("A", "B")),
       width = 9, height = 3.6)

## numbers
d <- o$diagnostics
out <- character(0)
say <- function(...) { l <- sprintf(...); out <<- c(out, l); cat(l, "\n") }
say("Primary fit: %d chains of %d iterations, burn-in %d, %d retained draws per chain",
    nch, o$config$n_iter, burn_in, nrow(fits[[1]]$z))
say("split-Rhat: %s", paste(sprintf("%s %.3f", d$param, d$rhat), collapse = "; "))
say("effective sample size: %s", paste(sprintf("%s %d", d$param, round(d$ess)), collapse = "; "))
say("minimum pairwise correlation of inclusion probabilities across chains: %.3f",
    min(cc[upper.tri(cc)]))
say("maximum across-chain standard deviation of inclusion probability: %.3f",
    max(o$taxa$P_incl_sd_across_chains))
say("sweeps changing at least one coordinate: %.1f%% to %.1f%% of iterations by chain",
    100 * min(sapply(fits, `[[`, "move_rate")), 100 * max(sapply(fits, `[[`, "move_rate")))
say("coordinates changed per sweep after burn-in: median %.0f, mean %.1f",
    median(mv$moved[mv$iteration > burn_in]), mean(mv$moved[mv$iteration > burn_in]))
say("sign relabelling triggered in %.2f%% to %.2f%% of iterations by chain",
    100 * min(sapply(fits, `[[`, "flip_rate")), 100 * max(sapply(fits, `[[`, "flip_rate")))
## all eight fits
all_files <- c(list.files(file.path(results_dir, "primary"), "^pfat_.*\\.rds$", full.names = TRUE),
               list.files(file.path(results_dir, "sensitivity"), "^pfat_.*\\.rds$", full.names = TRUE))
oo <- lapply(all_files, readRDS)
say("across all %d fits: largest split-Rhat %.3f, smallest effective sample size %d, smallest cross-chain inclusion correlation %.3f",
    length(oo), max(sapply(oo, function(x) max(x$diagnostics$rhat))),
    round(min(sapply(oo, function(x) min(x$diagnostics$ess)))),
    min(sapply(oo, `[[`, "pip_cross_cor")))
writeLines(out, file.path(results_dir, "convergence_numbers.txt"))
