## =============================================================================
## 04_figure6.R -- Figure 6: direct and indirect effects (left) and path
## coefficients (right) for the eight MEC-APS fits. Run from mec/ after
## 03_fit_mec.R (or submit_mec.sh) has produced all eight fits.
##
## Rows: the primary model (eta = 1/20, diabetes-adjusted, metformin users
## retained; filled circle) and seven fits that each change one thing (open
## triangles): eta = 1/3 or 1/5; adjustment for metformin or for neither; no
## adiposity covariate; metformin users excluded; participants with diabetes
## excluded (adjusting for neither).
## Output: <results_dir>/mec_forest_combined.pdf and .png
## =============================================================================
source("config.R")
suppressPackageStartupMessages({ library(ggplot2); library(cowplot) })
theme_set(theme_cowplot(font_size = 12))

pri_dir  <- file.path(results_dir, "primary")
sen_dir  <- file.path(results_dir, "sensitivity")
out_stem <- file.path(results_dir, "mec_forest_combined")

## One change at a time around the primary model.
pri_file <- function(a, e)
  sprintf("%s/pfat_ethall_metfRetained_adj%s_eta1over%d_gibbs_4chains_iter5000.rds",
          pri_dir, a, e)
spec <- rbind(
  data.frame(file  = pri_file("diabet", 20),
             label = "η = 1/20, diabetes",
             group = "primary"),
  data.frame(file  = pri_file("diabet", c(3, 5)),
             label = c("η = 1/3, diabetes", "η = 1/5, diabetes"),
             group = "sensitivity"),
  data.frame(file  = pri_file(c("metformin", "none"), 20),
             label = c("η = 1/20, metformin", "η = 1/20, neither"),
             group = "sensitivity"),
  data.frame(
    file = c(sprintf("%s/pfat_ethall_metfRetained_adjdiabet_covnoadip_eta1over20_gibbs_4chains_iter5000.rds", sen_dir),
             sprintf("%s/pfat_ethall_noMetf_adjdiabet_eta1over20_gibbs_4chains_iter5000.rds", sen_dir),
             sprintf("%s/pfat_ethall_noT2D_adjnone_eta1over20_gibbs_4chains_iter5000.rds", sen_dir)),
    label = c("no adiposity",
              "metformin excluded",
              "diabetes excluded"),
    group = "sensitivity"))
if (!all(file.exists(spec$file)))
  stop("missing fit(s) -- run submit_mec.sh (and its aggregate step) first:\n  ",
       paste(spec$file[!file.exists(spec$file)], collapse = "\n  "))

D <- do.call(rbind, lapply(seq_len(nrow(spec)), function(i) {
  d <- readRDS(spec$file[i])$diagnostics
  data.frame(setting = spec$label[i], group = spec$group[i],
             param = c("DE", "IE", "alpha", "beta"),
             mean = d[[4]], lo = d[[5]], hi = d[[6]])
}))
D$setting <- factor(D$setting, levels = rev(spec$label))   # first spec row on top
D$group   <- factor(D$group, levels = c("primary", "sensitivity"))
D$param   <- factor(D$param, levels = c("DE", "IE", "alpha", "beta"),
                    labels = c("DE (γ)", "IE (αβ)",
                               "α (B → LBP)", "β (LBP → y)"))

## zero reference only where zero is inside the plotted range (IE, beta);
## in the DE and alpha panels every CrI sits far above zero, and anchoring
## their axes at zero wastes most of the panel width
vline_df <- data.frame(param = factor(c("IE (αβ)", "β (LBP → y)"),
                                      levels = levels(D$param)))
forest <- function(sub, ttl) {
  vl <- vline_df[vline_df$param %in% unique(sub$param), , drop = FALSE]
  ggplot(sub, aes(mean, setting, shape = group)) +
  geom_vline(data = vl, aes(xintercept = 0), linetype = 2, colour = "grey60") +
  geom_pointrange(aes(xmin = lo, xmax = hi), orientation = "y", size = .4) +
  geom_point(data = function(d) d[d$group == "primary", ], size = 2.6) +
  scale_shape_manual(values = c(primary = 16, sensitivity = 2), name = NULL) +
  facet_wrap(~ param, scales = "free_x", nrow = 1) +
  labs(x = "posterior mean (95% CrI)", y = NULL, title = ttl) +
  theme(panel.spacing = unit(1.4, "lines"))
}

p_eff  <- forest(D[grepl("^DE|^IE", D$param), ], "Effect decomposition") +
  theme(legend.position = "none")
p_path <- forest(D[grepl("LBP", D$param), ], "Path coefficients") +
  theme(axis.text.y = element_blank(), legend.position = "none")
leg <- get_legend(forest(D, NULL) + theme(legend.position = "bottom"))

fig <- plot_grid(plot_grid(p_eff, p_path, ncol = 2, rel_widths = c(1.45, 1)),
                 leg, nrow = 2, rel_heights = c(1, .07))
## PDF with Greek glyphs: quartz on macOS (cairo is broken there when XQuartz
## is absent -- cairo_pdf "succeeds" but writes nothing); cairo_pdf elsewhere
if (capabilities("aqua")) {
  quartz(type = "pdf", file = paste0(out_stem, ".pdf"), width = 13, height = 4.1)
  print(fig); dev.off()
} else {
  ggsave(paste0(out_stem, ".pdf"), fig, width = 13, height = 4.1, device = cairo_pdf)
}
ggsave(paste0(out_stem, ".png"), fig, width = 13, height = 4.1, dpi = 150)
stopifnot(file.size(paste0(out_stem, ".pdf")) > 0)
cat("wrote", out_stem, ".pdf / .png\n")
