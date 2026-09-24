# =============================================================================
# sim_report.R -- tables and figures of Section 4 from the outputs of sim_run.R.
# Runs on whatever cells are present, so it can be used while the array finishes.
#
# Usage (from simulations/):
#   Rscript sim_report.R [in_dir] [out_dir]
#     in_dir   default results
#     out_dir  default report   (tables/, figures/ and paper/ subdirectories)
#
# Paper outputs, written to <out_dir>/paper/:
#   tab_sim_IE.tex              Table 1 (indirect-effect estimation)
#   sim_selection_tpr_fdr.pdf   Figure 3 (TPR versus FDR)
#   sim_sens_combined.pdf       Figure 4 (sensitivity and confounding)
#   sim_gibbs_vs_mh.pdf         Figure 5 (Gibbs versus Metropolis-Hastings)
#   section4_numbers.txt        numbers quoted in the text of Section 4.3
# Supporting tables and diagnostic figures go to <out_dir>/tables and figures.
# =============================================================================

args    <- commandArgs(trailingOnly = TRUE)
in_dir  <- if (length(args) >= 1) args[1] else "results"
out_dir <- if (length(args) >= 2) args[2] else "report"
tdir <- file.path(out_dir, "tables"); fdir <- file.path(out_dir, "figures")
dir.create(tdir, recursive = TRUE, showWarnings = FALSE)
dir.create(fdir, recursive = TRUE, showWarnings = FALSE)

## ---- load cells ------------------------------------------------------------
files <- list.files(in_dir, pattern = "\\.rds$", full.names = TRUE)
files <- files[!grepl("/report/", files)]
if (!length(files)) stop("no cell .rds found in ", in_dir)
cells <- lapply(files, readRDS)
cells <- Filter(function(o) all(c("row","summary","per_rep") %in% names(o)), cells)
cat("loaded", length(cells), "cells from", in_dir, "\n")

rparams <- c("study","n","d","scenario","beta","eta","a_plus","delta_m","delta_y","sampler","n_iter","burn_in","idx")
# attach only the row params NOT already present (per_rep carries idx/cell_id;
# summary carries study/cell_id) -> avoid duplicate columns, which ggplot rejects.
addcols <- function(base, r, np) {
  extra <- setdiff(rparams, names(base))
  cbind(base, r[rep(1, np), extra, drop = FALSE], row.names = NULL)
}
S <- do.call(rbind, lapply(cells, function(o) addcols(o$summary, o$row, nrow(o$summary))))
P <- do.call(rbind, lapply(cells, function(o) addcols(o$per_rep, o$row, nrow(o$per_rep))))
studies <- sort(unique(S$study))
cat("studies present:", paste(studies, collapse = ", "), "\n\n")

## replicate-averaged TPR at a controlled FDR (default 0.1) for BalMed, per cell:
## walk taxa in decreasing PIP order and stop at the target true FDR. Used as the
## selection summary in the sensitivity (S1-S3) and confounding (S0) figures.
tpr_at_fdr <- function(o, method = "BalMed", target = 0.1) {  # returns mean + Monte Carlo SE
  if (!"score_vec" %in% names(o$per_rep)) return(c(est = NA_real_, se = NA_real_))
  nact  <- o$row$a_plus + o$row$a_minus
  label <- c(rep(1L, nact), rep(0L, o$row$d - nact))
  sv <- o$per_rep$score_vec[o$per_rep$method == method]
  sv <- sv[!vapply(sv, is.null, logical(1))]
  if (!length(sv)) return(c(est = NA_real_, se = NA_real_))
  tprs <- vapply(sv, function(s) {
    ord <- order(s, decreasing = TRUE); lab <- label[ord]; k <- seq_along(lab)
    tp <- cumsum(lab); fdr <- (k - tp) / k; tpr <- tp / nact
    ok <- fdr <= target; if (any(ok)) max(tpr[ok]) else 0
  }, numeric(1))
  c(est = mean(tprs), se = sd(tprs) / sqrt(length(tprs)))
}
## Monte Carlo SE of RMSE (delta method) from the per-replicate squared errors.
rmse_se <- function(o, method = "BalMed") {
  pr <- o$per_rep[o$per_rep$method == method, ]
  e2 <- (pr$IE - pr$IE_true)^2; e2 <- e2[is.finite(e2)]
  if (length(e2) < 2) return(NA_real_)
  r <- sqrt(mean(e2)); if (r == 0) return(NA_real_)
  sd(e2) / sqrt(length(e2)) / (2 * r)
}
cid  <- vapply(cells, function(o) o$row$cell_id, character(1))
tt   <- vapply(cells, tpr_at_fdr, numeric(2))          # 2 x ncells (est, se)
S$tpr_fdr10    <- unname(setNames(tt["est", ], cid)[S$cell_id])
S$tpr_fdr10_se <- unname(setNames(tt["se",  ], cid)[S$cell_id])
S$RMSE_se      <- unname(setNames(vapply(cells, rmse_se, numeric(1)), cid)[S$cell_id])
{ ss <- S[S$study %in% c("S1","S2","S3","S0") & S$method == "BalMed", ]
  cat("TPR at FDR=0.1 (BalMed) -- sensitivity/confounding cells:\n")
  print(ss[order(ss$study, ss$eta, ss$d, ss$a_plus, ss$delta_m, ss$beta),
           c("study","n","beta","eta","d","a_plus","delta_m","tpr_fdr10")], row.names = FALSE) }

## ---- helpers ---------------------------------------------------------------
fmtpm  <- function(x, se) ifelse(is.na(x), "--",
                                 sprintf("%.3f (%.3f)", x, ifelse(is.na(se), 0, se)))
fmt3   <- function(x) ifelse(is.na(x), "--", sprintf("%.3f", x))
esc    <- function(s) gsub("_", "\\\\_", s)

# minimal booktabs LaTeX table from a character matrix
write_latex <- function(mat, header, file, caption, label, align = NULL) {
  if (is.null(align)) align <- paste(rep("l", ncol(mat)), collapse = "")
  con <- file(file, "w")
  writeLines(c("\\begin{table}[t]\\centering",
               sprintf("\\caption{%s}\\label{%s}", caption, label),
               sprintf("\\begin{tabular}{%s}\\toprule", align),
               paste(paste(esc(header), collapse = " & "), "\\\\ \\midrule"),
               apply(mat, 1, function(r) paste(paste(r, collapse = " & "), "\\\\")),
               "\\bottomrule\\end{tabular}\\end{table}"), con)
  close(con); cat("  wrote", file, "\n")
}

## ---- ggplot setup ----------------------------------------------------------
suppressPackageStartupMessages({ library(ggplot2); library(cowplot) })
theme_set(theme_cowplot())
save_gg <- function(p, file, w, h) {
  ggsave(file, p, width = w, height = h, device = "pdf")
  cat("  wrote", file, "\n")
}
pdir <- file.path(out_dir, "paper")
dir.create(pdir, recursive = TRUE, showWarnings = FALSE)
save_paper <- function(p, paper_name, w, h) save_gg(p, file.path(pdir, paper_name), w, h)
numbers <- character(0)                       # lines for section4_numbers.txt
note <- function(...) { line <- sprintf(...); numbers <<- c(numbers, line); cat(line, "\n") }

## ===========================================================================
## TABLE 1 -- main-grid IE estimation (BalMed vs Oracle)
## ===========================================================================
if ("main" %in% studies) {
  d <- S[S$study == "main" & S$method %in% c("BalMed","Oracle"), ]
  d <- d[order(d$scenario, d$n, d$beta, d$method), ]
  d$effect <- ifelse(d$beta == 0, "TypeI", "Power")
  mat <- cbind(d$scenario, d$n, sprintf("%g", d$beta), d$method,
               fmtpm(d$bias, d$bias_se), fmt3(d$RMSE), fmt3(d$cover95),
               fmt3(d$width), fmt3(d$rej0))
  hdr <- c("Scenario","n","beta","Method","Bias (SE)","RMSE","Cover95","Width","Reject")
  write.csv(d[, c("scenario","n","beta","method","bias","bias_se","RMSE",
                  "cover95","width","rej0","rej0_se")],
            file.path(tdir, "main_IE.csv"), row.names = FALSE)
  write_latex(mat, hdr, file.path(tdir, "main_IE.tex"),
              "Indirect-effect estimation on the main grid. Reject = Type-I error at $\\beta=0$ and power at $\\beta=2$.",
              "tab:sim-IE")
}


## ---- Table 1 in the layout of the paper (BalMed, Oracle, PrinBal) ----------
if ("main" %in% studies) {
  num <- function(x) {                          # $0.123$ or $\phantom{-}0.123$ / $-0.123$
    v <- sprintf("%.3f", x)
    v[v == "-0.000"] <- "0.000"
    ifelse(substr(v, 1, 1) == "-", sprintf("$%s$", v), sprintf("$\\phantom{-}%s$", v))
  }
  plain <- function(x) sprintf("$%.3f$", x)
  get_row <- function(sc, nn, b, mth) S[S$study == "main" & S$scenario == sc & S$n == nn &
                                          S$beta == b & S$method == mth, ][1, ]
  lines <- c("\\begin{table}[h!]\\centering",
    "\\caption{Indirect-effect estimation under scenarios I \\& II. The oracle method is given the true support and signs of $\\bz$---though not the scenario-II weights, whereas the BalMed method needs to infer $\\bz$ from the data. PrinBal fixes the balance from the leading log-ratio principal component, without reference to the mediator or outcome, and returns a point estimate but no credible interval (its coverage, width, and rejection are therefore omitted, shown as \\,--\\,).}\\label{tab:sim-IE}",
    "\\begin{tabular}{cccllrrrr}\\toprule",
    "Scenario & $n$ & $\\beta$ & Method & Bias (SE) & RMSE & Cover95 & Width & Reject \\\\ \\midrule")
  for (sc in c("I", "II")) {
    for (nn in c(100, 200)) {
      for (b in c(0, 2)) {
        for (mth in c("BalMed", "Oracle", "PrinBal")) {
          r <- get_row(sc, nn, b, mth)
          lead <- if (mth == "BalMed") {
            paste0(if (nn == 100 && b == 0) sprintf("\\multirow{12}{*}{%s}", sc) else "", " & ",
                   if (b == 0) sprintf("\\multirow{6}{*}{$%d$}", nn) else "", " & ",
                   sprintf("\\multirow{3}{*}{$%g$}", b), "\n")
          } else ""
          stats <- if (mth == "PrinBal")
            sprintf("%s (%s) & %s & -- & -- & --", num(r$bias), plain(r$bias_se), plain(r$RMSE))
          else
            sprintf("%s (%s) & %s & %s & %s & %s", num(r$bias), plain(r$bias_se), plain(r$RMSE),
                    plain(r$cover95), plain(r$width), plain(r$rej0))
          lines <- c(lines, sprintf("%s & %s & %-7s & %s \\\\", lead, if (mth == "BalMed") "" else " &", mth, stats))
        }
        if (b == 0) lines <- c(lines, "\\cmidrule(lr){3-9}")
      }
      if (nn == 100) lines <- c(lines, "\\cmidrule(lr){2-9}")
    }
    if (sc == "I") lines <- c(lines, "\\midrule")
  }
  lines <- c(lines, "\\bottomrule\\end{tabular}\\end{table}")
  writeLines(lines, file.path(pdir, "tab_sim_IE.tex"))
  cat("  wrote", file.path(pdir, "tab_sim_IE.tex"), "\n")
}

## ===========================================================================
## TABLE 2 -- main-grid taxon selection (AUPRC / F1) + balance recovery
## ===========================================================================
if ("main" %in% studies) {
  d <- S[S$study == "main" & S$method %in% c("BalMed","regDOC","PrinBal"), ]
  d <- d[order(d$scenario, d$n, d$beta, d$method), ]
  mat <- cbind(d$scenario, d$n, sprintf("%g", d$beta), d$method,
               fmtpm(d$AUPRC, d$AUPRC_se), fmt3(d$F1), fmt3(d$balrec))
  hdr <- c("Scenario","n","beta","Method","AUPRC (SE)","F1","Bal.recov")
  write.csv(d[, c("scenario","n","beta","method","AUPRC","AUPRC_se","F1","balrec")],
            file.path(tdir, "main_selection.csv"), row.names = FALSE)
  write_latex(mat, hdr, file.path(tdir, "main_selection.tex"),
              "Taxon selection (AUPRC, $F_1$) and balance recovery on the main grid.",
              "tab:sim-selection")
}

## ===========================================================================
## FIGURE 1 -- IE distribution (main grid, BalMed): null (panel A) and non-null
##             (panel B) as two ggplots (each faceted n x scenario), plot_grid'd.
## ===========================================================================
if ("main" %in% studies) {
  pm <- P[P$study == "main" & P$method == "BalMed", ]
  pm$n_f  <- factor(paste0("n = ", pm$n))
  pm$sc_f <- factor(paste0("Scenario ", pm$scenario))
  make_ie <- function(sub, ttl) {
    tr <- unique(sub[, c("n_f", "sc_f", "IE_true")])
    ggplot(sub, aes(x = "", y = IE)) +
      geom_boxplot(width = 0.5, fill = "grey90", outlier.size = 0.5) +
      geom_point(data = tr, aes(x = "", y = IE_true), colour = "red", shape = 18, size = 3) +
      facet_grid(n_f ~ sc_f) +
      labs(x = NULL, y = expression("estimated IE " * alpha * beta), title = ttl)
  }
  pA <- make_ie(pm[pm$beta == 0, ], "null (beta = 0)")
  pB <- make_ie(pm[pm$beta == 2, ], "non-null (beta = 2)")
  pg <- plot_grid(pA, pB, labels = c("A", "B"), ncol = 2)
  save_gg(pg, file.path(fdir, "IE_boxplot_main.pdf"), 11, 5)
}

## ===========================================================================
## FIGURE 2 -- taxon selection: true positive rate vs false discovery rate,
##   BalMed vs regDOC on the main grid. Each curve is traced by the PIP
##   threshold: per replicate we walk the taxa in decreasing PIP order (top-k
##   path) and compute the true FDR and TPR from the known active set, then
##   average over replicates and over beta (selection is ~beta-invariant).
##   The dotted line marks FDR = 0.1; the TPR there is the controlled-FDR
##   operating point (selecting at the PIP threshold that yields FDR = 0.1).
## ===========================================================================
if ("main" %in% studies) {
  fdr_grid <- seq(0, 0.5, by = 0.005)
  tpr_path <- function(score, label) {                 # top-k path -> TPR on fdr_grid
    ord <- order(score, decreasing = TRUE); lab <- label[ord]
    tp <- cumsum(lab); k <- seq_along(lab)
    fdr <- (k - tp) / k; tpr <- tp / sum(label)
    vapply(fdr_grid, function(f) { ok <- fdr <= f; if (any(ok)) max(tpr[ok]) else 0 }, numeric(1))
  }
  rows <- list()
  for (o in cells) {
    if (o$row$study != "main") next
    nact  <- o$row$a_plus + o$row$a_minus
    label <- c(rep(1L, nact), rep(0L, o$row$d - nact))
    for (mth in c("BalMed", "regDOC")) {
      sv <- o$per_rep$score_vec[o$per_rep$method == mth]
      sv <- sv[!vapply(sv, is.null, logical(1))]
      if (!length(sv)) next
      mt <- rowMeans(vapply(sv, tpr_path, numeric(length(fdr_grid)), label = label))
      rows[[length(rows) + 1]] <- data.frame(
        scenario = o$row$scenario, n = o$row$n, method = mth, fdr = fdr_grid, tpr = mt)
    }
  }
  if (length(rows)) {
    cv <- aggregate(tpr ~ scenario + n + method + fdr, do.call(rbind, rows), mean)  # avg over beta
    op <- cv[abs(cv$fdr - 0.1) < 1e-9, ]
    cat("TPR at controlled FDR = 0.1 (BalMed vs regDOC):\n")
    print(op[order(op$scenario, op$n, op$method), c("scenario","n","method","tpr")], row.names = FALSE)
    for (mth in c("BalMed", "regDOC"))
      note("TPR at FDR = 0.1, %s, main grid (averaged over beta): %.0f%% to %.0f%%",
           mth, 100 * min(op$tpr[op$method == mth]), 100 * max(op$tpr[op$method == mth]))
    cv$n_f  <- factor(paste0("n = ", cv$n))
    cv$sc_f <- factor(paste0("Scenario ", cv$scenario))
    ptf <- ggplot(cv, aes(fdr, tpr, colour = method, linetype = n_f,
                          group = interaction(method, n_f))) +
      geom_vline(xintercept = 0.1, linetype = 3, colour = "grey50") +
      geom_line(linewidth = 0.7) +
      facet_wrap(~ sc_f) +
      scale_linetype_manual(values = c("n = 100" = "dashed", "n = 200" = "solid")) +
      coord_cartesian(xlim = c(0, 0.5), ylim = c(0, 1)) +
      labs(x = "false discovery rate", y = "true positive rate", colour = NULL, linetype = NULL)
    save_paper(ptf, "sim_selection_tpr_fdr.pdf", 9, 3.2)
  }
}

## ===========================================================================
## FIGURES 3 -- one-factor sensitivity (BalMed): AUPRC & IE-RMSE vs factor.
##   sens_row() returns the ggplot; each design is saved on its own AND the
##   three are stacked into one 3-row figure (sens_combined.pdf) for the paper.
## ===========================================================================
sens_row <- function(study, xvar, xlab, ttl, frac = FALSE, keep_n = NULL) {
  d <- S[S$study == study & S$method == "BalMed", ]
  if (!is.null(keep_n)) d <- d[d$n == keep_n, ]
  if (!nrow(d)) return(NULL)
  d$beta_f <- factor(paste0("beta=", d$beta))          # beta by line type, black
  long <- rbind(
    data.frame(x = d[[xvar]], beta_f = d$beta_f, metric = "IE RMSE",          val = d$RMSE,      se = d$RMSE_se),
    data.frame(x = d[[xvar]], beta_f = d$beta_f, metric = "TPR at FDR = 0.1", val = d$tpr_fdr10, se = d$tpr_fdr10_se))
  long$metric <- factor(long$metric, levels = c("IE RMSE", "TPR at FDR = 0.1"))  # IE left, TPR right
  if (frac) {                                            # label x by fractions, e.g. 1/3, 1/5, 1/20
    lv <- sprintf("1/%d", round(1 / sort(unique(long$x), decreasing = TRUE)))
    long$xg <- factor(sprintf("1/%d", round(1 / long$x)), levels = lv)
    p <- ggplot(long, aes(xg, val, linetype = beta_f, group = beta_f))
  } else {
    p <- ggplot(long, aes(x, val, linetype = beta_f, group = beta_f))
  }
  p + geom_line() + geom_point(size = 2) +
    geom_errorbar(aes(ymin = val - se, ymax = val + se), width = 0.12, linetype = 1, na.rm = TRUE) +
    scale_linetype_manual(values = c("beta=0" = "solid", "beta=2" = "dashed"), name = NULL) +
    facet_wrap(~ metric, scales = "free_y") +
    labs(x = xlab, y = NULL, title = ttl) + 
    theme(plot.title = element_text(size = 12))
}
r1 <- if ("S1" %in% studies) sens_row("S1", "eta",     "prior on z  (eta)",    "(A) Prior on z", frac = TRUE)
r2 <- if ("S2" %in% studies) sens_row("S2", "d",       "number of taxa  d",    "(B) Dimension")
r3 <- if ("S3" %in% studies) sens_row("S3", "a_plus",  "active taxa per side", "(C) Signal density")
r4 <- if ("S0" %in% studies) sens_row("S0", "delta_m", expression(delta),
                                      "(D) Unmeasured confounding", keep_n = 100)
if (!is.null(r1)) save_gg(r1, file.path(fdir, "sens_S1_prior.pdf"),    9, 4.2)
if (!is.null(r2)) save_gg(r2, file.path(fdir, "sens_S2_dim.pdf"),      9, 4.2)
if (!is.null(r3)) save_gg(r3, file.path(fdir, "sens_S3_sparsity.pdf"), 9, 4.2)
rows <- Filter(Negate(is.null), list(r1, r2, r3, r4))
if (length(rows) >= 3) {          # one figure, one row per design (A-C sensitivity, D confounding)
  leg  <- get_legend(rows[[1]] + theme(legend.position = "right"))   # line-type legend on the right
  body <- plot_grid(plotlist = lapply(rows, function(p) p + theme(legend.position = "none")),
                    ncol = 2, align = "v")
  save_paper(plot_grid(body, leg, ncol = 2, rel_widths = c(1, 0.12)),
             "sim_sens_combined.pdf", 12, 1.5 * length(rows))
}

## (S0 unmeasured confounding is now row D of sens_combined.pdf, above.)

## ===========================================================================
## GIBBS vs MH -- paired per-replicate SCATTER (x=Gibbs, y=MH) + diagonal
##   (same seeds -> same data, so metrics should lie on y = x)
##   efficiency numbers (ESS, ESS/sec, move/accept) are written to the CSV.
## ===========================================================================
if ("GvMH" %in% studies) {
  d <- S[S$study == "GvMH" & S$method == "BalMed", ]
  eff <- intersect(c("ess_ie","fit_cpu","ess_per_sec","mix_rate"), names(d))
  write.csv(d[, c("sampler","IE_mean","bias","RMSE","cover95","width","AUPRC","balrec", eff)],
            file.path(tdir, "gibbs_vs_mh.csv"), row.names = FALSE)
  cat("  wrote", file.path(tdir, "gibbs_vs_mh.csv"), "\n")
  g <- P[P$study == "GvMH" & P$method == "BalMed" & P$sampler == "gibbs", ]
  m <- P[P$study == "GvMH" & P$method == "BalMed" & P$sampler == "mh", ]
  if (nrow(g) && nrow(m)) {
    ## per-replicate TPR at a controlled FDR = 0.1, from the per-taxon scores
    tpr_one <- function(s, label, target = 0.1) {
      ord <- order(s, decreasing = TRUE); lab <- label[ord]; k <- seq_along(lab)
      tp <- cumsum(lab); fdr <- (k - tp) / k; tpr <- tp / sum(label)
      ok <- fdr <= target; if (any(ok)) max(tpr[ok]) else 0
    }
    nact  <- 2 * g$a_plus[1]; label <- c(rep(1L, nact), rep(0L, g$d[1] - nact))
    g$tpr <- vapply(g$score_vec, tpr_one, numeric(1), label = label)
    m$tpr <- vapply(m$score_vec, tpr_one, numeric(1), label = label)
    mets <- c(IE = "IE", DE = "DE", TPR = "tpr")                # one row: IE, DE, TPR
    long <- do.call(rbind, lapply(names(mets), function(nm) {
      mm <- merge(g[, c("seed", mets[[nm]])], m[, c("seed", mets[[nm]])],
                  by = "seed", suffixes = c(".g", ".m"))
      data.frame(metric = nm, gibbs = mm[[paste0(mets[[nm]], ".g")]],
                 mh = mm[[paste0(mets[[nm]], ".m")]])
    }))
    long$metric <- factor(long$metric, levels = c("IE", "DE", "TPR"))
    strip_labs <- c(IE = "IE", DE = "DE", TPR = "TPR at FDR = 0.1")
    p <- ggplot(long, aes(gibbs, mh)) +
      geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey40") +
      geom_point(alpha = 0.5, size = 1) +
      facet_wrap(~ metric, scales = "free", nrow = 1,      # aligned in one row (by column)
                 labeller = as_labeller(strip_labs)) +
      labs(x = "Gibbs", y = "MH",
           title = "Gibbs vs MH agreement (paired by replicate; dashed = y = x)") +
      theme(aspect.ratio = 1)                # square panels: x and y visually equal length
    save_gg(p, file.path(fdir, "gibbs_vs_mh.pdf"), 12, 4.5)
    p_agree <- p                             # keep for the combined figure below
  }
  ## efficiency as Gibbs/MH ratios, normalized so run length (Gibbs 2e4 vs MH 1e5
  ## iterations) does not confound: effective samples per CPU-second, move rate
  ## (already per-iteration), and CPU per iteration.
  if (nrow(d) >= 2 && all(c("ess_per_sec", "fit_cpu", "mix_rate", "n_iter") %in% names(d))) {
    gv <- function(col) d[[col]][d$sampler == "gibbs"]
    mv <- function(col) d[[col]][d$sampler == "mh"]
    cpu_iter <- function(who) who("fit_cpu") / who("n_iter")
    eff_df <- data.frame(
      metric = factor(c("ESS/sec", "move rate", "CPU/iter"),
                      levels = c("ESS/sec", "move rate", "CPU/iter")),
      ratio  = c(gv("ess_per_sec") / mv("ess_per_sec"),
                 gv("mix_rate")     / mv("mix_rate"),
                 cpu_iter(gv)       / cpu_iter(mv)))
    if (all(is.finite(eff_df$ratio))) {
      pe <- ggplot(eff_df, aes(metric, ratio)) +
        geom_hline(yintercept = 1, linetype = 2, colour = "grey40") +
        geom_segment(aes(xend = metric, y = 1, yend = ratio), colour = "grey60") +
        geom_point(size = 3) +
        geom_text(aes(label = sprintf("%.2g", ratio)), vjust = -0.7, size = 3.5) +
        scale_y_log10(expand = expansion(mult = c(0.05, 0.12))) +
        ## let the top label draw into the margin (never clipped by the panel edge,
        ## robust to shrinking the figure) and reserve top margin so it clears the title.
        coord_cartesian(clip = "off") +
        theme(plot.margin = margin(t = 12, r = 6, b = 6, l = 6)) +
        labs(x = NULL, y = "Gibbs / MH ratio (log10; dashed = equal)",
             title = "Gibbs vs MH: efficiency-summary ratios")
      save_gg(pe, file.path(fdir, "gibbs_vs_mh_efficiency.pdf"), 8, 4.5)
      ## one figure: agreement (left) + efficiency (right)
      if (exists("p_agree")) {
        ## aspect.ratio = 1 makes Panel B a square panel; with rel_widths 3:1 its
        ## column is ~one A-facet wide, so its height matches each A subplot.
        save_paper(plot_grid(p_agree + labs(title = "Agreement (paired by replicate)"),
                             pe + labs(title = "Efficiency (Gibbs / MH)") +
                                theme(aspect.ratio = 1.2),
                             ncol = 2, rel_widths = c(3, 1), labels = c("A", "B")),
                   "sim_gibbs_vs_mh.pdf", 13, 4.5)
      }
    }
  }
}

## ---- numbers quoted in Section 4.3 -----------------------------------------
sens_tpr <- function(study, xvar) {
  d <- S[S$study == study & S$method == "BalMed", c(xvar, "beta", "tpr_fdr10")]
  a <- aggregate(tpr_fdr10 ~ get(xvar), data = d, FUN = mean)
  names(a) <- c(xvar, "tpr_mean_over_beta")
  a
}
if ("S1" %in% studies) {
  a <- sens_tpr("S1", "eta")
  note("S1 TPR at FDR = 0.1 by eta (mean over beta): %s",
       paste(sprintf("eta=1/%d: %.2f", round(1 / a$eta), a$tpr_mean_over_beta), collapse = "; "))
}
if ("S3" %in% studies) {
  a <- sens_tpr("S3", "a_plus")
  note("S3 TPR at FDR = 0.1 by active taxa per side (mean over beta): %s",
       paste(sprintf("a=%d: %.2f", a$a_plus, a$tpr_mean_over_beta), collapse = "; "))
}
if ("GvMH" %in% studies) {
  d <- S[S$study == "GvMH" & S$method == "BalMed", ]
  g <- d[d$sampler == "gibbs", ]; h <- d[d$sampler == "mh", ]
  if (nrow(g) && nrow(h)) {
    note("Gibbs/MH CPU per iteration: %.2f", (g$fit_cpu / g$n_iter) / (h$fit_cpu / h$n_iter))
    note("Gibbs move rate / MH acceptance rate: %.1f (MH acceptance %.3f)", g$mix_rate / h$mix_rate, h$mix_rate)
    note("Gibbs/MH effective samples per CPU-second (IE): %.1f", g$ess_per_sec / h$ess_per_sec)
    note("Gibbs CPU seconds per sweep: %.3f", g$fit_cpu / g$n_iter)
  }
}
writeLines(numbers, file.path(pdir, "section4_numbers.txt"))
cat("\n=== report written to", normalizePath(out_dir), "===\n")
