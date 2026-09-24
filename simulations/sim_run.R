# =============================================================================
# sim_run.R -- replication driver for the simulation study in Section 4.
# Runs one grid cell (R replicates in parallel) and saves per-replicate rows plus
# an aggregated summary with Monte Carlo standard errors. Run from simulations/.
# Results go to $SIM_RESULTS (default: results/).
#
# CLUSTER USAGE (recommended -- resilient to time-limit cancellations):
# The SLURM array enumerates (cell x replicate-CHUNK). Each task runs a chunk of
# replicates of one cell and writes ONE file PER replicate (with resume), so a
# task killed by the time limit loses nothing; just resubmit -- finished reps are
# skipped. build_tasks() must use the SAME chunk as the --array range.
#   Rscript sim_run.R ntasks [chunk]                 # #tasks -> set --array=0-(N-1)
#   Rscript sim_run.R task <task_id> [cores] [chunk] [n_iter] [burn_in] [comp]
#   Rscript sim_run.R aggregate                       # per-rep parts -> per-cell
#                                                     # summaries + completeness report
# Other:
#   Rscript sim_run.R count | grid | tasktable [chunk]
#   Rscript sim_run.R <cell_idx> [cores] [R] [n_iter] [burn_in] [comp]  # legacy whole-cell
#
# Output: results/parts/<cell_id>/rep_<seed>.rds   (one per replicate)
#         results/<cell_id>.rds / _summary.csv     (from `aggregate`)
#         results/ALL_summaries.csv, completeness.csv
# =============================================================================

sim_results_dir <- function() Sys.getenv("SIM_RESULTS", "results")

## ---- experiment grid (Section 4.1) ------------------------------------------
build_grid <- function() {
  base <- list(d = 50, rho = 0.2, eta = 1/5, alpha = 1/3, gamma = 1/3,  # default matches MEC moderate prior
               a_plus = 3, a_minus = 3, delta_m = 0, delta_y = 0,
               sampler = "gibbs", n_iter = 2e4, burn_in = 1e4, R = 200,
               competitors = FALSE)
  rows <- list()
  add  <- function(study, ...) rows[[length(rows) + 1]] <<-
    modifyList(base, c(list(study = study), list(...)))

  ## Main grid (2x2x2): competitors ON
  for (n in c(100, 200)) for (b in c(0, 2)) for (sc in c("I", "II"))
    add("main", n = n, beta = b, scenario = sc, competitors = TRUE)
  ## S0 unmeasured confounding (Scenario III)
  for (dl in list(c(0, 0), c(1, 1), c(2, 2))) for (b in c(0, 2)) for (n in c(100, 200))
    add("S0", n = n, beta = b, scenario = "III", delta_m = dl[1], delta_y = dl[2])
  ## S1 prior on z  (matches the MEC prior levels: 1/3, 1/5, 1/20)
  for (e in c(1/3, 1/5, 1/20)) for (b in c(0, 2))
    add("S1", n = 100, beta = b, scenario = "I", eta = e)
  ## S2 dimensionality
  for (dd in c(30, 50, 100)) for (b in c(0, 2))
    add("S2", n = 100, beta = b, scenario = "I", d = dd)
  ## S3 balance sparsity (number of active taxa per side)
  for (a in c(1, 3, 5, 10)) for (b in c(0, 2))
    add("S3", n = 100, beta = b, scenario = "I", a_plus = a, a_minus = a)
  ## Gibbs-vs-MH (same default datasets/seeds, both samplers; MH runs longer)
  add("GvMH", n = 100, beta = 2, scenario = "I", sampler = "gibbs")
  add("GvMH", n = 100, beta = 2, scenario = "I", sampler = "mh",
      n_iter = 1e5, burn_in = 1e4)

  g <- do.call(rbind, lapply(rows, function(r)
    as.data.frame(r, stringsAsFactors = FALSE)))
  g$idx     <- seq_len(nrow(g)) - 1L                       # 0-based (SLURM array)
  g$cell_id <- sprintf("%02d_%s_n%d_d%d_sc%s_b%g_eta%.3g_a%d_dm%g_%s",
                       g$idx, g$study, g$n, g$d, g$scenario, g$beta, g$eta,
                       g$a_plus, g$delta_m, g$sampler)
  g
}

## ---- aggregate R replicate rows -> one summary row per method --------------
summarise_cell <- function(res) {
  se <- function(x) { x <- x[!is.na(x)]; if (!length(x)) NA else sd(x)/sqrt(length(x)) }
  do.call(rbind, lapply(split(res, res$method), function(d) {
    IEt <- d$IE_true[1]
    data.frame(method = d$method[1], R = nrow(d), n_est = sum(!is.na(d$IE)),
               IE_mean = mean(d$IE, na.rm = TRUE),
               bias = mean(d$IE - IEt, na.rm = TRUE), bias_se = se(d$IE),
               RMSE = sqrt(mean((d$IE - IEt)^2, na.rm = TRUE)),
               cover95 = mean(d$cover, na.rm = TRUE), cover_se = se(d$cover),
               width = mean(d$width, na.rm = TRUE),
               rej0 = mean(d$reject0, na.rm = TRUE), rej0_se = se(d$reject0),
               DE_bias = if ("DE" %in% names(d)) mean(d$DE - d$DE_true, na.rm = TRUE) else NA,
               DE_cover95 = if ("DE_cover" %in% names(d)) mean(d$DE_cover, na.rm = TRUE) else NA,
               DE_width = if ("DE_width" %in% names(d)) mean(d$DE_width, na.rm = TRUE) else NA,
               AUPRC = mean(d$auprc, na.rm = TRUE), AUPRC_se = se(d$auprc),
               F1 = mean(d$f1, na.rm = TRUE),
               balrec = mean(d$brecov, na.rm = TRUE), balrec_se = se(d$brecov),
               ess_ie = if ("ess_ie" %in% names(d)) mean(d$ess_ie, na.rm = TRUE) else NA,
               fit_cpu = if ("fit_cpu" %in% names(d)) mean(d$fit_cpu, na.rm = TRUE) else NA,
               ess_per_sec = if ("ess_ie_per_sec" %in% names(d)) mean(d$ess_ie_per_sec, na.rm = TRUE) else NA,
               mix_rate = if ("mix_rate" %in% names(d)) mean(d$mix_rate, na.rm = TRUE) else NA,
               row.names = NULL)
  }))
}

## ---- run one cell ----------------------------------------------------------
run_cell <- function(task_id, cores, R_override = NA, ni_override = NA, bi_override = NA,
                     comp_override = NA) {
  library(parallel)
  g <- build_grid()
  if (!(task_id %in% g$idx)) stop("task_id ", task_id, " out of range [0, ", max(g$idx), "]")
  row <- g[g$idx == task_id, ]
  cfg <- list(n = row$n, d = row$d, rho = row$rho, scenario = row$scenario,
              beta = row$beta, eta = row$eta, alpha = row$alpha, gamma = row$gamma,
              a_plus = row$a_plus, a_minus = row$a_minus,
              delta = c(row$delta_m, row$delta_y), sampler = row$sampler,
              n_iter = if (!is.na(ni_override)) ni_override else row$n_iter,
              burn_in = if (!is.na(bi_override)) bi_override else row$burn_in,
              competitors = if (!is.na(comp_override)) as.logical(as.integer(comp_override))
                            else as.logical(row$competitors))
  R <- if (!is.na(R_override)) R_override else row$R

  cat(sprintf("cell %d (%s): n=%d d=%d scen=%s beta=%g eta=%.3g a=(%d,%d) delta=(%g,%g) %s | R=%d | cores=%d\n",
              task_id, row$study, cfg$n, cfg$d, cfg$scenario, cfg$beta, cfg$eta,
              cfg$a_plus, cfg$a_minus, cfg$delta[1], cfg$delta[2], cfg$sampler, R, cores))
  t0 <- Sys.time()
  reps <- mclapply(seq_len(R), function(s) {
    r <- tryCatch(run_replicate(cfg, seed = s), error = function(e) NULL)
    r
  }, mc.cores = cores)
  fail <- which(vapply(reps, is.null, logical(1)))
  if (length(fail)) cat("  WARNING: ", length(fail), " replicate(s) failed:",
                        paste(head(fail, 10), collapse = ","), "\n")
  res <- do.call(rbind, reps)
  elapsed_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  n_done <- length(unique(res$seed))
  cat(sprintf("  done in %.0f s (%d/%d replicates)\n", elapsed_sec, n_done, R))

  summ <- summarise_cell(res)
  summ <- cbind(cell_id = row$cell_id, study = row$study, summ,
                elapsed_sec = round(elapsed_sec, 1),          # cell wall time (all reps, on 1 node)
                cores = cores, sec_per_rep = round(elapsed_sec / max(1, n_done), 2),
                row.names = NULL)
  out_dir <- sim_results_dir(); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  saveRDS(list(config = cfg, row = row, R = R, per_rep = res, summary = summ,
               pip_by_method = avg_pip(res),          # replicate-averaged per-taxon PIP/score
               elapsed_sec = elapsed_sec, cores = cores, n_done = n_done,
               created = Sys.time()),
          file.path(out_dir, paste0(row$cell_id, ".rds")))
  write.csv(summ, file.path(out_dir, paste0(row$cell_id, "_summary.csv")), row.names = FALSE)
  cat("  saved: ", file.path(out_dir, paste0(row$cell_id, ".rds")), "\n")
  sm <- summ[, c("method","bias","RMSE","cover95","rej0","AUPRC","balrec")]
  num <- vapply(sm, is.numeric, logical(1)); sm[num] <- round(sm[num], 3)
  print(sm, row.names = FALSE)
}

## ---- 2-D task enumeration: (cell x replicate-chunk) for the SLURM array -----
cfg_from_row <- function(row, ni = NA, bi = NA, comp = NA) {
  list(n = row$n, d = row$d, rho = row$rho, scenario = row$scenario,
       beta = row$beta, eta = row$eta, alpha = row$alpha, gamma = row$gamma,
       a_plus = row$a_plus, a_minus = row$a_minus,
       delta = c(row$delta_m, row$delta_y), sampler = row$sampler,
       n_iter = if (!is.na(ni)) ni else row$n_iter,
       burn_in = if (!is.na(bi)) bi else row$burn_in,
       competitors = if (!is.na(comp)) as.logical(as.integer(comp)) else as.logical(row$competitors))
}

# one row per (cell, replicate-chunk); task_id is the flat 0-based SLURM index.
build_tasks <- function(chunk_size = 20) {
  g <- build_grid()
  tk <- do.call(rbind, lapply(seq_len(nrow(g)), function(i) {
    starts <- seq.int(1L, g$R[i], by = chunk_size)
    data.frame(cell_idx = g$idx[i], cell_id = g$cell_id[i],
               seed_from = starts, seed_to = pmin(starts + chunk_size - 1L, g$R[i]),
               stringsAsFactors = FALSE)
  }))
  tk$task_id <- seq_len(nrow(tk)) - 1L
  tk[, c("task_id", "cell_idx", "cell_id", "seed_from", "seed_to")]
}

# run one array task = a chunk of replicates of one cell; one file PER replicate,
# written atomically, skipped if already present (resume after a cancellation).
run_task <- function(task_id, cores, chunk_size = 20,
                     ni_override = NA, bi_override = NA, comp_override = NA) {
  library(parallel)
  g  <- build_grid(); tk <- build_tasks(chunk_size)
  if (!(task_id %in% tk$task_id))
    stop("task_id ", task_id, " out of range [0, ", max(tk$task_id), "] for chunk=", chunk_size)
  t   <- tk[tk$task_id == task_id, ]
  row <- g[g$idx == t$cell_idx, ]
  cfg <- cfg_from_row(row, ni_override, bi_override, comp_override)
  seeds <- t$seed_from:t$seed_to
  part_dir <- file.path(sim_results_dir(), "parts", row$cell_id)
  dir.create(part_dir, recursive = TRUE, showWarnings = FALSE)
  cat(sprintf("task %d | cell %d %s | seeds %d-%d | %s n=%d d=%d scen=%s b=%g eta=%.3g a=%d dm=%g | cores=%d\n",
              task_id, row$idx, row$cell_id, t$seed_from, t$seed_to, row$sampler,
              cfg$n, cfg$d, cfg$scenario, cfg$beta, cfg$eta, cfg$a_plus, cfg$delta_m, cores))
  t0 <- Sys.time()
  invisible(mclapply(seeds, function(s) {
    fp <- file.path(part_dir, sprintf("rep_%04d.rds", s))
    if (file.exists(fp)) return(invisible(NULL))            # resume: already done
    r <- tryCatch(run_replicate(cfg, seed = s), error = function(e) NULL)
    if (!is.null(r)) {
      r$cell_id <- row$cell_id; r$idx <- row$idx
      tmp <- paste0(fp, ".tmp-", Sys.getpid())              # atomic write: no half files
      saveRDS(r, tmp); file.rename(tmp, fp)
    }
    invisible(NULL)
  }, mc.cores = cores))
  done <- length(list.files(part_dir, pattern = "^rep_[0-9]+\\.rds$"))
  cat(sprintf("  %.0f s; cell now has %d/%d replicate files\n",
              as.numeric(difftime(Sys.time(), t0, units = "secs")), done, row$R))
}

## ---- aggregate per-replicate parts -> per-cell summary + combined table -----
aggregate_parts <- function() {
  out_dir <- sim_results_dir(); parts <- file.path(out_dir, "parts")
  g <- build_grid(); summaries <- list(); comp <- list()
  for (i in seq_len(nrow(g))) {
    row <- g[i, ]; pd <- file.path(parts, row$cell_id)
    fs <- list.files(pd, pattern = "^rep_[0-9]+\\.rds$", full.names = TRUE)
    if (!length(fs)) next
    res  <- do.call(rbind, lapply(fs, readRDS))
    nrep <- length(unique(res$seed))
    summ <- summarise_cell(res)
    summ <- cbind(cell_id = row$cell_id, study = row$study, summ,
                  n_rep = nrep, R_target = row$R, complete = (nrep >= row$R),
                  row.names = NULL)
    saveRDS(list(row = row, per_rep = res, summary = summ, n_rep = nrep,
                 pip_by_method = avg_pip(res),         # replicate-averaged per-taxon PIP/score
                 created = Sys.time()), file.path(out_dir, paste0(row$cell_id, ".rds")))
    write.csv(summ, file.path(out_dir, paste0(row$cell_id, "_summary.csv")), row.names = FALSE)
    summaries[[length(summaries) + 1]] <- summ
    comp[[length(comp) + 1]] <- data.frame(cell_idx = row$idx, cell_id = row$cell_id,
      study = row$study, n_rep = nrep, R_target = row$R, complete = (nrep >= row$R))
  }
  if (!length(summaries)) stop("no replicate parts under ", parts,
                               " -- run tasks first (Rscript sim_run.R task ...)")
  write.csv(do.call(rbind, summaries), file.path(out_dir, "ALL_summaries.csv"), row.names = FALSE)
  compdf <- do.call(rbind, comp)
  write.csv(compdf, file.path(out_dir, "completeness.csv"), row.names = FALSE)
  cat("aggregated", nrow(compdf), "cells with data:",
      sum(compdf$complete), "complete,", sum(!compdf$complete), "incomplete\n")
  inc <- compdf[!compdf$complete, , drop = FALSE]
  if (nrow(inc)) {
    cat("INCOMPLETE (resubmit the array to fill; finished reps are skipped):\n")
    for (j in seq_len(nrow(inc)))
      cat(sprintf("  cell %d %s : %d/%d reps\n", inc$cell_idx[j], inc$cell_id[j],
                  inc$n_rep[j], inc$R_target[j]))
  }
}

## ---- CLI dispatch ----------------------------------------------------------
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  cmd  <- if (length(args) >= 1) args[1] else "grid"
  if (cmd == "count") {                                   # number of CELLS
    cat(nrow(build_grid()), "\n")
  } else if (cmd == "ntasks") {                           # number of TASKS (needs chunk)
    chunk <- if (length(args) >= 2) as.integer(args[2]) else 20L
    cat(nrow(build_tasks(chunk)), "\n")
  } else if (cmd == "grid") {
    print(build_grid()[, c("idx","study","n","d","scenario","beta","eta","a_plus",
                           "delta_m","sampler","R","competitors")])
  } else if (cmd == "tasktable") {
    chunk <- if (length(args) >= 2) as.integer(args[2]) else 20L
    print(build_tasks(chunk))
  } else if (cmd == "aggregate") {
    source("sim_lib.R"); aggregate_parts()
  } else if (cmd == "task") {                             # run one chunk of replicates
    source("sim_lib.R")
    tid   <- as.integer(args[2])
    cores <- if (length(args) >= 3) as.integer(args[3]) else parallel::detectCores()
    chunk <- if (length(args) >= 4) as.integer(args[4]) else 20L
    ni    <- if (length(args) >= 5) as.numeric(args[5]) else NA
    bi    <- if (length(args) >= 6) as.numeric(args[6]) else NA
    comp  <- if (length(args) >= 7) as.integer(args[7]) else NA
    run_task(tid, cores, chunk, ni, bi, comp)
  } else {                                                # legacy: numeric cell_idx -> whole cell
    source("sim_lib.R")
    task_id <- as.integer(cmd)
    cores   <- if (length(args) >= 2) as.integer(args[2]) else parallel::detectCores()
    R_over  <- if (length(args) >= 3) as.integer(args[3]) else NA
    ni_over <- if (length(args) >= 4) as.numeric(args[4]) else NA
    bi_over <- if (length(args) >= 5) as.numeric(args[5]) else NA
    comp_over <- if (length(args) >= 6) as.integer(args[6]) else NA
    run_cell(task_id, cores, R_over, ni_over, bi_over, comp_over)
  }
}
