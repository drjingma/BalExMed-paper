# Reproducibility code: *Mediation Analysis with Compositional Exposures*

This folder contains the code that produced every table, figure and number in the
simulation study (Section 4) and the MEC-APS application (Section 5) of the paper. It
also contains the aggregated simulation results, so the Section 4 tables and figures can
be rebuilt in a few seconds without rerunning the simulations.

The MEC-APS data are not publicly available and are not included. The MEC-APS scripts
are provided so that the analysis can be inspected, and rerun by anyone with access to
the data.

For new analyses, use the R package **BalExMed**, which implements the method with
documentation, input checks and a faster sampler.

## Contents

```
R/
  balance_mediation_sampler.R   Gibbs and Metropolis-Hastings samplers used for the paper
  regDOC.R                      regDOC competitor (Wang et al., 2019)
simulations/
  sim_lib.R                     data generation, competing methods and metrics
  sim_run.R                     runs the simulation grid (one cell or one SLURM task)
  sim_report.R                  builds Table 1, Figures 3-5 and the Section 4.3 numbers
  submit_sim.sh                 SLURM array for the full grid
  results/                      aggregated results for all 42 cells (200 replicates each)
mec/
  config.R                      data and output paths
  01_read_asv_data.R            reads the ASV table and SILVA taxonomy
  02_build_mediation_inputs.R   builds the analysis object (genera, LBP, liver fat, covariates)
  03_fit_mec.R                  fits one model setting
  submit_mec.sh                 SLURM array for the eight fits in Figure 6
  04_figure6.R                  Figure 6
  05_table2.R                   Table 2
  06_section5_numbers.R         numbers quoted in the text of Section 5
  07_convergence_figures.R      convergence diagnostics for the supplement
sessionInfo.txt                 R and package versions used to check this folder
LICENSE                         GNU General Public License, version 3
```

Figures 1 and 2 are diagrams drawn directly in LaTeX.

## Requirements

R 4.4 was used; any recent version of R should work. Install the packages with

```r
install.packages(c("MASS", "coda", "PRROC", "mvtnorm", "statmod",
                   "ggplot2", "cowplot", "dplyr", "data.table"))
```

`parallel` ships with R. On macOS and Linux, chains and replicates run in parallel
through forking. `RhpcBLASctl` is optional and only limits BLAS threads.

## Where each result comes from

| Paper item | Script (run from its folder) | Output |
|---|---|---|
| Table 1 | `simulations/sim_report.R` | `simulations/report/paper/tab_sim_IE.tex` |
| Figure 3 | `simulations/sim_report.R` | `simulations/report/paper/sim_selection_tpr_fdr.pdf` |
| Figure 4 | `simulations/sim_report.R` | `simulations/report/paper/sim_sens_combined.pdf` |
| Figure 5 | `simulations/sim_report.R` | `simulations/report/paper/sim_gibbs_vs_mh.pdf` |
| Numbers in Section 4.3 | `simulations/sim_report.R` | `simulations/report/paper/section4_numbers.txt` |
| Figure 6 | `mec/04_figure6.R` | `mec/results/mec_forest_combined.pdf` |
| Table 2 | `mec/05_table2.R` | `mec/results/table2.tex` |
| Numbers in Section 5 | `mec/06_section5_numbers.R` | `mec/results/section5_numbers.txt` |
| Figures S1 and S2 | `mec/07_convergence_figures.R` | `mec/results/convergence_trace.pdf`, `convergence_mixing.pdf` |

## Simulations (Section 4)

### Rebuild the tables and figures from the included results

```sh
cd simulations
Rscript sim_report.R
```

This reads `results/` and writes `report/`. Running it on the included results
reproduces Table 1 and the numbers quoted in Section 4.3.


### Rerun the simulations

The design has 42 cells: the main grid (scenarios I and II, n = 100 and 200, beta = 0
and 2), unmeasured confounding (scenario III), sensitivity to the prior, to the number
of taxa and to the number of active taxa, and the Gibbs versus Metropolis-Hastings
comparison. Each cell has 200 replicates. Replicate r of every cell generates its data
after `set.seed(r)`. List the cells with

```sh
cd simulations
Rscript sim_run.R grid
```

A Gibbs fit with 2 x 10^4 iterations at n = 100 and d = 50 took roughly 75 CPU-minutes
with the earlier sampler (about 0.23 seconds per iteration on one core of an Intel Xeon
Gold 6254), so the full study needs a cluster. With SLURM:

```sh
cd simulations
sbatch submit_sim.sh             # 420 array tasks, 20 replicates each
Rscript sim_run.R aggregate      # per-replicate files -> per-cell results
Rscript sim_report.R
```

Tasks write one file per replicate and skip replicates already done, so a task that
times out can be resubmitted without losing work. To run a single task locally, for
example task 0 on 4 cores:

```sh
Rscript sim_run.R task 0 4 20
```

Arguments after the task id are the number of cores, the chunk size, and optionally the
number of iterations, burn-in and whether to run the competitors (1 or 0), which is
useful for a quick test:

```sh
SIM_RESULTS=test_results Rscript sim_run.R task 0 4 20 200 100 1
SIM_RESULTS=test_results Rscript sim_run.R aggregate
```

## MEC-APS application (Section 5)

### Data

Set the paths in `mec/config.R`, or through the environment variables named there. Step 1
expects `MEC_assemblage_9_2024_asv_table.txt` and `MEC_assemblage_9_2024_asv_taxonomy.csv`
in `raw_asv_dir`. Step 2 also expects `MEC_assemblage_9_2024_sample_info.csv` in
`raw_asv_dir`, and `mec-ids-data.csv`, `mediation_v3.csv` and
`03-Clean/mec-{bugs,ost,ocl}-data.csv` in `raw_meta_dir`. Steps 1 and 2 write
`MEC_asv_data.rds` and `MEC_genus_mediation.rds` to `data_dir`.

### Steps

Run every script from `mec/`.

```sh
cd mec
Rscript 01_read_asv_data.R
Rscript 02_build_mediation_inputs.R
```

Step 2 keeps one valid stool sample per participant and complete cases on the
covariates, mediator and outcome. It aggregates ASVs to genera and keeps genera with
relative abundance of at least 0.01% in at least 10% of samples, which leaves 135
genera. It then replaces zeros by 0.5 and log-transforms percent liver fat.

Step 3 fits one model setting. Without arguments it fits the primary model: eta = 1/20,
adjusted for diabetes, metformin users retained, 4 chains of 5000 iterations with 1000
burn-in.

```sh
Rscript 03_fit_mec.R                  # primary model, 4 chains in parallel
```

One chain takes about half an hour on a single core (roughly 0.4 seconds per iteration
on a recent desktop processor), so the eight fits in Figure 6 are best run as a SLURM
array of 32 single-chain tasks, or 32 parallel processes on a multicore machine:

```sh
sbatch submit_mec.sh
bash submit_mec.sh aggregate          # after all tasks finish
```

Then build the figure, the table and the numbers:

```sh
Rscript 04_figure6.R
Rscript 05_table2.R
Rscript 06_section5_numbers.R
Rscript 07_convergence_figures.R
```

`03_fit_mec.R` documents all arguments. Chain c is seeded with `set.seed(1000 + c)`, and
chain 1 starts from the leading principal component of the centered log-ratio transformed
genera.

## License

GNU General Public License, version 3. See `LICENSE`.
