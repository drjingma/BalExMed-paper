## =============================================================================
## 05_table2.R -- Table 2: genera selected (PIP > 0.5) in the primary MEC-APS fit,
## with their signed inclusion P(z = 1) - P(z = -1) and the number of the eight
## Figure 6 fits in which they are selected. Run from mec/ after all eight fits.
##   Usage: Rscript 05_table2.R [PIP_THRESHOLD]     (default 0.5)
## Output: <results_dir>/table2.csv  one row per selected genus: signed inclusion,
##                                   PIP, PIP in each fit, number of fits selected,
##                                   family, genus and number of ASVs
##         <results_dir>/table2.tex  tabular body in the layout of the paper
## Labels: underscores become spaces, and genus placeholders that are codes
## (containing digits) are prefixed with their SILVA family when the family is not
## already part of the name, e.g. DTU089 -> Ruminococcaceae DTU089.
## =============================================================================
source("config.R")
a <- commandArgs(trailingOnly = TRUE)
thr <- if (length(a) >= 1) as.numeric(a[1]) else 0.5

pri_dir <- file.path(results_dir, "primary")
sen_dir <- file.path(results_dir, "sensitivity")
pri_file <- function(a, e)
  sprintf("%s/pfat_ethall_metfRetained_adj%s_eta1over%d_gibbs_4chains_iter5000.rds", pri_dir, a, e)
## same eight fits, same order, as mec_fig6_rev.R; the first row is the primary
spec <- data.frame(
  file = c(pri_file("diabet", 20), pri_file("diabet", 3), pri_file("diabet", 5),
           pri_file("metformin", 20), pri_file("none", 20),
           sprintf("%s/pfat_ethall_metfRetained_adjdiabet_covnoadip_eta1over20_gibbs_4chains_iter5000.rds", sen_dir),
           sprintf("%s/pfat_ethall_noMetf_adjdiabet_eta1over20_gibbs_4chains_iter5000.rds", sen_dir),
           sprintf("%s/pfat_ethall_noT2D_adjnone_eta1over20_gibbs_4chains_iter5000.rds", sen_dir)),
  key  = c("P", "A1_eta3", "A2_eta5", "A3_metformin", "A4_neither",
           "S2_noadip", "S3_nometf", "S4_not2d"),
  stringsAsFactors = FALSE)
miss <- spec$file[!file.exists(spec$file)]
if (length(miss)) stop("missing fit(s):\n  ", paste(miss, collapse = "\n  "))

taxa <- lapply(spec$file, function(f) readRDS(f)$taxa)
tx0  <- taxa[[1]]                                     # primary
## PIP of every taxon in every fit (taxa x fits), matched by name
PIP <- sapply(taxa, function(t) t$P_incl[match(tx0$taxon, t$taxon)])
colnames(PIP) <- spec$key
stopifnot(!anyNA(PIP))                                # same feature set in all fits
n_fits <- rowSums(PIP > thr)

sel <- tx0$P_incl > thr
out <- data.frame(taxon  = tx0$taxon,
                  signed = round(tx0$P_plus - tx0$P_minus, 2),
                  PIP    = round(tx0$P_incl, 3),
                  n_fits = n_fits,
                  round(PIP, 3), check.names = FALSE, stringsAsFactors = FALSE)
for (v in c("family", "genus", "n_asv")) if (v %in% names(tx0)) out[[v]] <- tx0[[v]]
out <- out[sel, ]
## numerator block first (signed decreasing), then denominator (most negative first)
num <- out[out$signed > 0, ]; num <- num[order(-num$signed), ]
den <- out[out$signed < 0, ]; den <- den[order( den$signed), ]
out <- rbind(num, den)

cat(sprintf("primary: %d taxa at PIP > %g (%d numerator, %d denominator); %d selected in all %d fits\n",
            nrow(out), thr, nrow(num), nrow(den), sum(out$n_fits == nrow(spec)), nrow(spec)))
print(out[, c("taxon", "signed", "PIP", "n_fits")], row.names = FALSE)

csv <- file.path(results_dir, "table2.csv")
write.csv(out, csv, row.names = FALSE)

## LaTeX rows: label & $+0.99$ & 8 \\
label_of <- function(d) {
  lab <- gsub("_", " ", d$taxon)
  fam <- if ("family" %in% names(d)) gsub("_", " ", d$family) else rep(NA_character_, nrow(d))
  code <- grepl("[0-9]", lab) & !is.na(fam) & !mapply(grepl, fam, lab, fixed = TRUE)
  lab[code] <- paste(fam[code], lab[code])
  lab <- sub("^Candidatus ", "\\\\emph{Candidatus} ", lab)
  gsub("([&%#])", "\\\\\\1", lab)
}
row <- function(d) sprintf("%s & $%s%.2f$ & %d \\\\", label_of(d),
                           ifelse(d$signed > 0, "+", "-"), abs(d$signed), d$n_fits)
tex <- c("\\begin{tabular}{lcc}\\toprule",
         sprintf("Taxon & Signed inclusion & Models (of %d) \\\\ \\midrule", nrow(spec)),
         row(num), "\\midrule", row(den), "\\bottomrule\\end{tabular}")
texf <- file.path(results_dir, "table2.tex")
writeLines(tex, texf)
cat("wrote", csv, "and", texf, "\n")
