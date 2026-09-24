## =============================================================================
## 01_read_asv_data.R -- read the MEC-APS ASV count table and its SILVA 138
## taxonomy (Section 5). Run from mec/.
##
## Input  (raw_asv_dir in config.R):
##   MEC_assemblage_9_2024_asv_table.txt     ASV x sample count table
##   MEC_assemblage_9_2024_asv_taxonomy.csv  one row per ASV: Feature ID, Taxon, Confidence
## Output (data_dir): MEC_asv_data.rds, a list with
##   asv_counts  integer matrix, rows = ASVs, columns = samples
##   taxonomy    data frame with columns asv, taxon, confidence, domain, ..., species
## =============================================================================

source("config.R")
library(data.table)

asv_file <- file.path(raw_asv_dir, "MEC_assemblage_9_2024_asv_table.txt")
tax_file <- file.path(raw_asv_dir, "MEC_assemblage_9_2024_asv_taxonomy.csv")

## ---- ASV count table -------------------------------------------------
## Line 1 is a comment ("# Constructed from biom file"); line 2 is the
## header, whose first field is "#OTU ID" and which ends in a trailing tab.
## Counts are stored as floats ("4977.0") but are integers in value.
asv_dt <- fread(asv_file, sep = "\t", skip = 1, header = TRUE)
if (names(asv_dt)[ncol(asv_dt)] == "") asv_dt[, (ncol(asv_dt)) := NULL]  # drop empty trailing column
asv_ids <- asv_dt[[1]]
asv_counts <- as.matrix(asv_dt[, -1])
storage.mode(asv_counts) <- "integer"
rownames(asv_counts) <- asv_ids
stopifnot(!anyDuplicated(rownames(asv_counts)), !anyDuplicated(colnames(asv_counts)))

## ---- Taxonomy --------------------------------------------------------
## Columns: Feature ID, Taxon, Confidence. Taxon is a semicolon-separated
## string with rank prefixes, e.g. "d__Bacteria; p__Bacteroidota; ...".
taxonomy <- fread(tax_file)
setnames(taxonomy, c("asv", "taxon", "confidence"))
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
tax_split <- tstrsplit(taxonomy$taxon, ";\\s*", fixed = FALSE)
length(tax_split) <- length(ranks)                  # pad with NULL if fewer than 7 ranks
tax_split <- lapply(tax_split, function(x) if (is.null(x)) NA_character_ else x)
names(tax_split) <- ranks
taxonomy[, (ranks) := lapply(tax_split, function(x) sub("^[a-z]__", "", x))]
taxonomy[taxonomy == ""] <- NA                      # empty rank -> NA
taxonomy <- as.data.frame(taxonomy)

## ---- Align and check --------------------------------------------------
stopifnot(setequal(taxonomy$asv, rownames(asv_counts)))
taxonomy <- taxonomy[match(rownames(asv_counts), taxonomy$asv), ]
rownames(taxonomy) <- NULL

cat("ASVs:", nrow(asv_counts), " samples:", ncol(asv_counts), "\n")
cat("Total reads:", sum(asv_counts), "\n")
print(table(taxonomy$phylum, useNA = "ifany"))

dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(list(asv_counts = asv_counts, taxonomy = taxonomy),
        file = file.path(data_dir, "MEC_asv_data.rds"))
