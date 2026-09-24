## =============================================================================
## 02_build_mediation_inputs.R -- build the analysis object for Section 5 from
## the ASV counts and the MEC-APS phenotype files. Run from mec/.
##
## Steps
##   1. keep one valid stool sample per participant;
##   2. merge phenotypes: education recoded to three levels, Yes/No to 1/0,
##      season of stool collection, metformin use;
##   3. keep complete cases on the covariates, mediator and outcome;
##   4. aggregate ASVs to genera using the SILVA 138 labels (ASVs without a named
##      genus are labelled by their lowest named rank), keep genera with relative
##      abundance >= 0.01% in >= 10% of samples, and replace zeros by 0.5;
##   5. log-transform the outcome and save.
##
## Inputs
##   data_dir/MEC_asv_data.rds                          from 01_read_asv_data.R
##   raw_asv_dir/MEC_assemblage_9_2024_sample_info.csv
##   raw_meta_dir/mec-ids-data.csv, mediation_v3.csv,
##   raw_meta_dir/03-Clean/mec-bugs-data.csv, mec-ost-data.csv, mec-ocl-data.csv
## Output
##   data_dir/MEC_genus_mediation.rds: list(exposures, mediator, response,
##                                          covariates, metformin, taxonomy)
## =============================================================================

source("config.R")
suppressMessages({ library(data.table); library(dplyr) })

asv_dir  <- raw_asv_dir
meta_dir <- raw_meta_dir

## ---- 1. ASV counts and one valid sample per participant ---------------
asv <- readRDS(file.path(data_dir, "MEC_asv_data.rds"))
counts <- asv$asv_counts                         # ASV x sample

sample_info <- fread(file.path(asv_dir, "MEC_assemblage_9_2024_sample_info.csv"))
mec_ids <- fread(file.path(meta_dir, "mec-ids-data.csv"))[sample_def == "Valid"]
stopifnot(!anyDuplicated(mec_ids$p01_id))

## Valid stool_id is the ASV-table column name. Cross-check the participant
## ID against sample_info where both are present.
mec_ids <- mec_ids[stool_id %in% colnames(counts)]
chk <- merge(mec_ids[, .(stool_id, p01_id)],
             sample_info[, .(stool_id = sample_name, P01ID)], by = "stool_id")
stopifnot(all(chk$p01_id == chk$P01ID, na.rm = TRUE))

counts <- counts[, mec_ids$stool_id]
colnames(counts) <- mec_ids$p01_id
cat("Participants with a valid ASV sample:", ncol(counts), "\n")

## ---- 2. Phenotypes ------------------------------------------------------
fermented_vars <- c("OST_diet_kimchi", "OST_diet_miso", "OST_diet_natto",
                    "OST_diet_pickles", "OST_diet_tempeh")
factor_vars <- c("OST_diet_veget", "OQ3_diabet", "OMRI_NAFLD_55")
eth_codes <- c("African American" = "B", "Native Hawaiian" = "H",
               "Japanese American" = "J", "Latino" = "L", "White" = "W")

bugs <- read.csv(file.path(meta_dir, "03-Clean/mec-bugs-data.csv"))
bugs$asian_fermented_any <- as.integer(
  rowSums(bugs[, fermented_vars] == "Yes", na.rm = TRUE) > 0)

bugs <- bugs %>%
  mutate(edu_cat = case_when(
    Q1_YRSSCHL %in% c("Did not complete 6th Grade", "6th-8th Grade",
                      "9th-10th Grade", "11th-12th Grade") ~ "≤HighSchool",
    Q1_YRSSCHL %in% c("Vocational School", "Some College") ~ "SomeCollege",
    Q1_YRSSCHL %in% c("Graduated College",
                      "Graduate or Professional School") ~ "College+"),
    edu_cat = factor(edu_cat, levels = c("≤HighSchool", "SomeCollege", "College+")),
    Q1_eth = unname(eth_codes[Q1_eth]),           # one-letter ethnicity codes
    asian_fermented_any = factor(asian_fermented_any)) %>%
  mutate(across(where(is.character), as.factor)) %>%
  mutate(across(where(~ is.factor(.) &&
                        all(levels(.) %in% c("Yes", "No", "YES", "NO", "yes", "no"))),
                ~ ifelse(tolower(as.character(.)) == "yes", 1, 0))) %>%
  mutate(across(all_of(factor_vars), as.factor)) %>%
  select(-all_of(fermented_vars)) %>%
  inner_join(as.data.frame(mec_ids)[, c("p_id", "p01_id")], by = "p_id")
stopifnot(!anyDuplicated(bugs$p01_id))

season_of <- function(d) {
  m <- as.integer(format(d, "%m"))
  factor(c("Winter", "Spring", "Summer", "Fall")[
    ifelse(m %in% c(12, 1, 2), 1, ifelse(m <= 5, 2, ifelse(m <= 8, 3, 4)))],
    levels = c("Winter", "Spring", "Summer", "Fall"))
}
ost <- read.csv(file.path(meta_dir, "03-Clean/mec-ost-data.csv")) %>%
  transmute(p01_id, OST_sample_date = as.Date(OST_sample_date),
            season = season_of(OST_sample_date))
ocl <- read.csv(file.path(meta_dir, "03-Clean/mec-ocl-data.csv")) %>%
  transmute(p01_id, OCL_meds_diabetes_metformin =
              ifelse(OCL_meds_diabetes_metformin >= 1, 1, 0))
masld <- read.csv(file.path(meta_dir, "mediation_v3.csv"), check.names = FALSE) %>%
  transmute(p01_id = P01ID, OD_MASLD = factor(OD_MASLD)) %>%
  distinct(p01_id, .keep_all = TRUE)

meta <- bugs %>%
  left_join(ost, by = "p01_id") %>%
  left_join(ocl, by = "p01_id") %>%
  left_join(masld, by = "p01_id") %>%
  mutate(Q1_eth = as.factor(Q1_eth)) %>%
  filter(p01_id %in% colnames(counts))
cat("Participants with phenotypes and ASV sample:", nrow(meta), "\n")

## ---- 3. Sample exclusions --------------------------------------------------
vars_NA <- c("Q1_CORR_SEX", "Q1_eth", "Q1_POB", "edu_cat", "OQ3_diabet",
             "OST_sample_age", "season", "OCL_anthro_BMI",
             "OQ3_DP_AHEI2010_TOTSCORE", "ODXA_pfat_tot_corr",
             "OQ3_ethanol", "OQ3_packyrs", "OQ3_ac_hrsit",
             "O_LBP", "OMRI_pct_liver_fat_corr")
meta_clean <- meta[rowSums(is.na(meta[, vars_NA])) == 0, ]
cat("Complete cases on covariates/mediator/outcomes:", nrow(meta_clean),
    " (of which OD_MASLD missing:", sum(is.na(meta_clean$OD_MASLD)), ")\n")

## No trimming on alcohol or pack-years; pack-years enters as a covariate.
meta_filt <- meta_clean

## ---- 4. Aggregate ASVs to genus, then filter ---------------------------
## Aggregation key (no taxonomy is inferred; labels come from SILVA 138):
##   * a named genus                    -> the genus name
##   * genus "uncultured"/"Incertae_Sedis" -> "<lowest named rank> (<rank>; <genus>)",
##     because these SILVA labels recur under many unrelated families
##   * no genus assigned                -> "<lowest named rank> (<rank>)"
## "Named" excludes NA and "uncultured".
ranks <- c("domain", "phylum", "class", "order", "family", "genus")
placeholder <- c("uncultured", "Incertae_Sedis")
genus_key <- function(tx) {
  apply(tx[, ranks], 1, function(r) {
    if (!is.na(r["genus"]) && !r["genus"] %in% placeholder) return(unname(r["genus"]))
    named <- which(!is.na(r[1:5]) & r[1:5] != "uncultured")
    i <- max(named)
    if (is.na(r["genus"])) sprintf("%s (%s)", r[i], ranks[i])
    else sprintf("%s (%s; %s)", r[i], ranks[i], r["genus"])
  })
}
tx <- asv$taxonomy
tx$key <- genus_key(tx)
stopifnot(identical(tx$asv, rownames(counts)))

G <- rowsum(counts[, as.character(meta_filt$p01_id)], tx$key)  # genus x sample
G <- t(G)                                                       # sample x genus
cat("Genus-level features:", ncol(G), " (from", nrow(tx), "ASVs;",
    sum(is.na(tx$genus)), "ASVs without genus)\n")

min_rel  <- 1e-4     # 0.01% relative abundance
min_prev <- 0.10     # in at least 10% of samples
rel <- G / rowSums(G)
keep <- colMeans(rel >= min_rel) >= min_prev
X <- G[, keep]
cat("Genera retained:", ncol(X), "of", length(keep),
    "; share of reads retained:", round(sum(X) / sum(G), 4), "\n")
cat("Fraction of retained counts that are zero:", round(mean(X == 0), 3), "\n")
X_filt <- X
X_filt[X_filt == 0] <- 0.5

## genus-level taxonomy: ranks shared by all member ASVs, else NA
taxonomy <- do.call(rbind, lapply(colnames(X_filt), function(k) {
  m <- tx[tx$key == k, ranks]
  data.frame(feature = k, n_asv = nrow(m),
             lapply(m, function(v) if (length(unique(v)) == 1) v[1] else NA_character_),
             stringsAsFactors = FALSE)
}))
rownames(taxonomy) <- NULL

## ---- 5. Assemble ----------------------------------------------------------
vars.y <- "OMRI_pct_liver_fat_corr"
vars <- c("OST_sample_age", "OST_sample_date", "season", "Q1_CORR_SEX",
          "Q1_eth", "Q1_POB", "OCL_anthro_BMI", "OQ3_DP_AHEI2010_TOTSCORE",
          "ODXA_pfat_tot_corr", "OD_MASLD", 
          "OQ3_ethanol", "OQ3_packyrs", "OQ3_diabet", "edu_cat", "OQ3_ac_hrsit")

data.list <- list(
  exposures  = X_filt,
  mediator   = meta_filt$O_LBP,
  response   = as_tibble(meta_filt[, vars.y, drop = FALSE]) %>% mutate(across(everything(), log)),
  covariates = as_tibble(meta_filt[, vars]) %>% droplevels(),
  metformin  = meta_filt$OCL_meds_diabetes_metformin,
  taxonomy   = taxonomy)
stopifnot(identical(rownames(data.list$exposures), as.character(meta_filt$p01_id)))

cat("\nn =", nrow(data.list$response), " genera =", ncol(data.list$exposures),
    " metformin users =", sum(data.list$metformin == 1, na.rm = TRUE),
    " missing metformin =", sum(is.na(data.list$metformin)), "\n")
print(table(data.list$covariates$Q1_eth))
saveRDS(data.list, file = file.path(data_dir, "MEC_genus_mediation.rds"))
cat("saved", file.path(data_dir, "MEC_genus_mediation.rds"), "\n")
