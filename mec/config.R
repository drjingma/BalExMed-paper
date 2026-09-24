## Paths used by the MEC-APS scripts (run every script from mec/).
## The MEC-APS data are not publicly available; see the README for the files
## each step expects. Edit the defaults below or set the environment variables.
raw_asv_dir  <- Sys.getenv("MEC_RAW_ASV_DIR",  "path/to/MEC-ASV-based_microbiome")   # ASV table and taxonomy
raw_meta_dir <- Sys.getenv("MEC_RAW_META_DIR", "path/to/MEC-Master-Clean/Data")      # phenotype files
data_dir     <- Sys.getenv("MEC_DATA_DIR",     "data")      # derived analysis objects (steps 1-2)
results_dir  <- Sys.getenv("MEC_RESULTS_DIR",  "results")   # model fits, figure, tables (steps 3-6)
