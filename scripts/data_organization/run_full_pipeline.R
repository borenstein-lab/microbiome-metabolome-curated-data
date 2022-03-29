# Load datasets one by one and process into a common format
# Files will be saved at: data/prelim_data
source("data_organization/1_load_all.R")

# Process taxonomic data in a uniform manner
# Files will be saved at: data/processed_data
source("data_organization/2A_unify_genera.R")
source("data_organization/2B_unify_species.R")

# Reformat (transpose) tables into row-per-sample and column-per-feature
# Files will be overwritten at: data/processed_data
source("data_organization/3_transpose_feature_tables.R")

# Files will be overwritten at: data/processed_data
source("data_organization/4_zip_large_files.R")
