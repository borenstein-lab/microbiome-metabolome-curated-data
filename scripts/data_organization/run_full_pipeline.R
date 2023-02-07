# Load datasets one by one and process into a common format
# Files will be saved at: data/prelim_data
source("data_organization/1_load_all.R")
message('****************************************\n* Completed step 1 of 4\n****************************************')

# Process taxonomic data in a uniform manner
# Files will be saved at: data/processed_data
source("data_organization/2A_unify_genera.R")
source("data_organization/2B_unify_species.R")
message('****************************************\n* Completed step 2 of 4\n****************************************')

# Reformat (transpose) tables into row-per-sample and column-per-feature
# Files will be overwritten at: data/processed_data
source("data_organization/3_transpose_feature_tables.R")
message('****************************************\n* Completed step 3 of 4\n****************************************')

# Files will be overwritten at: data/processed_data
source("data_organization/4_zip_large_files.R")
message('****************************************\n* Completed step 4 of 4\n****************************************')
