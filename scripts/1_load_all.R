# ----------------------------------------------------------------
# Run all dataset-specific load_data scripts, found in the
#  "load_original_data" directory.
# After completing, each dataset should have a folder with all 
#  processed data in "data/processed_data" directory.
# ----------------------------------------------------------------

load.data.scripts <- list.files("load_original_data")
load.data.scripts <- grep("^load_data_.*\\.R", load.data.scripts, value = TRUE)
load.data.scripts <- file.path(getwd(),"load_original_data",load.data.scripts)

tmp <- sapply(load.data.scripts,
       FUN = function(x) {
         message(paste("Running script:", basename(x)))
         tmp.env <- new.env()
         source(x, echo = FALSE, local = tmp.env)
         rm(tmp.env)
         message(paste("Completed script:", basename(x)))
       })

message("Completed data loading")
