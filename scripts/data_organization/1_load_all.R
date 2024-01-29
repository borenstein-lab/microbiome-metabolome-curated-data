# ----------------------------------------------------------------
# Run all dataset-specific load_data scripts, found in the
#  "load_original_data" directory.
# After completing, each dataset should have a folder with all 
#  processed data in "data/prelim_data" directory.
# ----------------------------------------------------------------

library(gdata)
library(readr)
library(dplyr)
library(stringr)
library(MetaboAnalystR)
source("data_organization/utils.R")

load.data.scripts <- list.files("data_organization/load_original_data")
load.data.scripts <- grep("^load_data_.*\\.R", load.data.scripts, value = TRUE)
load.data.scripts <- file.path(getwd(),"data_organization/load_original_data",load.data.scripts)

tmp <- sapply(load.data.scripts,
       FUN = function(x) {
         message(paste("Running script:", basename(x)))
         tmp.env <- new.env()
         source(x, echo = FALSE, local = tmp.env)
         rm(tmp.env)
         message(paste("Completed script:", basename(x)))
       })

message("Completed data loading")

# --------------------------------
# Patches (involving >1 dataset)
# --------------------------------

# Mark duplicated samples in datasets ERAWIJANTARI_GASTRIC_CANCER_2020 & YACHIDA_CRC_2019
metadata.yach <- read_delim("../data/prelim_data/YACHIDA_CRC_2019/metadata.tsv", 
                            delim = "\t", escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE)
metadata.eraw <- read_delim("../data/prelim_data/ERAWIJANTARI_GASTRIC_CANCER_2020/metadata.tsv", 
                            delim = "\t", escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE)
shared.subjects <- intersect(gsub("\\..*","",metadata.eraw$Sample),
                             as.character(metadata.yach$Sample))
metadata.yach$Shared.w.ERAWIJANTARI_2020 <- as.character(metadata.yach$Sample) %in% shared.subjects
metadata.eraw$Shared.w.YACHIDA_2019 <- gsub("\\..*","",metadata.eraw$Sample) %in% shared.subjects

save.to.files("YACHIDA_CRC_2019", "prelim_data", metadata = metadata.yach)
save.to.rdata("YACHIDA_CRC_2019", "prelim_data", metadata = metadata.yach)
save.to.files("ERAWIJANTARI_GASTRIC_CANCER_2020", "prelim_data", metadata = metadata.eraw)
save.to.rdata("ERAWIJANTARI_GASTRIC_CANCER_2020", "prelim_data", metadata = metadata.eraw)

# Remove MetaboAnalyst unneeded files
for (f in c("master_compound_db.qs", "master_syn_nms.qs", "name_map.csv")) {
  if (file.exists(f)) file.remove(f)
}
