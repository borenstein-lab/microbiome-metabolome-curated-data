# ----------------------------------------------------------------
# Unify species names in species-feature-tables.
# Override relevant files (tsv + RData) with the unified version.
# ----------------------------------------------------------------

require(vegan)
require(dplyr)
require(cgwtools)
source("data_organization/utils.R")

all.data <- load.all.datasets("processed_data")
for(i in 1:length(all.data)) assign(names(all.data)[i], all.data[[i]])
rm(all.data)
datasets <- names(species)

# --------------------------------
# 1. Remove non-bacteria
# --------------------------------

# We remove non-bacteria entities 
#  (before re-normalizing each sample to 100%)

for (dataset in datasets) {
  tmp <- species[[dataset]]
  tmp.n.row <- nrow(tmp) # Recording number before we erase all non-bacteria
  tmp <- tmp[!grepl("Viruses", tmp$Species, ignore.case = T),]
  tmp <- tmp[!grepl("Archaea", tmp$Species, ignore.case = T),]
  tmp <- tmp[!grepl("Eukaryota",tmp$Species, ignore.case = T),]
  message(paste("Removed",tmp.n.row-nrow(tmp),"non-bacteria rows from dataset",dataset))
  species[[dataset]] <- tmp
}
rm(tmp.n.row)

# --------------------------------
# 2. Transform to relative abund'
# --------------------------------

species.new <- list()
for (dataset in datasets) {
  tmp <- decostand(species[[dataset]][,-1], method = "total", MARGIN = 2)
  tmp$Species <- species[[dataset]]$Species
  tmp <- tmp %>% relocate(Species)
  species.new[[dataset]] <- tmp
  # Sanity: print(apply(tmp[,-1], 2, sum))
}

# Code for a quick indication of whether genus names are consistent:
# x <- c(); for(dataset in datasets) {x <- c(x, species.new[[dataset]]$Species)}; 
# x <- data.frame(sp = x) %>% group_by(sp) %>% summarise(N = n()) %>% tidyr::separate(col = "sp", into = c("d","p","c","o","f","g","s"), sep = ";", remove = FALSE)

# --------------------------------
# 3. Save
# --------------------------------

# Override RData files and "genera" text tables

for (dataset in data.dirs[basename(data.dirs) %in% datasets]) {
  species <- species.new[[basename(dataset)]]
  save.to.files(basename(dataset), "processed_data", species = species)
  save.to.rdata(basename(dataset), "processed_data", species = species)
}


