# ----------------------------------------------------------------
# Unify species names in species-feature-tables.
# Override relevant files (tsv + RData) with the unified version.
# ----------------------------------------------------------------

require(vegan)
require(dplyr)
require(cgwtools)
source("utils.R")
all.data <- load.all.datasets("processed_data")
for(i in 1:length(all.data)) assign(names(all.data)[i], all.data[[i]])
rm(all.data)
datasets <- names(species)

# --------------------------------
# 1. Remove non-bacteria
# --------------------------------

# We remove non-bacteria entities 
#  (before re-normalizing each sample to 100%)
# Note: this has already been performed in some datasets, 
#  by authors/ depending on exact metagenomics processing.

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

# --------------------------------
# 3. Unify species names 
# --------------------------------

# 3.1. Format changes

# x <- c(); for(dataset in datasets) {x <- c(x, species.new[[dataset]]$Species)}; x <- unique(x)

for (dataset in datasets) {
  new.names <- species.new[[dataset]]$Species
  
  new.names <- gsub("\\|p__",";p__",new.names)
  new.names <- gsub("\\|c__",";c__",new.names)
  new.names <- gsub("\\|o__",";o__",new.names)
  new.names <- gsub("\\|f__",";f__",new.names)
  new.names <- gsub("\\|g__",";g__",new.names)
  new.names <- gsub("\\|s__",";s__",new.names)
  new.names <- gsub("\\|__",";__",new.names)
  
  new.names <- gsub(";p__;",";__;",new.names)
  new.names <- gsub(";c__;",";__;",new.names)
  new.names <- gsub(";o__;",";__;",new.names)
  new.names <- gsub(";f__;",";__;",new.names)
  new.names <- gsub(";g__;",";__;",new.names)
  new.names <- gsub(";s__$",";__",new.names)
  
  new.names <- gsub("[ocfg]__[0-9a-zA-Z_]*_noname;","__;",new.names)
  
  message(paste(dataset, "- reformatted",
                sum(new.names != species.new[[dataset]]$Species),
                "out of", length(new.names), "species entities"))
  
  species.new[[dataset]]$Species <- new.names
  # print(paste("After:",new.names[1]))
}
rm(new.names)

# --------------------------------
# 4. Save
# --------------------------------

# Override RData files and "genera" text tables
source("load_original_data/utils.R")

for (dataset in data.dirs[basename(data.dirs) %in% datasets]) {
  species <- species.new[[basename(dataset)]]
  save.to.files(basename(dataset), "processed_data", species = species)
  save.to.rdata(basename(dataset), "processed_data", species = species)
}


