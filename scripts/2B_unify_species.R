# ----------------------------------------------------------------
# Unify species names in species-feature-tables.
# Override relevant files (tsv + RData) with the unified version.
# ----------------------------------------------------------------

require(vegan)
require(dplyr)
require(cgwtools)
source("load_all_datasets.R")
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

for (dataset in datasets) {
  new.names <- species.new[[dataset]]$Species
  
  new.names <- gsub(";p__","\\|p__",new.names)
  new.names <- gsub(";c__","\\|c__",new.names)
  new.names <- gsub(";o__","\\|o__",new.names)
  new.names <- gsub(";f__","\\|f__",new.names)
  new.names <- gsub(";g__","\\|g__",new.names)
  new.names <- gsub(";s__","\\|s__",new.names)
  new.names <- gsub(";__","\\|__",new.names)
  
  new.names <- gsub("\\|p__\\|","\\|__\\|",new.names)
  new.names <- gsub("\\|c__\\|","\\|__\\|",new.names)
  new.names <- gsub("\\|o__\\|","\\|__\\|",new.names)
  new.names <- gsub("\\|f__\\|","\\|__\\|",new.names)
  new.names <- gsub("\\|g__\\|","\\|__\\|",new.names)
  new.names <- gsub("\\|s__$","\\|__",new.names)
  
  new.names <- gsub("[ocfg]__Bacteroidetes_noname\\|","__\\|",new.names)
  new.names <- gsub("[ocfg]__Bacillales_noname\\|","__\\|",new.names)
  new.names <- gsub("[ocfg]__Clostridiales_noname\\|","__\\|",new.names)
  new.names <- gsub("[ocfg]__Lachnospiraceae_noname\\|","__\\|",new.names)
  
  message(paste(dataset, "- reformatted",
                sum(new.names != species.new[[dataset]]$Species),
                "out of", length(new.names), "species entities"))
  
  species.new[[dataset]]$Species <- new.names
  # print(paste("After:",new.names[1]))
}
rm(new.names)

# 3.2. Fix reference database discrepancies

# First, extract only species names 
x <- c(); for(dataset in datasets) {x <- c(x, species.new[[dataset]]$Species)}; x <- unique(x)
sp.names.mapping <- data.frame(Sp = x)
sp.names.mapping <- sp.names.mapping %>%
  mutate(Sp.short = gsub(".*\\|s__", "", Sp))

duplicated.taxa <- sp.names.mapping$Sp.short[duplicated(sp.names.mapping$Sp.short)]

message(paste("A total of",nrow(sp.names.mapping),
              "unique species strings (full taxonomy) were found"))
message(paste(length(unique(sp.names.mapping$Sp.short)),
              "of which are actually unique taxa, and",
              length(duplicated.taxa),"are duplicates"))

## Create a mapping between current names and unified names
# First we create a table with mapping info
sp.names.mapping <- sp.names.mapping %>% 
  group_by(Sp.short) %>% 
  dplyr::mutate(Sp.new = first(Sp)) %>%
  ungroup()

# We now turn the mapping into a named vector
map.sp.names <- sp.names.mapping$Sp.new
names(map.sp.names) <- sp.names.mapping$Sp

# Finally, use this mapping vector to map current species names into unified ones
for (dataset in datasets) {
  tmp <- unname(map.sp.names[species.new[[dataset]]$Species])
  species.new[[dataset]]$Species <- tmp
  message(paste("Cleaned",dataset))
}
rm(sp.names.mapping, duplicated.taxa, map.sp.names, dataset, tmp)

# --------------------------------
# 4. Save
# --------------------------------

# Override RData files and "genera" text tables
source("load_original_data/utils.R")

for (dataset in data.dirs[basename(data.dirs) %in% datasets]) {
  species <- species.new[[basename(dataset)]]
  resave(species, file = file.path(dataset, ".RData"))
  save.to.files(basename(dataset), species = species)
}


