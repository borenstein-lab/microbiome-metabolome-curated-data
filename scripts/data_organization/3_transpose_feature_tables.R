# -----------------------------------------------------------------
# Transpose feature tables, from row-per-feature to row-per-sample.
# Override relevant files (tsv + RData) with the updated versions.
# -----------------------------------------------------------------

require(tibble)
require(dplyr)
source("data_organization/utils.R")

all.data <- load.all.datasets("processed_data")
for(i in 1:length(all.data)) assign(names(all.data)[i], all.data[[i]])
rm(all.data)
datasets <- basename(data.dirs)

curr_genera <- genera; rm(genera)
curr_species <- species; rm(species)
curr_genera.counts <- genera.counts; rm(genera.counts)
curr_species.counts <- species.counts; rm(species.counts)
curr_mtb <- mtb; rm(mtb)
has_species <- sapply(datasets, function(x) {x %in% names(curr_species)})

# Transpose
transpose.feat.table <- function(df, new.id.col.name = "Sample") {
  # Check if already transposed
  if (colnames(df)[1] == new.id.col.name) {
    message("Already transposed, leaving untouched")
    return(df)
  }
  new.df <- df %>% remove_rownames() %>% column_to_rownames(colnames(df)[1]) 
  new.df <- data.frame(t(new.df), check.names = FALSE)
  new.df <- new.df %>% rownames_to_column(new.id.col.name)
  return(new.df)
}

for(dataset in data.dirs) {
  message(paste("Transposing tables of:", basename(dataset)))
  
  genera <- transpose.feat.table(curr_genera[[basename(dataset)]])
  genera.counts <- transpose.feat.table(curr_genera.counts[[basename(dataset)]])
  mtb <- transpose.feat.table(curr_mtb[[basename(dataset)]])
  
  species = NULL
  species.counts = NULL
  if (has_species[basename(dataset)]) {
    species <- transpose.feat.table(curr_species[[basename(dataset)]])
    species.counts <- transpose.feat.table(curr_species.counts[[basename(dataset)]])
  } 
  
  save.to.files(basename(dataset), "processed_data", genera = genera, genera.counts = genera.counts, mtb = mtb, species = species, species.counts = species.counts)
  save.to.rdata(basename(dataset), "processed_data", genera = genera, genera.counts = genera.counts, mtb = mtb, species = species, species.counts = species.counts)
  message("Saved transposed genera, [species,] and mtb tables")
}
