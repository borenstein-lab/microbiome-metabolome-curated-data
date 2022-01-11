# -----------------------------------------------------------------
# Transpose feature tables, from row-per-feature to row-per-sample.
# Override relevant files (tsv + RData) with the updated versions.
# -----------------------------------------------------------------

require(tibble)
require(dplyr)
require(cgwtools)
source("load_all_datasets.R")
source("load_original_data/utils.R")
datasets <- basename(data.dirs)

curr_genera <- genera; rm(genera)
curr_species <- species; rm(species)
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
  mtb <- transpose.feat.table(curr_mtb[[basename(dataset)]])
  
  if (has_species[basename(dataset)]) {
    species <- transpose.feat.table(curr_species[[basename(dataset)]])
    resave(genera, mtb, species, file = file.path(dataset, ".RData"))
    save.to.files(basename(dataset), genera = genera, mtb = mtb, species = species)
    message("Saved transposed genera, species and mtb tables")
  } else {
    resave(genera, mtb, file = file.path(dataset, ".RData"))
    save.to.files(basename(dataset), genera = genera, mtb = mtb)
    message("Saved transposed genera and mtb tables")
  }
}
