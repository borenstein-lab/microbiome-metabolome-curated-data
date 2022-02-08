# -----------------------------------------------------------------
# Transpose feature tables, from row-per-feature to row-per-sample.
# Override relevant files (tsv + RData) with the updated versions.
# -----------------------------------------------------------------

require(tibble)
require(dplyr)
source("utils.R")
source("data_organization/utils.R")

all.data <- load.all.datasets("processed_data")
for(i in 1:length(all.data)) assign(names(all.data)[i], all.data[[i]])
rm(all.data)
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
    save.to.files(basename(dataset), "processed_data", genera = genera, mtb = mtb, species = species)
    save.to.rdata(basename(dataset), "processed_data", genera = genera, mtb = mtb, species = species)
    message("Saved transposed genera, species and mtb tables")
  } else {
    save.to.files(basename(dataset), "processed_data", genera = genera, mtb = mtb)
    save.to.rdata(basename(dataset), "processed_data", genera = genera, mtb = mtb)
    message("Saved transposed genera and mtb tables")
  }
}
