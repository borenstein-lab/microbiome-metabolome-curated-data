# ----------------------------------------------------------------
# The following script loads all processed datasets 
# ----------------------------------------------------------------

# Get all processed datasets
data.dirs <- list.dirs("../data/processed_data")[-1]

# Initialize table lists
metadata <- list()
mtb <- list()
mtb.map <- list()
genera <- list()
species <- list()

for (x in data.dirs) {
  # Create a temporary environment to hold all processed tables
  tmp.env <- new.env()
  dataset.name <- basename(x)
  
  # Load and save tables
  load(file.path(x, ".RData"), tmp.env)
  mtb[[dataset.name]] <- get('mtb', tmp.env) 
  mtb.map[[dataset.name]] <- get('mtb.map', tmp.env) 
  genera[[dataset.name]] <- get('genera', tmp.env) 
  metadata[[dataset.name]] <- get('metadata', tmp.env)
  if ("species" %in% ls(tmp.env)) species[[dataset.name]] <- get('species', tmp.env) 
  
  # Clean up
  rm(tmp.env)
}

message("Datasets loaded successfully")
rm(x, dataset.name)
