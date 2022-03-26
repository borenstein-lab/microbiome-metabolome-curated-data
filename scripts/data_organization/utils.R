require(readr)

# ----------------------------------------------------------------
# Utility functions for data loading
# ----------------------------------------------------------------

map.compound.names.MetaboAnalyst <- function(cmpds.to.search, search.by = "name") {
  require(MetaboAnalystR)
  
  # Create a MetaboAnalyst object
  MA.obj <- InitDataObjects(data.type = "NA", 
                            anal.type = "utils")
  MA.obj <- Setup.MapData(mSetObj = MA.obj, 
                          qvec = cmpds.to.search)
  
  # Map compound names to references
  MA.obj <- CrossReferencing(mSetObj = MA.obj, 
                             q.type = search.by, 
                             hmdb = T,
                             metlin = F, 
                             pubchem = F, 
                             chebi = F)
  MA.obj <- CreateMappingResultTable(MA.obj)
  
  # Organize table of MetaboAnalyst matches
  MA.matches <- data.frame(MA.obj$dataSet$map.table)
  MA.matches <- MA.matches %>% 
    # Fix missing values
    mutate(across(where(is.character), ~na_if(., "NA"))) %>%
    mutate(across(where(is.character), ~na_if(., ""))) %>%
    select(Query, Match, HMDB, KEGG) %>%
    rename(MA.Name.Match = Match)
  
  return(MA.matches)
}

# ----------------------------------------------------------------
# Utility functions for saving processed data
# ----------------------------------------------------------------

save.to.files <- function(new.folder, 
                          parent.folder,
                          metadata = NULL, 
                          mtb = NULL, 
                          mtb.map = NULL, 
                          genera = NULL, 
                          species = NULL) {
  full.folder.path <- file.path("../data", parent.folder, new.folder)
  dir.create(full.folder.path, showWarnings = FALSE, recursive = TRUE)
  
  if (!is.null(metadata)) write_delim(metadata, 
                                      file.path(full.folder.path, "metadata.tsv"),
                                      delim = "\t")
  if (!is.null(mtb))      write_delim(mtb, 
                                      file.path(full.folder.path, "mtb.tsv"),
                                      delim = "\t")
  if (!is.null(mtb.map))  write_delim(mtb.map, 
                                      file.path(full.folder.path, "mtb.map.tsv"),
                                      delim = "\t")
  if (!is.null(genera))   write_delim(genera, 
                                      file.path(full.folder.path, "genera.tsv"),
                                      delim = "\t")
  if (!is.null(species))  write_delim(species, 
                                      file.path(full.folder.path, "species.tsv"),
                                      delim = "\t")
  message("Wrote data to text files")
}

save.to.rdata <- function(new.folder,
                          parent.folder,
                          metadata = NULL, 
                          mtb = NULL, 
                          mtb.map = NULL, 
                          genera = NULL, 
                          species = NULL,
                          override.all = FALSE) {
  full.folder.path <- file.path("../data", parent.folder, new.folder)
  dir.create(full.folder.path, showWarnings = FALSE, recursive = TRUE)
  rdata.file.path <- file.path(full.folder.path, ".RData")
  
  if (override.all | ! file.exists(rdata.file.path)) {
    files.to.save <- c("metadata","mtb","mtb.map","genera")
    if (!is.null(species)) files.to.save <- c(files.to.save, "species")
    save(list = files.to.save, file = rdata.file.path)
    message("Wrote data to rdata file")
  } else {
    require(cgwtools)
    if (!is.null(metadata)) resave(metadata, file = rdata.file.path)
    if (!is.null(mtb))      resave(mtb, file = rdata.file.path)
    if (!is.null(mtb.map))  resave(mtb.map, file = rdata.file.path)
    if (!is.null(genera))   resave(genera, file = rdata.file.path)
    if (!is.null(species))  resave(species, file = rdata.file.path)
    message("Updated rdata file")
  }
}

# ----------------------------------------------------------------
# Utility functions for data processing and data analysis
# ----------------------------------------------------------------

load.all.datasets <- function(parent.folder = "processed_data") {
  # Get all processed datasets
  data.dirs <- list.dirs(file.path("../data", parent.folder))[-1]
  
  # Initialize table lists
  all.data <- list()
  all.data$data.dirs <- data.dirs
  all.data$metadata <- list()
  all.data$mtb <- list()
  all.data$mtb.map <- list()
  all.data$genera <- list()
  all.data$species <- list()
  
  for (x in data.dirs) {
    # Create a temporary environment to hold all processed tables
    tmp.env <- new.env()
    dataset.name <- basename(x)
    
    # Load and save tables
    load(file.path(x, ".RData"), tmp.env)
    all.data$mtb[[dataset.name]] <- get('mtb', tmp.env) 
    all.data$mtb.map[[dataset.name]] <- get('mtb.map', tmp.env) 
    all.data$genera[[dataset.name]] <- get('genera', tmp.env) 
    all.data$metadata[[dataset.name]] <- get('metadata', tmp.env)
    if ("species" %in% ls(tmp.env)) all.data$species[[dataset.name]] <- get('species', tmp.env) 
    
    # Clean up
    rm(tmp.env)
  }
  
  message("Datasets loaded successfully")
  return(all.data)
}

get.genera.dataset.stats <- function(genera, datasets) {
  # Initialize table for genera statistics 
  genera.dataset.stats <- 
    data.frame(Taxon = character(0),
               Dataset = character(0),
               Dataset.N = integer(0),
               Taxon.Mean.Abundance = numeric(0),
               Taxon.Var.Abundance = numeric(0),
               Taxon.Perc.of.Non.Zeros = numeric(0))
  
  # We iterate over each dataset and record a few basic stats 
  #  per each genus in each dataset
  for (dataset in datasets) {
    tmp <- genera[[dataset]] %>% select(-Sample)
    tmp.genera <- colnames(tmp)
    tmp.means <- unname(apply(tmp, MARGIN = 2, mean))
    tmp.vars <- unname(apply(tmp, MARGIN = 2, var))
    tmp.non.zero.perc <- unname(apply(tmp, MARGIN = 2, 
                                      function(v) {100*sum(v>0)/length(v)}))
    genera.dataset.stats <- 
      bind_rows(genera.dataset.stats, 
                data.frame(Taxon = tmp.genera, 
                           Dataset = dataset,
                           Dataset.N = nrow(tmp),
                           Taxon.Mean.Abundance = tmp.means,
                           Taxon.Var.Abundance = tmp.vars,
                           Taxon.Perc.of.Non.Zeros = tmp.non.zero.perc))
  }
  
  return(genera.dataset.stats)
}

get.metab.dataset.stats <- function(mtb.map, datasets) {
  metabolites.per.dataset <- data.frame(stringsAsFactors = F)
  
  for (dataset in datasets) {
    tmp <- mtb.map[[dataset]] %>%
      select(Compound, KEGG, HMDB) %>%
      rename(Orig.Compound = Compound) %>%
      tidyr::pivot_longer(cols = c("KEGG","HMDB"), 
                          names_to = "Type", 
                          values_to = "Compound", 
                          values_drop_na = TRUE) %>%
      mutate(Dataset = dataset) 
    metabolites.per.dataset <- bind_rows(metabolites.per.dataset, tmp)
  }
  
  return(metabolites.per.dataset)
}

# Get significant marks given p values, for plotting
get.signif.marks <- function(p.vals) {
  signif.marks <- sapply(p.vals,
                         function(x) {
                           ifelse(x <= 0.001, "***", 
                                  ifelse(x <= 0.01, "**", 
                                         ifelse(x <= 0.05, "*", 
                                                ifelse(x <= 0.1, "^", ""))))
                         })
  signif.marks[is.na(signif.marks)] <- ""
  return(signif.marks)
}

add.var.to.formula <- function(form, new.var) {
  require(formula.tools)
  new.form.str <- paste(as.character(form), "+", new.var)
  return(as.formula(new.form.str))
}