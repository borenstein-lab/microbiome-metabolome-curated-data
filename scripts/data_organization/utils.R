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

# --------------------------------------------------------
# Prepare a map between metaphlan species name to full taxonomy
#  based on the file in: 
#  https://bitbucket.org/biobakery/metaphlan2/downloads/mpa_v20_m200_marker_info.txt.bz2
# The new map will be saved to the file: 
#  species_to_full_taxonomy_map.tsv
# --------------------------------------------------------
make.metaphlan.species.mapper <- function() {
  
  # Load a table with full metaphlan species names 
  species.mapping <- read_delim("../references/metaphlan/parsed_marker_gene_info.tab", 
                                "\t", escape_double = FALSE, 
                                trim_ws = TRUE, 
                                comment = "#")
  # We only need the list of taxa, discarding other columns
  species.mapping <- species.mapping[,c("taxon")]
  
  # Remove strains, noted with t__ prefix
  species.mapping$taxon <- gsub("\\|t__.*","",species.mapping$taxon) 
  
  # Remove duplicates (metaphlan may have several marker genes per species/strain)
  species.mapping <- species.mapping %>% distinct(taxon)
  
  # Add a column for taxonomy at genus level (will have duplicates)
  species.mapping$taxon.genus <- gsub("\\|s__.*","",species.mapping$taxon) 
  
  # Extract the "short names" to be able to merge with this data
  short.names <- unname(sapply(species.mapping$taxon, 
                               FUN = function(x) {
                                 tmp <- tail(strsplit(x, split = '__', fixed = TRUE)[[1]],1)
                                 return(tmp)}))
  species.mapping$species.short <- short.names
  
  write_delim(x = species.mapping, 
              file = "../references/metaphlan/species_to_full_taxonomy_map.tsv",
              delim = "\t")
}

get.metaphlan.species.mapper <- function() {
  read_delim("../references/metaphlan/species_to_full_taxonomy_map.tsv", 
             delim = "\t", escape_double = FALSE, 
             trim_ws = TRUE)
}

# --------------------------------------------------------
# Map original metaphlan species labels to full 
#  genus-level taxonomy.
# --------------------------------------------------------
get.genus.level <- function(species, species.mapping) {
  #  The "tmp" table will hold all these intermediate mappings. 
  #  Due to the messiness of the raw data, 
  #  and mapping file limitations, we use a 
  #  few different mapping "patches".
  
  # 1. We map the raw data table to metaphlan genera labels 
  #  (ideally everything should have been mapped here)
  tmp <- merge(species, species.mapping, 
               by.x = "OTU", 
               by.y = "species.short", 
               all.x = T, all.y = FALSE) # For classified species
  
  # 2. Now create a mapping to genus level, 
  #  ONLY FOR ENTITIES IN THE METAPHLAN TABLE WHERE THE 
  #  GENUS IS CLASSIFIED BUT THE SPECIES IS NOT.
  #  (so the metaphlan entities are in the form of 
  #  "<genus_name>_unclassified")
  
  # First create a table that maps original OTU strings to full genus name
  genus.mapping <- species.mapping[,"taxon.genus"]
  
  # Remove duplicates (from multiple species in same genera)
  genus.mapping <- genus.mapping[!duplicated(genus.mapping),] 
  
  # Remove taxa that do not have a classified genus (family level or above only)
  genus.mapping <- genus.mapping[grepl("g__",genus.mapping$taxon.genus),] 
  
  # We add this dummy "unclassified" in order to match the "unclassified" entities in raw data 
  genus.mapping$genus.short <- paste0(gsub(".*g__","",genus.mapping$taxon.genus),"_unclassified") 
  
  # Match
  tmp <- merge(tmp, genus.mapping, 
               by.x = "OTU", 
               by.y = "genus.short", 
               all.x = T, all.f = FALSE) 
  
  # We now take the mapping from #1 or #2 above (whatever found)
  tmp$Genus <- ifelse(is.na(tmp$taxon.genus.y), tmp$taxon.genus.x, tmp$taxon.genus.y) 
  
  # 3. Some species are not in our metaphlan mapping 
  #  file but their genera are. 
  #  (verison issue? need to re-process with metaphlan3)  
  # In these cases the species are of form <known_genus>_<unknown_species>, 
  #  and there's only a ...g__<known_genus> record in the species.mapping file. 
  # The code below is a patch for these cases.
  for (unmapped.otu in tmp[is.na(tmp$Genus),"OTU"]) {
    genus.only <- strsplit(unmapped.otu, split = "_")[[1]][1]
    candidate.full.taxonomy <- grep(genus.only, 
                                    species.mapping$taxon, 
                                    value = TRUE)
    if (length(candidate.full.taxonomy)==0) { message(paste("Couldn't find a mapping for",unmapped.otu)); next; }
    if (length(candidate.full.taxonomy)>1) { message(paste("Found more than one mapping for",unmapped.otu)); next; }
    # print(paste("Mapping",unmapped.otu,"to",candidate.full.taxonomy))
    tmp$Genus[tmp$OTU == unmapped.otu] <- candidate.full.taxonomy
  }
  
  # Lastly, we add a few manual mappings
  # View(tmp[is.na(tmp$Genus),c("OTU","taxon","taxon.genus.x","taxon.genus.y","Genus")])
  # Note: some of these are unclassified genera, given by names such as "Bacteroidales_bacterium_ph8". 
  #  Some are viruses/phages, so we skip their mappings.
  tmp$Genus[tmp$OTU == "Eubacterium_infirmum"] <- "k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Eubacteriaceae|g__Eubacterium"
  tmp$Genus[tmp$OTU == "Morganella_morganii"] <- "k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacteriales|f__Enterobacteriaceae|g__Morganella"
  tmp$Genus[tmp$OTU == "Tropheryma_whipplei"] <- "k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Microbacteriaceae|g__Tropheryma"
  
  # Unclassified genera (bacteria only) are labeled "unclassified"
  tmp$Genus[is.na(tmp$Genus)] <- 'Unclassified'
  
  # Lastly, we group into genera. 
  genera <- tmp %>%
    select(-one_of(c("taxon","taxon.genus.y","taxon.genus.x","OTU"))) %>%
    filter(!is.na(Genus)) %>%
    group_by(Genus) %>%
    summarise_all(sum)
  # Sanity: hist(apply(genera[,-1], 2, sum), breaks = 20)
  
  # Here we add a species table, with full names instead of short names. We aggregate entities where species-level was unclassified 
  species <- tmp %>%
    rename(Species = taxon) %>%
    relocate(Species) %>%
    select(-one_of(c("OTU","taxon.genus.y","taxon.genus.x","Genus"))) %>%
    filter(!is.na(Species)) 
  # Sanity: hist(apply(species[,-1], 2, sum))
  
  return(list(species = species, genera = genera))
}

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

# load.gtdb.taxonomy <- function() {
#   # Load GTDB table, obtained from https://data.gtdb.ecogenomic.org/releases/latest/
#   bac120_taxonomy <- read_delim("../references/gtdb/bac120_taxonomy.tsv", 
#                                 delim = "\t", escape_double = FALSE, 
#                                 col_names = FALSE, trim_ws = TRUE)
#   tax.mapper <- bac120_taxonomy %>%
#     select(-1) %>%
#     distinct() %>%
#     rename(full.name = 1) %>%
#     mutate(species = gsub(".*;s__", "", full.name)) %>%
#     mutate(full.name.no.species = gsub(";s__.*", "", full.name)) %>%
#     mutate(genus = gsub(".*;g__", "", full.name.no.species))
#   
#   return(list(tax.mapper.sp = tax.mapper %>% 
#                 select(full.name, species) %>% 
#                 distinct(), 
#               tax.mapper.ge = tax.mapper %>% 
#                 select(full.name.no.species, genus) %>% 
#                 distinct()))
# }


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

enrichment <- function(correlated.genera, all.genera, genus.groups, adj = "fdr") {
  # Calculate one-sided Fisher exact test per group
  fish.per.group <- lapply(genus.groups, function(gg) {
    
    # Build 2X2 confusion matrix for fisher.test
    non.cor.genera <- all.genera[!all.genera %in% correlated.genera]
    N.ncor.in.group  <- sum(non.cor.genera %in% gg)
    N.ncor.nin.group <- length(non.cor.genera) - N.ncor.in.group
    N.cor.in.group   <- sum(correlated.genera %in% gg)
    N.cor.nin.group  <- length(correlated.genera) - N.cor.in.group
    fmat <- matrix(c(N.cor.in.group, N.ncor.in.group, 
                     N.cor.nin.group, N.ncor.nin.group), nrow = 2, 
                   ncol = 2, byrow = F)
    colnames(fmat) <- c("in.group", "not.in.group")
    rownames(fmat) <- c("correlated", "not.correlated")
    
    # Calculate Fisher
    fish <- fisher.test(fmat, alternative = "greater")
    
    # Organize results
    p.value <- fish$p.value
    N.in.group <- N.ncor.in.group + N.cor.in.group
    data.frame("N.genera.correlated.in.phylum" = N.cor.in.group, 
               "N.genera.in.phylum" = N.in.group, 
               "enrichment.p.value" = p.value)
  })
  
  fish.per.group <- bind_rows(fish.per.group, .id = "phylum") %>%
    arrange(enrichment.p.value) %>%
    mutate(enrichment.fdr = p.adjust(enrichment.p.value, method = adj))
  
  return(fish.per.group)
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
