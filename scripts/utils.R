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

get.gtdb.mapper <- function() {
  require(stringr)
  require(readr)
  
  add.prefix.to.silva <- function(x) {
    if (is.na(x)) return(NA)
    if (str_count(x, ";") != 6) return(NA)
    return(gsub("(.*);(.*);(.*);(.*);(.*);(.*);(.*)",
                "d__\\1;p__\\2;c__\\3;o__\\4;f__\\5;g__\\6;s__\\7", x))
  }
  
  bac120_metadata_r202 <- read_delim("../references/gtdb/bac120_metadata_r202.tsv", 
                                     delim = "\t", 
                                     escape_double = FALSE, 
                                     trim_ws = TRUE, 
                                     na = c(""," ","NA","none"))
  bac120_metadata_r202$ssu_silva_taxonomy <- sapply(bac120_metadata_r202$ssu_silva_taxonomy, 
                                                    add.prefix.to.silva, 
                                                    USE.NAMES = FALSE)
  bac120_metadata_r202 <- bac120_metadata_r202 %>%
    select(gtdb_taxonomy, 
           ssu_silva_taxonomy,
           ssu_gg_taxonomy,
           ncbi_taxonomy) %>%
    distinct() %>%
    tidyr::pivot_longer(cols = c("ssu_silva_taxonomy","ssu_gg_taxonomy","ncbi_taxonomy"),
                        names_to = "ref_db", values_to = "ref_taxonomy",
                        values_drop_na = TRUE) %>%
    mutate(gtdb_genus = gsub(";s__.*$","",gtdb_taxonomy)) %>%
    mutate(ref_genus = gsub(";s__.*$","",ref_taxonomy)) 
  
  return(bac120_metadata_r202)
}


# get.gtdb.genus <- function(genus.str, gtdb.map, use.ncbi = T, use.gg = F, use.silva = T) {
#   # message(paste("Mapping:", genus.str))
#   genus.only <- gsub(".*;g__","g__",genus.str)
#   
#   # Filter mapping table to relevant DBs only
#   gtdb.map.tmp <- gtdb.map
#   if (! use.ncbi) gtdb.map.tmp <- gtdb.map.tmp %>% filter(ref_db != "ncbi_taxonomy")
#   if (! use.gg) gtdb.map.tmp <- gtdb.map.tmp %>% filter(ref_db != "ssu_gg_taxonomy")
#   if (! use.silva) gtdb.map.tmp <- gtdb.map.tmp %>% filter(ref_db != "ssu_silva_taxonomy")
#   
#   # Check if taxon already matches GTDB taxonomy
#   if (genus.str %in% gtdb.map.tmp$gtdb_genus) return(c(genus = genus.str, comment = ""))
#   
#   # Otherwise, check if taxon matches reference taxonomies
#   tmp <- gtdb.map.tmp %>% filter(ref_genus == genus.str) %>% group_by(gtdb_genus) %>% summarise(N=n())
#   
#   # If the entire string is not found, we use only the genus part
#   if (nrow(tmp) == 0) {
#     # message(paste0("No matches for full name. Searching '", genus.only, "'"))
#     tmp <- gtdb.map.tmp %>% filter(grepl(genus.only, ref_genus)) %>% group_by(gtdb_genus) %>% summarise(N=n())
#   }
#   
#   # Single match
#   if (nrow(tmp) == 1) return(c(genus = tmp$gtdb_genus, comment = ""))
#   
#   # Multiple matches
#   if (nrow(tmp) > 1) {
#     tmp <- tmp %>% filter(N == max(N)) %>% pull(gtdb_genus)
#     # message("Found more than one match.")
#     if (length(tmp) == 1) return(c(genus = tmp, comment = "Multiple matches found. Choosing genus with max matches."))
#     
#     # In this rare case, a genus is mapped to >1 GTDB genera with no single dominating mapping. We therefore choose the one closer to current genus name.
#     best <- tmp[grepl(genus.only, tmp)]
#     if (length(best) == 0) best <- tmp
#     return(c(genus = best, comment = paste("Multiple matches found, with a tie in top matches. Choosing one of:", paste(tmp, collapse = ", ")))) 
#   }
#   # message(paste0("No match at all found for: ", genus.str, " --> keeping as is"))
#   return(c(genus = genus.str, comment = "No matches found. Keeping original name."))
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
