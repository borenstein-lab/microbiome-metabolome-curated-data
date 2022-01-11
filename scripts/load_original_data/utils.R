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

# Prepare a map between metaphlan species name to full taxonomy
#  based on the file in: 
#  https://bitbucket.org/biobakery/metaphlan2/downloads/mpa_v20_m200_marker_info.txt.bz2
# The new map will be saved to the file: 
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

get.genus.level <- function(species, species.mapping) {
  # We now want to map the original metaphlan labels to 
  #  a full genus-level taxonomy. 
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

save.to.files <- function(save_to_folder, 
                          metadata = NULL, 
                          mtb = NULL, 
                          mtb.map = NULL, 
                          genera = NULL, 
                          species = NULL) {
  require(readr)
  
  full.folder.path <- file.path("../data/processed_data", save_to_folder)
  dir.create(full.folder.path, showWarnings = FALSE)
  
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

save.to.rdata <- function(save_to_folder, 
                          metadata, 
                          mtb, 
                          mtb.map, 
                          genera, 
                          species = NULL) {
  full.folder.path <- file.path("../data/processed_data", save_to_folder)
  dir.create(full.folder.path, showWarnings = FALSE)
  
  files.to.save <- c("metadata","mtb","mtb.map","genera")
  if (!is.null(species)) files.to.save <- c(files.to.save, "species")
  
  save(list = files.to.save, 
       file = file.path(full.folder.path, ".RData"))
  
  message("Wrote data to rdata files")
}
