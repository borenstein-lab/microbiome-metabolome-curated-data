# ----------------------------------------------------------------
# The following script loads and processes data from:
# 
# Lloyd-Price, Jason, et al. "Multi-omics of the gut microbial 
#  ecosystem in inflammatory bowel diseases." 
#  Nature 569.7758 (2019): 655-662.
# 
# The output tables are:
# 1) metadata - subject metadata
# 2) genera - taxonomic profiles per subject (genus-level)
# 3) species - taxonomic profiles per subject (species-level)
# 4) mtb - metabolomic profiles per subject       
# 5) mtb.map - additional details per metabolite 
# ----------------------------------------------------------------

require(readr)
require(dplyr)
require(MetaboAnalystR)

source("data_organization/utils.R")
source("data_organization/gtdb_utils.R")

# --------------------------------
# Required files & info
# --------------------------------

# For details about the source of each file below see: <COMPLETE>
METADATA_FILE <- "../data/original_data/iHMP_IBDMDB_2019/hmp2_metadata.csv"
TAXONOMY_FILE_SP <- '../data/original_data/iHMP_IBDMDB_2019/kraken/kraken_species_level_taxonomy.tsv'
TAXONOMY_FILE_GE <- '../data/original_data/iHMP_IBDMDB_2019/kraken/kraken_genus_level_taxonomy.tsv'
TAXONOMY_SAMPLE_MAP <- '../data/original_data/iHMP_IBDMDB_2019/kraken/ena_download.txt'
METABOLOMICS_FILE <- "../data/original_data/iHMP_IBDMDB_2019/iHMP_metabolomics.csv"

PUBLICATION_NAME <- 'Multi-omics of the gut microbial ecosystem in inflammatory bowel diseases'
DOI <- '10.1038/s41586-019-1237-9'
DATASET_NAME <- 'iHMP_IBDMDB_2019'

# --------------------------------
# Load metadata
# --------------------------------

metadata <- read_csv(METADATA_FILE, 
                     guess_max = 7000, 
                     show_col_types = FALSE)

# Organize column names & order
metadata <- metadata %>%
  # Keep only metagenomics+metabolomics samples
  filter(data_type %in% c('metabolomics','metagenomics')) %>%
  select(`External ID`,`Participant ID`,
         site_sub_coll, ProjectSpecificID,
         week_num, date_of_receipt,
         interval_days, visit_num, site_name,
         `Age at diagnosis`, consent_age,
         diagnosis, Antibiotics, Weight, BMI,
         Height, Weight.1, is_inflamed, race,
         sex, `smoking status`, stool_id, 
         fecalcal) %>%
  distinct() %>%
  rename(Sample = `External ID`) %>%
  rename(Subject = `Participant ID`) %>%
  rename(Study.Group = diagnosis) %>%
  mutate(Age.Units = "Years") %>%
  rename(Gender = sex) %>%
  # Add publication info
  mutate(Dataset = DATASET_NAME) %>%
  mutate(DOI = DOI) %>%
  mutate(Publication.Name = PUBLICATION_NAME) %>%
  # Reorder
  relocate(c(Dataset, Sample, Subject, Study.Group, Gender, BMI, DOI, Publication.Name, consent_age, Age.Units)) %>%
  # Deal with some cases where a record is duplicated due to the fecalcal value 
  mutate(fecalcal_dummy = ifelse(is.na(fecalcal),0,fecalcal)) %>%
  group_by(Dataset, Sample, Subject) %>%
  slice_min(fecalcal_dummy) %>%
  ungroup() %>%
  select(-fecalcal_dummy)

# --------------------------------
# Load taxonomic profiles 
# --------------------------------

# Load kraken files
species <- read_delim(TAXONOMY_FILE_SP, "\t", 
                      escape_double = FALSE, 
                      trim_ws = TRUE,
                      show_col_types = FALSE)
names(species)[1] <- 'Species'

genera <- read_delim(TAXONOMY_FILE_GE, "\t", 
                     escape_double = FALSE, 
                     trim_ws = TRUE,
                     show_col_types = FALSE)
names(genera)[1] <- 'Genus'

tax.map <- read_delim(TAXONOMY_SAMPLE_MAP, 
                       delim = "\t", 
                       escape_double = FALSE, 
                       trim_ws = TRUE, 
                       col_select = c("run_accession", "experiment_title"),
                      show_col_types = FALSE)
tax.map <- tax.map %>%
  filter(grepl("WGS", experiment_title)) %>%
  mutate(Sample = gsub("Illumina HiSeq 2000 sequencing; ", "", experiment_title)) %>%
  mutate(Sample = gsub(" stool WGS", "", Sample)) %>%
  filter(run_accession %in% names(genera)) %>%
  select(-experiment_title)

tax.map.vec <- c("Species", "Genus", tax.map$Sample)
names(tax.map.vec) <- c("Species", "Genus", tax.map$run_accession)

# Map file names to sample id's
names(species) <- unname(tax.map.vec[names(species)])
names(genera) <- unname(tax.map.vec[names(genera)])

# We remove samples with less than 50K reads 
low.depth.samples <- names(which(colSums(genera %>% select(-Genus)) < 50000))
species <- species[,! names(species) %in% low.depth.samples]
genera <- genera[,! names(genera) %in% low.depth.samples]

# Map species/genus short names to long versions
species$Species <- map.gtdb.short.to.long(species$Species, level = "species")
genera$Genus <- map.gtdb.short.to.long(genera$Genus, level = "genera")

# --------------------------------
# Load metabolomic profiles
# --------------------------------

mtb <- read_csv(METABOLOMICS_FILE, show_col_types = FALSE)
mtb$Compound <- paste(mtb$Compound, mtb$Metabolite, sep = "__")

mtb.map <- data.frame(Compound = mtb$Compound,
                      HMDB = mtb$`HMDB (*Representative ID)`,
                      m.z = mtb$`m/z`,
                      Method = mtb$Method,
                      Compound.Name = mtb$Metabolite,
                      High.Confidence.Annotation = TRUE,
                      stringsAsFactors = FALSE)

# Use representative HMDB ID's as final ID's, but mark them with low confidence... 
#  (erase stars and non-HMDB ID's)
mtb.map$High.Confidence.Annotation[grepl("\\*",mtb.map$HMDB)] <- FALSE
mtb.map$HMDB <- gsub("\\*","",mtb.map$HMDB)
mtb.map$HMDB <- gsub("HMDB","HMDB00",mtb.map$HMDB)
mtb.map$HMDB[!is.na(mtb.map$HMDB) & mtb.map$HMDB == "redundant ion"] <- NA

# Keep only compound name and data in main mtb table
mtb <- mtb[,7:ncol(mtb)]

# We'll now use MetaboAnalyst to also get kegg id's. 
# (we will search by hmdb id)
cmpds.to.search <- mtb.map$HMDB
cmpds.to.search <- cmpds.to.search[!is.na(cmpds.to.search)]
cmpds.to.search <- unique(cmpds.to.search)

MA.matches <- map.compound.names.MetaboAnalyst(cmpds.to.search, search.by = "hmdb")

# Some HMDB ID's weren't matched by MetaboAnalyst for some reason: View(MA.matches[is.na(MA.matches$MA.Name.Match),])
# For some, this is because their HMDB ID was updated. We add them manually:
MA.matches[MA.matches$Query == "HMDB0062548","HMDB"] <- "HMDB0031057"
MA.matches[MA.matches$Query == "HMDB0062548","MA.Name.Match"] <- "2-Hydroxyhexadecanoic acid"
MA.matches[MA.matches$Query == "HMDB0062781","HMDB"] <- "HMDB0000208"
MA.matches[MA.matches$Query == "HMDB0059655","KEGG"] <- "C00026"
MA.matches[MA.matches$Query == "HMDB0062781","MA.Name.Match"] <- "alpha-Ketoglutaric acid"
MA.matches[MA.matches$Query == "HMDB0041876","HMDB"] <- "HMDB0002172" 
MA.matches[MA.matches$Query == "HMDB0041876","KEGG"] <- "C03413"
MA.matches[MA.matches$Query == "HMDB0041876","MA.Name.Match"] <- "N1,N12-Diacetylspermine"
MA.matches[MA.matches$Query == "HMDB0061710","HMDB"] <- "HMDB0037397"
MA.matches[MA.matches$Query == "HMDB0061710","MA.Name.Match"] <- "17-Methyloctadecanoic acid"

MA.matches[MA.matches$Query == "HMDB0004159","HMDB"] <- "HMDB0004160" 
MA.matches[MA.matches$Query == "HMDB0004159","KEGG"] <- "C05794"
MA.matches[MA.matches$Query == "HMDB0004159","MA.Name.Match"] <- "Urobilin"
MA.matches[MA.matches$Query == "HMDB0004161","HMDB"] <- "HMDB0004160" 
MA.matches[MA.matches$Query == "HMDB0004161","KEGG"] <- "C05794"
MA.matches[MA.matches$Query == "HMDB0004161","MA.Name.Match"] <- "Urobilin"

# Additional manual mappings:
MA.matches[MA.matches$Query == "HMDB0000208","KEGG"] <- "C00026"
MA.matches[MA.matches$Query == "HMDB0000651","KEGG"] <- "C03299"
MA.matches[MA.matches$Query == "HMDB0000688","KEGG"] <- "C20826"
MA.matches[MA.matches$Query == "HMDB0000212","KEGG"] <- "C01074"

# We now merge the MetaboAnalyst mappings with the main mtb.map table
# (the final HMDB will be taken from MA.matches$HMDB and not original data, as some HMDB ID's are not updated in the original table)
names(mtb.map)[2] <- "Query" 
mtb.map <- merge(mtb.map, MA.matches, by = "Query", all = T)
mtb.map <- mtb.map %>%
  select(-Query, -Compound.Name) %>%
  rename(Compound.Name = MA.Name.Match)

# Mark cases of duplicated HMDN/KEGG ID as lower confidence
mtb.map$KEGG <- trimws(mtb.map$KEGG)
mtb.map$HMDB <- trimws(mtb.map$HMDB)
kegg.dups <- names(table(mtb.map$KEGG)[table(mtb.map$KEGG) > 1])
hmdb.dups <- names(table(mtb.map$HMDB)[table(mtb.map$HMDB) > 1])
mtb.map$High.Confidence.Annotation[mtb.map$KEGG %in% kegg.dups] <- FALSE
mtb.map$High.Confidence.Annotation[mtb.map$HMDB %in% hmdb.dups] <- FALSE

# --------------------------------
# Keep only samples with all data
# --------------------------------

sample.intersect <- Reduce(intersect, list(names(mtb)[-1], names(species)[-1], names(genera)[-1], metadata$Sample))
message(paste(length(sample.intersect),"samples have all data types"))

mtb <- mtb[,c("Compound",sample.intersect)]
species <- species[,c("Species",sample.intersect)]
genera <- genera[,c("Genus",sample.intersect)]
metadata <- metadata[metadata$Sample %in% sample.intersect,]

# --------------------------------
# Save to files + R objects
# --------------------------------

save.to.files(DATASET_NAME, "prelim_data", metadata = metadata, mtb = mtb, mtb.map = mtb.map, genera = genera, species = species)
save.to.rdata(DATASET_NAME, "prelim_data", metadata = metadata, mtb = mtb, mtb.map = mtb.map, genera = genera, species = species, override.all = T)

rm(list = ls())

