# ----------------------------------------------------------------
# The following script loads and processes data from:
# 
# Yachida, Shinichi, et al. "Metagenomic and metabolomic analyses 
#  reveal distinct stage-specific phenotypes of the gut microbiota 
#  in colorectal cancer." Nature medicine 25.6 (2019): 968-976.
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

source("load_original_data/utils.R")

# --------------------------------
# Required files & info
# --------------------------------

# For details about the source of each file below see: <COMPLETE>
METADATA_FILE <- "~/MICROBIOME_METABOLOME/YACHIDA_CRC_2019/RAW_DATA/Table_S2-1_Metadata.tsv" # I read the metadata table from S2-1 which is the metadata of subjects for which a metagenome sample was taken. 
TAXONOMY_FILE <- "~/MICROBIOME_METABOLOME/YACHIDA_CRC_2019/RAW_DATA/Table_S9_Metaphlan_Profiles.tsv"
METABOLOMICS_FILE <- "~/MICROBIOME_METABOLOME/YACHIDA_CRC_2019/RAW_DATA/Table_S13_MS_Profiles.tsv"

PUBLICATION_NAME <- 'Metagenomic and metabolomic analyses reveal distinct stage-specific phenotypes of the gut microbiota in colorectal cancer'
DOI <- '10.1038/s41591-019-0458-7'
DATASET_NAME <- 'YACHIDA_CRC_2019'

# --------------------------------
# Load metadata
# --------------------------------


metadata <- read_delim(METADATA_FILE, 
                       "\t", 
                       escape_double = FALSE, 
                       trim_ws = TRUE)

# Organize column names & order
metadata <- metadata %>%
  rename(Sample = 1) %>%
  mutate(Subject = Sample) %>%
  rename(Study.Group = Group) %>%
  mutate(Age.Units = "Years") %>%
  mutate(Gender = ifelse(Gender == "F", "Female", "Male")) %>%
  # Add publication info
  mutate(Dataset = DATASET_NAME) %>%
  mutate(DOI = DOI) %>%
  mutate(Publication.Name = PUBLICATION_NAME) %>%
  # Reorder
  relocate(c(Dataset, Sample, Subject, Study.Group, Age, Age.Units, Gender, BMI, DOI, Publication.Name)) 


# --------------------------------
# Load taxonomic profiles (species)
# --------------------------------

mtg <- read_delim(TAXONOMY_FILE, "\t", 
                  escape_double = FALSE, 
                  trim_ws = TRUE)
names(mtg)[1] <- 'OTU'
mtg$OTU <- gsub(" ","_",mtg$OTU)

# Sanity (abundances sum to ~100): table(cut(colSums(mtg[,2:ncol(mtg)]), breaks = c(0,96,97,98,99,101)))

# --------------------------------
# Load metabolomic profiles
# --------------------------------

mtb <- read_delim(METABOLOMICS_FILE, 
                  "\t", escape_double = FALSE, 
                  trim_ws = TRUE)
names(mtb)[1] <- "Compound"

# Create mapping file to KEGG ID's, by parsing metab names
mtb.map <- data.frame(Compound = mtb$Compound)

# Extract KEGG compound IDs when available. 
mtb.map$KEGG <- gsub("_.*","",mtb.map$Compound) 

# If there's no KEGG ID, the previous command may return 
#  parts of the metabolite name in some cases, 
#  so we further fix this:
mtb.map$KEGG[!grepl("C\\d{5}",mtb.map$KEGG)] <- NA
mtb.map$Compound.Name <- gsub("(^C\\d{5}_)|(^_)","",mtb.map$Compound)

# We'll now use MetaboAnalyst to also get hmdb id's. 
#  We will search by name but validate against the 
#  already available KEGG ID
cmpds.to.search <- mtb.map$Compound.Name
MA.matches <- map.compound.names.MetaboAnalyst(cmpds.to.search)
names(MA.matches)[4] <- "KEGG.new" 

# Merge with the names in our data
mtb.map <- merge(mtb.map, MA.matches, 
                 by.x = "Compound.Name", 
                 by.y = "Query", all = T)
mtb.map$High.Confidence.Annotation <- TRUE

# Any inconsistencies?
# View(mtb.map[!is.na(mtb.map$KEGG) & !is.na(mtb.map$KEGG.new) & mtb.map$KEGG != mtb.map$KEGG.new,])

# Fix inconsistencies, take authors mapping + 
#  mark them with low confidence
mtb.map[mtb.map$Compound.Name == "2,3-Diaminopropionate","High.Confidence.Annotation"] <- FALSE
mtb.map[mtb.map$Compound.Name == "3-Phenyllactate","High.Confidence.Annotation"] <- FALSE
mtb.map[mtb.map$Compound.Name == "Citramalate","High.Confidence.Annotation"] <- FALSE
mtb.map[mtb.map$Compound.Name == "G3P","High.Confidence.Annotation"] <- FALSE
mtb.map[mtb.map$Compound.Name == "Malate","High.Confidence.Annotation"] <- FALSE
mtb.map[mtb.map$Compound.Name == "Sorbitol 6-phosphate","High.Confidence.Annotation"] <- FALSE

# Additional mappings (based on searching in kegg/hmdb)
mtb.map[mtb.map$Compound.Name == "5-Hydroxytryptophan","KEGG"] <- 'C00643'
mtb.map[mtb.map$Compound.Name == "5-Hydroxytryptophan","High.Confidence.Annotation"] <- FALSE
mtb.map[mtb.map$Compound.Name == "Hydroxyproline","KEGG"] <- 'C01015'
mtb.map[mtb.map$Compound.Name == "Hydroxyproline","High.Confidence.Annotation"] <- FALSE
mtb.map[mtb.map$Compound.Name == "N-Acetylneuraminate","KEGG"] <- 'C00270'
mtb.map[mtb.map$Compound.Name == "N-Acetylneuraminate","High.Confidence.Annotation"] <- FALSE
mtb.map[mtb.map$Compound.Name == "Nicotine","KEGG"] <- 'C00745'
mtb.map[mtb.map$Compound.Name == "Nicotine","High.Confidence.Annotation"] <- FALSE
mtb.map[mtb.map$Compound.Name == "Quinate","KEGG"] <- 'C00296'
mtb.map[mtb.map$Compound.Name == "Quinate","High.Confidence.Annotation"] <- FALSE
mtb.map[mtb.map$Compound.Name == "1-Aminocyclopropane-1-carboxylate","HMDB"] <- 'HMDB0036458'
mtb.map[mtb.map$Compound.Name == "2',3'-cCMP","HMDB"] <- 'HMDB0011691'
mtb.map[!is.na(mtb.map$KEGG) & mtb.map$KEGG == 'C00141','HMDB'] <- 'HMDB0000019'
mtb.map[!is.na(mtb.map$KEGG) & mtb.map$KEGG == 'C00099','HMDB'] <- 'HMDB0000056'
mtb.map[!is.na(mtb.map$KEGG) & mtb.map$KEGG == 'C00085','HMDB'] <- 'HMDB0000124'
mtb.map[!is.na(mtb.map$KEGG) & mtb.map$KEGG == 'C00051','HMDB'] <- 'HMDB0000125'
mtb.map[!is.na(mtb.map$KEGG) & mtb.map$KEGG == 'C00148','HMDB'] <- 'HMDB0000162'
mtb.map[!is.na(mtb.map$KEGG) & mtb.map$KEGG == 'C02356','HMDB'] <- 'HMDB0000452'
mtb.map[!is.na(mtb.map$KEGG) & mtb.map$KEGG == 'C00199','HMDB'] <- 'HMDB0000618'
mtb.map[!is.na(mtb.map$KEGG) & mtb.map$KEGG == 'C00064','HMDB'] <- 'HMDB0000641'
mtb.map[!is.na(mtb.map$KEGG) & mtb.map$KEGG == 'C05824','HMDB'] <- 'HMDB0000731'
mtb.map[!is.na(mtb.map$KEGG) & mtb.map$KEGG == 'C00438','HMDB'] <- 'HMDB0000828'
mtb.map[!is.na(mtb.map$KEGG) & mtb.map$KEGG == 'C03145','HMDB'] <- 'HMDB0001015'
mtb.map[!is.na(mtb.map$KEGG) & mtb.map$KEGG == 'C00354','HMDB'] <- 'HMDB0001058'
mtb.map[!is.na(mtb.map$KEGG) & mtb.map$KEGG == 'C05382','HMDB'] <- 'HMDB0001068'
mtb.map[!is.na(mtb.map$KEGG) & mtb.map$KEGG == 'C01239','HMDB'] <- 'HMDB0001104'
mtb.map[!is.na(mtb.map$KEGG) & mtb.map$KEGG == 'C00019','HMDB'] <- 'HMDB0001185'
mtb.map[!is.na(mtb.map$KEGG) & mtb.map$KEGG == 'C02567','HMDB'] <- 'HMDB0001186'
mtb.map[!is.na(mtb.map$KEGG) & mtb.map$KEGG == 'C00630','HMDB'] <- 'HMDB0001243'
mtb.map[!is.na(mtb.map$KEGG) & mtb.map$KEGG == 'C03793','HMDB'] <- 'HMDB0001325'
mtb.map[!is.na(mtb.map$KEGG) & mtb.map$KEGG == 'C00103','HMDB'] <- 'HMDB0001586'
mtb.map[!is.na(mtb.map$KEGG) & mtb.map$KEGG == 'C00127','HMDB'] <- 'HMDB0003337'
mtb.map[!is.na(mtb.map$KEGG) & mtb.map$KEGG == 'C00631','HMDB'] <- 'HMDB0003391'
mtb.map[!is.na(mtb.map$KEGG) & mtb.map$KEGG == 'C01771','HMDB'] <- 'HMDB0010720'
mtb.map[!is.na(mtb.map$KEGG) & mtb.map$KEGG == 'C00371','HMDB'] <- 'HMDB0012204'
mtb.map[!is.na(mtb.map$KEGG) & mtb.map$KEGG == 'C02595','HMDB'] <- 'HMDB0029718'
mtb.map[!is.na(mtb.map$KEGG) & mtb.map$KEGG == 'C03173','HMDB'] <- 'HMDB0035646'
mtb.map[!is.na(mtb.map$KEGG) & mtb.map$KEGG == 'C08304','HMDB'] <- 'HMDB0035762'
mtb.map[!is.na(mtb.map$KEGG) & mtb.map$KEGG == 'C03969','HMDB'] <- 'HMDB0062225'
mtb.map[!is.na(mtb.map$KEGG) & mtb.map$KEGG == 'C02630','HMDB'] <- 'HMDB0059655'
mtb.map[!is.na(mtb.map$KEGG) & mtb.map$KEGG == 'C02129','HMDB'] <- 'HMDB0061881'
mtb.map[!is.na(mtb.map$KEGG) & mtb.map$KEGG == 'C03678','HMDB'] <- 'HMDB0061877'
mtb.map[!is.na(mtb.map$KEGG) & mtb.map$KEGG == 'C05341','HMDB'] <- 'HMDB0060442'
mtb.map[!is.na(mtb.map$KEGG) & mtb.map$KEGG == 'C02273','HMDB'] <- 'HMDB0039721'
mtb.map[!is.na(mtb.map$KEGG) & mtb.map$KEGG == 'C00092','HMDB'] <- 'HMDB0001401'
mtb.map[!is.na(mtb.map$KEGG) & mtb.map$KEGG == 'C08737','HMDB'] <- 'HMDB0041922'
mtb.map[!is.na(mtb.map$KEGG) & mtb.map$KEGG == 'C01073','HMDB'] <- 'HMDB0061880'
mtb.map[!is.na(mtb.map$KEGG) & mtb.map$KEGG == 'C03958','HMDB'] <- 'HMDB0060080'
mtb.map[!is.na(mtb.map$KEGG) & mtb.map$KEGG == 'C03172','HMDB'] <- 'HMDB0038670'
mtb.map[!is.na(mtb.map$KEGG) & mtb.map$KEGG == 'C04501','HMDB'] <- 'HMDB0001367'
mtb.map[!is.na(mtb.map$KEGG) & mtb.map$KEGG == 'C04501',"High.Confidence.Annotation"] <- FALSE
mtb.map[!is.na(mtb.map$KEGG) & mtb.map$KEGG == 'C02713','HMDB'] <- 'HMDB0060493'
mtb.map[!is.na(mtb.map$KEGG) & mtb.map$KEGG == 'C01044','HMDB'] <- 'HMDB0060495'
mtb.map[!is.na(mtb.map$KEGG) & mtb.map$KEGG == 'C01047','HMDB'] <- 'HMDB0034365'
mtb.map[!is.na(mtb.map$KEGG) & mtb.map$KEGG == 'C02721','HMDB'] <- 'HMDB0094692'
mtb.map[!is.na(mtb.map$KEGG) & mtb.map$KEGG == 'C01046','HMDB'] <- 'HMDB0062660'
mtb.map[mtb.map$Compound.Name == "N1,N8-Diacetylspermidine",'HMDB'] <- 'HMDB0041947'
mtb.map[mtb.map$Compound.Name == "Ophthalmate",'HMDB'] <- 'HMDB0005765'

# Clean mapping table
mtb.map$Compound.Name <- NULL
mtb.map$KEGG.new <- NULL

# Mark duplicated kegg/hmdb identifiers
#  table(mtb.map$KEGG)[table(mtb.map$KEGG) > 1]
#  table(mtb.map$HMDB)[table(mtb.map$HMDB) > 1]
mtb.map <- mtb.map %>%
  group_by(KEGG) %>%
  mutate(High.Confidence.Annotation = ifelse(!is.na(KEGG) & n() > 1,
                                             FALSE,
                                             High.Confidence.Annotation)) %>%
  ungroup() %>%
  group_by(HMDB) %>%
  mutate(High.Confidence.Annotation = ifelse(!is.na(HMDB) & n() > 1,
                                             FALSE,
                                             High.Confidence.Annotation)) %>%
  ungroup()
  
# --------------------------------
# Keep only samples with all data
# --------------------------------

sample.intersect <- Reduce(intersect, list(names(mtb)[-1], names(mtg)[-1], metadata$Sample))
message(paste(length(sample.intersect),"samples have all data types"))

mtb <- mtb[,c("Compound",sample.intersect)]
mtg <- mtg[,c("OTU",sample.intersect)]
metadata <- metadata[metadata$Sample %in% sample.intersect,]

# --------------------------------
# Collapse species to genera
# --------------------------------

# Get metaphlan species name mapper
species.mapping <- get.metaphlan.species.mapper()

# Convert current table to a species table + genus-level table
l <- get.genus.level(mtg, species.mapping)
species <- l$species
genera <- l$genera

species <- species %>% 
  filter(!grepl("k__Archaea", Species)) 

# Sanity: hist(apply(species[,-1], 2, sum))

# --------------------------------
# Save to files + R objects
# --------------------------------

save.to.files(DATASET_NAME, metadata, mtb, mtb.map, genera, species)
save.to.rdata(DATASET_NAME, metadata, mtb, mtb.map, genera, species)
rm(list = ls())
