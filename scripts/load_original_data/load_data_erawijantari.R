# ----------------------------------------------------------------
# The following script loads and processes data from:
# 
# Erawijantari et al. Influence of gastrectomy for gastric cancer 
#  treatment on faecal microbiome and metabolome profiles. 
#  Gut. 2020 Aug;69(8):1404-1415. doi: 10.1136/gutjnl-2019-319188. 
#  Epub 2020 Jan 16. PMID: 31953253; PMCID: PMC7398469.
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
require(MetaboAnalystR) # See https://github.com/xia-lab/MetaboAnalystR

source("load_original_data/utils.R")

# --------------------------------
# Required files & info
# --------------------------------

# For details about the source of each file below see: <COMPLETE>
METADATA_FILE <- '../data/original_data/ERAWIJANTARI_GASTRIC_CANCER_2020/metadata.txt'  
TAXONOMY_FILE <- '../data/original_data/ERAWIJANTARI_GASTRIC_CANCER_2020/Raw_MetaPlhan2.tsv'
METABOLOMICS_FILE <- '../data/original_data/ERAWIJANTARI_GASTRIC_CANCER_2020/Raw_Metabolome.tsv'

PUBLICATION_NAME <- 'Influence of gastrectomy for gastric cancer treatment on faecal microbiome and metabolome profiles.'
DOI <- '10.1136/gutjnl-2019-319188'
DATASET_NAME <- 'ERAWIJANTARI_GASTRIC_CANCER_2020'

# --------------------------------
# Load metadata
# --------------------------------

metadata <- read_delim(METADATA_FILE, "\t", escape_double = FALSE, trim_ws = TRUE)

# Organize column names & order
metadata <- metadata %>%
  rename(Sample = 1) %>%
  mutate(Subject = Sample) %>%
  rename(Study.Group = Status) %>%
  mutate(Age.Units = "Years") %>%
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

# --------------------------------
# Load metabolomic profiles
# --------------------------------

mtb <- read_delim(METABOLOMICS_FILE, "\t", 
                  escape_double = FALSE, 
                  trim_ws = TRUE)
names(mtb)[1] <- "Compound"

# Create a mapping file to KEGG ID's, 
#  by parsing metabolite names in original table
mtb.map <- data.frame(Compound = mtb$Compound)

# Extract KEGG compound IDs when available
mtb.map$KEGG <- gsub("_.*", "", mtb.map$Compound) 

#  If missing a KEGG ID, put NA
mtb.map$KEGG[!grepl("C\\d{5}", mtb.map$KEGG)] <- NA

# Extract compound name
mtb.map$Compound.Name <- gsub("(^C\\d{5}_)|(^_)|(^\\-_)","",
                              mtb.map$Compound)

# Use MetaboAnalyst to retrieve HMDB ID's. 
# We will search by name but validate against  
#  the already available KEGG ID's
cmpds.to.search <- mtb.map$Compound.Name
MA.matches <- map.compound.names.MetaboAnalyst(cmpds.to.search)
MA.matches <- MA.matches %>% rename(KEGG.from.MA = KEGG) 

# Merge with the names in our data
mtb.map <- merge(mtb.map, MA.matches, 
                 by.x = "Compound.Name", 
                 by.y = "Query", all = T)

# Mark "lower confidence annotation" in cases where provided 
#  KEGG ID and MetaboAnalyst-idenfied KEGG ID do not match.
mtb.map <- mtb.map %>%
  mutate(High.Confidence.Annotation = 
           (is.na(KEGG) | is.na(KEGG.from.MA) | KEGG == KEGG.from.MA))

# Additional manual mappings based on HMDB searches
# (either by provided KEGG ID or by provided compound name)
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
mtb.map[!is.na(mtb.map$KEGG) & mtb.map$KEGG == 'C07091','HMDB'] <- 'HMDB0060461'
mtb.map[!is.na(mtb.map$KEGG) & mtb.map$KEGG == 'C00813','HMDB'] <- 'HMDB0041833'
mtb.map[!is.na(mtb.map$KEGG) & mtb.map$KEGG == 'C04483','HMDB'] <- 'HMDB0000626'
mtb.map[!is.na(mtb.map$KEGG) & mtb.map$KEGG == 'C00842','HMDB'] <- 'HMDB0001328'
mtb.map[!is.na(mtb.map$KEGG) & mtb.map$KEGG == 'C00134','HMDB'] <- 'HMDB0001414'
mtb.map[!is.na(mtb.map$KEGG) & mtb.map$KEGG == 'C00134',"High.Confidence.Annotation"] <- FALSE
mtb.map[mtb.map$Compound.Name == "1-Aminocyclopropane-1-carboxylate","HMDB"] <- 'HMDB0036458'
mtb.map[mtb.map$Compound.Name == "2',3'-cCMP","HMDB"] <- 'HMDB0011691'
mtb.map[mtb.map$Compound.Name == "3PG",'HMDB'] <- 'HMDB0000807'
mtb.map[mtb.map$Compound.Name == "3PG","High.Confidence.Annotation"] <- FALSE
mtb.map[!is.na(mtb.map$KEGG) & mtb.map$KEGG == 'C00117','HMDB'] <- 'HMDB0001548'
mtb.map[mtb.map$Compound.Name == "N1,N8-Diacetylspermidine",'HMDB'] <- 'HMDB0041947'
mtb.map[mtb.map$Compound.Name == "Ophthalmate",'HMDB'] <- 'HMDB0005765'

# Remove unneeded columns
mtb.map <- mtb.map %>%
  mutate(KEGG = coalesce(KEGG, KEGG.from.MA)) %>%
  select(-Compound.Name, -KEGG.from.MA, -MA.Name.Match)

# Mark cases of duplicated HMDN/KEGG ID as lower confidence
kegg.dups <- names(table(mtb.map$KEGG)[table(mtb.map$KEGG) > 1])
hmdb.dups <- names(table(mtb.map$HMDB)[table(mtb.map$HMDB) > 1])
mtb.map$High.Confidence.Annotation[mtb.map$KEGG %in% kegg.dups] <- FALSE
mtb.map$High.Confidence.Annotation[mtb.map$HMDB %in% hmdb.dups] <- FALSE

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

# --------------------------------
# Save to files + R objects
# --------------------------------

save.to.files(DATASET_NAME, metadata, mtb, mtb.map, genera, species)
save.to.rdata(DATASET_NAME, metadata, mtb, mtb.map, genera, species)
rm(list = ls())
