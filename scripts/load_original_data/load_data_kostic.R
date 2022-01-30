# ----------------------------------------------------------------
# The following script loads and processes data from:
# 
# Kostic, Aleksandar D., et al. "The dynamics of the human 
#  infant gut microbiome in development and in progression toward 
#  type 1 diabetes." Cell host & microbe 17.2 (2015): 260-273.
# 
# The output tables are:
# 1) metadata - subject metadata
# 2) genera - taxonomic profiles per subject (genus-level)
# 3) mtb - metabolomic profiles per subject       
# 4) mtb.map - additional details per metabolite 
# ----------------------------------------------------------------

require(gdata)
require(readr)
require(dplyr)
require(MetaboAnalystR)
source("load_original_data/utils.R")

# --------------------------------
# Required files & info
# --------------------------------

# For details about the source of each file below see: <COMPLETE>
METADATA_FILE <- "../data/original_data/KOSTIC_INFANTS_DIABETES_2015/diabimmune_t1d_16s_metadata.xlsx"
METADATA_FILE2 <-"../data/original_data/KOSTIC_INFANTS_DIABETES_2015/diabimmune_sample_id_mapping.txt"
TAXONOMY_FILE <- "../data/original_data/KOSTIC_INFANTS_DIABETES_2015/feature-table.tsv"
METABOLOMICS_FILE <- "../data/original_data/KOSTIC_INFANTS_DIABETES_2015/MxP_DIABIMMUNE.xlsx"

PUBLICATION_NAME <- 'The Dynamics of the Human Infant Gut Microbiome in Development and in Progression towards Type 1 Diabetes'
DOI <- '10.1016/j.chom.2015.01.001'
DATASET_NAME <- 'KOSTIC_INFANTS_DIABETES_2015'

# --------------------------------
# Load metadata
# --------------------------------

metadata <- read.xls(METADATA_FILE, 
                     sheet = 1,
                     header = TRUE, 
                     stringsAsFactors = FALSE)

# Organize column names & order
metadata <- metadata %>%
  rename(Sample = 1) %>%
  rename(Subject = Subject_ID) %>%
  rename(Study.Group = Case_Control) %>%
  rename(Age = Age_at_Collection) %>%
  mutate(Age.Units = "Days") %>%
  mutate(Gender = ifelse(Gender == "male", "Male", "Female")) %>%
  # Add publication info
  mutate(Dataset = DATASET_NAME) %>%
  mutate(DOI = DOI) %>%
  mutate(Publication.Name = PUBLICATION_NAME) %>%
  # Reorder
  relocate(c(Dataset, Sample, Subject, Study.Group, Age, 
             Age.Units, Gender, DOI, Publication.Name)) 


# --------------------------------
# Load taxonomic profiles (genera)
# --------------------------------

# Read genus-level abundances from qiime analysis
genera <- read_delim(TAXONOMY_FILE, 
                     "\t", escape_double = FALSE, 
                     trim_ws = TRUE, 
                     skip = 1)
names(genera)[1] <- 'Genus'

# --------------------------------
# Load metabolomic profiles
# --------------------------------

mtb1 <- read.xls(METABOLOMICS_FILE, 
                sheet = 2, 
                header = TRUE, 
                stringsAsFactors = FALSE)
mtb2 <- read.xls(METABOLOMICS_FILE, 
                sheet = 3, 
                header = TRUE, 
                stringsAsFactors = FALSE)
mtb2 <- mtb2 %>% 
  rename(Metabolite = Accepted.Compound.ID) %>%
  mutate(Compound = as.character(Compound))
mtb3 <- read.xls(METABOLOMICS_FILE, 
                sheet = 4, 
                header = TRUE, 
                stringsAsFactors = FALSE)
mtb <- bind_rows(mtb1,mtb2,mtb3)
rm(mtb1,mtb2,mtb3)

# Start building a mapping file with info per metabolite
mtb.map <- mtb %>%
  select(Metabolite, Compound, m.z, RT..min.) %>%
  rename(Compound.Original.ID = Compound) %>%
  rename(Compound = Metabolite) %>%
  mutate(High.Confidence.Annotation = TRUE)

mtb <- mtb %>%
  select(-Compound, -m.z, -RT..min.) %>%
  rename(Compound = Metabolite) 

# Remove QC samples
mtb[,grep("Pool", names(mtb))] <- NULL

# Mapping ID's
diabimmune_sample_id_mapping <- read.delim(METADATA_FILE2)
patient.IDs <- as.character(diabimmune_sample_id_mapping$G_Project)
names(patient.IDs) <- as.character(gsub("-","\\.",diabimmune_sample_id_mapping$Root_Samples))
mtb <- plyr::rename(mtb, patient.IDs, warn_missing = FALSE)
rm(patient.IDs)

# Mark ambiguous metabolites
mtb.map$High.Confidence.Annotation[grepl("[a-zA-Z]/", mtb.map$Compound)] <- FALSE

# Get mapping to KEGG/other compound ID's, using MetaboAnalyst 
cmpds.to.search <- mtb.map$Compound
MA.matches <- map.compound.names.MetaboAnalyst(cmpds.to.search)

# Merge with the names in our data
mtb.map <- merge(mtb.map, MA.matches, 
                 by.x = "Compound", 
                 by.y = "Query", all = T)
mtb.map$MA.Name.Match <- NULL

# Few manual mappings
# View(mtb.map %>% filter(is.na(HMDB) | is.na(KEGG)) %>% filter(!grepl("[a-zA-Z]/",Compound)))
mtb.map[mtb.map$Compound == "alpha-glycerophosphocholine","KEGG"] <- 'C00670'
mtb.map[mtb.map$Compound == "alpha-glycerophosphocholine","HMDB"] <- 'HMDB0000086'
mtb.map[mtb.map$Compound == "NMMA","KEGG"] <- 'C03884'
mtb.map[mtb.map$Compound == "NMMA","HMDB"] <- 'HMDB0029416'
mtb.map[mtb.map$Compound == "8.11.14-Eicosatrienoic acid","KEGG"] <- 'C03242'
mtb.map[mtb.map$Compound == "8.11.14-Eicosatrienoic acid","HMDB"] <- 'HMDB0002925'

# Mark cases of duplicated HMDN/KEGG ID as lower confidence
kegg.dups <- names(table(mtb.map$KEGG)[table(mtb.map$KEGG) > 1])
hmdb.dups <- names(table(mtb.map$HMDB)[table(mtb.map$HMDB) > 1])
mtb.map$High.Confidence.Annotation[mtb.map$KEGG %in% kegg.dups] <- FALSE
mtb.map$High.Confidence.Annotation[mtb.map$HMDB %in% hmdb.dups] <- FALSE

# --------------------------------
# Keep only samples with all data
# --------------------------------

sample.intersect <- Reduce(intersect, list(names(mtb)[-1], names(genera)[-1], metadata$Sample))
message(paste(length(sample.intersect),"samples have all data types"))

mtb <- mtb[,c("Compound",sample.intersect)]
genera <- genera[,c("Genus",sample.intersect)]
metadata <- metadata[metadata$Sample %in% sample.intersect,]

# --------------------------------
# Save to files + R objects
# --------------------------------

save.to.files(DATASET_NAME, "prelim_data", metadata, mtb, mtb.map, genera)
save.to.rdata(DATASET_NAME, "prelim_data", metadata, mtb, mtb.map, genera)
rm(list = ls())

