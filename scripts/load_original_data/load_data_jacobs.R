# ----------------------------------------------------------------
# The following script loads and processes data from:
# 
# Jacobs, Jonathan P., et al. "A disease-associated microbial 
#  and metabolomics state in relatives of pediatric inflammatory 
#  bowel disease patients." Cellular and molecular 
#  gastroenterology and hepatology 2.6 (2016): 750-766.
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
METADATA_FILE <- "~/MICROBIOME_METABOLOME/JACOBS_IBD_FAMILIES_2016/RAW_DATA/Pediatric_family_cohort_map.txt"
TAXONOMY_FILE <- "~/MICROBIOME_METABOLOME/JACOBS_IBD_FAMILIES_2016/JACOBS_IBD_FAMILIES_2016_QIIME/6_Taxonomy/taxa_genus_table_exported/feature-table.tsv"
METABOLOMICS_FILE1 <- "~/MICROBIOME_METABOLOME/JACOBS_IBD_FAMILIES_2016/RAW_DATA/Metabolites_normalized_POS.txt"
METABOLOMICS_FILE2 <- "~/MICROBIOME_METABOLOME/JACOBS_IBD_FAMILIES_2016/RAW_DATA/Metabolites_normalized_NEG.txt"
METABOLOMICS_INFO <- "C:/Users/efrat/Documents/MICROBIOME_METABOLOME/JACOBS_IBD_FAMILIES_2016/RAW_DATA/Data file S1.xlsx"

PUBLICATION_NAME <- 'A Disease-Associated Microbial and Metabolomics State in Relatives of Pediatric Inflammatory Bowel Disease Patients'
DOI <- '10.1016/j.jcmgh.2016.06.004'
DATASET_NAME <- 'JACOBS_IBD_FAMILIES_2016'

# --------------------------------
# Load metadata
# --------------------------------

metadata <- read_delim(METADATA_FILE, 
                       "\t", escape_double = FALSE, 
                       trim_ws = TRUE)

# Organize column names & order
metadata <- metadata %>%
  rename(Sample = 1) %>%
  mutate(Subject = Sample) %>%
  rename(Study.Group = IBD_Status) %>%
  mutate(Age.Units = "Years") %>%
  mutate(Gender = ifelse(Gender=="F", "Female", "Male")) %>%
  # Add publication info
  mutate(Dataset = DATASET_NAME) %>%
  mutate(DOI = DOI) %>%
  mutate(Publication.Name = PUBLICATION_NAME) %>%
  # Reorder
  relocate(c(Dataset, Sample, Subject, Study.Group, Age, Age.Units, Gender, DOI, Publication.Name)) 

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

mtb_POS <- read_delim(METABOLOMICS_FILE1, 
                      "\t", escape_double = FALSE, 
                      trim_ws = TRUE)
mtb_NEG <- read_delim(METABOLOMICS_FILE2, 
                      "\t", escape_double = FALSE, 
                      trim_ws = TRUE)
mtb <- bind_rows("Positive" = mtb_POS, "Negative" = mtb_NEG, .id = "POS.NEG")
rm(mtb_NEG, mtb_POS)
names(mtb)[2] <- "m.z..r.t"

# Read mapping table from excel
mtb.map <- read.xls(METABOLOMICS_INFO, 
                    sheet = 1,
                    header = TRUE, 
                    stringsAsFactors = FALSE, 
                    na.strings = c(""))

# We create a single metabolite identifier "Compound": <mode>_<m/z>_<retention time>, in both tables
mtb.map$Compound <- paste(mtb.map$ESI.mode, mtb.map$m.z, mtb.map$Retention.time, sep = "_")
mtb$Compound <- paste(mtb$POS.NEG, mtb$m.z..r.t, sep = "_")

# Clean mtb table
mtb <- mtb %>%
  relocate(Compound) %>%
  select(-POS.NEG, -m.z..r.t)

# Reorganize mtb.map table
mtb.map <- mtb.map[,c("Compound","Validated.ID","m.z","KEGG.identification","HMDB.identification")]
names(mtb.map)[2] <- "Compound.Name"
mtb.map <- mtb.map[match(mtb$Compound, mtb.map$Compound),] # Reorder mapping file according to mtb - just for convenience




High.Confidence.Annotation

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

save.to.files(DATASET_NAME, metadata, mtb, mtb.map, genera///, species)
save.to.rdata(DATASET_NAME, metadata, mtb, mtb.map, genera///, species)
rm(list = ls())








# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Previous code


# Map the few (~50) validated IDs to their KEGG/HMDB IDs - using MetaboAnalyst
mSet <- InitDataObjects("NA", "utils", FALSE) # ignore warning
cmpds <- mtb.map[,"Compound.Name"]
cmpds <- cmpds[!is.na(cmpds)]
cmpds <- unique(cmpds) # The validated metabolites are non-unique...
mSet <- Setup.MapData(mSet, cmpds)
mSet <- CrossReferencing(mSet, "name", metlin = F, pubchem = F, chebi = F)
mSet <- CreateMappingResultTable(mSet)
match.table <- data.table(mSet$dataSet$map.table)
rm(mSet, cmpds)
match.table <- match.table[,c("Query","HMDB","KEGG")]
match.table[match.table== "NA"] <- NA
match.table[match.table== ""] <- NA

# Add mappings to mtb.map
mtb.map <- merge(mtb.map, match.table, by.x = "Compound.Name", by.y = "Query", all = T)
mtb.map$High.Confidence <- TRUE

# Manual completions (mark those with different masses/retention times in same pos/neg mode as lower confidence, or those with unclear match to hmdb/kegg)
# View(mtb.map[,c("Compound.Name","Compound","KEGG.identification","KEGG")])
mtb.map[!is.na(mtb.map$Compound.Name) & mtb.map$Compound.Name == "5-aminosalicylic acid","KEGG"] <- 'C07138'
mtb.map[!is.na(mtb.map$Compound.Name) & mtb.map$Compound.Name == "5-aminosalicylic acid","High.Confidence"] <- FALSE
mtb.map[!is.na(mtb.map$Compound.Name) & mtb.map$Compound.Name == "Acetyl-glutamic acid","KEGG"] <- 'C00624'
mtb.map[!is.na(mtb.map$Compound.Name) & mtb.map$Compound.Name == "Acetyl-glutamic acid","HMDB"] <- 'HMDB0001138'
mtb.map[!is.na(mtb.map$Compound.Name) & mtb.map$Compound.Name == "Azaleic acid","KEGG"] <- 'C08261'
mtb.map[!is.na(mtb.map$Compound.Name) & mtb.map$Compound.Name == "Azaleic acid","HMDB"] <- 'HMDB0000784'
mtb.map[!is.na(mtb.map$Compound.Name) & mtb.map$Compound.Name == "Gamma-Glu-isoLeu","HMDB"] <- 'HMDB0011170'
mtb.map[!is.na(mtb.map$Compound.Name) & mtb.map$Compound.Name == "Gamma-Glu-isoLeu","High.Confidence"] <- FALSE
mtb.map[!is.na(mtb.map$Compound.Name) & mtb.map$Compound.Name == "N-acetyl-ornithine","KEGG"] <- 'C00437'
mtb.map[!is.na(mtb.map$Compound.Name) & mtb.map$Compound.Name == "N-acetyl-ornithine","HMDB"] <- 'HMDB0003357'
mtb.map[!is.na(mtb.map$Compound.Name) & mtb.map$Compound.Name == "Serinyl tryptophan","HMDB"] <- 'HMDB0029050'
mtb.map[!is.na(mtb.map$Compound.Name) & mtb.map$Compound.Name == "Serinyl tryptophan","High.Confidence"] <- FALSE # Mass doesn't match, name does
# Same mode different retention times
mtb.map[!is.na(mtb.map$Compound.Name) & mtb.map$Compound.Name == "7-Ketodeoxycholic acid","High.Confidence"] <- FALSE
mtb.map[!is.na(mtb.map$Compound.Name) & mtb.map$Compound.Name == "Chenodeoxycholic acid sulfate","High.Confidence"] <- FALSE
mtb.map[!is.na(mtb.map$Compound.Name) & mtb.map$Compound.Name == "N-Acetylcadaverine","High.Confidence"] <- FALSE
mtb.map[!is.na(mtb.map$Compound.Name) & mtb.map$Compound.Name == "Taurochenodeoxycholate","High.Confidence"] <- FALSE
mtb.map[!is.na(mtb.map$Compound.Name) & mtb.map$Compound.Name == "Taurine","High.Confidence"] <- FALSE


