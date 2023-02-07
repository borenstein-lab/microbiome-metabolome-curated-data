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
source("data_organization/utils.R")

# --------------------------------
# Required files & info
# --------------------------------

# For details about the source of each file below see: <COMPLETE>
METADATA_FILE <- "../data/original_data/JACOBS_IBD_FAMILIES_2016/Pediatric_family_cohort_map.txt"
TAXONOMY_FILE <- "../data/original_data/JACOBS_IBD_FAMILIES_2016/feature_table_gtdb_207.tsv"
METABOLOMICS_FILE1 <- "../data/original_data/JACOBS_IBD_FAMILIES_2016/Metabolites_normalized_POS.txt"
METABOLOMICS_FILE2 <- "../data/original_data/JACOBS_IBD_FAMILIES_2016/Metabolites_normalized_NEG.txt"
METABOLOMICS_INFO <- "../data/original_data/JACOBS_IBD_FAMILIES_2016/Data file S1.xlsx"

PUBLICATION_NAME <- 'A Disease-Associated Microbial and Metabolomics State in Relatives of Pediatric Inflammatory Bowel Disease Patients'
DOI <- '10.1016/j.jcmgh.2016.06.004'
DATASET_NAME <- 'JACOBS_IBD_FAMILIES_2016'

# --------------------------------
# Load metadata
# --------------------------------

metadata <- read_delim(METADATA_FILE, 
                       "\t", escape_double = FALSE, 
                       show_col_types = FALSE,
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
                     show_col_types = FALSE,
                     trim_ws = TRUE)
names(genera)[1] <- 'Genus'

# --------------------------------
# Load metabolomic profiles
# --------------------------------

mtb_POS <- read_delim(METABOLOMICS_FILE1, 
                      "\t", escape_double = FALSE, 
                      show_col_types = FALSE,
                      trim_ws = TRUE)
mtb_NEG <- read_delim(METABOLOMICS_FILE2, 
                      "\t", escape_double = FALSE, 
                      show_col_types = FALSE,
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
mtb.map <- mtb.map %>%
  relocate(Compound) %>%
  rename(Compound.Name = Validated.ID)
mtb.map <- mtb.map[match(mtb$Compound, mtb.map$Compound),] # Reorder mapping file according to mtb - just for convenience

# Map the few (~50) validated IDs to their KEGG/HMDB IDs - using MetaboAnalyst
cmpds <- mtb.map$Compound.Name
cmpds <- cmpds[!is.na(cmpds)]
cmpds <- unique(cmpds)
MA.matches <- map.compound.names.MetaboAnalyst(cmpds)

# Merge with the names in our data
mtb.map <- merge(mtb.map, MA.matches, 
                 by.x = "Compound.Name", 
                 by.y = "Query", all = T)
mtb.map$MA.Name.Match <- NULL
mtb.map <- mtb.map %>% relocate(Compound)
mtb.map$High.Confidence.Annotation <- TRUE

# Manual completions 
# (mark those with different masses/retention times 
#  in same pos/neg mode as lower confidence, 
#  or those with ambiguous match to hmdb/kegg)
# View(mtb.map %>% select(Compound.Name, MA.Name.Match, HMDB.identification, KEGG.identification, HMDB, KEGG, m.z))
mtb.map$HMDB[which(mtb.map$Compound.Name == "11-Oxo-androsterone [glucuronide]")] <- "HMDB0010338"
mtb.map$KEGG[which(mtb.map$Compound.Name == "5-aminosalicylic acid")] <- "C07138"
mtb.map$High.Confidence.Annotation[which(mtb.map$Compound.Name == "5-aminosalicylic acid")] <- FALSE
mtb.map$KEGG[which(mtb.map$Compound.Name == "Acetyl-glutamic acid")] <- "C00624"
mtb.map$HMDB[which(mtb.map$Compound.Name == "Acetyl-glutamic acid")] <- "HMDB0001138"
mtb.map$High.Confidence.Annotation[which(mtb.map$Compound.Name == "Gamma-Glu-isoLeu")] <- FALSE
mtb.map$HMDB[which(mtb.map$Compound.Name == "Gamma-Glu-isoLeu")] <- "HMDB0011170"
mtb.map$KEGG[which(mtb.map$Compound.Name == "Gamma-Glu-isoLeu")] <- "C18135"
mtb.map$KEGG[which(mtb.map$Compound.Name == "Azaleic acid")] <- "C08261"
mtb.map$HMDB[which(mtb.map$Compound.Name == "Azaleic acid")] <- "HMDB0000784"
mtb.map$High.Confidence.Annotation[which(mtb.map$Compound.Name == "Serinyl tryptophan")] <- FALSE
mtb.map$HMDB[which(mtb.map$Compound.Name == "Serinyl tryptophan")] <- "HMDB0029050"
mtb.map$HMDB[which(mtb.map$Compound.Name == "N-acetyl-ornithine")] <- "HMDB0003357"
mtb.map$KEGG[which(mtb.map$Compound.Name == "N-acetyl-ornithine")] <- "C00437"

# Some discrepancy in retention time/m.z...
mtb.map$High.Confidence.Annotation[which(mtb.map$Compound.Name == "7-Ketodeoxycholic acid")] <- FALSE
mtb.map$High.Confidence.Annotation[which(mtb.map$Compound.Name == "Chenodeoxycholic acid sulfate")] <- FALSE
mtb.map$High.Confidence.Annotation[which(mtb.map$Compound.Name == "N-Acetylcadaverine")] <- FALSE
mtb.map$High.Confidence.Annotation[which(mtb.map$Compound.Name == "Taurochenodeoxycholate")] <- FALSE
mtb.map$High.Confidence.Annotation[which(mtb.map$Compound.Name == "Taurine")] <- FALSE

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

