# ----------------------------------------------------------------
# The following script loads and processes data from:
# 
# Kang, Dae-Wook, et al. "Differences in fecal microbial 
#  metabolites and microbiota of children with autism spectrum 
#  disorders." Anaerobe 49 (2018): 121-131.
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
source("load_original_data/utils.R")

# --------------------------------
# Required files & info
# --------------------------------

# For details about the source of each file below see: <COMPLETE>
METADATA_FILE <- "~/MICROBIOME_METABOLOME/KANG_AUTISM_2017/KANG_AUTISM_2017_QIIME/combined.metadata.tsv"
TAXONOMY_FILE <- "~/MICROBIOME_METABOLOME/KANG_AUTISM_2017/KANG_AUTISM_2017_QIIME/6_Taxonomy/taxa_genus_table_exported/feature-table.tsv"
METABOLOMICS_FILE <- "C:/Users/efrat/Documents/MICROBIOME_METABOLOME/KANG_AUTISM_2017/RAW_DATA/1-s2.0-S1075996417302305-mmc2.xlsx"

PUBLICATION_NAME <- 'Differences in fecal microbial metabolites and microbiota of children with autism spectrum disorders'
DOI <- '10.1016/j.anaerobe.2017.12.007'
DATASET_NAME <- 'KANG_AUTISM_2017'

# --------------------------------
# Load metadata
# --------------------------------

metadata <- read.delim(METADATA_FILE, 
                       comment.char="#", 
                       stringsAsFactors=FALSE)

# Organize column names & order
metadata <- metadata %>%
  rename(RawSampleID = 1) %>%
  rename(Sample = 2) %>%
  mutate(Subject = Sample) %>%
  rename(Study.Group = Group) %>%
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

# Mapping ID's to those used in paper
patient.IDs <- metadata$Sample
names(patient.IDs) <- metadata$RawSampleID
genera <- plyr::rename(genera, patient.IDs, warn_missing = FALSE)
rm(patient.IDs)

# --------------------------------
# Load metabolomic profiles
# --------------------------------

# Read metabolic profiles
mtb <- read.xls(METABOLOMICS_FILE, 
                sheet = 'metabolites', 
                header = TRUE,
                pattern = 'Sample',
                stringsAsFactors = FALSE)
names(mtb)[1] <- 'Compound'

# Manually fix metabolite with special chars
mtb$Compound[54] <- "beta-Alanine"

mtb.map <- mtb %>% select(Compound)

# Get mapping to KEGG/other compound ID's, using MetaboAnalyst 
cmpds.to.search <- mtb.map$Compound
cmpds.to.search <- unique(cmpds.to.search)
MA.matches <- map.compound.names.MetaboAnalyst(cmpds.to.search)

# Merge with the names in our data
mtb.map <- merge(mtb.map, MA.matches, 
                 by.x = "Compound", 
                 by.y = "Query", all = T)
mtb.map$MA.Name.Match <- NULL
mtb.map$High.Confidence.Annotation <- TRUE

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

save.to.files(DATASET_NAME, metadata, mtb, mtb.map, genera)
save.to.rdata(DATASET_NAME, metadata, mtb, mtb.map, genera)
rm(list = ls())
