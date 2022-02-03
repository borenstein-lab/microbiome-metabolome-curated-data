# ----------------------------------------------------------------
# The following script loads and processes data from:
# 
# He, Xuan, et al. "Fecal microbiome and metabolome of infants 
#  fed bovine MFGM supplemented formula or standard formula with 
#  breast-fed infants as reference: a randomized controlled trial." 
#  Scientific reports 9.1 (2019): 1-14.
# 
# The output tables are:
# 1) metadata - subject metadata
# 2) genera - taxonomic profiles per subject (genus-level)
# 3) mtb - metabolomic profiles per subject       
# 4) mtb.map - additional details per metabolite 
# ----------------------------------------------------------------

require(readr)
require(dplyr)
source("load_original_data/utils.R")

# --------------------------------
# Required files & info
# --------------------------------

# For details about the source of each file below see: <COMPLETE>
METADATA_FILE <- "../data/original_data/HE_INFANTS_MFGM_2019/63189_mapping_file__small.txt"
TAXONOMY_FILE <- "../data/original_data/HE_INFANTS_MFGM_2019/feature_table_gtdb.tsv"
METABOLOMICS_FILE <- "../data/original_data/HE_INFANTS_MFGM_2019/Tumme stool concentration 4.11.2017 final uM.csv"

PUBLICATION_NAME <- 'Fecal microbiome and metabolome of infants fed bovine MFGM supplemented formula or standard formula with breast-fed infants as reference: a randomized controlled trial'
DOI <- '10.1038/s41598-019-47953-4'
DATASET_NAME <- 'HE_INFANTS_MFGM_2019'

# --------------------------------
# Load metadata
# --------------------------------

metadata <- read_delim(METADATA_FILE, 
                       "\t", escape_double = FALSE, 
                       col_types = cols(gender = col_character(), 
                                        host_age = col_integer()), 
                       trim_ws = TRUE)

# Organize column names & order
metadata <- metadata %>%
  rename(Sample = 1) %>%
  mutate(Sample = gsub("12021\\.","",Sample)) %>% # Reformat sample names, to match other tables
  rename(Subject = host_subject_id) %>%
  rename(Study.Group = group_name) %>%
  mutate(Age.Units = "Months") %>%
  select(-host_age_units) %>%
  rename(Age = host_age) %>%
  rename(Gender = gender) %>%
  mutate(Gender = ifelse(Gender=="M","Male","Female")) %>%
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
                     trim_ws = TRUE)
names(genera)[1] <- 'Genus'

# Map ID's to uniform format
names(genera) <- gsub("12021\\.","",names(genera)) 

# --------------------------------
# Load metabolomic profiles
# --------------------------------

mtb <- read_csv(METABOLOMICS_FILE, 
                col_types = cols(`Bristol score` = col_skip(), 
                                 Diet = col_skip(), 
                                 Gender = col_skip(), 
                                 H2O.Perc = col_skip(),  
                                 Month = col_skip(), 
                                 Treatment = col_skip(), 
                                 cnx = col_skip()))

# The following processing steps are taken from the R notebook sent by authors (slightly adjusted, all are mentioned in paper)
## Rename some metabolites:
names(mtb)[names(mtb) == "Ethyl-b-D-glucuronide"]<- "Ethyl-beta-D-glucuronide"
names(mtb)[names(mtb) == "N-Carbamoyl-b-alanine"]<- "N-Carbamoyl-beta-alanine"
names(mtb)[names(mtb) == "b-Alanine"]<- 'beta-Alanine'
names(mtb)[names(mtb) == "b-Pseudouridine"]<- 'beta-Pseudouridine'

## Remove outliers ("BF 35-2 dry weight lost in the freeze drying process")
outlier.list<-c("BF 35-2")

## Exclude these samples/subjects
## From paper: "After data generation, a few samples were further excluded, including infants who stopped consumption of study formula (n=2), infants who had no record of a food diary (n=1), infants from the BF group that were heavily mixed-fed (n=3), infants from the FF group who consumed another formula (n=4), and infants who had a record of antibiotics use (n=3). For fecal metabolome data, samples that contained urea were suspected of being contaminated by urine since normal fecal samples do not contain urea due to bacterial urease activity. Therefore, fecal samples containing urea were excluded from analysis."
exclude.sample <- c("EP 52-2", "EP 52-3", "EP 52-4", "EP 65-2", "EP 65-3", "EP 65-4", "EF 79-2", "EF 79-3", "EF 79-4")
exclude.subject <- c("EF 52", "BP 15", "BP 26", "BP 27", "EF 23", "EF 38", "EF 39", "EP 7", "EP 15", "BP 1")

mtb<-mtb[!mtb$StudyID %in% outlier.list,]
mtb<-mtb[!mtb$StudyID %in% exclude.sample,]
mtb<-mtb[!mtb$InfantID %in% exclude.subject,]

## Remove those with urea contaimination (stool collected from diapers)
mtb<-mtb[mtb$Urea==0,] 

## From author script: "Fix those FF samples that have HMOs, but are really negligible when double check the cnx file." (?)
mtb$`3' fucosyl lactose`[mtb$StudyID == "EF 7-1"] <- 0
mtb$`3'-sialyllactose`[mtb$StudyID == "EF 7-1"] <- 0
mtb$`3' fucosyl lactose`[mtb$StudyID == "EF 41-2"] <- 0
mtb$`2' fucosyl lactose`[mtb$StudyID == "EP 49-2"] <- 0
mtb$`3' fucosyl lactose`[mtb$StudyID == "EF 34-3"] <- 0
mtb$`3' fucosyl lactose`[mtb$StudyID == "EP 51-3"] <- 0
mtb$`3' fucosyl lactose`[mtb$StudyID == "EF 12-4"] <- 0

mtb$InfantID <- NULL

# Transpose
x <- mtb %>% tibble::column_to_rownames("StudyID") 
x <- data.frame(t(x))
x <- x %>% tibble::rownames_to_column("Compound")
mtb <- x

# Create metabolite mapping table
mtb.map <- mtb %>%
  select(Compound) %>%
  mutate(High.Confidence.Annotation = TRUE)

# Get mapping to KEGG/other compound ID's, using MetaboAnalyst 
cmpds.to.search <- mtb.map$Compound
MA.matches <- map.compound.names.MetaboAnalyst(cmpds.to.search)

# Merge with the names in our data
mtb.map <- merge(mtb.map, MA.matches, 
                 by.x = "Compound", 
                 by.y = "Query", all = T)
mtb.map$MA.Name.Match <- NULL

# Few manual corrections/completions
# View(mtb.map %>% filter(is.na(HMDB) | is.na(KEGG)))
mtb.map[mtb.map$Compound == "DL-methionine-sulfoxide","KEGG"] <- 'C02989'
mtb.map[mtb.map$Compound == "Ethyl-beta-D-glucuronide","HMDB"] <- "HMDB0010325" 
mtb.map[mtb.map$Compound == "Methylsuccinate","KEGG"] <- "C02614" # The ID given by MetaboAnalyst does not appear in KEGG DB anymore (C08645)
mtb.map[mtb.map$Compound == "2' fucosyl lactose","HMDB"] <- "HMDB0002098" 
mtb.map[mtb.map$Compound == "2-Hydroxyisobutyrate","KEGG"] <- "C21297" 
mtb.map[mtb.map$Compound == "3' fucosyl lactose","HMDB"] <- "HMDB0002094" 
mtb.map[mtb.map$Compound == "6' sialyllactose","HMDB"] <- "HMDB0006569" 
mtb.map[mtb.map$Compound == "N-Acetylneuraminate","HMDB"] <- "HMDB0000230" 

# Add some manual fixes of sample names (TODO: verify with author + ask why these samples are missing from metadata: "BP.4.1","EF.15.1","EF.51.1","EP.77.4" )
names(mtb)[names(mtb) == "EF.3.1"] <- "EF.3.1a"
names(mtb)[names(mtb) == "BF.9.1"] <- "BF.9.1b"
names(mtb)[names(mtb) == "EP.51.2"] <- "EP.51.2a"
names(mtb)[names(mtb) == "BP.4.3"] <- "BP.4.3b"
names(mtb)[names(mtb) == "BP.40.3"] <- "BF.40.3b"

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
