# ----------------------------------------------------------------
# The following script loads and processes data from:
# 
# Poyet, M., et al. "A library of human gut bacterial isolates 
#  paired with longitudinal multiomics data enables mechanistic 
#  microbiome research." Nature medicine 25.9 (2019): 1442-1452.
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
METADATA_FILE <- "../data/original_data/POYET_BIO_ML_2019/41591_2019_559_MOESM3_ESM.xlsx"
METADATA_FILE2 <- "../data/original_data/POYET_BIO_ML_2019/SraRunTable.txt"
TAXONOMY_FILE <- "../data/original_data/POYET_BIO_ML_2019/feature_table_gtdb.tsv"
METABOLOMICS_FILE1 <- "../data/original_data/POYET_BIO_ML_2019/HILIC_NEG_ION_MODE.txt"
METABOLOMICS_FILE2 <- "../data/original_data/POYET_BIO_ML_2019/HILIC_POS_ION_MODE.txt"
METABOLOMICS_FILE3 <- "../data/original_data/POYET_BIO_ML_2019/Rev_Phase_NEG_ION_MODE.txt"
METABOLOMICS_FILE4 <- "../data/original_data/POYET_BIO_ML_2019/Rev_Phase_POS_ION_MODE.txt"

PUBLICATION_NAME <- 'A library of human gut bacterial isolates paired with longitudinal multiomics data enables mechanistic microbiome research'
DOI <- '10.1038/s41591-019-0559-3'
DATASET_NAME <- 'POYET_BIO_ML_2019'

# --------------------------------
# Load metadata
# --------------------------------

# Taken from supp info. - metadata per donor (subject)
metadata <- read.xls(METADATA_FILE, 
                     sheet = 'SupTable2', 
                     header = TRUE,
                     skip = 3)
# Remove empty columns
metadata <- metadata[,colSums(is.na(metadata)) < nrow(metadata)]
metadata <- metadata %>% rename(Subject = Donor)

# We merge this with metadata per sample (there are multiple samples per donor), taken from SRA Run Selector.
SraRunTable <- read_csv(METADATA_FILE2, guess_max = 5000)

# Get only samples from hosts (not isolates), and only WGS ones
SraRunTable <- SraRunTable %>% filter(!is.na(host_subject_ID)) %>% filter(`Assay Type` == "AMPLICON")

# Take only relevant columns
SraRunTable <- SraRunTable %>% select(`Sample Name`, host_subject_ID, Run, Collection_Date)
names(SraRunTable)[1] <- "Sample"
names(SraRunTable)[2] <- "Subject"

# Merge
metadata <- merge(SraRunTable, metadata, by = "Subject", all.x = T, all.y = F)

# Organize column names & order
metadata <- metadata %>%
  rename(Gender = Sex) %>%
  mutate(Age.Units = "Years") %>%
  mutate(Sample = gsub("_16S","",Sample)) %>%
  # Add publication info
  mutate(Dataset = DATASET_NAME) %>%
  mutate(DOI = DOI) %>%
  mutate(Publication.Name = PUBLICATION_NAME) %>%
  # Reorder
  relocate(c(Dataset, Sample, Subject, Age, Age.Units, Gender, BMI, DOI, Publication.Name)) 

# --------------------------------
# Load taxonomic profiles (genera)
# --------------------------------

# Read genus-level abundances from qiime analysis
genera <- read_delim(TAXONOMY_FILE, 
                     "\t", escape_double = FALSE, trim_ws = TRUE)
names(genera)[1] <- 'Genus'

# --------------------------------
# Load metabolomic profiles
# --------------------------------

mtb1 <- read_delim(METABOLOMICS_FILE1, show_col_types = FALSE,
                   "\t", escape_double = FALSE, trim_ws = TRUE, comment = "Factors")
names(mtb1)[1] <- "Compound"
mtb1$Method <- "HILIC_NEG"
mtb2 <- read_delim(METABOLOMICS_FILE2, show_col_types = FALSE,
                   "\t", escape_double = FALSE, trim_ws = TRUE, comment = "Factors")
names(mtb2)[1] <- "Compound"
mtb2$Method <- "HILIC_POS"
mtb3 <- read_delim(METABOLOMICS_FILE3, 
                   "\t", escape_double = FALSE, trim_ws = TRUE, comment = "Factors")
names(mtb3)[1] <- "Compound"
mtb3$Method <- "REV_NEG"
mtb4 <- read_delim(METABOLOMICS_FILE4, 
                   "\t", escape_double = FALSE, trim_ws = TRUE, comment = "Factors")
names(mtb4)[1] <- "Compound"
mtb4$Method <- "REV_POS"

mtb <- bind_rows(mtb1, mtb2, mtb3, mtb4)
mtb$Compound.Name <- mtb$Compound
mtb$Compound <- paste(mtb$Method, mtb$Compound, sep = "_")
mtb$Method <- NULL
rm(mtb1, mtb2, mtb3, mtb4)

# Create metabolite mapping file
mtb.map <- mtb %>% select(Compound, Compound.Name)
mtb$Compound.Name <- NULL
mtb.map$High.Confidence.Annotation <- TRUE

# Mark ambiguous metabolites
mtb.map$High.Confidence.Annotation[grepl("[a-zA-Z]/", mtb.map$Compound.Name)] <- FALSE

# Slightly reformat metabolite names to get rid of characters that fail MetaboAnalyst 
mtb.map$Compound.Name <- gsub('_', ',', mtb.map$Compound.Name)

# Get mapping to KEGG/other compound ID's, using MetaboAnalyst 
cmpds.to.search <- unique(mtb.map$Compound.Name)
MA.matches <- map.compound.names.MetaboAnalyst(cmpds.to.search)

# Merge with the names in our data
mtb.map <- merge(mtb.map, MA.matches, 
                 by.x = "Compound.Name", 
                 by.y = "Query", 
                 all.x = T, all.y = F)
mtb.map$MA.Name.Match <- NULL
mtb.map <- mtb.map %>% relocate(Compound)

# Additional mappings
# View(mtb.map %>% filter(is.na(HMDB)))
mtb.map[mtb.map$Compound == 'HILIC_NEG_3-4-dihydroxybenzylamine','HMDB'] <- 'HMDB0012153'
mtb.map[mtb.map$Compound == 'HILIC_NEG_3-4-dihydroxyphenylglycol','HMDB'] <- 'HMDB0000318'
mtb.map[mtb.map$Compound == 'HILIC_NEG_3-4-dihydroxyphenylglycol','KEGG'] <- 'C05576'                      
mtb.map[mtb.map$Compound == 'HILIC_NEG_3-hydroxykynurenate','HMDB'] <- 'HMDB0000732'
mtb.map[mtb.map$Compound == 'HILIC_NEG_3-hydroxykynurenate','KEGG'] <- 'C02794'                           
mtb.map[mtb.map$Compound == 'REV_NEG_acetytyrosine','HMDB'] <- 'HMDB0000866'
mtb.map[mtb.map$Compound == 'REV_NEG_acetytyrosine','KEGG'] <- 'C01657'
mtb.map[mtb.map$Compound == 'HILIC_NEG_ascorbate','HMDB'] <- 'HMDB0000044'                                  
mtb.map[mtb.map$Compound == 'HILIC_POS_alanyalanine','HMDB'] <- 'HMDB0028680' 
mtb.map[mtb.map$Compound == 'HILIC_NEG_malonate','HMDB'] <- 'HMDB0000691'                                        
mtb.map[mtb.map$Compound == 'HILIC_POS_N-a-acetyarginine','HMDB'] <- 'HMDB0004620'                           
mtb.map[mtb.map$Compound == 'HILIC_POS_N-acetyphenylalanine','HMDB'] <- 'HMDB0000512'
mtb.map[mtb.map$Compound == 'HILIC_POS_N-acetyphenylalanine','KEGG'] <- 'C03519'                          
mtb.map[mtb.map$Compound == 'HILIC_POS_NMMA','HMDB'] <- 'HMDB0029416'
mtb.map[mtb.map$Compound == 'HILIC_POS_NMMA','KEGG'] <- 'C03884'                                       
mtb.map[mtb.map$Compound == 'REV_POS_palmithoylethanolamide','HMDB'] <- 'HMDB0002100'
mtb.map[mtb.map$Compound == 'REV_POS_palmithoylethanolamide','KEGG'] <- 'C16512'
mtb.map[mtb.map$Compound == 'HILIC_NEG_tetrahydro-1-methyl-beta-carboline-3-carboxylate','HMDB'] <- 'HMDB0037942'
mtb.map[mtb.map$Compound == 'HILIC_NEG_tetrahydro-1-methyl-beta-carboline-3-carboxylate','High.Confidence.Annotation'] <- FALSE
mtb.map[(!is.na(mtb.map$HMDB)) & mtb.map$HMDB == 'HMDB0062263','HMDB'] <- 'HMDB0000187'

# Mark cases of duplicated HMDN/KEGG ID as lower confidence 
# (These are mostly caused by overlap between MS methods)
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
