# ----------------------------------------------------------------
# The following script loads and processes data from:
# 
# Wang, Xifan, et al. "Aberrant gut microbiota alters host 
#  metabolome and impacts renal failure in humans and rodents." 
#  Gut 69.12 (2020): 2131-2142.
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
require(gdata)

source("data_organization/utils.R")
source("data_organization/gtdb_utils.R")

# --------------------------------
# Required files & info
# --------------------------------

# For details about the source of each file below see: <COMPLETE>
METADATA_FILE <- "../data/original_data/WANG_ESRD_2020/MetaboLights/s_MTBLS700.txt" # I read the metadata table from S2-1 which is the metadata of subjects for which a metagenome sample was taken. 
TAXONOMY_FILE_SP <- '../data/original_data/WANG_ESRD_2020/kraken/species_level_taxonomy.tsv'
TAXONOMY_FILE_GE <- '../data/original_data/WANG_ESRD_2020/kraken/genus_level_taxonomy.tsv'
TAXONOMY_SAMPLE_MAP <- '../data/original_data/WANG_ESRD_2020/kraken/SraRunTable.txt'
METABOLOMICS_FILE <- "../data/original_data/WANG_ESRD_2020/MetaboLights/Fecal metabolite analysis.xlsx"
PUBLICATION_NAME <- 'Aberrant gut microbiota alters host metabolome and impacts renal failure in humans and rodents'
DOI <- '10.1136/gutjnl-2019-319766'
DATASET_NAME <- 'WANG_ESRD_2020'

# --------------------------------
# Load metadata
# --------------------------------

metadata <- read_delim(METADATA_FILE, 
                       delim = "\t", 
                       escape_double = FALSE, 
                       trim_ws = TRUE,
                       show_col_types = FALSE)

# Organize column names & order
metadata <- metadata %>%
  filter(`Characteristics[Organism part]` == "feces") %>%
  select(`Sample Name`,
         `Factor Value[ESRD status]`,
         `Factor Value[Age]`,
         `Factor Value[Gender]`,
         `Factor Value[BMI]`,
         `Source Name`,
         `Factor Value[Creatinine]`,
         `Factor Value[Dialysis frequency]`,
         `Factor Value[Urea]`,
         `Factor Value[eGFR]`) %>%
  rename(Sample = `Sample Name`) %>%
  mutate(Subject = gsub("\\-(F|Fb)","",Sample)) %>%
  rename(Study.Group = `Factor Value[ESRD status]`) %>%
  mutate(Study.Group = ifelse(Study.Group == 0, "Control", "ESRD")) %>%
  rename(Age = `Factor Value[Age]`) %>%
  mutate(Age.Units = "Years") %>%
  rename(Gender = `Factor Value[Gender]`) %>%
  mutate(Gender = ifelse(Gender == "female", "Female", "Male")) %>%
  # Add publication info
  mutate(Dataset = DATASET_NAME) %>%
  mutate(DOI = DOI) %>%
  mutate(Publication.Name = PUBLICATION_NAME) %>%
  rename(BMI = `Factor Value[BMI]`) %>%
  rename(eGFR = `Factor Value[eGFR]`) %>%
  rename(Urea = `Factor Value[Urea]`) %>%
  rename(Creatinine = `Factor Value[Creatinine]`) %>%
  rename(Dialysis.Frequency = `Factor Value[Dialysis frequency]`) %>%
  # Reorder
  relocate(c(Dataset, Sample, Subject, Study.Group, Age, Age.Units, Gender, BMI, DOI, Publication.Name)) 

# --------------------------------
# Load taxonomic profiles (species)
# --------------------------------

species <- read_delim(TAXONOMY_FILE_SP, "\t", 
                      escape_double = FALSE, 
                      trim_ws = TRUE,
                      show_col_types = FALSE)
names(species)[1] <- 'Species'
names(species) <- gsub("_Gut_microbiota.*","",names(species))

genera <- read_delim(TAXONOMY_FILE_GE, "\t", 
                     escape_double = FALSE, 
                     trim_ws = TRUE,
                     show_col_types = FALSE)
names(genera)[1] <- 'Genus'
names(genera) <- gsub("_Gut_microbiota.*","",names(genera))

# Prepare mapping from file names to sample IDs (as in metadata)
tax.map <- read_csv(TAXONOMY_SAMPLE_MAP, 
                    col_select = c("Run", "Sample Name"),
                    show_col_types = FALSE) %>% 
  rename(Sample = `Sample Name`) %>%
  filter(Sample %in% metadata$Subject) # The 3 samples that drop here do not have metabolomic data anyway

# Sanity: table(names(species)[-1] %in% tax.map$Run)
tax.map.vec <- c("Species", "Genus", tax.map$Sample)
names(tax.map.vec) <- c("Species", "Genus", tax.map$Run)

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

mtb.raw <- read.xls(METABOLOMICS_FILE,
               header = TRUE,
               na.strings = c("NA",""), 
               check.names = FALSE,
               stringsAsFactors = FALSE)
# Transpose
mtb <- mtb.raw %>%
  tibble::column_to_rownames("#SampleID") %>%
  t() %>%
  data.frame() %>%
  tibble::rownames_to_column("Compound")

# Fix sample names
names(mtb)[-1] <- mtb.raw$`#SampleID`

# We'll now use MetaboAnalyst to get kegg + hmdb id's from metabolite names
cmpds.to.search <- mtb$Compound
MA.matches <- map.compound.names.MetaboAnalyst(cmpds.to.search)

# Create mapping file
mtb.map <- MA.matches %>%
  rename(Compound = Query) %>% 
  mutate(High.Confidence.Annotation = TRUE)

# Additional mappings (based on searching in kegg/hmdb)
# View(mtb.map %>% filter(!grepl("unknown", Compound)) %>% filter(is.na(KEGG) | is.na(HMDB)))
mtb.map[mtb.map$Compound == "Disiloxane, hexamethyl-","HMDB"] <- 'HMDB0033532'
mtb.map[mtb.map$Compound == "Hexadecanal","KEGG"] <- 'C00517'
mtb.map[mtb.map$Compound == "Hexadecanal","HMDB"] <- 'HMDB0001551'
mtb.map[mtb.map$Compound == "Nonanal","HMDB"] <- 'HMDB0059835'
mtb.map[mtb.map$Compound == "Furan, 2-pentyl-","HMDB"] <- 'HMDB0013824'
mtb.map[mtb.map$Compound == "Ethanone, 1-(2-aminophenyl)-","HMDB"] <- 'HMDB0032630'
mtb.map[mtb.map$Compound == "Benzyl methyl ketone","KEGG"] <- 'C15512'
mtb.map[mtb.map$Compound == "Butanal, 2-methyl-","KEGG"] <- 'C02223'
mtb.map[mtb.map$Compound == "Butanal, 2-methyl-","HMDB"] <- 'HMDB0031526'
mtb.map[mtb.map$Compound == "Butanoic acid, 2-methyl-","KEGG"] <- 'C18319'
mtb.map[mtb.map$Compound == "Butanoic acid, 3-methyl-","HMDB"] <- 'HMDB0000718'
mtb.map[mtb.map$Compound == "Butanoic acid, 3-methyl-","KEGG"] <- 'C08262'
mtb.map[mtb.map$Compound == "Anethole","HMDB"] <- 'HMDB0030837'
mtb.map[mtb.map$Compound == "3-Pentanol","HMDB"] <- 'HMDB0303831'
mtb.map[mtb.map$Compound == "Isobutyl isovalerate","HMDB"] <- 'HMDB0038040'
mtb.map[mtb.map$Compound == "Methylamine, N,N-dimethyl-","HMDB"] <- 'HMDB0000906'
mtb.map[mtb.map$Compound == "Methylamine, N,N-dimethyl-","KEGG"] <- 'C00565'
mtb.map[mtb.map$Compound == "5,9-Undecadien-2-one, 6,10-dimethyl-","HMDB"] <- 'HMDB0031846'
mtb.map[mtb.map$Compound == "5,9-Undecadien-2-one, 6,10-dimethyl-","KEGG"] <- 'C13297'
mtb.map[mtb.map$Compound == " Propyl 2-methylpropanoate","HMDB"] <- 'HMDB0036209'
mtb.map[mtb.map$Compound == "1-Propanol, 2-methyl-","HMDB"] <- 'HMDB0006006'
mtb.map[mtb.map$Compound == "1-Propanol, 2-methyl-","KEGG"] <- 'C14710'
mtb.map[mtb.map$Compound == "2,4-Imidazolidinedione, 1-methyl-","HMDB"] <- 'HMDB0003646'
mtb.map[mtb.map$Compound == "2,4-Imidazolidinedione, 1-methyl-","KEGG"] <- 'C02565'
mtb.map[mtb.map$Compound == "2-Octenal, (E)-","HMDB"] <- 'HMDB0013809'
mtb.map[mtb.map$Compound == "2-Octenal, (E)-","KEGG"] <- 'C21138'
mtb.map[mtb.map$Compound == "Thiophene, 2-pentyl-","HMDB"] <- 'HMDB0040240'
mtb.map[mtb.map$Compound == "Butanoic acid, 3-methylbutyl ester","HMDB"] <- 'HMDB0040221'
mtb.map[mtb.map$Compound == "Butanoic acid, butyl ester","HMDB"] <- 'HMDB0039620'
mtb.map[mtb.map$Compound == "Cyclohexene, 1-methyl-4-(1-methylethylidene)-","KEGG"] <- 'C06075'
mtb.map[mtb.map$Compound == "Cyclohexene, 1-methyl-4-(1-methylethylidene)-","HMDB"] <- 'HMDB0036994'
mtb.map[mtb.map$Compound == "(+)-2-Bornanone","HMDB"] <- 'HMDB0059838'
mtb.map[mtb.map$Compound == "Propanoic acid, 2-methyl-, anhydride","HMDB"] <- 'HMDB0062788'
mtb.map[mtb.map$Compound == "Ethanol, 2-ethoxy-","HMDB"] <- 'HMDB0031213'
mtb.map[mtb.map$Compound == "Ethanol, 2-ethoxy-","KEGG"] <- 'C14687'
mtb.map[mtb.map$Compound == "1-Butene, 4-isothiocyanato-","HMDB"] <- 'HMDB0033867'
mtb.map[mtb.map$Compound == "Benzene, (2-isothiocyanatoethyl)-","HMDB"] <- 'HMDB0038445'
mtb.map[mtb.map$Compound == "Heptacosane","HMDB"] <- 'HMDB0302068'
mtb.map[mtb.map$Compound == "Pentadecane, 2,6,10,14-tetramethyl-","HMDB"] <- 'HMDB0034497'
mtb.map[mtb.map$Compound == "4-Hydroxy-3-hexanone","KEGG"] <- 'C02948'
mtb.map[mtb.map$Compound == "3-Pentanone, 2-methyl-","HMDB"] <- 'HMDB0005846'
mtb.map[mtb.map$Compound == "Benzaldehyde, 4-hydroxy-","KEGG"] <- 'C00633'
mtb.map[mtb.map$Compound == "Benzaldehyde, 4-hydroxy-","HMDB"] <- 'HMDB0011718'
mtb.map[mtb.map$Compound == "Tributyl phosphate","KEGG"] <- 'C14439'
mtb.map[mtb.map$Compound == "Tributyl phosphate","HMDB"] <- 'HMDB0259164'
mtb.map[mtb.map$Compound == "cis-9-Hexadecenal","HMDB"] <- 'HMDB0253124'
mtb.map[mtb.map$Compound == "Hexanoic acid, butyl ester","HMDB"] <- 'HMDB0040211'
mtb.map[mtb.map$Compound == "cis-11-Hexadecenal","HMDB"] <- 'HMDB0243751'
mtb.map[mtb.map$Compound == "Cyclohexanone, 2,2,6-trimethyl-","HMDB"] <- 'HMDB0033794'

# Clean mapping table
mtb.map$MA.Name.Match <- NULL

# Mark duplicated kegg/hmdb identifiers
#  table(mtb.map$KEGG)[table(mtb.map$KEGG) > 1]
#  table(mtb.map$HMDB)[table(mtb.map$HMDB) > 1]
  
# --------------------------------
# Keep only samples with all data
# --------------------------------

sample.intersect <- Reduce(intersect, list(names(mtb)[-1], names(genera)[-1], names(species)[-1], metadata$Sample))
message(paste(length(sample.intersect),"samples have all data types"))

mtb <- mtb[,c("Compound",sample.intersect)]
genera <- genera[,c("Genus",sample.intersect)]
species <- species[,c("Species",sample.intersect)]
metadata <- metadata[metadata$Sample %in% sample.intersect,]

# --------------------------------
# Save to files + R objects
# --------------------------------

save.to.files(DATASET_NAME, "prelim_data", metadata, mtb, mtb.map, genera, species)
save.to.rdata(DATASET_NAME, "prelim_data", metadata, mtb, mtb.map, genera, species)
rm(list = ls())
