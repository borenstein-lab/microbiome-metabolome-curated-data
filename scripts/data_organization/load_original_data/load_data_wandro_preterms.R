# ----------------------------------------------------------------
# The following script loads and processes data from:
# 
# Wandro, Stephen, et al. "The microbiome and metabolome of 
#  preterm infant stool are personalized and not driven by health 
#  outcomes, including necrotizing enterocolitis and late-onset 
#  sepsis." Msphere 3.3 (2018): e00104-18.
# 
# The output tables are:
# 1) metadata - subject metadata
# 2) genera - taxonomic profiles per subject (genus-level)
# 3) mtb - metabolomic profiles per subject       
# 4) mtb.map - additional details per metabolite 
# ----------------------------------------------------------------

require(readr)
require(dplyr)
source("data_organization/utils.R")

# --------------------------------
# Required files & info
# --------------------------------

# For details about the source of each file below see: <COMPLETE>
RDATA_FILE <- "../data/original_data/WANDRO_PRETERMS_2018/.RData"
TAXONOMY_FILE <- "../data/original_data/WANDRO_PRETERMS_2018/feature_table_gtdb_207.tsv"

PUBLICATION_NAME <- 'The Microbiome and Metabolome of Preterm Infant Stool Are Personalized and Not Driven by Health Outcomes, Including Necrotizing Enterocolitis and Late-Onset Sepsis'
DOI <- '10.1128/mSphere.00104-18'
DATASET_NAME <- 'WANDRO_PRETERMS_2018'

e.tmp <- new.env() # ls(e.tmp)
load(RDATA_FILE, e.tmp)

# --------------------------------
# Load metadata
# --------------------------------

metadata <- get('infant.file',e.tmp)

# Organize column names & order
metadata <- metadata %>%
  rename(Sample = 1) %>%
  rename(Subject = Patient) %>%
  rename(Study.Group = Group) %>%
  rename(Age = Age_at_sample_collection) %>%
  mutate(Age.Units = "Days") %>%
  # Add publication info
  mutate(Dataset = DATASET_NAME) %>%
  mutate(DOI = DOI) %>%
  mutate(Publication.Name = PUBLICATION_NAME) %>%
  # Reorder
  relocate(c(Dataset, Sample, Subject, Study.Group, Age, Age.Units, DOI, Publication.Name)) 

# --------------------------------
# Load taxonomic profiles (genera)
# --------------------------------

# Read genus-level abundances from qiime analysis
genera <- read_delim(TAXONOMY_FILE, show_col_types = FALSE,
                     "\t", escape_double = FALSE, trim_ws = TRUE)
names(genera)[1] <- 'Genus'
# Sanity: apply(genera[,-1], 2, sum)

# --------------------------------
# Load metabolomic profiles
# --------------------------------

# Load metabolite data
met.file <- get('met.file', e.tmp)

# Transpose
mtb <- setNames(data.frame(t(met.file[,-1])), met.file[,1])
mtb <- cbind(Compound = as.character(rownames(mtb)), mtb, stringsAsFactors = FALSE)
mtb$Compound <- trimws(mtb$Compound)
rownames(mtb) <- NULL

# Prepare table of mappings of metabolites to known IDs
met.info <- get('met.info', e.tmp)
met.info$BinBase.name <- trimws(met.info$BinBase.name)
mtb.map <- merge(data.frame(Compound = mtb$Compound),
                      met.info[,c('BinBase.name','KEGG','PubChem')],
                      by.x='Compound',by.y='BinBase.name',
                      all.x = TRUE,
                      all.y = FALSE)
mtb.map$KEGG[mtb.map$KEGG == ""] <- NA
mtb.map$PubChem[mtb.map$PubChem == ""] <- NA
rm(met.file, met.info)

# Use MetaboAnalyst to retrieve HMDB ID's. 
# We will search by name but validate against  
#  the already available KEGG ID's
cmpds.to.search <- mtb.map$Compound
MA.matches <- map.compound.names.MetaboAnalyst(cmpds.to.search)
MA.matches <- MA.matches %>% rename(KEGG.from.MA = KEGG) 

# Merge with the names in our data
mtb.map <- merge(mtb.map, MA.matches, 
                 by.x = "Compound", 
                 by.y = "Query", all = T)

# Mark "lower confidence annotation" in cases where provided 
#  KEGG ID and MetaboAnalyst-idenfied KEGG ID do not match.
mtb.map <- mtb.map %>%
  mutate(High.Confidence.Annotation = 
           (is.na(KEGG) | is.na(KEGG.from.MA) | KEGG == KEGG.from.MA))

# For a few specific cases out of the above, we do switch the kegg ID 
# (for these cases we assume the original kegg annotation was an error)
# View(mtb.map %>% filter(!(is.na(KEGG) | is.na(KEGG.from.MA) | KEGG == KEGG.from.MA)))
mtb.map$KEGG[mtb.map$Compound == "cellobiose"] <- "C06422"
mtb.map$KEGG[mtb.map$Compound == "cholic acid"] <- "C00695"
mtb.map$KEGG[mtb.map$Compound == "hexitol"] <- "C01697"
mtb.map$KEGG[mtb.map$Compound == "putrescine"] <- "C00134"

mtb.map <- mtb.map %>%
  mutate(KEGG = coalesce(KEGG, KEGG.from.MA)) %>%
  select(-KEGG.from.MA, -MA.Name.Match)

# Additional mappings for unmapped named metabolites 
# View(mtb.map %>% filter((is.na(mtb.map$HMDB) | is.na(KEGG)) & grepl("[a-zA-Z]", mtb.map$Compound)))
mtb.map[mtb.map$Compound == "2-hydroxyglutaric acid","HMDB"] <- 'HMDB0059655' 
mtb.map[mtb.map$Compound == "4-hydroxymandelic acid","HMDB"] <- 'HMDB0000822' 
mtb.map[mtb.map$Compound == "4-hydroxymandelic acid","KEGG"] <- 'C11527' 
mtb.map[mtb.map$Compound == "adenosine-5-monophosphate","HMDB"] <- 'HMDB0000045' 
mtb.map[mtb.map$Compound == "alanine","HMDB"] <- 'HMDB0000161' 
mtb.map[mtb.map$Compound == "cystine","HMDB"] <- 'HMDB0000192' 
mtb.map[mtb.map$Compound == "cystine","KEGG"] <- 'C00491' 
mtb.map[mtb.map$Compound == "indole-3-lactate","HMDB"] <- 'HMDB0000671' 




mtb.map[mtb.map$Compound == "tocopherol beta NIST","HMDB"] <- 'HMDB0006335'; 
mtb.map[mtb.map$Compound == "cerotinic aci","KEGG"] <- 'C21931' 
mtb.map[mtb.map$Compound == "cerotinic aci","HMDB"] <- 'HMDB0002356' 
mtb.map[mtb.map$Compound == "cerotinic aci","High.Confidence.Annotation"] <- FALSE; 
mtb.map[mtb.map$Compound == "5-hydroxymethyl-2-furoic acid NIST","HMDB"] <- 'HMDB0002432'
mtb.map[mtb.map$Compound == "acetophenone NIST","HMDB"] <- 'HMDB0033910'
mtb.map[mtb.map$Compound == "7-methylguanine NIST","KEGG"] <- 'C02242'
mtb.map[mtb.map$Compound == "7-methylguanine NIST","HMDB"] <- 'HMDB0000897'
mtb.map[mtb.map$Compound == "beta-mannosylglycerate","KEGG"] <- 'C16699'
mtb.map[mtb.map$Compound == "beta-mannosylglycerate","HMDB"] <- 'HMDB0012152'
mtb.map[mtb.map$Compound == "n-acetyl-d-hexosamine","KEGG"] <- 'C03136'
mtb.map[mtb.map$Compound == "butyrolactam NIST","KEGG"] <- 'C11118'
mtb.map[mtb.map$Compound == "butyrolactam NIST","HMDB"] <- 'HMDB0002039'
mtb.map[mtb.map$Compound == "parabanic acid NIST","HMDB"] <- 'HMDB0062802'
mtb.map[mtb.map$Compound == "tocopherol alpha-","KEGG"] <- 'C02477'
mtb.map[mtb.map$Compound == "tocopherol alpha-","HMDB"] <- 'HMDB0001893'
mtb.map$High.Confidence.Annotation[mtb.map$Compound == "tocopherol alpha-"] <- FALSE # Contradicts author annotation
mtb.map[mtb.map$Compound == "tocopherol gamma-","HMDB"] <- 'HMDB0001492' 
mtb.map[mtb.map$Compound == "xylonolactone NIST","HMDB"] <- 'HMDB0011676' 
mtb.map[mtb.map$Compound == "xylulose NIST","HMDB"] <- 'HMDB0000751' 
mtb.map[mtb.map$Compound == "xylonic acid","KEGG"] <- 'C05411' 
mtb.map[mtb.map$Compound == "xylonic acid","HMDB"] <- 'HMDB0060256' 
mtb.map$High.Confidence.Annotation[mtb.map$Compound == "xylonic acid"] <- FALSE # Contradicts author annotation
mtb.map[mtb.map$Compound == "delta-tocopherol NIST","KEGG"] <- 'C14151' 
mtb.map[mtb.map$Compound == "delta-tocopherol NIST","HMDB"] <- 'HMDB0002902'
mtb.map[mtb.map$Compound == "glycerol-3-galactoside","HMDB"] <- 'HMDB0006790' 
mtb.map[mtb.map$Compound == "inositol-4-monophosphate","HMDB"] <- 'HMDB0001313' 
mtb.map$High.Confidence.Annotation[mtb.map$Compound == "isothreitol"] <- FALSE
mtb.map[mtb.map$Compound == "lyxitol","HMDB"] <- 'HMDB0001851' 
mtb.map$High.Confidence.Annotation[mtb.map$Compound == "lyxitol"] <- FALSE 
mtb.map[mtb.map$Compound == "conduritol-beta-epoxide","HMDB"] <- 'HMDB0247281' 
mtb.map$High.Confidence.Annotation[mtb.map$Compound == "tagatose"] <- FALSE
mtb.map[(!is.na(mtb.map$HMDB)) & mtb.map$HMDB == 'HMDB0062263','HMDB'] <- 'HMDB0000187'

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

sample.intersect <- Reduce(intersect, list(names(mtb)[-1], names(genera)[-1], metadata$Sample))
message(paste(length(sample.intersect),"samples have all data types"))

mtb <- mtb[,c("Compound",sample.intersect)]
genera <- genera[,c("Genus",sample.intersect)]
metadata <- metadata[metadata$Sample %in% sample.intersect,]

# --------------------------------
# Save to files + R objects
# --------------------------------

save.to.files(DATASET_NAME, "prelim_data", metadata = metadata, mtb = mtb, mtb.map = mtb.map, genera = genera)
save.to.rdata(DATASET_NAME, "prelim_data", metadata = metadata, mtb = mtb, mtb.map = mtb.map, genera = genera, override.all = T)

rm(list = ls())

