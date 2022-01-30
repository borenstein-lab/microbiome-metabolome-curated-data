# ----------------------------------------------------------------
# The following script loads and processes data from:
# 
# Sinha, Rashmi, et al. "Fecal microbiota, fecal metabolome, and 
#  colorectal cancer interrelations." 
#  PloS one 11.3 (2016): e0152126.
# 
# The output tables are:
# 1) metadata - subject metadata
# 2) genera - taxonomic profiles per subject (genus-level)
# 3) mtb - metabolomic profiles per subject       
# 4) mtb.map - additional details per metabolite 
# ----------------------------------------------------------------

require(MetaboAnalystR)
require(readr)
require(dplyr)
source("load_original_data/utils.R")

# --------------------------------
# Required files & info
# --------------------------------

# For details about the source of each file below see: <COMPLETE>
METADATA_FILE <- "../data/original_data/SINHA_CRC_2016/journal.pone.0152126.s004.CSV"
TAXONOMY_FILE <- "../data/original_data/SINHA_CRC_2016/journal.pone.0152126.s003.CSV"
METABOLOMICS_FILE <- "../data/original_data/SINHA_CRC_2016/journal.pone.0152126.s002.CSV"

PUBLICATION_NAME <- 'Fecal Microbiota, Fecal Metabolome, and Colorectal Cancer Interrelations'
DOI <- '10.1371/journal.pone.0152126'
DATASET_NAME <- 'SINHA_CRC_2016'

# --------------------------------
# Load metadata
# --------------------------------

metadata <- read_csv(METADATA_FILE, 
                     col_types = 
                       cols_only(ID = col_character(), 
                                 age = col_integer(), 
                                 bmi = col_double(), 
                                 case = col_integer(), 
                                 hosp1 = col_integer(), 
                                 hosp2 = col_integer(), 
                                 race = col_integer(), 
                                 sex = col_integer()))
# Organize column names & order
metadata <- metadata %>%
  rename(Sample = 1) %>%
  mutate(Subject = Sample) %>%
  rename(Study.Group = case) %>%
  rename(Age = age) %>%
  mutate(Age.Units = "Years") %>%
  rename(Gender = sex) %>%
  mutate(Gender = ifelse(Gender == 0, "Male", "Female")) %>%
  rename(BMI = bmi) %>%
  # Add publication info
  mutate(Dataset = DATASET_NAME) %>%
  mutate(DOI = DOI) %>%
  mutate(Publication.Name = PUBLICATION_NAME) %>%
  # Reorder
  relocate(c(Dataset, Sample, Subject, Study.Group, Age, Age.Units, Gender, BMI, DOI, Publication.Name)) 

# --------------------------------
# Load taxonomic profiles (genera)
# --------------------------------

# This table provided in the paper supp' includes 
#  taxonomic relative abundance of all taxonomic levels.
# From paper: "The current analysis was restricted to the 
#  220 microbes (across taxonomic levels, including 
#  91 Firmicutes, 33 Bacteroidetes, 45 Proteobacteria, 
#  11 Actinobacteria, 5 Fusobacteria, and 35 in other phyla) 
#  that were detected in at least 13 (10%) of the subjects."
# We will first reformat the table and then extract only 
#  genus-level abundances.
genera <- read_csv(TAXONOMY_FILE)

# Transpose
genera <- genera %>% tibble::column_to_rownames("ID")
genera <- data.frame(t(genera), check.names = FALSE)
genera <- genera %>% tibble::rownames_to_column(var = "Genus")

# Get genus-level only
genera <- genera[grepl(";g__", genera$Genus),]
genera$Genus <- gsub("Root;","",genera$Genus)
# Verify (most samples should sum to 1): hist(apply(genera[,2:ncol(genera)],2,sum))

# Fix typos
genera$Genus <- gsub("c__Erysipelotrichi;","c__Erysipelotrichia;",genera$Genus)

# --------------------------------
# Load metabolomic profiles
# --------------------------------

# HPLC-GC/MS-MS, Metabolon
# Note from paper: "The current analysis was restricted 
#  to the 530 metabolites that were detected in at least 
#  118 (90%) of the subjects."
mtb <- read_csv(METABOLOMICS_FILE)
mtb.map <- data.frame(Compound = names(mtb)[-1], 
                      stringsAsFactors = FALSE)

# Transpose original mtb file
mtb <- mtb %>% tibble::column_to_rownames("ID")
mtb <- data.frame(t(mtb), check.names = FALSE)
mtb <- mtb %>% tibble::rownames_to_column(var = "Compound")

# Slightly reformat metabolite names to get rid of characters that fail MetaboAnalyst
mtb.map$Compound.Name <- gsub("^_","",mtb.map$Compound)
mtb.map$Compound.Name <- gsub('([0-9])(_)([0-9])', '\\1,\\3', mtb.map$Compound.Name)
mtb.map$Compound.Name <- gsub('_', ' ', mtb.map$Compound.Name)
mtb.map$High.Confidence.Annotation <- TRUE

# Get mapping to KEGG/other compound ID's, using MetaboAnalyst 
cmpds.to.search <- mtb.map$Compound.Name
MA.matches <- map.compound.names.MetaboAnalyst(cmpds.to.search)

# Merge with the names in our data
mtb.map <- merge(mtb.map, MA.matches, 
                 by.x = "Compound.Name", 
                 by.y = "Query", all = T)
mtb.map$MA.Name.Match <- NULL

# Manual mappings
mtb.map[mtb.map$Compound == "_13_METHYLMYRISTIC_ACID","HMDB"] <- 'HMDB0061707'
mtb.map[mtb.map$Compound == "_15_METHYLPALMITATE","HMDB"] <- 'HMDB0061709'
mtb.map[mtb.map$Compound == "_2_HYDROXYGLUTARATE","KEGG"] <- 'C02630'
mtb.map[mtb.map$Compound == "_2_HYDROXYGLUTARATE","HMDB"] <- 'HMDB0059655'
mtb.map[mtb.map$Compound == "_2_OXINDOLE_3_ACETATE","HMDB"] <- 'HMDB0094657'
mtb.map[mtb.map$Compound == "_3_4_DIHYDROXYBENZOATE","KEGG"] <- 'C00230'
mtb.map[mtb.map$Compound == "_3_4_DIHYDROXYBENZOATE","HMDB"] <- 'HMDB0001856'
mtb.map[mtb.map$Compound == "_3_7_DIMETHYLURATE","KEGG"] <- 'C16360'
mtb.map[mtb.map$Compound == "_3_7_DIMETHYLURATE","HMDB"] <- 'HMDB0001982'
mtb.map[mtb.map$Compound == "_6_OXOPIPERIDINE_2_CARBOXYLIC","HMDB"] <- 'HMDB0061705'
mtb.map[mtb.map$Compound == "_8_HYDROXYOCTANOATE","HMDB"] <- 'HMDB0061914'
mtb.map[mtb.map$Compound == "ALPHA_AMYRIN","KEGG"] <- 'C08615'
mtb.map[mtb.map$Compound == "ASPARAGYLLEUCINE","HMDB"] <- 'HMDB0028735'
mtb.map[mtb.map$Compound == "PALMITOYL_SPHINGOMYELIN","HMDB"] <- 'HMDB0061712'
mtb.map[mtb.map$Compound == "PYROGLUTAMYLVALINE","HMDB"] <- 'HMDB0094651'
mtb.map[mtb.map$Compound == "VAL_VAL_VAL","HMDB"] <- 'HMDB0094676'
mtb.map[mtb.map$Compound == "N_ACETYLPROLINE","HMDB"] <- 'HMDB0094701'
mtb.map[mtb.map$Compound == "N_ACETYLMURAMATE","KEGG"] <- 'C02713'
mtb.map[mtb.map$Compound == "N_ACETYLMURAMATE","HMDB"] <- 'HMDB0060493'
mtb.map[mtb.map$Compound == "METHYLPHOSPHATE","HMDB"] <- 'HMDB0061711'
mtb.map[mtb.map$Compound == "LINOLEAMIDE__18_2N6","HMDB"] <- 'HMDB0062656'
mtb.map[mtb.map$Compound == "DIHYDROFERULIC_ACID","HMDB"] <- 'HMDB0062121'
mtb.map[mtb.map$Compound == "DIHOMO_LINOLEATE__20_2N6","HMDB"] <- 'HMDB0061864'

mtb.map[mtb.map$Compound == '_3_HYDROXYISOBUTYRATE','HMDB'] <- 'HMDB0000023' 
mtb.map[mtb.map$Compound == '_3_HYDROXYISOBUTYRATE','KEGG'] <- 'C06001'
mtb.map[mtb.map$Compound == '_3_HYDROXYISOBUTYRATE','High.Confidence.Annotation'] <- FALSE
mtb.map[mtb.map$Compound == '_1_3_7_TRIMETHYLURATE','HMDB'] <- 'HMDB0002123'	
mtb.map[mtb.map$Compound == '_1_HEXADECANOL','HMDB'] <- 'HMDB0003424' 
mtb.map[mtb.map$Compound == '_1_HEXADECANOL','KEGG'] <- 'C00823'
mtb.map[mtb.map$Compound == 'N1_METHYLADENOSINE','HMDB'] <- 'HMDB0003331' 
mtb.map[mtb.map$Compound == 'N1_METHYLADENOSINE','KEGG'] <- 'C02494'
mtb.map[mtb.map$Compound == '_1_METHYLXANTHINE','HMDB'] <- 'HMDB0010738'	
mtb.map[mtb.map$Compound == '_2_ETHYLHEXANOATE','HMDB'] <- 'HMDB0031230'	
mtb.map[mtb.map$Compound == '_2_HYDROXYMYRISTATE','HMDB'] <- 'HMDB0002261'	
mtb.map[mtb.map$Compound == '_3_DEHYDROCARNITINE','HMDB'] <- 'HMDB0012154'	
mtb.map[mtb.map$Compound == '_3_DEHYDROCARNITINE','KEGG'] <- 'C02636'
mtb.map[mtb.map$Compound == '_3_HYDROXYBENZOATE','HMDB'] <- 'HMDB0002466'	
mtb.map[mtb.map$Compound == '_3_HYDROXYBENZOATE','KEGG'] <- 'C00587'
mtb.map[mtb.map$Compound == '_3_HYDROXYDECANOATE','HMDB'] <- 'HMDB0002203'	
mtb.map[mtb.map$Compound == 'BETA_HYDROXYISOVALERATE','HMDB'] <- 'HMDB0000754'	
mtb.map[mtb.map$Compound == '_3_HYDROXYOCTANOATE','HMDB'] <- 'HMDB0001954'	
mtb.map[mtb.map$Compound == '_3_METHYL_2_OXOVALERATE','HMDB'] <- 'HMDB0000491'	
mtb.map[mtb.map$Compound == '_3_METHYL_2_OXOVALERATE','KEGG'] <- 'C03465'
mtb.map[mtb.map$Compound == '_4_ACETAMIDOBUTANOATE','HMDB'] <- 'HMDB0003681'	
mtb.map[mtb.map$Compound == '_4_ACETAMIDOBUTANOATE','KEGG'] <- 'C02946'
mtb.map[mtb.map$Compound == '_4_GUANIDINOBUTANOATE','HMDB'] <- 'HMDB0003464'	
mtb.map[mtb.map$Compound == '_4_GUANIDINOBUTANOATE','KEGG'] <- 'C01035'
mtb.map[mtb.map$Compound == 'P_HYDROXYBENZALDEHYDE','HMDB'] <- 'HMDB0011718'	
mtb.map[mtb.map$Compound == 'P_HYDROXYBENZALDEHYDE','KEGG'] <- 'C00633'
mtb.map[mtb.map$Compound == '_4_HYDROXYBENZOATE','HMDB'] <- 'HMDB0000500'	
mtb.map[mtb.map$Compound == '_4_HYDROXYBENZOATE','KEGG'] <- 'C00156'
mtb.map[mtb.map$Compound == '_4_HYDROXYCINNAMATE','HMDB'] <- 'HMDB0002035'	
mtb.map[mtb.map$Compound == '_4_HYDROXYCINNAMATE','KEGG'] <- 'C00811'
mtb.map[mtb.map$Compound == '_5_AMINOVALERATE','HMDB'] <- 'HMDB0003355'	
mtb.map[mtb.map$Compound == '_5_AMINOVALERATE','KEGG'] <- 'C00431'
mtb.map[mtb.map$Compound == '_6_HYDROXYNICOTINATE','HMDB'] <- 'HMDB0002658'	
mtb.map[mtb.map$Compound == '_6_HYDROXYNICOTINATE','KEGG'] <- 'C01020'
mtb.map[mtb.map$Compound == '_7_METHYLGUANINE','HMDB'] <- 'HMDB0000897'	
mtb.map[mtb.map$Compound == '_7_METHYLGUANINE','KEGG'] <- 'C02242'
mtb.map[mtb.map$Compound == '_4_ACETAMIDOPHENOL','HMDB'] <- 'HMDB0001859'	
mtb.map[mtb.map$Compound == 'N_ACETYLGLYCINE','HMDB'] <- 'HMDB0000532'	
mtb.map[mtb.map$Compound == '_2_HYDROXYISOBUTYRATE','HMDB'] <- 'HMDB0000729'	
mtb.map[mtb.map$Compound == 'ALPHA_TOCOPHEROL','HMDB'] <- 'HMDB0001893'	
mtb.map[mtb.map$Compound == 'ALPHA_TOCOPHEROL','KEGG'] <- 'C02477'
mtb.map[mtb.map$Compound == '_2_AMINOADIPATE','HMDB'] <- 'HMDB0000510'	
mtb.map[mtb.map$Compound == '_2_AMINOADIPATE','KEGG'] <- 'C00956'
mtb.map[mtb.map$Compound == 'BETA_SITOSTEROL','HMDB'] <- 'HMDB0000852'	
mtb.map[mtb.map$Compound == 'BETA_SITOSTEROL','KEGG'] <- 'C01753'
mtb.map[mtb.map$Compound == 'DHEA_S','HMDB'] <- 'HMDB0001032'	
mtb.map[mtb.map$Compound == 'DHEA_S','KEGG'] <- 'C04555'
mtb.map[mtb.map$Compound == '_5_HYDROXYHEXANOATE','HMDB'] <- 'HMDB0000453'	
mtb.map[mtb.map$Compound == 'DELTA_TOCOPHEROL','HMDB'] <- 'HMDB0002902'	
mtb.map[mtb.map$Compound == 'DELTA_TOCOPHEROL','KEGG'] <- 'C14151'
mtb.map[mtb.map$Compound == 'ALPHA_HYDROXYISOCAPROATE','HMDB'] <- 'HMDB0000624'	
mtb.map[mtb.map$Compound == 'ALPHA_HYDROXYISOCAPROATE','KEGG'] <- 'C03264'
mtb.map[mtb.map$Compound == 'GAMMA_GLUTAMYLISOLEUCINE','HMDB'] <- 'HMDB0011170'	
mtb.map[mtb.map$Compound == 'GAMMA_GLUTAMYLLEUCINE','HMDB'] <- 'HMDB0011171'	
mtb.map[mtb.map$Compound == 'GAMMA_GLUTAMYLPHENYLALANINE','HMDB'] <- 'HMDB0000594'	
mtb.map[mtb.map$Compound == 'GAMMA_GLUTAMYLVALINE','HMDB'] <- 'HMDB0011172'	
mtb.map[mtb.map$Compound == 'GAMMA_TOCOPHEROL','HMDB'] <- 'HMDB0001492'	
mtb.map[mtb.map$Compound == 'ALPHA_GLUTAMYLGLUTAMATE','HMDB'] <- 'HMDB0028818'	
mtb.map[mtb.map$Compound == 'ALPHA_GLUTAMYLTRYPTOPHAN','HMDB'] <- 'HMDB0028830'	
mtb.map[mtb.map$Compound == 'ALPHA_GLUTAMYLTYROSINE','HMDB'] <- 'HMDB0028831'	
mtb.map[mtb.map$Compound == 'ALPHA_GLUTAMYLVALINE','HMDB'] <- 'HMDB0028832'	
mtb.map[mtb.map$Compound == '_1_PALMITOYLGLYCEROL','HMDB'] <- 'HMDB0031074'	
mtb.map[mtb.map$Compound == '_3_PHENYLPROPIONATE','HMDB'] <- 'HMDB0000764'	
mtb.map[mtb.map$Compound == '_3_PHENYLPROPIONATE','KEGG'] <- 'C05629'
mtb.map[mtb.map$Compound == '_3_HYDROXYPROPANOATE','HMDB'] <- 'HMDB0000700'	
mtb.map[mtb.map$Compound == '_3_HYDROXYPROPANOATE','KEGG'] <- 'C01013'
mtb.map[mtb.map$Compound == '_4_METHYL_2_OXOPENTANOATE','HMDB'] <- 'HMDB0000695'	
mtb.map[mtb.map$Compound == '_4_METHYL_2_OXOPENTANOATE','KEGG'] <- 'C00233'
mtb.map[mtb.map$Compound == 'ALLO_THREONINE','HMDB'] <- 'HMDB0004041'	
mtb.map[mtb.map$Compound == 'ALLO_THREONINE','KEGG'] <- 'C05519'
mtb.map[mtb.map$Compound == '_2_AMINOBUTYRATE','HMDB'] <- 'HMDB0000452'	
mtb.map[mtb.map$Compound == '_2_AMINOBUTYRATE','KEGG'] <- 'C02356'
mtb.map[mtb.map$Compound == 'N6_ACETYLLYSINE','HMDB'] <- 'HMDB0000206'	
mtb.map[mtb.map$Compound == 'N6_ACETYLLYSINE','KEGG'] <- 'C02727'
mtb.map[mtb.map$Compound == 'N_ACETYLGLUCOSAMINE','HMDB'] <- 'HMDB0000215'	
mtb.map[mtb.map$Compound == 'N_ACETYLGLUCOSAMINE','KEGG'] <- 'C00140'
mtb.map[mtb.map$Compound == 'N_ACETYLGLUTAMATE','HMDB'] <- 'HMDB0001138'	
mtb.map[mtb.map$Compound == 'N_ACETYLGLUTAMATE','KEGG'] <- 'C00624'
mtb.map[mtb.map$Compound == 'N_ACETYLGLUTAMINE','HMDB'] <- 'HMDB0006029'	
mtb.map[mtb.map$Compound == 'N_ACETYLHISTIDINE','HMDB'] <- 'HMDB0032055'	
mtb.map[mtb.map$Compound == 'N_ACETYLHISTIDINE','KEGG'] <- 'C02997'
mtb.map[mtb.map$Compound == 'N_ACETYLALANINE','HMDB'] <- 'HMDB0000766'	
mtb.map[mtb.map$Compound == 'N_ACETYLLEUCINE','HMDB'] <- 'HMDB0011756'	
mtb.map[mtb.map$Compound == 'N_ACETYLLEUCINE','KEGG'] <- 'C02710'
mtb.map[mtb.map$Compound == 'N_ACETYLMETHIONINE','HMDB'] <- 'HMDB0011745'	
mtb.map[mtb.map$Compound == 'N_ACETYLMETHIONINE','KEGG'] <- 'C02712'
mtb.map[mtb.map$Compound == 'N_ACETYLPHENYLALANINE','HMDB'] <- 'HMDB0000512'	
mtb.map[mtb.map$Compound == 'N_ACETYLPHENYLALANINE','KEGG'] <- 'C03519'
mtb.map[mtb.map$Compound == 'N_ACETYLTYROSINE','HMDB'] <- 'HMDB0000866'	
mtb.map[mtb.map$Compound == 'N_ACETYLTYROSINE','KEGG'] <- 'C01657'
mtb.map[mtb.map$Compound == 'N_ACETYLNEURAMINATE','HMDB'] <- 'HMDB0000230'	
mtb.map[mtb.map$Compound == 'N_ACETYLNEURAMINATE','KEGG'] <- 'C19910'
mtb.map[mtb.map$Compound == 'N_ACETYLPUTRESCINE','HMDB'] <- 'HMDB0002064'	
mtb.map[mtb.map$Compound == 'N_ACETYLPUTRESCINE','KEGG'] <- 'C02714'
mtb.map[mtb.map$Compound == 'N_ACETYLTRYPTOPHAN','HMDB'] <- 'HMDB0013713'	
mtb.map[mtb.map$Compound == 'N_ACETYLVALINE','HMDB'] <- 'HMDB0011757'	
mtb.map[mtb.map$Compound == 'N2_ACETYLLYSINE','HMDB'] <- 'HMDB0000446'	
mtb.map[mtb.map$Compound == 'N2_ACETYLLYSINE','KEGG'] <- 'C12989'
mtb.map[mtb.map$Compound == 'N_METHYLPHENYLALANINE','HMDB'] <- 'HMDB0029224'	
mtb.map[mtb.map$Compound == '_1_OCTADECANOL','HMDB'] <- 'HMDB0002350'	
mtb.map[mtb.map$Compound == '_4_HYDROXYPHENYLACETATE','HMDB'] <- 'HMDB0000020'	
mtb.map[mtb.map$Compound == '_4_HYDROXYPHENYLACETATE','KEGG'] <- 'C00642'
mtb.map[mtb.map$Compound == '_5_OXOPROLINE','HMDB'] <- 'HMDB0000267'	
mtb.map[mtb.map$Compound == '_5_OXOPROLINE','KEGG'] <- 'C01879'
mtb.map[mtb.map$Compound == 'SERYLPHENYALANINE','HMDB'] <- 'HMDB0029046'	
mtb.map[mtb.map$Compound == '_4_HYDROXYBUTYRATE__GHB','HMDB'] <- 'HMDB0000710'	
mtb.map[mtb.map$Compound == '_4_HYDROXYBUTYRATE__GHB','KEGG'] <- 'C00989'
mtb.map[mtb.map$Compound == '_7_BETA_HYDROXYCHOLESTEROL','HMDB'] <- 'HMDB0006119'	
mtb.map[mtb.map$Compound == 'ARACHIDATE__20_0','HMDB'] <- 'HMDB0002212'	
mtb.map[mtb.map$Compound == 'ARACHIDATE__20_0','KEGG'] <- 'C06425'
mtb.map[mtb.map$Compound == 'ARACHIDONATE__20_4N6','HMDB'] <- 'HMDB0001043'	
mtb.map[mtb.map$Compound == 'ARACHIDONATE__20_4N6','KEGG'] <- 'C00219'
mtb.map[mtb.map$Compound == 'DIMETHYLARGININE__SDMA___ADMA','HMDB'] <- 'HMDB0001539'	
mtb.map[mtb.map$Compound == 'DIMETHYLARGININE__SDMA___ADMA','KEGG'] <- 'C03626'
mtb.map[mtb.map$Compound == 'AZELATE__NONANEDIOATE','HMDB'] <- 'HMDB0000784'	
mtb.map[mtb.map$Compound == 'AZELATE__NONANEDIOATE','KEGG'] <- 'C08261'
mtb.map[mtb.map$Compound == 'BEHENATE__22_0','HMDB'] <- 'HMDB0000944'	
mtb.map[mtb.map$Compound == 'BEHENATE__22_0','KEGG'] <- 'C08281'
mtb.map[mtb.map$Compound == 'CAPRATE__10_0','HMDB'] <- 'HMDB0000511'	
mtb.map[mtb.map$Compound == 'CAPRATE__10_0','KEGG'] <- 'C01571'
mtb.map[mtb.map$Compound == 'CAPROATE__6_0','HMDB'] <- 'HMDB0000535'	
mtb.map[mtb.map$Compound == 'CAPROATE__6_0','KEGG'] <- 'C01585'
mtb.map[mtb.map$Compound == 'CAPRYLATE__8_0','HMDB'] <- 'HMDB0000482'	
mtb.map[mtb.map$Compound == 'CAPRYLATE__8_0','KEGG'] <- 'C06423'
mtb.map[mtb.map$Compound == '_2__DEOXYGUANOSINE','HMDB'] <- 'HMDB0000085'	
mtb.map[mtb.map$Compound == '_2__DEOXYGUANOSINE','KEGG'] <- 'C00330'
mtb.map[mtb.map$Compound == '_2__DEOXYINOSINE','HMDB'] <- 'HMDB0000071'	
mtb.map[mtb.map$Compound == '_2__DEOXYINOSINE','KEGG'] <- 'C05512'
mtb.map[mtb.map$Compound == '_2__DEOXYURIDINE','HMDB'] <- 'HMDB0000012'	
mtb.map[mtb.map$Compound == '_2__DEOXYURIDINE','KEGG'] <- 'C00526'
mtb.map[mtb.map$Compound == '_1_STEAROYLGLYCERO','HMDB'] <- 'HMDB0011133'	
mtb.map[mtb.map$Compound == '_1_STEAROYLGLYCERO','KEGG'] <- 'C03805'
mtb.map[mtb.map$Compound == 'GAMMA_AMINOBUTYRATE__GABA','HMDB'] <- 'HMDB0000112'	
mtb.map[mtb.map$Compound == 'GAMMA_AMINOBUTYRATE__GABA','KEGG'] <- 'C00334'
mtb.map[mtb.map$Compound == 'GLUTAMINE_GLUTAMINE','HMDB'] <- 'HMDB0028795'	
mtb.map[mtb.map$Compound == 'GLUTAMINE_ISOLEUCINE','HMDB'] <- 'HMDB0028800'	
mtb.map[mtb.map$Compound == 'GLYCEROL_3_PHOSPHATE__G3P','HMDB'] <- 'HMDB0000126'	
mtb.map[mtb.map$Compound == 'GLYCEROL_3_PHOSPHATE__G3P','KEGG'] <- 'C00093'
mtb.map[mtb.map$Compound == 'LEVULINATE__4_OXOVALERATE','HMDB'] <- 'HMDB0000720'	
mtb.map[mtb.map$Compound == 'METHYL_4_HYDROXYBENZOATE','HMDB'] <- 'HMDB0032572'	
mtb.map[mtb.map$Compound == 'N_ACETYLARGININE','HMDB'] <- 'HMDB0004620'	
mtb.map[mtb.map$Compound == 'N_ACETYLASPARTATE__NAA','HMDB'] <- 'HMDB0000812'	
mtb.map[mtb.map$Compound == 'N_ACETYLASPARTATE__NAA','KEGG'] <- 'C01042'
mtb.map[mtb.map$Compound == 'P_CRESOL_SULFATE','HMDB'] <- 'HMDB0011635'	
mtb.map[mtb.map$Compound == 'PHENYLLACTATE__PLA','HMDB'] <- 'HMDB0000779'	
mtb.map[mtb.map$Compound == 'PHENYLLACTATE__PLA','KEGG'] <- 'C01479'
mtb.map[mtb.map$Compound == 'RIBOFLAVIN__VITAMIN_B2','HMDB'] <- 'HMDB0000244'	
mtb.map[mtb.map$Compound == 'RIBOFLAVIN__VITAMIN_B2','KEGG'] <- 'C00255'
mtb.map[mtb.map$Compound == 'STEARATE__18_0','HMDB'] <- 'HMDB0000827'	
mtb.map[mtb.map$Compound == 'STEARATE__18_0','KEGG'] <- 'C01530'
mtb.map[mtb.map$Compound == 'LIGNOCERATE__24_0','HMDB'] <- 'HMDB0002003'	
mtb.map[mtb.map$Compound == 'LIGNOCERATE__24_0','KEGG'] <- 'C08320'
mtb.map[mtb.map$Compound == 'THIAMIN__VITAMIN_B1','HMDB'] <- 'HMDB0000235'	
mtb.map[mtb.map$Compound == 'THIAMIN__VITAMIN_B1','KEGG'] <- 'C00378'
mtb.map[mtb.map$Compound == 'THREONYLALANINE','HMDB'] <- 'HMDB0029054'	
mtb.map[mtb.map$Compound == 'THREONYLGLUTAMATE','HMDB'] <- 'HMDB0029060'	
mtb.map[mtb.map$Compound == 'THREONYLISOLEUCINE','HMDB'] <- 'HMDB0029064'	
mtb.map[mtb.map$Compound == 'THREONYLLEUCINE','HMDB'] <- 'HMDB0029065'	
mtb.map[mtb.map$Compound == 'THREONYLPHENYLALANINE','HMDB'] <- 'HMDB0029068'	
mtb.map[mtb.map$Compound == 'THREONYLSERINE','HMDB'] <- 'HMDB0029070'	
mtb.map[mtb.map$Compound == 'THREONYLVALINE','HMDB'] <- 'HMDB0029074'	
mtb.map[mtb.map$Compound == 'MARGARATE__17_0','HMDB'] <- 'HMDB0002259'	
mtb.map[mtb.map$Compound == 'MYRISTATE__14_0','HMDB'] <- 'HMDB0000806'	
mtb.map[mtb.map$Compound == 'MYRISTATE__14_0','KEGG'] <- 'C06424'
mtb.map[mtb.map$Compound == 'MYRISTOLEATE__14_1N5','HMDB'] <- 'HMDB0002000'	
mtb.map[mtb.map$Compound == 'MYRISTOLEATE__14_1N5','KEGG'] <- 'C08322'
mtb.map[mtb.map$Compound == 'N_6_TRIMETHYLLYSINE','HMDB'] <- 'HMDB0001325'	
mtb.map[mtb.map$Compound == 'N_6_TRIMETHYLLYSINE','KEGG'] <- 'C03793'
mtb.map[mtb.map$Compound == 'PALMITATE__16_0','HMDB'] <- 'HMDB0000220'	
mtb.map[mtb.map$Compound == 'PALMITATE__16_0','KEGG'] <- 'C00249'
mtb.map[mtb.map$Compound == 'PALMITOLEATE__16_1N7','HMDB'] <- 'HMDB0003229'	
mtb.map[mtb.map$Compound == 'PALMITOLEATE__16_1N7','KEGG'] <- 'C08362'
mtb.map[mtb.map$Compound == 'PELARGONATE__9_0','HMDB'] <- 'HMDB0000847'	
mtb.map[mtb.map$Compound == 'PELARGONATE__9_0','KEGG'] <- 'C01601'
mtb.map[mtb.map$Compound == 'PENTADECANOATE__15_0','HMDB'] <- 'HMDB0000826'	
mtb.map[mtb.map$Compound == 'PENTADECANOATE__15_0','KEGG'] <- 'C16537'
mtb.map[mtb.map$Compound == 'STEARIDONATE__18_4N3','HMDB'] <- 'HMDB0006547'	
mtb.map[mtb.map$Compound == 'STEARIDONATE__18_4N3','KEGG'] <- 'C16300'
mtb.map[mtb.map$Compound == 'UNDECANOATE__11_0','HMDB'] <- 'HMDB0000947'	
mtb.map[mtb.map$Compound == 'UNDECANOATE__11_0','KEGG'] <- 'C17715'
mtb.map[mtb.map$Compound == '_4_HYDROXYBENZYL_ALCOHOL','HMDB'] <- 'HMDB0011724'	
mtb.map[mtb.map$Compound == '_4_HYDROXYBENZYL_ALCOHOL','KEGG'] <- 'C17467'
mtb.map[mtb.map$Compound == '_5_DODECENOATE__12_1N7','HMDB'] <- 'HMDB0000529'	
mtb.map[mtb.map$Compound == 'LAURATE__12_0','HMDB'] <- 'HMDB0000638'	
mtb.map[mtb.map$Compound == 'LAURATE__12_0','KEGG'] <- 'C02679'
mtb.map[mtb.map$Compound == 'EICOSENOATE__20_1N9_OR_11','HMDB'] <- 'HMDB0002231'	
mtb.map[mtb.map$Compound == 'EICOSENOATE__20_1N9_OR_11','KEGG'] <- 'C16526'
mtb.map[mtb.map$Compound == 'HEPTANOATE__7_0','HMDB'] <- 'HMDB0000666'	
mtb.map[mtb.map$Compound == 'HEPTANOATE__7_0','KEGG'] <- 'C17714'
mtb.map[mtb.map$Compound == 'LINOLEATE__18_2N6','HMDB'] <- 'HMDB0000673'	
mtb.map[mtb.map$Compound == 'LINOLEATE__18_2N6','KEGG'] <- 'C01595'
mtb.map[mtb.map$Compound == '_1_11_UNDECANEDICARBOXYLATE','HMDB'] <- 'HMDB0002327'	
mtb.map[mtb.map$Compound == '_12_METHYLTRIDECANOIC_ACID','HMDB'] <- 'HMDB0031072'	
mtb.map[mtb.map$Compound == '_3_HYDROXYBUTYRATE__BHBA','HMDB'] <- 'HMDB0000357'	
mtb.map[mtb.map$Compound == '_3_HYDROXYBUTYRATE__BHBA','KEGG'] <- 'C01089'
mtb.map[mtb.map$Compound == 'GLYCYLGLYCINE','HMDB'] <- 'HMDB0011733'	
mtb.map[mtb.map$Compound == 'GLYCYLGLYCINE','KEGG'] <- 'C02037'
mtb.map[mtb.map$Compound == '_1_OLEOYLGLYCEROPHOSPHOCHOLIN','HMDB'] <- 'HMDB0002815'	
mtb.map[mtb.map$Compound == '_1_OLEOYLGLYCEROPHOSPHOCHOLIN','KEGG'] <- 'C04230'
mtb.map[mtb.map$Compound == 'N__2_FUROYL_GLYCINE','HMDB'] <- 'HMDB0000439'	
mtb.map[mtb.map$Compound == 'ASPARAGYLVALINE','HMDB'] <- 'HMDB0028744'	
mtb.map[mtb.map$Compound == 'HISTIDYLISOLEUCINE','HMDB'] <- 'HMDB0028888'	
mtb.map[mtb.map$Compound == 'HISTIDYLVALINE','HMDB'] <- 'HMDB0028898'	
mtb.map[mtb.map$Compound == 'INOSITOL_1_PHOSPHATE__I1P','HMDB'] <- 'HMDB0000213'	
mtb.map[mtb.map$Compound == 'INOSITOL_1_PHOSPHATE__I1P','KEGG'] <- 'C04006'
mtb.map[mtb.map$Compound == 'NERVONATE__24_1N9','HMDB'] <- 'HMDB0002368'	
mtb.map[mtb.map$Compound == 'NERVONATE__24_1N9','KEGG'] <- 'C08323'
mtb.map[mtb.map$Compound == 'N_FORMYLMETHIONINE','HMDB'] <- 'HMDB0001015'	
mtb.map[mtb.map$Compound == 'N_FORMYLMETHIONINE','KEGG'] <- 'C03145'
mtb.map[mtb.map$Compound == 'NONADECANOATE__19_0','HMDB'] <- 'HMDB0000772'	
mtb.map[mtb.map$Compound == 'NONADECANOATE__19_0','KEGG'] <- 'C16535'
mtb.map[mtb.map$Compound == 'OLEATE__18_1N9','HMDB'] <- 'HMDB0000207'	
mtb.map[mtb.map$Compound == 'OLEATE__18_1N9','KEGG'] <- 'C00712'
mtb.map[mtb.map$Compound == '_13_HOTE__6Z_9Z_11E_OR_9Z_11E','HMDB'] <- 'HMDB0010203'	
mtb.map[mtb.map$Compound == '_2_HYDROXYBUTYRATE__AHB','HMDB'] <- 'HMDB0000008'	
mtb.map[mtb.map$Compound == '_2_HYDROXYBUTYRATE__AHB','KEGG'] <- 'C05984'
mtb.map[mtb.map$Compound == '_5_6_DIHYDROURACIL','HMDB'] <- 'HMDB0000076'	
mtb.map[mtb.map$Compound == '_5_6_DIHYDROURACIL','KEGG'] <- 'C00429'

# Fix bad mappings (discovered due to duplicated mappings)
mtb.map[mtb.map$Compound == 'LEUCYLGLYCINE','HMDB'] <- 'HMDB0028929'
mtb.map[mtb.map$Compound == 'LEUCYLGLYCINE','KEGG'] <- NA
mtb.map[mtb.map$Compound == 'HOMOSTACHYDRINE','HMDB'] <- 'HMDB0033433'
mtb.map[mtb.map$Compound == 'HOMOSTACHYDRINE','KEGG'] <- "C08283"

# Slightly reorganize mapping table
mtb.map <- mtb.map %>% select(-Compound.Name) 
names(mtb.map)[names(mtb.map)=="Match"] <- "Compound.Name"

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

