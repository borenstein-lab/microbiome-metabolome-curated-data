# ----------------------------------------------------------------
# The following script loads and processes data from:
# 
# Mars, Ruben AT, et al. "Longitudinal multi-omics reveals 
#  subset-specific mechanisms underlying irritable bowel syndrome." 
#  Cell 182.6 (2020): 1460-1473.
# 
# The output tables are:
# 1) metadata - subject metadata
# 2) genera - taxonomic profiles per subject (genus-level)
# 3) species - taxonomic profiles per subject (species-level)
# 4) mtb - metabolomic profiles per subject       
# 5) mtb.map - additional details per metabolite 
# ----------------------------------------------------------------

require(gdata)
require(dplyr)
options("scipen"=100, "digits"=4)

source("data_organization/utils.R")
source("data_organization/gtdb_utils.R")

# --------------------------------
# Required files & info
# --------------------------------

# For details about the source of each file below see: <COMPLETE>
METADATA_FILE <- "../data/original_data/MARS_IBS_2020/1-s2.0-S0092867420309983-mmc1.xlsx"
TAXONOMY_FILE_SP <- '../data/original_data/MARS_IBS_2020/kraken/kraken_species_level_taxonomy.tsv'
TAXONOMY_FILE_GE <- '../data/original_data/MARS_IBS_2020/kraken/kraken_genus_level_taxonomy.tsv'
TAXONOMY_SAMPLE_MAP <- '../data/original_data/MARS_IBS_2020/kraken/PRJEB37924.txt'
METABOLOMICS_FILE_1 <- "../data/original_data/MARS_IBS_2020/Supplementary data II source data NMR metabolomics.xlsx"
METABOLOMICS_FILE_2 <- "../data/original_data/MARS_IBS_2020/Supplementary data IV source data Tryptophan metabolomics.xlsx"
METABOLOMICS_FILE_3 <- "../data/original_data/MARS_IBS_2020/Supplementary data V source data bile acid metabolomics.xlsx"

PUBLICATION_NAME <- 'Longitudinal Multi-omics Reveals Subset-Specific Mechanisms Underlying Irritable Bowel Syndrome'
DOI <- '10.1016/j.cell.2020.08.007'
DATASET_NAME <- 'MARS_IBS_2020'

# --------------------------------
# Load metadata
# --------------------------------

metadata <- read.xls(METADATA_FILE , 
                     sheet = 2, header = TRUE, 
                     na.strings=c("NA","#DIV/0!",""," "),
                     stringsAsFactors = FALSE)
metadata <- metadata %>%
  select(SampleID, study_id, Timepoint,
         ID_on_tube, Cohort, Age, BMI, Gender,
         survey_event, survey_event_date,
         IBS_symptom_severity_1_500, Effect_on_life_in_general_0_100,
         Flare, Flare_timepoint.y,
         Antibiotics, which_abx,
         Antibiotic_in_the_past_month._1_no_2_yes, Antibiotic_used,
         X_probiotic_in_the_past_2_weeks_1_no_2_yes, Probiotic_used,
         unique_id) %>%
  rename(Sample = SampleID) %>%
  rename(Subject = study_id) %>%
  rename(Study.Group = Cohort) %>%
  mutate(Gender = ifelse(Gender == "F", "Female", "Male")) %>%
  mutate(Age.Units = "Years") %>%
  # Add publication info
  mutate(Dataset = DATASET_NAME) %>%
  mutate(DOI = DOI) %>%
  mutate(Publication.Name = PUBLICATION_NAME) %>%
  # Reorder
  relocate(c(Dataset, Sample, Subject, Study.Group, 
             Age, Age.Units, Gender, BMI, DOI, Publication.Name)) 

# --------------------------------
# Load taxonomic profiles 
# --------------------------------

species <- read_delim(TAXONOMY_FILE_SP, "\t", 
                      escape_double = FALSE, 
                      show_col_types = FALSE,
                      trim_ws = TRUE)
names(species)[1] <- 'Species'

genera <- read_delim(TAXONOMY_FILE_GE, "\t", 
                     escape_double = FALSE, 
                     show_col_types = FALSE,
                     trim_ws = TRUE)
names(genera)[1] <- 'Genus'

tax.map <- read_delim(TAXONOMY_SAMPLE_MAP, 
                      delim = "\t", 
                      show_col_types = FALSE,
                      col_select = c("run_accession", "submitted_ftp"), 
                      escape_double = FALSE, 
                      trim_ws = TRUE) %>%
  filter(grepl('Study\\.ID', submitted_ftp)) %>%
  mutate(Sample = gsub(".*Study\\.ID\\.", "", submitted_ftp)) %>%
  mutate(Sample = gsub("\\.S[0-9]*\\.R1\\.001\\.fq\\.gz$","",Sample)) %>%
  mutate(run_accession = paste0(run_accession, '_Illumina_HiSeq_4000_sequencing')) %>%
  filter(run_accession %in% names(genera))

tax.map.vec <- c("Species", "Genus", tax.map$Sample)
names(tax.map.vec) <- c("Species", "Genus", tax.map$run_accession)

# Map file names to sample id's
species <- species %>% select(any_of(names(tax.map.vec)))
genera <- genera %>% select(any_of(names(tax.map.vec)))
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

# 1
mtb1 <- read.xls(METABOLOMICS_FILE_1, 
                 sheet = 1, header = TRUE, 
                 check.names = FALSE, 
                 stringsAsFactors = FALSE)
mtb1 <- mtb1 %>%
  mutate(Sample = paste0(ID_on_tube, ".T.",
                         ifelse(`Time point` == 7, 
                                "Flare", 
                                `Time point`))) %>%
  select(-Matlab_no, -ID_on_tube, -`Time point`,
         -Gender, -`Gender code`, 
         -Diagnosis, -Group) 
  
# Sanity: sum(mtb1$Sample %in% metadata$Sample) == nrow(mtb1)
mtb1_t <- data.frame(t(tibble::column_to_rownames(mtb1, "Sample")), check.names = FALSE) %>%
  tibble::rownames_to_column("Compound") %>%
  mutate(Method = "NMR")

# 2
mtb2 <- read.xls(METABOLOMICS_FILE_2, 
                 sheet = 1, header = TRUE, 
                 check.names = FALSE, 
                 na.strings = c("ND","","NA"), 
                 stringsAsFactors = FALSE)
mtb2 <- mtb2 %>%
  rename(Sample = SampleID) %>%
  select(-BIOME_ID_nr,-Cohort,-SubjectID,-time_point,-alt_name,-ID_on_tube)

# Sanity: sum(mtb2$Sample %in% metadata$Sample) == nrow(mtb2)
mtb2_t <- data.frame(t(tibble::column_to_rownames(mtb2, "Sample")), check.names = FALSE) %>%
  tibble::rownames_to_column("Compound") %>%
  mutate(Method = "Tryptophan_Metabolites")

# 3
mtb3 <- read.xls(METABOLOMICS_FILE_3, 
                 sheet = 1, header = TRUE, 
                 check.names = FALSE, 
                 na.strings = c("","NA"), 
                 stringsAsFactors = FALSE)
mtb3 <- mtb3 %>%
  rename(Sample = SampleID) %>%
  select(-ID_on_tube,-`Time point`,-Flare,-Group,-Age,-BMI,-Gender)

# Sanity: sum(mtb3$Sample %in% metadata$Sample) == nrow(mtb3)
mtb3_t <- data.frame(t(tibble::column_to_rownames(mtb3, "Sample")), check.names = FALSE) %>%
  tibble::rownames_to_column("Compound") %>%
  mutate(Method = "BA_Metabolites")

# Combine all metabolites
mtb <- bind_rows(mtb1_t, mtb2_t, mtb3_t)
rm(mtb1_t, mtb2_t, mtb3_t, mtb1, mtb2, mtb3)

mtb$Compound.Name <- mtb$Compound
mtb$Compound <- paste(mtb$Method, mtb$Compound, sep = "_")
mtb$Method <- NULL

# Create metabolite mapping file
mtb.map <- mtb %>% select(Compound, Compound.Name)
mtb$Compound.Name <- NULL

# We will mark annotations with low-confidence as FALSE
mtb.map$High.Confidence.Annotation <- TRUE
mtb.map$HMDB <- NA
mtb.map$KEGG <- NA

# Manual mappings
mtb.map[mtb.map$Compound.Name == '2-Methylbutyrate','HMDB'] <- 'HMDB0002176'
mtb.map[mtb.map$Compound.Name == '2-Methylbutyrate','KEGG'] <- 'C18319'
mtb.map[mtb.map$Compound.Name == 'Lactate','HMDB'] <- 'HMDB0000190'
mtb.map[mtb.map$Compound.Name == 'Lactate','KEGG'] <- 'C00256'
mtb.map[mtb.map$Compound.Name == 'Alanine','HMDB'] <- 'HMDB0000161'
mtb.map[mtb.map$Compound.Name == 'Alanine','KEGG'] <- 'C00041'
mtb.map[mtb.map$Compound.Name == 'Tyrosine','HMDB'] <- 'HMDB0000158'
mtb.map[mtb.map$Compound.Name == 'Tyrosine','KEGG'] <- 'C00082'
mtb.map[mtb.map$Compound.Name == 'Isoleucine','HMDB'] <- 'HMDB0000172'
mtb.map[mtb.map$Compound.Name == 'Isoleucine','KEGG'] <- 'C00407'
mtb.map[mtb.map$Compound.Name == 'Leucine','HMDB'] <- 'HMDB0000687'
mtb.map[mtb.map$Compound.Name == 'Leucine','KEGG'] <- 'C00123'
mtb.map[mtb.map$Compound.Name == 'Valine','HMDB'] <- 'HMDB0000883'
mtb.map[mtb.map$Compound.Name == 'Valine','KEGG'] <- 'C00183'
mtb.map[mtb.map$Compound.Name == 'Lysine','HMDB'] <- 'HMDB0000182'
mtb.map[mtb.map$Compound.Name == 'Lysine','KEGG'] <- 'C00047'
mtb.map[mtb.map$Compound.Name == 'Succinate','HMDB'] <- 'HMDB0000254'
mtb.map[mtb.map$Compound.Name == 'Succinate','KEGG'] <- 'C00042'
mtb.map[mtb.map$Compound.Name == 'Glycine','HMDB'] <- 'HMDB0000123'
mtb.map[mtb.map$Compound.Name == 'Glycine','KEGG'] <- 'C00037'
mtb.map[mtb.map$Compound.Name == 'b_arabinose','HMDB'] <- 'HMDB0000646'
mtb.map[mtb.map$Compound.Name == 'b_arabinose','KEGG'] <- 'C02479'
#mtb.map[mtb.map$Compound.Name == 'b_xylose','HMDB'] <- 'HMDB0000000'
#mtb.map[mtb.map$Compound.Name == 'b_xylose','KEGG'] <- 'C00000'
mtb.map[mtb.map$Compound.Name == 'Acetate','HMDB'] <- 'HMDB0000042'
mtb.map[mtb.map$Compound.Name == 'Acetate','KEGG'] <- 'C00033'
mtb.map[mtb.map$Compound.Name == 'Propionate','HMDB'] <- 'HMDB0000237'
mtb.map[mtb.map$Compound.Name == 'Propionate','KEGG'] <- 'C00163'
mtb.map[mtb.map$Compound.Name == 'Butyrate','HMDB'] <- 'HMDB0000039'
mtb.map[mtb.map$Compound.Name == 'Butyrate','KEGG'] <- 'C00246'
mtb.map[mtb.map$Compound.Name == 'Glucose','HMDB'] <- 'HMDB0000122'
mtb.map[mtb.map$Compound.Name == 'Glucose','KEGG'] <- 'C00031'
mtb.map[mtb.map$Compound.Name == 'Isovalerate','HMDB'] <- 'HMDB0000718'
mtb.map[mtb.map$Compound.Name == 'Isovalerate','KEGG'] <- 'C08262'
mtb.map[mtb.map$Compound.Name == 'Uracil','HMDB'] <- 'HMDB0000300'
mtb.map[mtb.map$Compound.Name == 'Uracil','KEGG'] <- 'C00106'
mtb.map[mtb.map$Compound.Name == 'Hypoxanthine','HMDB'] <- 'HMDB0000157'
mtb.map[mtb.map$Compound.Name == 'Hypoxanthine','KEGG'] <- 'C00262'
mtb.map[mtb.map$Compound.Name == "5-hydroxyindole acetate",'HMDB'] <- 'HMDB0000763'
mtb.map[mtb.map$Compound.Name == "5-hydroxyindole acetate",'KEGG'] <- 'C05635'
mtb.map[mtb.map$Compound.Name == "melatonin (pg/mg)",'HMDB'] <- 'HMDB0001389'
mtb.map[mtb.map$Compound.Name == "melatonin (pg/mg)",'KEGG'] <- 'C01598'
mtb.map[mtb.map$Compound.Name == "anthranilate (pg/mg)",'HMDB'] <- 'HMDB0001123'
mtb.map[mtb.map$Compound.Name == "anthranilate (pg/mg)",'KEGG'] <- 'C00108'
mtb.map[mtb.map$Compound.Name == "3-methylindole (relative; raw peak area only)",'HMDB'] <- 'HMDB0000466'
mtb.map[mtb.map$Compound.Name == "3-methylindole (relative; raw peak area only)",'KEGG'] <- 'C08313'
mtb.map[mtb.map$Compound.Name == "Kynurenine (no IS)",'HMDB'] <- 'HMDB0000684'
mtb.map[mtb.map$Compound.Name == "Kynurenine (no IS)",'KEGG'] <- 'C00328'
mtb.map[mtb.map$Compound.Name == "indole-3-propionic acid",'HMDB'] <- 'HMDB0002302'
mtb.map[mtb.map$Compound.Name == "Kynurenate (no IS)",'HMDB'] <- 'HMDB0000715'
mtb.map[mtb.map$Compound.Name == "Kynurenate (no IS)",'KEGG'] <- 'C01717'
mtb.map[mtb.map$Compound.Name == "tryptamine",'HMDB'] <- 'HMDB0000303'
mtb.map[mtb.map$Compound.Name == "tryptamine",'KEGG'] <- 'C00398'
mtb.map[mtb.map$Compound.Name == "Deoxycholic acid",'HMDB'] <- 'HMDB0000626'
mtb.map[mtb.map$Compound.Name == "Deoxycholic acid",'KEGG'] <- 'C04483'
mtb.map[mtb.map$Compound.Name == "Lithocholic acid",'HMDB'] <- 'HMDB0000761'
mtb.map[mtb.map$Compound.Name == "indole-3-acetic acid (no IS)",'HMDB'] <- 'HMDB0000197'
mtb.map[mtb.map$Compound.Name == "indole-3-acetic acid (no IS)",'KEGG'] <- 'C00954'
mtb.map[mtb.map$Compound.Name == "Cholic acid",'HMDB'] <- 'HMDB0000619'
mtb.map[mtb.map$Compound.Name == "Cholic acid",'KEGG'] <- 'C00695'
mtb.map[mtb.map$Compound.Name == "Chenodeoxycholic acid",'HMDB'] <- 'HMDB0000518'
mtb.map[mtb.map$Compound.Name == "Chenodeoxycholic acid",'KEGG'] <- 'C02528'
mtb.map[mtb.map$Compound.Name == "Tryptophan",'HMDB'] <- 'HMDB0000929'
mtb.map[mtb.map$Compound.Name == "Tryptophan",'KEGG'] <- 'C00078'
mtb.map[mtb.map$Compound.Name == "Glycocholic acid",'HMDB'] <- 'HMDB0000138'
mtb.map[mtb.map$Compound.Name == "Glycodeoxycholic acid",'HMDB'] <- 'HMDB0000631'
mtb.map[mtb.map$Compound.Name == "Taurocholic acid",'HMDB'] <- 'HMDB0000036'
mtb.map[mtb.map$Compound.Name == "Taurocholic acid",'KEGG'] <- 'C05122'
mtb.map[mtb.map$Compound.Name == "Taurochenodeoxycholic acid",'HMDB'] <- 'HMDB0000951'
mtb.map[mtb.map$Compound.Name == "Taurochenodeoxycholic acid",'KEGG'] <- 'C05465'
mtb.map[mtb.map$Compound.Name == 'Serotonin','HMDB'] <- 'HMDB0000259'
mtb.map[mtb.map$Compound.Name == 'Serotonin','KEGG'] <- 'C00780'
mtb.map[mtb.map$Compound.Name == 'indole','HMDB'] <- 'HMDB0000738'
mtb.map[mtb.map$Compound.Name == 'indole','KEGG'] <- 'C00463'
mtb.map[mtb.map$Compound.Name == "Ursodeoxycholic acid",'HMDB'] <- 'HMDB0000946'
mtb.map[mtb.map$Compound.Name == "Ursodeoxycholic acid",'KEGG'] <- 'C07880'
mtb.map[mtb.map$Compound.Name == "Taurodeoxycholic acid",'HMDB'] <- 'HMDB0000896'
mtb.map[mtb.map$Compound.Name == "Taurodeoxycholic acid",'KEGG'] <- 'C05463'

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

sample.intersect <- Reduce(intersect, list(names(mtb)[-1], names(species)[-1], names(genera)[-1], metadata$Sample))
message(paste(length(sample.intersect),"samples have all data types"))

mtb <- mtb[,c("Compound",sample.intersect)]
species <- species[,c("Species",sample.intersect)]
genera <- genera[,c("Genus",sample.intersect)]
metadata <- metadata[metadata$Sample %in% sample.intersect,]

# --------------------------------
# Save to files + R objects
# --------------------------------

save.to.files(DATASET_NAME, "prelim_data", metadata = metadata, mtb = mtb, mtb.map = mtb.map, genera = genera, species = species)
save.to.rdata(DATASET_NAME, "prelim_data", metadata = metadata, mtb = mtb, mtb.map = mtb.map, genera = genera, species = species, override.all = T)

rm(list = ls())
