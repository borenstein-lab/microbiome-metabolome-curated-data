# ----------------------------------------------------------------
# The following script loads and processes data from:
# 
# Franzosa, Eric A., et al. "Gut microbiome structure 
#  and metabolic activity in inflammatory bowel disease." 
#  Nature microbiology 4.2 (2019): 293-305.
# 
# The output tables are:
# 1) metadata - subject metadata
# 2) genera - taxonomic profiles per subject (genus-level)
# 3) species - taxonomic profiles per subject (species-level)
# 4) mtb - metabolomic profiles per subject       
# 5) mtb.map - additional details per metabolite 
# ----------------------------------------------------------------

require(gdata)
require(readr)
require(dplyr)

source("load_original_data/utils.R")

# --------------------------------
# Required files & info
# --------------------------------

# For details about the source of each file below see: <COMPLETE>
METADATA_FILE <- '../data/original_data/FRANZOSA_IBD_2019/Metobolites.xlsx'  
TAXONOMY_FILE <- '../data/original_data/FRANZOSA_IBD_2019/Microbial_species_relative_abundances.xlsx'
METABOLOMICS_INFO <- '../data/original_data/FRANZOSA_IBD_2019/Metabolites_metadata.xlsx'

PUBLICATION_NAME <- 'Gut microbiome structure and metabolic activity in inflammatory bowel disease'
DOI <- '10.1038/s41564-018-0306-4'
DATASET_NAME <- 'FRANZOSA_IBD_2019'

# --------------------------------
# Load metadata
# --------------------------------

supp.tmp <- read.xls(METADATA_FILE,
                     header = TRUE,
                     skip = 1,
                     na.strings = c("0","#N/A"), 
                     stringsAsFactors = FALSE)
# Split metadata rows
metadata <- supp.tmp[1:7,]

# Transpose
rownames(metadata) <- metadata$X..Feature...Sample
metadata$X..Feature...Sample <- NULL
metadata <- as.data.frame(t(metadata), stringsAsFactors = FALSE)
metadata$Sample <- rownames(metadata)
rownames(metadata) <- NULL

# Organize column names & order
metadata <- metadata %>%
  mutate(Subject = Sample) %>%
  rename(Study.Group = Diagnosis) %>%
  mutate(Age = as.numeric(Age)) %>%
  mutate(Age.Units = "Years") %>%
  mutate(Fecal.Calprotectin = as.numeric(Fecal.Calprotectin)) %>%
  # Add publication info
  mutate(Dataset = DATASET_NAME) %>%
  mutate(DOI = DOI) %>%
  mutate(Publication.Name = PUBLICATION_NAME) %>%
  # Reorder
  relocate(c(Dataset, Sample, Subject, Study.Group, Age, Age.Units, DOI, Publication.Name)) 

# --------------------------------
# Load taxonomic profiles (species)
# --------------------------------

# Created by Metaphlan
mtg <- read.xls(TAXONOMY_FILE,
                header = TRUE,
                skip = 1,
                na.strings = c("#N/A"), 
                stringsAsFactors=FALSE)
names(mtg)[1] <- 'OTU'
mtg <- mtg[-(1:8),]

# Convert to numeric
for(i in c(2:ncol(mtg))) mtg[,i] <- as.numeric(as.character(mtg[,i]))

# Sanity (abundances sum to ~1): table(cut(colSums(mtg[,2:ncol(mtg)]), breaks = c(0,0.96,0.97,0.98,0.99,1.01)))

# --------------------------------
# Load metabolomic profiles
# --------------------------------

mtb <- supp.tmp[-(1:7),]
names(mtb)[1] <- 'Compound'

# Convert to numeric
for(i in c(2:ncol(mtb))) mtb[,i] <- as.numeric(as.character(mtb[,i]))

# Read metabolite info 
mtb.map <- read.xls(METABOLOMICS_INFO,
                     header = TRUE,
                     skip = 1,
                     na.strings = c("#N/A"), 
                     stringsAsFactors=FALSE)
names(mtb.map)[1] <- 'Compound'
names(mtb.map)[6] <- 'Compound.Name'
mtb <- merge(mtb, mtb.map[,c("Compound","Compound.Name")], by='Compound', all.x=TRUE, all.y=FALSE)

# We change the unique identifier to also include annotated name if exists
mtb$Compound <- paste0(mtb$Compound,": ",mtb$Compound.Name)
mtb.map$Compound <- paste0(mtb.map$Compound,": ",
                           mtb.map$Compound.Name)
mtb$Compound.Name <- NULL

# Build a mapping table which we will gradually fill with mappings. 
mtb.map$Compound.Name = trimws(mtb.map$Compound.Name)
mtb.map$High.Confidence.Annotation <- TRUE

# Mark ambiguous annotations as low confidence
ambig.rows <- grepl("/", mtb.map$Compound) | grepl("\\*", mtb.map$Compound)
mtb.map$High.Confidence.Annotation[ambig.rows] <- FALSE

# Get mapping to KEGG/other compound ID's, using MetaboAnalyst 
cmpds.to.search <- mtb.map$Compound.Name[!is.na(mtb.map$Compound.Name)]
cmpds.to.search <- unique(cmpds.to.search)
MA.matches <- map.compound.names.MetaboAnalyst(cmpds.to.search)

# Merge with the names in our data
mtb.map <- merge(mtb.map, MA.matches, 
                 by.x = "Compound.Name", 
                 by.y = "Query", all = T)

# Sanity - should not return anything: 
#  MA.matches$Query[!MA.matches$Query %in% mtb.map$Compound.Name]

# Reorder columns
mtb.map <- mtb.map %>%
  relocate(Compound, HMDB, KEGG, High.Confidence.Annotation) %>%
  select(-MA.Name.Match)

# Add manual mappings - exact matches in KEGG + verified by mass
# View unmapped ones: 
#  View(mtb.map[!is.na(mtb.map$Compound.Name) & (is.na(mtb.map$HMDB) | is.na(mtb.map$KEGG)),c("Compound","m.z","HMDB","KEGG")])
# Verify drop in na annotations after running the below: 
#  sum(is.na(mtb.map$KEGG)) + sum(is.na(mtb.map$HMDB))
mtb.map[mtb.map$Compound == "HILIC-neg_Cluster_0169: 2-hydroxyglutarate","KEGG"] <- 'C02630'
mtb.map[mtb.map$Compound == "C18-neg_Cluster_0927: alpha-muricholate","KEGG"] <- 'C17647'
mtb.map[mtb.map$Compound == 'HILIC-neg_Cluster_0480: 1-3-7-trimethylurate','HMDB'] <- 'HMDB0002123'	
mtb.map[mtb.map$Compound == 'HILIC-neg_Cluster_0480: 1-3-7-trimethylurate','KEGG'] <- 'C16361'	
mtb.map[mtb.map$Compound == 'HILIC-neg_Cluster_0622: 1,2,3,4-tetrahydro-1-methyl-beta-carboline-3-carboxylic acid','HMDB'] <- 'HMDB0037942'
mtb.map[mtb.map$Compound == 'C18-neg_Cluster_0393: 12.13-diHOME','HMDB'] <- 'HMDB0004705'	
mtb.map[mtb.map$Compound == 'C18-neg_Cluster_0393: 12.13-diHOME','KEGG'] <- 'C14829'	
mtb.map[mtb.map$Compound == 'HILIC-neg_Cluster_0169: 2-hydroxyglutarate','HMDB'] <- 'HMDB0059655'		
mtb.map[mtb.map$Compound == 'HILIC-neg_Cluster_0220: 3-methyladipate-pimelate','HMDB'] <- 'HMDB0000555'		
mtb.map[mtb.map$Compound == 'HILIC-neg_Cluster_0220: 3-methyladipate-pimelate','High.Confidence.Annotation'] <- FALSE
mtb.map[mtb.map$Compound == 'HILIC-pos_Cluster_0177: 4-hydroxy-3-methylacetophenone','HMDB'] <- 'HMDB0059824'		
mtb.map[mtb.map$Compound == 'C8-pos_Cluster_0100: 7-hexadecenoic acid methyl ester','HMDB'] <- 'HMDB0061858'		
mtb.map[mtb.map$Compound == 'C18-neg_Cluster_0394: 9.10-diHOME','HMDB'] <- 'HMDB0004704'		
mtb.map[mtb.map$Compound == 'HILIC-neg_Cluster_0103: acetylalanine','KEGG'] <- 'C01073'	
mtb.map[mtb.map$Compound == 'HILIC-neg_Cluster_0103: acetylalanine','High.Confidence.Annotation'] <- FALSE
mtb.map[mtb.map$Compound == 'C18-neg_Cluster_0090: acetytyrosine','HMDB'] <- 'HMDB0000866'	
mtb.map[mtb.map$Compound == 'C18-neg_Cluster_0090: acetytyrosine','KEGG'] <- 'C01657'	
mtb.map[mtb.map$Compound == 'HILIC-pos_Cluster_0298: dimethyllysine','HMDB'] <- 'HMDB0013287'	
mtb.map[mtb.map$Compound == 'HILIC-pos_Cluster_0298: dimethyllysine','KEGG'] <- 'C05545'
mtb.map[mtb.map$Compound == 'HILIC-pos_Cluster_0298: dimethyllysine','High.Confidence.Annotation'] <- FALSE
mtb.map[mtb.map$Compound == 'HILIC-pos_Cluster_0913: ethyl 9-hexadecenoate','HMDB'] <- 'HMDB0059871'		
mtb.map[mtb.map$Compound == 'HILIC-pos_Cluster_0657: geranyl acetoacetate*','HMDB'] <- 'HMDB0038256'		
mtb.map[mtb.map$Compound == 'HILIC-pos_Cluster_0657: geranyl acetoacetate*','High.Confidence.Annotation'] <- FALSE
mtb.map[mtb.map$Compound == 'HILIC-pos_Cluster_0658: geranyl acetoacetate*','HMDB'] <- 'HMDB0038256'		
mtb.map[mtb.map$Compound == 'HILIC-pos_Cluster_0658: geranyl acetoacetate*','High.Confidence.Annotation'] <- FALSE
mtb.map[mtb.map$Compound == 'HILIC-neg_Cluster_0401: glucurote','HMDB'] <- 'HMDB0000127'	
mtb.map[mtb.map$Compound == 'HILIC-neg_Cluster_0401: glucurote','KEGG'] <- 'C00191'	
mtb.map[mtb.map$Compound == 'C18-neg_Cluster_0900: ketodeoxycholate','HMDB'] <- 'HMDB0000391'		
mtb.map[mtb.map$Compound == 'C18-neg_Cluster_0900: ketodeoxycholate','High.Confidence.Annotation'] <- FALSE
mtb.map[mtb.map$Compound == 'HILIC-pos_Cluster_0521: L-1,2,3,4-tetrahydro-beta-carboline-3-carboxylic acid*','HMDB'] <- 'HMDB0035665'		
mtb.map[mtb.map$Compound == 'HILIC-neg_Cluster_0247: L-phenyllactic acid','HMDB'] <- 'HMDB0000748'	
mtb.map[mtb.map$Compound == 'HILIC-neg_Cluster_0247: L-phenyllactic acid','KEGG'] <- 'C05607'	
mtb.map[mtb.map$Compound == 'HILIC-neg_Cluster_0247: L-phenyllactic acid','High.Confidence.Annotation'] <- FALSE
mtb.map[mtb.map$Compound == 'HILIC-pos_Cluster_0082: N-methylproline','HMDB'] <- 'HMDB0094696'		
mtb.map[mtb.map$Compound == 'HILIC-neg_Cluster_0552: pantothete','HMDB'] <- 'HMDB0000210'	
mtb.map[mtb.map$Compound == 'HILIC-neg_Cluster_0552: pantothete','KEGG'] <- 'C00864'	
mtb.map[mtb.map$Compound == 'HILIC-neg_Cluster_0069: succite','HMDB'] <- 'HMDB0000254'	
mtb.map[mtb.map$Compound == 'HILIC-neg_Cluster_0069: succite','KEGG'] <- 'C00042'	
mtb.map[mtb.map$Compound == 'HILIC-neg_Cluster_0069: succite','High.Confidence.Annotation'] <- FALSE
mtb.map[mtb.map$Compound == 'HILIC-pos_Cluster_2093: urobilin*','HMDB'] <- 'HMDB0004160'	
mtb.map[mtb.map$Compound == 'HILIC-pos_Cluster_2093: urobilin*','KEGG'] <- 'C05793'	

# Mark cases of duplicated HMDN/KEGG ID as lower confidence
kegg.dups <- names(table(mtb.map$KEGG)[table(mtb.map$KEGG) > 1])
hmdb.dups <- names(table(mtb.map$HMDB)[table(mtb.map$HMDB) > 1])
mtb.map$High.Confidence.Annotation[mtb.map$KEGG %in% kegg.dups] <- FALSE
mtb.map$High.Confidence.Annotation[mtb.map$HMDB %in% hmdb.dups] <- FALSE

# --------------------------------
# Keep only samples with all data
# --------------------------------

sample.intersect <- Reduce(intersect, list(names(mtb)[-1], names(mtg)[-1], metadata$Sample))
message(paste(length(sample.intersect),"samples have all data types"))

mtb <- mtb[,c("Compound",sample.intersect)]
mtg <- mtg[,c("OTU",sample.intersect)]
metadata <- metadata[metadata$Sample %in% sample.intersect,]

# --------------------------------
# Collapse species to genera
# --------------------------------

# Get metaphlan species name mapper
species.mapping <- get.metaphlan.species.mapper()

# Convert current table to a species table + genus-level table
l <- get.genus.level(mtg, species.mapping)
species <- l$species
genera <- l$genera

# Remove Archea
genera <- genera[!grepl("k__Archaea",genera$Genus),]
species <- species[!grepl("k__Archaea",species$Species),]

# --------------------------------
# Save to files + R objects
# --------------------------------

save.to.files(DATASET_NAME, metadata, mtb, mtb.map, genera, species)
save.to.rdata(DATASET_NAME, metadata, mtb, mtb.map, genera, species)
rm(list = ls())

