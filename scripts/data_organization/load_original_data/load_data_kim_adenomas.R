# ----------------------------------------------------------------
# The following script loads and processes data from:
# 
# Kim, Minsuk, et al. "Fecal metabolomic signatures in colorectal 
#  adenoma patients are associated with gut microbiota and early 
#  events of colorectal cancer pathogenesis." 
#  MBio 11.1 (2020): e03186-19.
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
METADATA_FILE <- "../data/original_data/KIM_ADENOMAS_2020/mBio_Metadata_16S_mapping.csv"
TAXONOMY_FILE <- "../data/original_data/KIM_ADENOMAS_2020/feature_table_gtdb.tsv"
METABOLOMICS_FILE <- "../data/original_data/KIM_ADENOMAS_2020/inline-supplementary-material-6.xlsx"

PUBLICATION_NAME <- 'Fecal Metabolomic Signatures in Colorectal Adenoma Patients Are Associated with Gut Microbiota and Early Events of Colorectal Cancer Pathogenesis'
DOI <- '10.1128/mBio.03186-19'
DATASET_NAME <- 'KIM_ADENOMAS_2020'

# --------------------------------
# Load metadata
# --------------------------------

metadata <- read_csv(METADATA_FILE)

# Organize column names & order
metadata <- metadata %>%
  rename(Sample = 1) %>%
  rename(Fastq.Sample.ID = 2) %>%
  mutate(Subject = Sample) %>%
  rename(Study.Group = Group) %>%
  rename(Age = `Age group`) %>%
  mutate(Age.Units = "Years") %>%
  rename(Gender = Sex) %>%
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

genera$Genus[genera$Genus %in% c("Ambiguous_taxa;Ambiguous_taxa;Ambiguous_taxa;Ambiguous_taxa;Ambiguous_taxa;Ambiguous_taxa",
                                 "D_0__Bacteria;__;__;__;__;__",
                                 "Unassigned;__;__;__;__;__")] <- "Unclassified"

# --------------------------------
# Load metabolomic profiles
# --------------------------------

# Read metab' table from excel
mtb <- read.xls(METABOLOMICS_FILE, 
                sheet = 1,
                header = TRUE, 
                stringsAsFactors = FALSE)
names(mtb)[1] <- "Compound"

mtb.map <- mtb %>% select(1:4)
mtb <- mtb %>% select(-2:-4)
mtb.map$High.Confidence.Annotation <- TRUE

# Sanity: all mtbs are unique: sum(duplicated(mtb$Compound)) == 0

# Mark ambiguous metabolites
mtb.map$High.Confidence.Annotation[grepl("[a-zA-Z]/", mtb.map$Compound)] <- FALSE

# Mark starred metabolites with low conf (?)
mtb.map$High.Confidence.Annotation[grepl("\\*$", mtb.map$Compound)] <- FALSE

# Get mapping to KEGG/other compound ID's, using MetaboAnalyst 
cmpds.to.search <- mtb.map$Compound
MA.matches <- map.compound.names.MetaboAnalyst(cmpds.to.search)

# Merge with the names in our data
mtb.map <- merge(mtb.map, MA.matches, 
                 by.x = "Compound", 
                 by.y = "Query", all = T)
mtb.map$MA.Name.Match <- NULL

# Manually add more mappings
## View(mtb.map[is.na(mtb.map$HMDB),c("Compound","High.Confidence.Annotation","HMDB","KEGG")])
## sum(is.na(mtb.map$HMDB)); sum(is.na(mtb.map$KEGG)); 
mtb.map$HMDB[which(mtb.map$Compound == "1-oleoyl-GPE (18:1)")] <- "HMDB0011506"
mtb.map$HMDB[which(mtb.map$Compound == "1-oleoyl-GPC (18:1)")] <- "HMDB0240602"
mtb.map$HMDB[which(mtb.map$Compound == "1,2-dilinoleoyl-GPC (18:2/18:2)")] <- "HMDB0008138"
mtb.map$HMDB[which(mtb.map$Compound == "1,2-dipalmitoyl-GPC (16:0/16:0)")] <- "HMDB0000564"
mtb.map$HMDB[which(mtb.map$Compound == "1,7-dimethylurate")] <- "HMDB0011103"
mtb.map$HMDB[which(mtb.map$Compound == "1-(1-enyl-oleoyl)-GPE (P-18:1)*")] <- "HMDB0240599"
mtb.map$HMDB[which(mtb.map$Compound == "1-(1-enyl-stearoyl)-GPE (P-18:0)*")] <- "HMDB0240598"
mtb.map$HMDB[which(mtb.map$Compound == "1-linoleoyl-GPE (18:2)*")] <- "HMDB0011507"
mtb.map$HMDB[which(mtb.map$Compound == "1-oleoylglycerol (18:1)")] <- "HMDB0011567"
mtb.map$HMDB[which(mtb.map$Compound == "1-palmitoyl-GPG (16:0)*")] <- "HMDB0240601"
mtb.map$HMDB[which(mtb.map$Compound == "1-stearoyl-GPC (18:0)")] <- "HMDB0010384"
mtb.map$HMDB[which(mtb.map$Compound == "1-stearoyl-GPE (18:0)")] <- "HMDB0011130"
mtb.map$HMDB[which(mtb.map$Compound == "10-hydroxystearate")] <- "HMDB0037396"
mtb.map$HMDB[which(mtb.map$Compound == "17-methylstearate")] <- "HMDB0037397"
mtb.map$HMDB[which(mtb.map$Compound == "2-keto-3-deoxy-gluconate")] <- "HMDB0001353"
mtb.map$HMDB[which(mtb.map$Compound == "3-hydroxybutyrate (BHBA)")] <- "HMDB0000357"
mtb.map$HMDB[which(mtb.map$Compound == "3-hydroxypalmitate")] <- "HMDB0061658"
mtb.map$HMDB[which(mtb.map$Compound == "7-methylguanine")] <- "HMDB0000897"
mtb.map$HMDB[which(mtb.map$Compound == "8-hydroxyguanine")] <- "HMDB0002032"
mtb.map$HMDB[which(mtb.map$Compound == "9-HOTrE")] <- "HMDB0010224"
mtb.map$HMDB[which(mtb.map$Compound == "argininate*")] <- "HMDB0003148"
mtb.map$HMDB[which(mtb.map$Compound == "azelate (nonanedioate)")] <- "HMDB0000784"
mtb.map$HMDB[which(mtb.map$Compound == "beta-cryptoxanthin")] <- "HMDB0033844"
mtb.map$HMDB[which(mtb.map$Compound == "bilirubin (Z,Z)")] <- "HMDB0000054"
mtb.map$HMDB[which(mtb.map$Compound == "biliverdin")] <- "HMDB0001008"
mtb.map$HMDB[which(mtb.map$Compound == "glutamate, gamma-methyl ester")] <- "HMDB0061715"
mtb.map$HMDB[which(mtb.map$Compound == "carnitine")] <- "HMDB0000062"
mtb.map$HMDB[which(mtb.map$Compound == "choline phosphate")] <- "HMDB0001565"
mtb.map$HMDB[which(mtb.map$Compound == "creatinine")] <- "HMDB0000562"
mtb.map$HMDB[which(mtb.map$Compound == "cytidine")] <- "HMDB0000089"
mtb.map$HMDB[which(mtb.map$Compound == "dihomo-linoleate (20:2n6)")] <- "HMDB0061864"
mtb.map$HMDB[which(mtb.map$Compound == "dihydroferulic acid")] <- "HMDB0062121"
mtb.map$HMDB[which(mtb.map$Compound == "docosadienoate (22:2n6)")] <- "HMDB0061714"
mtb.map$HMDB[which(mtb.map$Compound == "enterolactone")] <- "HMDB0006101"
mtb.map$HMDB[which(mtb.map$Compound == "fucose")] <- "HMDB0000174"
mtb.map$HMDB[which(mtb.map$Compound == "glycerophosphoglycerol")] <- "HMDB0240316"
mtb.map$HMDB[which(mtb.map$Compound == "guanine")] <- "HMDB0000132"
mtb.map$HMDB[which(mtb.map$Compound == "guanosine")] <- "HMDB0000133"
mtb.map$HMDB[which(mtb.map$Compound == "hyocholate")] <- "HMDB0000760"
mtb.map$HMDB[which(mtb.map$Compound == "hypoxanthine")] <- "HMDB0000157"
mtb.map$HMDB[which(mtb.map$Compound == "imidazole propionate")] <- "HMDB0002271"
mtb.map$HMDB[which(mtb.map$Compound == "inosine")] <- "HMDB0000195"
mtb.map$HMDB[which(mtb.map$Compound == "lactosyl-N-palmitoyl-sphingosine (d18:1/16:0)")] <- "HMDB0006750"
mtb.map$HMDB[which(mtb.map$Compound == "pyrraline")] <- "HMDB0033143"
mtb.map$HMDB[which(mtb.map$Compound == "pterin")] <- "HMDB0000802"
mtb.map$HMDB[which(mtb.map$Compound == "maltose")] <- "HMDB0000163"
mtb.map$HMDB[which(mtb.map$Compound == "methionine sulfoxide")] <- "HMDB0002005"
mtb.map$HMDB[which(mtb.map$Compound == "N('1)-acetylspermidine")] <- "HMDB0001276"
mtb.map$HMDB[which(mtb.map$Compound == "N-alpha-acetylornithine")] <- "HMDB0003357"
mtb.map$HMDB[which(mtb.map$Compound == "sphingosine")] <- "HMDB0000252"
mtb.map$HMDB[which(mtb.map$Compound == "stachydrine")] <- "HMDB0004827"
mtb.map$HMDB[which(mtb.map$Compound == "N-acetylhistidine")] <- "HMDB0032055"
mtb.map$HMDB[which(mtb.map$Compound == "N-acetylmethionine")] <- "HMDB0011745"
mtb.map$HMDB[which(mtb.map$Compound == "N-formylmethionine")] <- "HMDB0001015"
mtb.map$HMDB[which(mtb.map$Compound == "arabonate")] <- "HMDB0000539"
mtb.map$HMDB[which(mtb.map$Compound == "xylonate")] <- "HMDB0060256"
mtb.map$HMDB[which(mtb.map$Compound == "xylulose")] <- "HMDB0000751"
mtb.map$HMDB[which(mtb.map$Compound == "gamma-tocopherol")] <- "HMDB0001492"
mtb.map$HMDB[which(mtb.map$Compound == "dihomo-a-linolenate")] <- "HMDB0060039"
mtb.map$HMDB[which(mtb.map$Compound == "N6-carboxymethyllysine")] <- "HMDB0240347"
mtb.map$HMDB[which(mtb.map$Compound == "oleoyl ethanolamide")] <- "HMDB0002088"
mtb.map$HMDB[which(mtb.map$Compound == "N-acetylarginine")] <- "HMDB0004620"
mtb.map$HMDB[which(mtb.map$Compound == "N-acetylcitrulline")] <- "HMDB0000856"
mtb.map$HMDB[which(mtb.map$Compound == "N-acetylmethionine sulfoxide")] <- "HMDB0240343"
mtb.map$HMDB[which(mtb.map$Compound == "thioproline")] <- "HMDB0062164"
mtb.map$HMDB[which(mtb.map$Compound == "threonylphenylalanine")] <- "HMDB0029068"
mtb.map$HMDB[which(mtb.map$Compound == "N-delta-acetylornithine")] <- "HMDB0240589"
mtb.map$HMDB[which(mtb.map$Compound == "N-formylphenylalanine")] <- "HMDB0240317"
mtb.map$HMDB[which(mtb.map$Compound == "N-palmitoyl-sphinganine (d18:0/16:0)")] <- "HMDB0011760"
mtb.map$HMDB[which(mtb.map$Compound == "1-oleoyl-2-linoleoyl-GPC (18:1/18:2)*")] <- "HMDB0008105"
mtb.map$HMDB[which(mtb.map$Compound == "2-hydroxypalmitate")] <- "HMDB0031057"
mtb.map$HMDB[which(mtb.map$Compound == "2-hydroxystearate")] <- "HMDB0062549" 
mtb.map$HMDB[which(mtb.map$Compound == "2-linoleoylglycerol (18:2)")] <- "HMDB0011538"
mtb.map$HMDB[which(mtb.map$Compound == "2,3-dimethylsuccinate")] <- "HMDB0245405" 
mtb.map$HMDB[which(mtb.map$Compound == "acisoga")] <- "HMDB0061384"
mtb.map$HMDB[which(mtb.map$Compound == "arachidoylcarnitine (C20)*")] <- "HMDB0006460"
mtb.map$HMDB[which(mtb.map$Compound == "N-acetyl-3-methylhistidine*")] <- "HMDB0240341"
mtb.map$HMDB[which(mtb.map$Compound == "gamma-glutamylisoleucine*")] <- "HMDB0011170" 
mtb.map$HMDB[which(mtb.map$Compound == "glycerophosphoserine*")] <- "HMDB0252850"
mtb.map$HMDB[which(mtb.map$Compound == "leucylglutamine*")] <- "HMDB0028927"
mtb.map$HMDB[which(mtb.map$Compound == "pyroglutamine*")] <- "HMDB0062558" 
mtb.map$HMDB[which(mtb.map$Compound == "riboflavin (Vitamin B2)")] <- "HMDB0000244" 
mtb.map$HMDB[which(mtb.map$Compound == "succinylcarnitine (C4-DC)")] <- "HMDB0061717" 
mtb.map$HMDB[which(mtb.map$Compound == "2-oxoarginine*")] <- "HMDB0004225" 
mtb.map$HMDB[which(mtb.map$Compound == "5-dodecenoate (12:1n7)")] <- "HMDB0000529" 
mtb.map$HMDB[which(mtb.map$Compound == "8-hydroxyguanine ")] <- "HMDB0002032" 
mtb.map$HMDB[which(mtb.map$Compound == "trigonelline (N'-methylnicotinate)")] <- "HMDB0000875" 
mtb.map$HMDB[which(mtb.map$Compound == "phenyllactate (PLA)")] <- "HMDB0000779" 
mtb.map$HMDB[which(mtb.map$Compound == "pimelate (heptanedioate)")] <- "HMDB0000857" 
mtb.map$HMDB[which(mtb.map$Compound == "thiamin (Vitamin B1)")] <- "HMDB0000235" 
mtb.map$HMDB[which(mtb.map$Compound == "nonadecanoate (19:0)")] <- "HMDB0000772"
mtb.map$HMDB[which(mtb.map$Compound == "arachidate (20:0)")] <- "HMDB0002212" 
mtb.map$HMDB[which(mtb.map$Compound == "glycerophosphorylcholine (GPC)")] <- "HMDB0000086" 
mtb.map$HMDB[which(mtb.map$Compound == "N-acetylglucosaminylasparagine ")] <- "HMDB0000489" 
mtb.map$HMDB[which(mtb.map$Compound == "N-acetylaspartate (NAA)")] <- "HMDB0000812" 
mtb.map$HMDB[which(mtb.map$Compound == "leucylglycine")] <- NA # MA bug

mtb.map$KEGG[which(mtb.map$Compound == "N-acetylaspartate (NAA)")] <- "C01042"
mtb.map$KEGG[which(mtb.map$Compound == "methylmalonate (MMA)")] <- "C02170"
mtb.map$KEGG[which(mtb.map$Compound == "homostachydrine*")] <- "C08283"
mtb.map$KEGG[which(mtb.map$Compound == "caprylate (8:0)")] <- "C06423"
mtb.map$KEGG[which(mtb.map$Compound == "docosahexaenoate (DHA; 22:6n3)")] <- "C06429"
mtb.map$KEGG[which(mtb.map$Compound == "N-acetylglucosaminylasparagine ")] <- "C04540"
mtb.map$KEGG[which(mtb.map$Compound == "nonadecanoate (19:0)")] <- "C16535"
mtb.map$KEGG[which(mtb.map$Compound == "palmitate (16:0)")] <- "C00249"
mtb.map$KEGG[which(mtb.map$Compound == "oxalate (ethanedioate)")] <- "C00209"
mtb.map$KEGG[which(mtb.map$Compound == "palmitoleate (16:1n7)")] <- "C08362"
mtb.map$KEGG[which(mtb.map$Compound == "arachidonate (20:4n6)")] <- "C00219"
mtb.map$KEGG[which(mtb.map$Compound == "phenyllactate (PLA)")] <- "C05607"
mtb.map$KEGG[which(mtb.map$Compound == "pimelate (heptanedioate)")] <- "C02656"
mtb.map$KEGG[which(mtb.map$Compound == "thiamin (Vitamin B1)")] <- "C00378"
mtb.map$KEGG[which(mtb.map$Compound == "arachidate (20:0)")] <- "C06425"
mtb.map$KEGG[which(mtb.map$Compound == "glycerophosphorylcholine (GPC)")] <- "C00670"
mtb.map$KEGG[which(mtb.map$Compound == "8-hydroxyguanine ")] <- "C20155"
mtb.map$KEGG[which(mtb.map$Compound == "trigonelline (N'-methylnicotinate)")] <- "C01004"
mtb.map$KEGG[which(mtb.map$Compound == "riboflavin (Vitamin B2)")] <- "C00255"
mtb.map$KEGG[which(mtb.map$Compound == "2-oxoarginine*")] <- "C03771"
mtb.map$KEGG[which(mtb.map$Compound == "1,7-dimethylurate")] <- "C16356"
mtb.map$KEGG[which(mtb.map$Compound == "2,3-dihydroxyisovalerate")] <- "C04039"
mtb.map$KEGG[which(mtb.map$Compound == "3-hydroxybutyrate (BHBA)")] <- "C01089"
mtb.map$KEGG[which(mtb.map$Compound == "azelate (nonanedioate)")] <- "C08261"
mtb.map$KEGG[which(mtb.map$Compound == "bilirubin (Z,Z)")] <- "C00486"
mtb.map$KEGG[which(mtb.map$Compound == "docosadienoate (22:2n6)")] <- "C16533"
mtb.map$KEGG[which(mtb.map$Compound == "inosine")] <- "C00294"
mtb.map$KEGG[which(mtb.map$Compound == "lactosyl-N-palmitoyl-sphingosine (d18:1/16:0)")] <- "C01290"
mtb.map$KEGG[which(mtb.map$Compound == "methionine sulfoxide")] <- "C02989"
mtb.map$KEGG[which(mtb.map$Compound == "N('1)-acetylspermidine")] <- "C00612"
mtb.map$KEGG[which(mtb.map$Compound == "N-alpha-acetylornithine")] <- "C00437"
mtb.map$KEGG[which(mtb.map$Compound == "N-acetylhistidine")] <- "C02997"
mtb.map$KEGG[which(mtb.map$Compound == "N-acetylmethionine")] <- "C02712"
mtb.map$KEGG[which(mtb.map$Compound == "N-formylmethionine")] <- "C03145"
mtb.map$KEGG[which(mtb.map$Compound == "nicotianamine")] <- "C05324"
mtb.map$KEGG[which(mtb.map$Compound == "2-hydroxystearate")] <- "C03045"
mtb.map$KEGG[which(mtb.map$Compound == "leucylglycine")] <- NA # MA bug

mtb.map$High.Confidence.Annotation[which(mtb.map$Compound == "palmitate (16:0)")] <- FALSE
mtb.map$High.Confidence.Annotation[which(mtb.map$Compound == "2-hydroxystearate")] <- FALSE 
mtb.map$High.Confidence.Annotation[which(mtb.map$Compound == "palmitoleate (16:1n7)")] <- FALSE
mtb.map$High.Confidence.Annotation[which(mtb.map$Compound == "5-dodecenoate (12:1n7)")] <- FALSE
mtb.map$High.Confidence.Annotation[which(mtb.map$Compound == "succinylcarnitine (C4-DC)")] <- FALSE
mtb.map$High.Confidence.Annotation[which(mtb.map$Compound == "9-HOTrE")] <- FALSE
mtb.map$High.Confidence.Annotation[which(mtb.map$Compound == "fucose")] <- FALSE
mtb.map$High.Confidence.Annotation[which(mtb.map$Compound == "N-acetylarginine")] <- FALSE
mtb.map$High.Confidence.Annotation[which(mtb.map$Compound == "N-acetylcitrulline")] <- FALSE

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

