library(dada2)
library(dplyr)
library(readr)

gtdb.for.dada <- "../references/gtdb/sbdi_16s_gtdb_database/gtdb-sbdi-sativa.r07rs207.1genome.assignTaxonomy.fna.gz"
datasets <- c("JACOBS_IBD_FAMILIES_2016","HE_INFANTS_MFGM_2019","KIM_ADENOMAS_2020","KOSTIC_INFANTS_DIABETES_2015","POYET_BIO_ML_2019","KANG_AUTISM_2017","WANDRO_PRETERMS_2018")

for (dataset in datasets) {
  message(paste("Processing", dataset))
  
  # File names
  seqs.file <- file.path("../data/original_data", dataset, "dna-sequences.fasta")
  asv.file <- file.path("../data/original_data", dataset, "feature_table_asv.tsv")
  
  # Read ASVs (previously processed using QIIME2)
  seqs <- getSequences(seqs.file)
  seqs.df <- data.frame(seqs) %>%
    tibble::rownames_to_column(var = "asv.id") %>%
    rename(asv = 2)
  
  # Get taxonomy using GTDB reference file
  tax <- assignTaxonomy(seqs, gtdb.for.dada, taxLevels = c("Kingdom", "Domain", "Phylum", "Class", "Order", "Family", "Genus"))
  
  # Organize taxonomy mappings
  tax <- tax %>% 
    replace(is.na(.), "") %>%
    data.frame() %>%
    tibble::rownames_to_column(var = "asv") %>%
    mutate(taxon = paste0('d__', Domain,
                          ';p__', Phylum,
                          ';c__', Class,
                          ';o__', Order,
                          ';f__', Family,
                          ';g__', Genus)) %>%
    left_join(seqs.df, by = "asv") %>%
    select(asv.id, taxon)
  
  # Read ASV feature table and map to assign taxonomy
  asvs <- read_delim(asv.file, skip = 1, delim = '\t', show_col_types = FALSE)
  message(paste("Sanity:",nrow(asvs)))
  asvs <- asvs %>%
    rename(asv.id = 1) %>%
    left_join(tax, by = "asv.id") %>%
    relocate(taxon) %>%
    select(-asv.id) 
  message(paste("Sanity (should be the same as above):",nrow(asvs)))
  
  # Collapse to genus abundances
  genera <- asvs %>%
    group_by(taxon) %>% 
    summarise(across(everything(), sum))
  message(paste("Sanity (should be less than above):",nrow(genera)))
  
  # Save
  write_delim(genera, file.path("../data/original_data", dataset, "feature_table_gtdb_207.tsv"), delim = "\t")
}

message("Done")

