# ----------------------------------------------------------------
# Unify genus names in genera-feature-tables.
# Override relevant files (tsv + RData) with the unified version.
# ----------------------------------------------------------------

require(vegan)
require(dplyr)
require(readr)
source("utils.R")

# Load data
all.data <- load.all.datasets("prelim_data")
for(i in 1:length(all.data)) assign(names(all.data)[i], all.data[[i]])
rm(all.data)
datasets <- basename(data.dirs)
# gtdb.map <- get.gtdb.mapper()

# --------------------------------
# 1. Remove non-bacteria
# --------------------------------

# We remove non-bacteria entities 
#  (before re-normalizing each sample to 100%)
# Note: this has already been performed in some datasets, 
#  by authors/ depending on exact metagenomics processing.

for (dataset in datasets) {
  tmp <- genera[[dataset]]
  tmp.n.row <- nrow(tmp) # Recording number before we erase all non-bacteria
  tmp <- tmp[!grepl("Viruses", tmp$Genus, ignore.case = T),]
  tmp <- tmp[!grepl("Archaea", tmp$Genus, ignore.case = T),]
  tmp <- tmp[!grepl("Eukaryota",tmp$Genus, ignore.case = T),]
  message(paste("Removed",tmp.n.row-nrow(tmp),"non-bacteria rows from dataset",dataset))
  genera[[dataset]] <- tmp
}
rm(tmp.n.row)

# --------------------------------
# 2. Transform to relative abund'
# --------------------------------

genera.new <- list()
for (dataset in datasets) {
  tmp <- decostand(genera[[dataset]][,-1], method = "total", MARGIN = 2)
  tmp$Genus <- genera[[dataset]]$Genus
  tmp <- tmp %>% relocate(Genus)
  genera.new[[dataset]] <- tmp
  # Sanity: print(apply(tmp, 2, sum))
}

# --------------------------------
# 3. Unify genera names 
# --------------------------------

# 3.1. Reformat genus taxonomy

# For the qiime-processed datasets, we change from silva-like 
#  convention to metaplhan-like convention.
# E.g. from "D_0__Bacteria;D_1__Actinobacteria;D_2__Actinobacteria;D_3__Actinomycetales;D_4__Actinomycetaceae;D_5__Actinotignum" to: 
#  "d__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Actinomycetales;f__Actinomycetaceae;g__Actinotignum".

for (dataset in datasets) {
  # Reformat the genera names
  new.names <- genera.new[[dataset]]$Genus
  #print(paste("Before:",new.names[1:5]))
  
  new.names <- gsub("\\|",";",new.names)
  new.names <- gsub("D_0__","d__",new.names)
  new.names <- gsub("D_1__","p__",new.names)
  new.names <- gsub("D_2__","c__",new.names)
  new.names <- gsub("D_3__","o__",new.names)
  new.names <- gsub("D_4__","f__",new.names)
  new.names <- gsub("D_5__","g__",new.names)
  
  # GTDB
  new.names <- gsub("k__Bacteria","d__Bacteria",new.names)
  
  # # Other adaptations related to taxa identified at a level higher than genus
  # new.names <- gsub(";Ambiguous_taxa",";__",new.names)
  # new.names <- gsub(";[pcof]__;",";__;",new.names)
  # new.names <- gsub(";g__$",";__",new.names)
  # new.names <- gsub(";f__gut metagenome;",";__;",new.names)
  # new.names <- gsub(";g__gut metagenome$",";__",new.names)
  # new.names <- gsub(";f__metagenome;",";__;",new.names)
  # new.names <- gsub(";g__metagenome$",";__",new.names)
  # new.names <- gsub(";f__uncultured bacterium;",";__;",new.names)
  # new.names <- gsub(";f__uncultured prokaryote;",";__;",new.names)
  # new.names <- gsub(";f__uncultured soil bacterium;",";__;",new.names)
  # new.names <- gsub(";f__uncultured;",";__;",new.names)
  # new.names <- gsub(";g__uncultured bacterium$",";__",new.names)
  # new.names <- gsub(";g__uncultured prokaryote$",";__",new.names)
  # new.names <- gsub(";g__uncultured soil bacterium$",";__",new.names)
  # new.names <- gsub(";g__uncultured$",";__",new.names)
  # new.names <- gsub(";f__uncultured organism;",";__;",new.names)
  # new.names <- gsub(";g__uncultured organism$",";__",new.names)
  # new.names <- gsub(";f__uncultured rumen bacterium;",";__;",new.names)
  # new.names <- gsub(";g__uncultured rumen bacterium$",";__",new.names)
  # new.names <- gsub(";g__.*_noname$",";__",new.names)
  # new.names <- gsub(";f__.*_noname;__$",";__;__",new.names)
  # new.names <- gsub(";o__.*_noname;__;__$",";__;__;__",new.names)
  # new.names <- gsub(";c__.*_noname;__;__;__$",";__;__;__;__",new.names)
  # 
  # # For a few cases where some phylogeny levels are consistently missing:
  # new.names <- gsub(";__;o__Nostocales;",";c__Cyanophyceae;o__Nostocales;",new.names)
  # new.names <- gsub(";__;o__Oscillatoriales;",";c__Cyanophyceae;o__Oscillatoriales;",new.names)
  # new.names <- gsub(";__;o__Pleurocapsales;",";c__Cyanophyceae;o__Pleurocapsales;",new.names)
  # new.names <- gsub(";__;o__Synechococcales;",";c__Cyanophyceae;o__Synechococcales;",new.names)
  # new.names <- gsub(";__;f__Leptospiraceae;",";o__Leptospirales;f__Leptospiraceae;",new.names)
  
  message(paste(dataset, "- reformatted",
                sum(new.names != genera.new[[dataset]]$Genus),
                "out of", length(new.names), "genus entities"))
  
  genera.new[[dataset]]$Genus <- new.names
  #print(paste("After:",new.names[1:5]))
}
rm(new.names)


# 3.2. Regroup unclassified genera
  
# We label unclassified entities as "Unclassified".
# Note: in cases where the entity is classified to a higher-level taxonomy 
#  (e.g. class-level, order-level, etc.), we leave it as is.

# Examples of entities we rename as "Unclassified":
#  "d__Bacteria;__;__;__;__;__"
#  "Unassigned;__;__;__;__;__"

# Save some statistics about this step, for sanity
debug.unclass.genera.stats <- data.frame(Dataset=character(0),
                                        Num_Rows_Unclassified=integer(0),
                                        Rel_abundance__Min=numeric(0),
                                        Rel_abundance__Max=numeric(0),
                                        Rel_abundance__Median=numeric(0))

# NO RUN
# x <- c(); for(dataset in datasets) {x <- c(x, genera.new[[dataset]]$Genus)}; x <- unique(x)
# print(c(grep(";__$", x,value = T), 
#         grep(";g__$", x,value = T), 
#         grep("_noname$", x,value = T), 
#         grep("Unclassified", x,value = T),
#         grep("g__uncultured", x,value = T),
#         grep("Ambiguous_taxa$", x,value = T),
#         grep("g__metagenome$", x,value = T),
#         grep("g__unidentified rumen bacterium", x,value = T),
#         grep("g__gut metagenome$", x,value = T)))

unclass.tax.strings <- c("d__Bacteria;p__;c__;o__;f__;g__",
                         "d__Bacteria;__;__;__;__;__",
                         "Unassigned;__;__;__;__;__",
                         "Ambiguous_taxa;__;__;__;__;__",
                         "Unclassified")

for (dataset in datasets) {
  tmp <- genera.new[[dataset]]
  unclassified.rows <- tmp$Genus[tmp$Genus %in% unclass.tax.strings]
  
  # After identifying all unclassified entities, 
  #  we regroup/merge them into a new entity named "Unclassified" (summing over values). 
  # Lastly we replace the original rows with the new "Unclassified" row.
  if (length(unclassified.rows) > 0) {
    tmp <- tmp %>%
      mutate(Genus = ifelse(Genus %in% unclassified.rows, 
                            "Unclassified",
                            Genus)) %>%
      group_by(Genus) %>% 
      summarise(across(everything(), sum))
    genera.new[[dataset]] <- tmp
    
    # Record statistics for later sanity checks
    tmp <- unlist(tmp[tmp$Genus == "Unclassified",-1])
    debug.unclass.genera.stats <- 
      bind_rows(debug.unclass.genera.stats,
                data.frame(Dataset = dataset,
                           Num_Rows_Unclassified = length(unclassified.rows),
                           Rel_abundance__Min = round(min(tmp),2),
                           Rel_abundance__Max = round(max(tmp),2),
                           Rel_abundance__Median = round(median(tmp),2)))
    
  } else {
    # Add stats row (empty)
    debug.unclass.genera.stats <- 
      bind_rows(debug.unclass.genera.stats, 
                data.frame(Dataset = dataset, 
                           Num_Rows_Unclassified = 0))
  }
}

# Check out statistics: debug.unclass.genera.stats 
rm(unclassified.rows, tmp)

  
# 3.3. 16s/WGS gaps
  
# We perform several manual curations of taxonomic annotations, for cases such as:  
#  * In the metaphlan datasets, we merge Escherichia and Shigella entities, 
#    which cannot be differentiated by 16s using SILVA database. Same for Hafnia/Obesumbacterium.
#  * Wherever silva OTU's are in the species/sub-genus level, we merge the entities to 
#    a genus-level entity for consistency (e.g. Eubacterium, Ruminococcus, ...)  

# A list of pairs, where in each pair there's a "new.string" element defining the unified genus name,
#  and an "old.strings" element including all genera that should be mapped to the "new.string".
strings.to.correct <- list()

strings.to.correct$Escherichia_Shigella <- 
  list(new.string = "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacteriales;f__Enterobacteriaceae;g__Escherichia-Shigella",
       old.strings = c("d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacteriales;f__Enterobacteriaceae;g__Escherichia",
                       "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacteriales;f__Enterobacteriaceae;g__Shigella",
                       "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacteriales;f__Enterobacteriaceae;g__Escherichia-Shigella"))

strings.to.correct$Hafnia <- 
  list(new.string = "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacteriales;f__Enterobacteriaceae;g__Hafnia-Obesumbacterium",
       old.strings = c("d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacteriales;f__Enterobacteriaceae;g__Hafnia",
                       "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacteriales;f__Enterobacteriaceae;g__Hafnia-Obesumbacterium",
                       "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacteriales;f__Enterobacteriaceae;g__Obesumbacterium"))

strings.to.correct$Coprococcus <- 
  list(new.string = "d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__Coprococcus",
       old.strings = c("d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__Coprococcus",
                       "d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__Coprococcus 1",
                       "d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__Coprococcus 2",
                       "d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__Coprococcus 3"))

strings.to.correct$Corynebacterium <- 
  list(new.string = "d__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Corynebacteriales;f__Corynebacteriaceae;g__Corynebacterium",
       old.strings = c("d__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Corynebacteriales;f__Corynebacteriaceae;g__Corynebacterium",
                       "d__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Corynebacteriales;f__Corynebacteriaceae;g__Corynebacterium 1"))

strings.to.correct$Prevotella <- 
  list(new.string = "d__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Prevotellaceae;g__Prevotella",
       old.strings = c("d__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Prevotellaceae;g__Prevotella",
                       "d__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Prevotellaceae;g__Prevotella 1",
                       "d__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Prevotellaceae;g__Prevotella 2",
                       "d__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Prevotellaceae;g__Prevotella 6",
                       "d__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Prevotellaceae;g__Prevotella 7",
                       "d__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Prevotellaceae;g__Prevotella 9",
                       "d__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__[Paraprevotellaceae];g__[Prevotella]"))

strings.to.correct$Ruminococcus <- 
  list(new.string = "d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Ruminococcaceae;g__Ruminococcus",
       old.strings = c("d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Ruminococcaceae;g__Ruminococcus",
                       "d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Ruminococcaceae;g__Ruminococcus 1",
                       "d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Ruminococcaceae;g__Ruminococcus 2",
                       "d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__[Ruminococcus]",
                       "d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__[Ruminococcus] gauvreauii group",
                       "d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__[Ruminococcus] gnavus group",
                       "d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__[Ruminococcus] torques group"))

strings.to.correct$Ruminiclostridium <- 
  list(new.string = "d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Ruminococcaceae;g__Ruminiclostridium",
       old.strings = c("d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Ruminococcaceae;g__Ruminiclostridium",
                       "d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Ruminococcaceae;g__Ruminiclostridium 1",
                       "d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Ruminococcaceae;g__Ruminiclostridium 5",
                       "d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Ruminococcaceae;g__Ruminiclostridium 6",
                       "d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Ruminococcaceae;g__Ruminiclostridium 9"))

strings.to.correct$Eubacterium <- 
  list(new.string = "d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Eubacteriaceae;g__Eubacterium",
       old.strings = c("d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Eubacteriaceae;g__Eubacterium",
                       "d__Bacteria;p__Firmicutes;c__Erysipelotrichia;o__Erysipelotrichales;f__Erysipelotrichaceae;g__[Eubacterium]",
                       "d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Ruminococcaceae;g__[Eubacterium] coprostanoligenes group",
                       "d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__[Eubacterium] eligens group",
                       "d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__[Eubacterium] fissicatena group",
                       "d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__[Eubacterium] hallii group",
                       "d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__[Eubacterium] ruminantium group",
                       "d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__[Eubacterium] ventriosum group",
                       "d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__[Eubacterium] xylanophilum group",
                       "d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Family XIII;g__[Eubacterium] brachy group",
                       "d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Family XIII;g__[Eubacterium] nodatum group",
                       "d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__[Eubacterium] oxidoreducens group"))

strings.to.correct$Clostridium <- 
  list(new.string = "d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiaceae;g__Clostridium",
       old.strings = c("d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiaceae 1;g__Clostridium sensu stricto 1",
                       "d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiaceae 1;g__Clostridium sensu stricto 2",
                       "d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiaceae 1;g__Clostridium sensu stricto 3",
                       "d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiaceae 1;g__Clostridium sensu stricto 6",
                       "d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiaceae 1;g__Clostridium sensu stricto 13",
                       "d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiaceae 1;g__Clostridium sensu stricto 18",
                       "d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiaceae;g__Clostridium",
                       "d__Bacteria;p__Firmicutes;c__Erysipelotrichia;o__Erysipelotrichales;f__Erysipelotrichaceae;g__[Clostridium] innocuum group",
                       "d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Ruminococcaceae;g__CAG-352",
                       "d__Bacteria;p__Cyanobacteria;c__Melainabacteria;o__Gastranaerophilales;f__Clostridium sp. K4410.MGS-306;g__Clostridium sp. K4410.MGS-306"))

strings.to.correct$Tyzzerella <- 
  list(new.string = "d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__Tyzzerella",
       old.strings = c("d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__Tyzzerella",
                       "d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__Tyzzerella 3",
                       "d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__Tyzzerella 4"))

strings.to.correct$Azospirillum <- 
  list(new.string = "d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhodospirillales;f__Rhodospirillaceae;g__Azospirillum",
       old.strings = c("d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhodospirillales;f__Rhodospirillaceae;g__Azospirillum",
                       "d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhodospirillales;__;g__Azospirillum sp. 47_25",
                       "d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhodospirillales;f__uncultured;g__Azospirillum sp. 47_25"))



strings.to.correct$Lachnoclostridium <- 
  list(new.string = "d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__Lachnoclostridium",
       old.strings = c("d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__Lachnoclostridium",
                       "d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__Lachnoclostridium 10",
                       "d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__Lachnoclostridium 12",
                       "d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__Lachnoclostridium 5"))

strings.to.correct$Selenomonas <- 
  list(new.string = "d__Bacteria;p__Firmicutes;c__Negativicutes;o__Selenomonadales;f__Veillonellaceae;g__Selenomonas",
       old.strings = c("d__Bacteria;p__Firmicutes;c__Negativicutes;o__Selenomonadales;f__Veillonellaceae;g__Selenomonas",
                       "d__Bacteria;p__Firmicutes;c__Negativicutes;o__Selenomonadales;f__Veillonellaceae;g__Selenomonas 3",
                       "d__Bacteria;p__Firmicutes;c__Negativicutes;o__Selenomonadales;f__Veillonellaceae;g__Selenomonas 4"))

strings.to.correct$Treponema <- 
  list(new.string = "d__Bacteria;p__Spirochaetes;c__Spirochaetia;o__Spirochaetales;f__Spirochaetaceae;g__Treponema",
       old.strings = c("d__Bacteria;p__Spirochaetes;c__Spirochaetia;o__Spirochaetales;f__Spirochaetaceae;g__Treponema",
                       "d__Bacteria;p__Spirochaetes;c__Spirochaetia;o__Spirochaetales;f__Spirochaetaceae;g__Treponema 2"))

strings.to.correct$Faecalibacterium <- 
  list(new.string = "d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Ruminococcaceae;g__Faecalibacterium",
       old.strings = c("d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Ruminococcaceae;g__Faecalibacterium",
                       "d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Ruminococcaceae;g__UBA1819"))

strings.to.correct$Candidatus_Soleaferrea <- 
  list(new.string = "d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Ruminococcaceae;g__Candidatus_Soleaferrea",
       old.strings = c("d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Ruminococcaceae;g__Candidatus_Soleaferrea",
                       "d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Ruminococcaceae;g__Candidatus Soleaferrea"))

strings.to.correct$Candidatus_Stoquefichus <- 
  list(new.string = "d__Bacteria;p__Firmicutes;c__Erysipelotrichia;o__Erysipelotrichales;f__Erysipelotrichaceae;g__Candidatus_Stoquefichus",
       old.strings = c("d__Bacteria;p__Firmicutes;c__Erysipelotrichia;o__Erysipelotrichales;f__Erysipelotrichaceae;g__Candidatus_Stoquefichus",
                       "d__Bacteria;p__Firmicutes;c__Erysipelotrichia;o__Erysipelotrichales;f__Erysipelotrichaceae;g__Candidatus Stoquefichus"))

strings.to.correct$Barnesiellaceae__ <-
  list(new.string = "d__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Barnesiellaceae;__",
       old.strings = c("d__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__[Barnesiellaceae];__",
                       "d__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Barnesiellaceae;__"))

strings.to.correct$Clostridiaceae__ <-
  list(new.string = "d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiaceae;__",
       old.strings = c("d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiaceae 1;__",
                       "d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiaceae;__"))

strings.to.correct$Acidaminococcaceae__ <- 
  list(new.string = "d__Bacteria;p__Firmicutes;c__Negativicutes;o__Selenomonadales;f__Acidaminococcaceae;__",
       old.strings = c("d__Bacteria;p__Firmicutes;c__Negativicutes;o__Selenomonadales;f__Acidaminococcaceae;g__Acidaminococcaceae_unclassified",
                       "d__Bacteria;p__Firmicutes;c__Negativicutes;o__Selenomonadales;f__Acidaminococcaceae;__"))

strings.to.correct$Clostridiales_Family_XIII_Incertae_Sedis <-
  list(new.string = "d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiales_Family_XIII_Incertae_Sedis;__",
       old.strings = c("d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiales_Family_XIII_Incertae_Sedis;__",
                       "d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiales_Family_XIII._Incertae_Sedis;__",
                       "d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiales_Family_XIII_Incertae_Sedis;g__Clostridiales_Family_XIII_Incertae_Sedis_unclassified"))

strings.to.correct$f__Leptotrichiaceae__ <-
  list(new.string = "d__Bacteria;p__Fusobacteria;c__Fusobacteriia;o__Fusobacteriales;f__Leptotrichiaceae;__",
       old.strings = c("d__Bacteria;p__Fusobacteria;c__Fusobacteriia;o__Fusobacteriales;f__Leptotrichiaceae;__",
                       "d__Bacteria;p__Fusobacteria;c__Fusobacteriia;o__Fusobacteriales;f__Leptotrichiaceae;g__Leptotrichiaceae_unclassified"))

strings.to.correct$Muribaculaceae__ <-
  list(new.string = "d__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Muribaculaceae;__",
       old.strings = c("d__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Muribaculaceae;__",
                       "d__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Muribaculaceae;g__uncultured Bacteroidales bacterium"))

strings.to.correct$Gastranaerophilales____ <-
  list(new.string = "d__Bacteria;p__Cyanobacteria;c__Melainabacteria;o__Gastranaerophilales;__;__",
       old.strings = c("d__Bacteria;p__Cyanobacteria;c__Melainabacteria;o__Gastranaerophilales;__;__",
                       "d__Bacteria;p__Cyanobacteria;c__Melainabacteria;o__Gastranaerophilales;f__uncultured cyanobacterium;g__uncultured cyanobacterium"))

strings.to.correct$Mollicutes_RF39____ <-
  list(new.string = "d__Bacteria;p__Tenericutes;c__Mollicutes;o__Mollicutes RF39;__;__",
       old.strings = c("d__Bacteria;p__Tenericutes;c__Mollicutes;o__Mollicutes RF39;__;__",
                       "d__Bacteria;p__Tenericutes;c__Mollicutes;o__Mollicutes RF39;f__wallaby gut metagenome;g__wallaby gut metagenome",
                       "d__Bacteria;p__Tenericutes;c__Mollicutes;o__Mollicutes RF39;f__unidentified rumen bacterium RF39;g__unidentified rumen bacterium RF39"))

strings.to.correct$Rikenella <-
  list(new.string = "d__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Rikenellaceae;g__Rikenella",
       old.strings = c("d__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Rikenellaceae;g__Rikenella",
                       "d__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Rikenellaceae;g__uncultured Rikenella sp."))

strings.to.correct$Thermobaculum <- 
  list(new.string = "d__Bacteria;p__Chloroflexi;__;__;__;g__Thermobaculum",
      old.strings = c("d__Bacteria;__;__;__;__;g__Thermobaculum"))

strings.to.correct$Haloplasma <- 
  list(new.string = "d__Bacteria;p__Tenericutes;c__Mollicutes;o__Haloplasmatales;f__Haloplasmataceae;g__Haloplasma",
       old.strings = c("d__Bacteria;__;__;o__Haloplasmatales;f__Haloplasmataceae;g__Haloplasma"))

strings.to.correct$Pyrinomonas <-
  list(new.string = "d__Bacteria;p__Acidobacteria;c__Blastocatellia;o__Blastocatellales;f__Pyrinomonadaceae;g__Pyrinomonas",
       old.strings = c("d__Bacteria;p__Acidobacteria;c__Blastocatellia;__;__;g__Pyrinomonas"))

strings.to.correct$Tropheryma <-
  list(new.string = "d__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Actinomycetales;f__Microbacteriaceae;g__Tropheryma",
       old.strings = c("d__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Micrococcales;__;g__Tropheryma",
                       "d__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Actinomycetales;f__Microbacteriaceae;g__Tropheryma"))

strings.to.correct$Phocaeicola <-
  list(new.string = "d__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Phocaeicola",
       old.strings = c("d__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;__;g__Phocaeicola"))

strings.to.correct$Dehalogenimonas <-
  list(new.string = "d__Bacteria;p__Chloroflexi;c__Dehalococcoidia;o__Dehalococcoidales;f__Dehalococcoidaceae;g__Dehalogenimonas",
       old.strings = c("d__Bacteria;p__Chloroflexi;c__Dehalococcoidia;__;__;g__Dehalogenimonas"))

strings.to.correct$Microcystis <-
  list(new.string = "d__Bacteria;p__Cyanobacteria;c__Cyanophyceae;o__Chroococcales;f__Microcystaceae;g__Microcystis",
       old.strings = c("d__Bacteria;p__Cyanobacteria;__;o__Chroococcales;f__Microcystaceae;g__Microcystis"))

strings.to.correct$Chroococcidiopsis <-
  list(new.string = "d__Bacteria;p__Cyanobacteria;c__Cyanophyceae;o__Chroococcidiopsidales;f__Chroococcidiopsidaceae;g__Chroococcidiopsis",
       old.strings = c("d__Bacteria;p__Cyanobacteria;__;o__Chroococcidiopsidales;f__Chroococcidiopsidaceae;g__Chroococcidiopsis"))

strings.to.correct$Gottschalkia <-
  list(new.string = "d__Bacteria;p__Firmicutes;c__Tissierellia;o__Tissierellales;f__Gottschalkiaceae;g__Gottschalkia",
       old.strings = c("d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;__;g__Gottschalkia"))

strings.to.correct$Methyloceanibacter <-
  list(new.string = "d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Hyphomicrobiales;f__Hyphomicrobiaceae;g__Methyloceanibacter",
       old.strings = c("d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhizobiales;__;g__Methyloceanibacter"))

strings.to.correct$Ideonella <-
  list(new.string = "d__Bacteria;p__Proteobacteria;c__Betaproteobacteria;o__Burkholderiales;f__Comamonadaceae;g__Ideonella",
       old.strings = c("d__Bacteria;p__Proteobacteria;c__Betaproteobacteria;o__Burkholderiales;__;g__Ideonella"))

strings.to.correct$Paucibacter <-
  list(new.string = "d__Bacteria;p__Proteobacteria;c__Betaproteobacteria;o__Burkholderiales;f__Comamonadaceae;g__Paucibacter",
       old.strings = c("d__Bacteria;p__Proteobacteria;c__Betaproteobacteria;o__Burkholderiales;__;g__Paucibacter"))

strings.to.correct$Roseateles <-
  list(new.string = "d__Bacteria;p__Proteobacteria;c__Betaproteobacteria;o__Burkholderiales;f__Comamonadaceae;g__Roseateles",
       old.strings = c("d__Bacteria;p__Proteobacteria;c__Betaproteobacteria;o__Burkholderiales;__;g__Roseateles"))

strings.to.correct$Rhizobacter <-
  list(new.string = "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudomonadaceae;g__Rhizobacter",
       old.strings = c("d__Bacteria;p__Proteobacteria;c__Betaproteobacteria;o__Burkholderiales;__;g__Rhizobacter"))

strings.to.correct$Rubrivivax <-
  list(new.string = "d__Bacteria;p__Proteobacteria;c__Betaproteobacteria;o__Burkholderiales;f__Comamonadaceae;g__Rubrivivax",
       old.strings = c("d__Bacteria;p__Proteobacteria;c__Betaproteobacteria;o__Burkholderiales;__;g__Rubrivivax"))

strings.to.correct$Sphaerotilus <-
  list(new.string = "d__Bacteria;p__Proteobacteria;c__Betaproteobacteria;o__Burkholderiales;f__Comamonadaceae;g__Sphaerotilus",
       old.strings = c("d__Bacteria;p__Proteobacteria;c__Betaproteobacteria;o__Burkholderiales;__;g__Sphaerotilus"))

strings.to.correct$Tepidimonas <-
  list(new.string = "d__Bacteria;p__Proteobacteria;c__Betaproteobacteria;o__Burkholderiales;f__Comamonadaceae;g__Tepidimonas",
       old.strings = c("d__Bacteria;p__Proteobacteria;c__Betaproteobacteria;o__Burkholderiales;__;g__Tepidimonas"))

strings.to.correct$Thiomonas <-
  list(new.string = "d__Bacteria;p__Proteobacteria;c__Betaproteobacteria;o__Burkholderiales;f__Comamonadaceae;g__Thiomonas",
       old.strings = c("d__Bacteria;p__Proteobacteria;c__Betaproteobacteria;o__Burkholderiales;f__Burkholderiales_noname;g__Thiomonas",
                       "d__Bacteria;p__Proteobacteria;c__Betaproteobacteria;o__Burkholderiales;__;g__Thiomonas"))

strings.to.correct$Plesiomonas <-
  list(new.string = "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Plesiomonas",
       old.strings = c("d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;__;g__Plesiomonas"))

# For debugging
debug.unify.genera <- data.frame()

for (dataset in datasets) {
  tmp <- genera.new[[dataset]]
  
  # Iterate over genera requiring unification
  for (fix.genus in names(strings.to.correct)) {
    old.strings <- strings.to.correct[[fix.genus]]$old.strings
    new.string <- strings.to.correct[[fix.genus]]$new.string
    tmp <- tmp %>%
      mutate(Genus = ifelse(Genus %in% old.strings, 
                            new.string, Genus))
  }
  
  message(paste(dataset, "- renamed",
                sum(tmp$Genus != genera.new[[dataset]]$Genus),
                "out of", nrow(tmp), "genus entities"))
  
  debug.unify.genera <- 
    bind_rows(debug.unify.genera,
              data.frame(dataset = dataset,
                         old = genera.new[[dataset]]$Genus,
                         new = tmp$Genus) %>%
                filter(old != new))
  
  # Sum over new genus assignments
  tmp <- tmp %>%
    group_by(Genus) %>% 
    summarise(across(everything(), sum))
  
  genera.new[[dataset]] <- tmp
}
# Check out statistics: debug.unify.genera 
rm(old.strings, new.string, fix.genus)

# NO RUN
# require(stringdist)
# x <- c(); for(dataset in datasets) {x <- c(x, genera.new[[dataset]]$Genus)}; x <- unique(x)
# xx <- data.frame(stringsimmatrix(a = x, b = x))
# colnames(xx) <- x; xx$G1 <- x; 
# xx <- xx %>% tidyr::pivot_longer(!G1, names_to = "G2", values_to = "similarity") %>% filter(G1 > G2)
# xx <- xx %>% filter(similarity > 0.5)
# xx <- xx %>% filter(gsub(".*;g__","",G1) != gsub(".*;g__","",G2))
# write.table(xx, "tmp.tsv", sep="\t", row.names = F)


# --------------------------------
# 4. Save
# --------------------------------

# We copy the new unified tables to the "data/processed_data" folder in which final tables will be stored.

source("load_original_data/utils.R")
file.copy(data.dirs, "../data/processed_data", recursive = TRUE)

for(dataset in data.dirs) {
  genera <- genera.new[[basename(dataset)]]
  save.to.files(basename(dataset), "processed_data", genera = genera)
  save.to.rdata(basename(dataset), "processed_data", genera = genera)
}
