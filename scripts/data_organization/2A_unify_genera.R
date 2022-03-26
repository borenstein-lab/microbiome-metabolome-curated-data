# ----------------------------------------------------------------
# Unify genus names in genera-feature-tables.
# Override relevant files (tsv + RData) with the unified version.
# ----------------------------------------------------------------

require(vegan)
require(dplyr)
require(readr)
source("data_organization/utils.R")

# Load data
all.data <- load.all.datasets("prelim_data")
for(i in 1:length(all.data)) assign(names(all.data)[i], all.data[[i]])
rm(all.data)
datasets <- basename(data.dirs)

# --------------------------------
# 1. Remove non-bacteria
# --------------------------------

# We remove non-bacteria entities 
#  (before re-normalizing each sample to 100%)

for (dataset in datasets) {
  tmp <- genera[[dataset]]
  tmp.n.row <- nrow(tmp) # Recording number before we erase all non-bacteria
  tmp <- tmp[!grepl("Viruses", tmp$Genus, ignore.case = T),]
  tmp <- tmp[!grepl("Archaea", tmp$Genus, ignore.case = T),]
  tmp <- tmp[!grepl("Eukaryota",tmp$Genus, ignore.case = T),]
  tmp <- tmp[!grepl("d__;",tmp$Genus, ignore.case = T),]
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
# 3. Relabel unclassified genera 
# --------------------------------

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

unclass.tax.strings <- c("d__Bacteria;p__;c__;o__;f__;g__",
                         "d__Bacteria;p__.;c__;o__;f__;g__")

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

# Check out statistics: View(debug.unclass.genera.stats)
rm(unclassified.rows, tmp)

# Code for a quick indication of whether genus names are consistent:
# x <- c(); for(dataset in datasets) {x <- c(x, genera.new[[dataset]]$Genus)}; 
# x <- data.frame(genus = x) %>% group_by(genus) %>% summarise(N = n()) %>% tidyr::separate(col = "genus", into = c("d","p","c","o","f","g"), sep = ";", remove = FALSE)
# Sanity #1 (should return empty): View(x %>% filter(g != "g__") %>% group_by(g) %>% filter(n()>1))
# Sanity #2 (should return empty): View(x %>% filter(g == "g__" & f != "f__") %>% group_by(f) %>% filter(n()>1))

# --------------------------------
# 4. Save
# --------------------------------

# We copy the new unified tables to the "data/processed_data" 
#  folder in which final tables will be stored.

file.copy(data.dirs, 
          "../data/processed_data", 
          recursive = TRUE, 
          overwrite = TRUE)

for(dataset in data.dirs) {
  message(paste("Saving:", basename(dataset)))
  genera <- genera.new[[basename(dataset)]]
  save.to.files(basename(dataset), 
                "processed_data", 
                genera = genera)
  save.to.rdata(basename(dataset), 
                "processed_data", 
                genera = genera)
}
