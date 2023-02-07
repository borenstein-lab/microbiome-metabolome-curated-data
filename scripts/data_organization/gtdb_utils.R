get.gtdb.mapper <- function() {
  require(stringr)
  require(readr)
  require(dplyr)
  
  add.prefix.to.silva <- function(x) {
    if (is.na(x)) return(NA)
    if (str_count(x, ";") != 6) return(NA)
    return(gsub("(.*);(.*);(.*);(.*);(.*);(.*);(.*)",
                "d__\\1;p__\\2;c__\\3;o__\\4;f__\\5;g__\\6;s__\\7", x))
  }
  
  bac_metadata <- read_delim("../references/gtdb/bac120_metadata_r207.tsv", 
                                     delim = "\t", 
                                     escape_double = FALSE, 
                                     trim_ws = TRUE, 
                                     na = c(""," ","NA","none"),
                                     show_col_types = FALSE)
  
  ar_metadata <- read_delim("../references/gtdb/ar53_metadata_r207.tsv", 
                                     delim = "\t", 
                                     escape_double = FALSE, 
                                     trim_ws = TRUE, 
                                     na = c(""," ","NA","none"),
                                    show_col_types = FALSE)
  
  bac_metadata <- bind_rows(bac_metadata, ar_metadata)
  rm(ar_metadata)
  
  # Add taxonomy prefixes to silva taxonomy
  bac_metadata$ssu_silva_taxonomy <- sapply(bac_metadata$ssu_silva_taxonomy, 
                                                    add.prefix.to.silva, 
                                                    USE.NAMES = FALSE)
  
  # Organize table as a map from gtdb to the other databases
  bac_metadata <- bac_metadata %>%
    select(gtdb_taxonomy, 
           ssu_silva_taxonomy,
           ssu_gg_taxonomy,
           ncbi_taxonomy,
           ncbi_taxid) %>%
    mutate(ncbi_taxid = as.character(ncbi_taxid)) %>%
    distinct() %>%
    tidyr::pivot_longer(cols = c("ssu_silva_taxonomy","ssu_gg_taxonomy","ncbi_taxonomy","ncbi_taxid"),
                        names_to = "ref_db", values_to = "ref_taxonomy",
                        values_drop_na = TRUE) %>%
    mutate(gtdb_genus = gsub(";s__.*$","",gtdb_taxonomy)) %>%
    mutate(ref_genus = gsub(";s__.*$","",ref_taxonomy)) 
  
  return(bac_metadata)
}

map.gtdb.short.to.long <- function(names.to.map, level = "species") {
  gtdb.map <- get.gtdb.mapper()
  
  # Create mapping vector from short to long names
  if (level == "species") {
    gtdb.map <- gtdb.map %>% 
      select(gtdb_taxonomy) %>% 
      rename(long = gtdb_taxonomy) %>%
      distinct() %>%
      mutate(short = gsub(".*;s__", "", long))
  } else if (level == "genera") {
    gtdb.map <- gtdb.map %>% 
      select(gtdb_genus) %>% 
      rename(long = gtdb_genus) %>%
      distinct() %>%
      mutate(short = gsub(".*;g__", "", long))
  } else return(NULL)
  
  unmappable <- sum(! names.to.map %in% gtdb.map$short)
  if (unmappable > 0) {
    print("Some names cannot be mapped")
  }
    
  map.vec <- gtdb.map$long
  names(map.vec) <- gtdb.map$short
  
  new.names <- unname(map.vec[names.to.map])
  return(new.names)
}

# --------------------------------------------------------
# Get the lowest level classified
# 1 = domain, 2 = phylum, 3 = class, 4 = order, 5 = family, 6 = genus
# --------------------------------------------------------
get.classified.level <- function(genus.str) {
  level <- 6
  if (grepl("g__$", genus.str)) level <- 5
  if (grepl("f__;g__$", genus.str)) level <- 4
  if (grepl("o__;f__;g__$", genus.str)) level <- 3
  if (grepl("c__;o__;f__;g__$", genus.str)) level <- 2
  if (grepl("p__;c__;o__;f__;g__$", genus.str)) level <- 1
  return(level)
}



get.last.level <- function(genus.str, level) {
  if (level == 6) tmp <- gsub(".*;g__", "g__", genus.str)
  if (level == 5) tmp <- gsub(".*;f__", "f__", gsub(";g__.*", "", genus.str))
  if (level == 4) tmp <- gsub(".*;o__", "o__", gsub(";f__.*", "", genus.str))
  if (level == 3) tmp <- gsub(".*;c__", "c__", gsub(";o__.*", "", genus.str))
  if (level == 2) tmp <- gsub(".*;p__", "p__", gsub(";c__.*", "", genus.str))
  if (level == 1) tmp <- gsub(";p__.*", "", genus.str)
  tmp <- gsub("\\[", "\\\\\\[", tmp)
  tmp <- gsub("\\]", "\\\\\\]", tmp)
  return(tmp)
}

make.unclassified.util <- function(level) {
  if (level == 6) return(c("", ""))
  if (level == 5) return(c("g__.*", "g__"))
  if (level == 4) return(c("f__.*;g__.*", "f__;g__"))
  if (level == 3) return(c("o__.*;f__.*;g__.*", "o__;f__;g__"))
  if (level == 2) return(c("c__.*;o__.*;f__.*;g__.*", "c__;o__;f__;g__"))
  if (level == 1) return(c("p__.*;c__.*;o__.*;f__.*;g__.*", "p__.;c__;o__;f__;g__"))
}

# --------------------------------------------------------
# Returns the (closest) GTDB taxonomy corresponding to a given 
#  genus name (based on another DB  such as GG/Silva/NCBI).
# For transparency regarding non-trivial mappings, a "comment" 
#  string is returned as well, i.e. both the mapping itself 
#  and a comment are returned as a string vector.
# Only to be used when it is not possible to directly use
#  GTDB when assigning taxonomy to raw data!
# --------------------------------------------------------
map.other.ref.to.gtdb.genus <- function(genus.str, gtdb.map, use.ncbi = T, use.ncbi.id = F, use.gg = F, use.silva = T) {
  level <- get.classified.level(genus.str)
  
  # Filter mapping table to relevant DBs only
  gtdb.map.tmp <- gtdb.map
  if (! use.ncbi) gtdb.map.tmp <- gtdb.map.tmp %>% filter(ref_db != "ncbi_taxonomy")
  if (! use.ncbi.id) gtdb.map.tmp <- gtdb.map.tmp %>% filter(ref_db != "ncbi_taxid")
  if (! use.gg) gtdb.map.tmp <- gtdb.map.tmp %>% filter(ref_db != "ssu_gg_taxonomy")
  if (! use.silva) gtdb.map.tmp <- gtdb.map.tmp %>% filter(ref_db != "ssu_silva_taxonomy")
  
  # Check if taxon already matches GTDB taxonomy as is
  if (gsub("^k__", "d__", genus.str) %in% gtdb.map.tmp$gtdb_genus)
    return(c(genus = gsub("^k__", "d__", genus.str), comment = ""))
  
  # Check if taxon matches reference taxonomies
  tmp <- gtdb.map.tmp %>% 
    filter(ref_genus == genus.str) %>% 
    group_by(gtdb_genus) %>% 
    summarise(N=n()) %>%
    arrange(-N)
  
  # Any single "dominating" match?
  #  (We add a small buffer (1/2) so that we don't get supposedly confident mappings if some ambiguity exists)
  if (nrow(tmp) > 1) tmp <- tmp %>% filter(N > max(N)*0.3) 
  
  # No match - try a "looser" search
  taxon.short <- get.last.level(genus.str, level)
  if (nrow(tmp) == 0) {
    tmp <- gtdb.map.tmp %>% 
      filter(grepl(taxon.short, ref_genus)) %>% 
      group_by(gtdb_genus) %>% 
      summarise(N=n()) %>%
      arrange(-N)
    if (nrow(tmp) > 1) tmp <- tmp %>% filter(N > max(N)*0.3)
  }
  
  # Single match
  if (nrow(tmp) == 1) 
    return(c(genus = tmp$gtdb_genus, comment = ""))
  
  # Multiple matches
  if (nrow(tmp) > 1 & level == 6) {
    message(paste0("The mapping for \'", genus.str, "\' should be verified manually"))
    return(c(genus = tmp$gtdb_genus[1], 
             comment = "Multiple matches found. Validate mapping manually."))
  }
  
  # Here we deal with the option of features unclassified at the genus or higher level.  
  #  We return a string based on GTDB taxonomy but classified only up to a reliable level, 
  #  where there's no ambiguity in mappings from reference DB to GTDB.
  if (nrow(tmp) > 1 & level < 6) {
    # Leave unclassified at relevant levels
    tmps <- make.unclassified.util(level)
    tmp <- tmp %>%
      mutate(gtdb_genus = gsub(tmps[1],tmps[2],gtdb_genus)) %>%
      group_by(gtdb_genus) %>% 
      summarise(N=sum(N)) %>%
      arrange(-N)
    
    message(paste0("The mapping for \'", genus.str, "\' should be verified manually"))
    return(c(genus = tmp$gtdb_genus[1], 
             comment = "Multiple matches found. Validate mapping manually."))
  }
  
  message(paste0("No match at all found for: \'", genus.str, "\'"))
  return(c(genus = gsub("^k__", "d__", genus.str), 
           comment = "No matches found. Keeping original name."))
}
