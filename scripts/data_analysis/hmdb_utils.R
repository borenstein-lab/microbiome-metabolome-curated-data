# -------------- Function for parsing xml HMDB object --------------
# The "store" slot within any hmdb object is not very convenient,
#  it's an xml converted to a list. Here I turn it into a readable 
#  string, which can be later explored with simple regular expressions.
hmdb.store.to.string <- function(root) {
  if (!("term" %in% names(root))) {return()}
  term <- root$term
  if (!("descendants" %in% names(root))) {return(paste0("[",term,"]"))}
  children <- root$descendants
  children.string <- ""
  for (k in 1:length(children)) {
    new.root <- children[[k]]
    new.string <- hmdb.store.to.string(new.root)
    if (k > 1) {children.string <- paste0(children.string, ", ")}
    children.string <- paste0(children.string, new.string)
  }
  res.string <- paste0("[",term,": ",children.string,"]")
  return(res.string)
}


get.hmdb.data.by.id <- function(hmdb.id) {
  require(httr)
  require(XML)
  hmdb.xml <- xmlParse(rawToChar(GET(paste0("https://www.hmdb.ca/metabolites/", 
                                            hmdb.id, ".xml"))$content))
  hmdb.xml.l <- xmlToList(hmdb.xml, simplify = F)
  
  # Get metabolite name in HMDB
  hmdb.name <- hmdb.xml.l$name
  
  # Get metabolite class (if available)
  if ("taxonomy" %in% names(hmdb.xml.l) &
      is.list(hmdb.xml.l$taxonomy)) {
    tmp <- hmdb.xml.l$taxonomy$sub_class
    hmdb.sub.class <- ifelse(!is.null(tmp), tmp, NA)
    tmp <- hmdb.xml.l$taxonomy$class
    hmdb.class <- ifelse(!is.null(tmp), tmp, NA)
    tmp <- hmdb.xml.l$taxonomy$super_class
    hmdb.super.class <- ifelse(!is.null(tmp), tmp, NA)
    tmp <- hmdb.xml.l$taxonomy$kingdom
    hmdb.kingdom <- ifelse(!is.null(tmp), tmp, NA)
  } else {
    hmdb.sub.class <- NA
    hmdb.class <- NA
    hmdb.super.class <- NA
    hmdb.kingdom <- NA
  }
  
  # Record sources if available
  hmdb.source.endogenous <- NA
  hmdb.source.food <- NA
  hmdb.source.microbes <- NA
  
  if ("ontology" %in% names(hmdb.xml.l)) {
    for (j in 1:length(hmdb.xml.l$ontology)) {
      ont.string <- hmdb.store.to.string(hmdb.xml.l$ontology[[j]])
      if (is.null(ont.string) | !grepl("^\\[Disposition:", ont.string)) { next }
      else { break }
    }
    
    hmdb.source.endogenous <- grepl("\\[Source.*\\[Endogenous", ont.string)
    hmdb.source.food <- grepl("\\[Source.*\\[Food", ont.string)
    hmdb.source.microbes <- grepl("\\[Source.*\\[Microbe", ont.string)
    # hmdb.source.microbes.List <- ifelse(hmdb.source.microbes, gsub("^.*\\[Microbe", "[Microbe", ont.string), NA)
    # hmdb.source.microbes.List <- gsub("(, \\[Plant).*|(\\]\\], \\[Biological location).*|(, \\[Fungi).*|(, \\[Animal).*|(]]]$)", 
    #                                               "", hmdb.source.microbes.List)
  }
  
  return(c("HMDB" = hmdb.id,
           "HMDB.Name" = hmdb.name,
           "HMDB.Sub.Class" = hmdb.sub.class,
           "HMDB.Class" = hmdb.class,
           "HMDB.Super.Class" = hmdb.super.class,
           "HMDB.Kingdom" = hmdb.kingdom,
           "HMDB.Source.Endogenous" = hmdb.source.endogenous,
           "HMDB.Source.Food" = hmdb.source.food,
           "HMDB.Source.Microbial" = hmdb.source.microbes))
}

get.hmdb.data.by.ids <- function(hmdb.ids) {
  require(dplyr, quietly = TRUE)
  return(data.frame(t(sapply(hmdb.ids, get.hmdb.data.by.id))) %>% tibble::remove_rownames())
}
