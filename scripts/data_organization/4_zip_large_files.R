# -----------------------------------------------------------------
# Zip data files that are to large for GitHub
# -----------------------------------------------------------------

all.data.files <- list.files("../data/processed_data", recursive = T, full.names = T)
file.sizes <- file.size(all.data.files)/1000

# Which files require zipping?
files.to.zip <- all.data.files[file.sizes > 99000]
message(paste("Zipping",length(files.to.zip),"files"))

# Zip
for (f in files.to.zip) {
  zip(zipfile = paste0(f, ".zip"), 
      files = f)
}
