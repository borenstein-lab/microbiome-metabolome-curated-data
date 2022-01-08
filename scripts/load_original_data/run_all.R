load.data.scripts <- list.files("load_original_data")
load.data.scripts <- grep("load_data_.*\\.R", load.data.scripts, value = TRUE)

sapply(load.data.scripts,
       FUN = function(x) {
         source(x)
         message(paste("Completed script:",x))
       })

message("Completed data loading")