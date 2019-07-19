# use gslist and gs.names from crossmeta

gslist <- crossmeta::gslist
gs.names <- crossmeta::gs.names

usethis::use_data(gslist, gs.names, internal = TRUE, overwrite = TRUE)
