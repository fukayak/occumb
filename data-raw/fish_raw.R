## code to prepare `fish_raw` dataset goes here
setwd("<occumb directory>")
devtools::load_all()

load("./data-raw/data_fish.Rdata")
riverbank_new <- riverbank
riverbank_new[riverbank_new == 0] <- "with_vegetation"
riverbank_new[riverbank_new == 1] <- "without_vegetation"
riverbank <- as.factor(riverbank_new)

fish_raw <- list(y = y, mismatch = mismatch, riverbank = riverbank)

usethis::use_data(fish_raw, overwrite = TRUE)
