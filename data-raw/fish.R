## code to prepare `fish` dataset goes here
setwd("<occumb directory>")
devtools::load_all()

load("./data-raw/data_fish.Rdata")
riverbank_new <- riverbank
riverbank_new[riverbank_new == 0] <- "with_vegetation"
riverbank_new[riverbank_new == 1] <- "without_vegetation"
riverbank <- as.factor(riverbank_new)

fish <- occumbData(
    y = y,
    spec_cov = list(mismatch = mismatch),
    site_cov = list(riverbank = riverbank),
    repl_cov = NULL
)

usethis::use_data(fish, overwrite = TRUE)
