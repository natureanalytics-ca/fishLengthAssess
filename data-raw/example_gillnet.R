#load the csv
raw_data_gillnet <- read.csv("gillnet_data_example.csv", header = FALSE)

# Save as .rda
usethis::use_data(raw_data_gillnet, overwrite = TRUE)
