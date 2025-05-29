#load the csv
raw_data_gillnet <- read.csv("gillnet_data_example.csv", header = FALSE)
colnames(raw_data_gillnet) <- c("MidLength", paste0("Mesh", 1:8))
# Save as .rda
usethis::use_data(raw_data_gillnet, overwrite = TRUE)
