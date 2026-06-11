# load required libraries
require(data.table)
require(stringr)

# load gbif data and filter by recent years
gbif_db <- fread("gbifData27102025.csv")
gbif_recent <- gbif_db[gbif_db$year > 2014, ]
gbif_clean <- gbif_recent[, c("species", "decimalLatitude", "decimalLongitude")]

# load secondary dataset
gracia_db <- fread("gracia.csv")
gracia_db$y <- as.numeric(str_replace(gracia_db$Latdd, ",", "."))
gracia_db$x <- as.numeric(str_replace(gracia_db$Londd, ",", "."))

# filter secondary dataset by year 2015-2024
gracia_recent <- gracia_db[gracia_db$YEAR >= 2015, ]
gracia_filtered <- gracia_recent[gracia_recent$YEAR < 2025, ]

gracia_clean <- gracia_filtered[, c("Taxon name", "y", "x")]
names(gracia_clean) <- c("species", "decimalLatitude", "decimalLongitude")

# merge datasets
final_bat_db <- rbind(gbif_clean, gracia_clean)

# save final compiled dataset
write.csv(final_bat_db, "bat_records_europe.csv", row.names = FALSE)