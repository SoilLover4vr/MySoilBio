# Script to calculate bacterial and fungal biomass, along with counts of protozoa and nematodes

# Load necessary libraries
library(dplyr)

# -------------------------------------------------------------------
# Instructions:
# 1. Place your data files in a folder named "data" within your project directory.
# 2. File requirements:
#    - metadata.csv: Metadata file with columns for dilution factors, bacterial counts, nematode counts, etc.
#    - fungal_data.csv: Fungal fragment data with columns for lengths and widths.
# 3. Output files will be saved in a folder named "output" in your project directory.
# -------------------------------------------------------------------

# Load the metadata and fungal data
metadata <- read.csv("data/metadata.csv", stringsAsFactors = FALSE)
fungal_data <- read.csv("data/fungal_data.csv", stringsAsFactors = FALSE)

# Convert 'ID' to character and 'Date' columns to Date format for consistency
metadata$ID <- as.character(metadata$ID)
metadata$Date <- as.Date(metadata$Date, format="%Y-%m-%d")
fungal_data$ID <- as.character(fungal_data$ID)
fungal_data$Fungal_Date <- as.Date(fungal_data$Fungal_Date, format="%Y-%m-%d")

# Convert dilution factors to numeric by removing "1:"
metadata <- metadata %>%
  mutate(
    Main.Dilution = as.numeric(gsub("1:", "", Main.Dilution)),
    Bacterial.Dilution = as.numeric(gsub("1:", "", Bacterial.Dilution)),
    Drops.per.mL = as.numeric(Drops.per.mL)
  )

# Translate 'Bacterial FOV' to fractional values
metadata <- metadata %>%
  mutate(
    FOV_Fraction = case_when(
      Bacterial.FOV == "Quarter" ~ 0.25,
      Bacterial.FOV == "Half" ~ 0.5,
      Bacterial.FOV == "Whole" ~ 1,
      TRUE ~ NA_real_
    )
  )

# Merge fungal_data and metadata on 'ID' and 'Date'
data <- merge(fungal_data, metadata, by.x = c("ID", "Fungal_Date"), by.y = c("ID", "Date"), all.x = TRUE)

# Constants for scaling
fov_diameter_um <- 450  # Field of view diameter in µm
fov_per_drop <- 2038    # Approximate number of FOVs per drop
density_constant_fungal <- 0.41  # Fungal density in pg/µm³
density_constant_bacterial <- 0.33  # Bacterial density in pg/µm³

# Define fungal biomass calculation function
calculate_fungal_biomass <- function(length_proportions, diameters, dilution, drops_per_ml) {
  biomass_per_fragment <- mapply(function(length_prop, diameter) {
    length_um <- length_prop * fov_diameter_um
    radius <- diameter / 2
    volume <- pi * radius^2 * length_um
    biomass <- volume * density_constant_fungal
    return(biomass)
  }, length_proportions, diameters)
  
  total_biomass_pg <- sum(biomass_per_fragment, na.rm = TRUE) / 25 * fov_per_drop * drops_per_ml * dilution
  total_biomass_ug <- total_biomass_pg / 1e6
  return(total_biomass_ug)
}

# Define bacterial biomass calculation function
calculate_bacterial_biomass <- function(counts, dilution, fov_fraction, drops_per_ml) {
  avg_count <- mean(counts, na.rm = TRUE)
  full_fov_count <- avg_count * (1 / fov_fraction)
  total_count_per_drop <- full_fov_count * fov_per_drop
  biomass_pg <- total_count_per_drop * density_constant_bacterial
  biomass_ug_g <- (biomass_pg * dilution * drops_per_ml) / 1e6
  return(biomass_ug_g)
}

# Apply calculations across each unique ID-Date
data <- data %>%
  group_by(ID, Fungal_Date) %>%
  summarize(
    BacBio = calculate_bacterial_biomass(c(Bac1, Bac2, Bac3, Bac4, Bac5), unique(Bacterial.Dilution), unique(FOV_Fraction), unique(Drops.per.mL)),
    FunBio = calculate_fungal_biomass(FunL, FunW, unique(Main.Dilution), unique(Drops.per.mL)),
    `F:B` = FunBio / BacBio,
    .groups = "drop"
  )

# Define output file name
output_file <- "output/Final_Summary_Results.csv"

# Save the results to CSV
write.csv(data, output_file, row.names = FALSE)

# Print results
print("Final Summary Results:")
print(data)
