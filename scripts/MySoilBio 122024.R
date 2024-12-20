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
    Bacterial.Dilution = as.numeric(gsub("1:", "", Bacterial.Dilution))
  )

# Merge fungal_data and metadata on 'ID' and 'Date'
data <- merge(fungal_data, metadata, by.x = c("ID", "Fungal_Date"), by.y = c("ID", "Date"), all.x = TRUE)

# Constants for scaling
field_number <- 18            # Eyepiece field number
objective_magnification <- 40  # Objective lens magnification
fov_diameter_mm <- field_number / objective_magnification  # FoV diameter in mm
fov_diameter_um <- fov_diameter_mm * 1000  # Convert to µm

# Define fungal biomass calculation function
# This function calculates fungal biomass based on length and diameter of fragments.
calculate_fungal_biomass <- function(length_proportions, diameters) {
  biomass_per_fragment <- mapply(function(length_prop, diameter) {
    length_um <- length_prop * fov_diameter_um  # Convert length proportion to absolute length in µm
    radius <- diameter / 2
    volume <- pi * radius^2 * length_um  # Volume in µm³
    biomass <- volume * 0.41  # Convert volume to biomass using fungal density constant (pg/µm³)
    return(biomass)
  }, length_proportions, diameters)
  
  total_biomass <- sum(biomass_per_fragment, na.rm = TRUE)  # Total biomass per sample
  return(total_biomass)  # Returns biomass in picograms (pg)
}

# Apply fungal biomass calculation across each unique ID and date
data <- data %>%
  group_by(ID, Fungal_Date) %>%
  mutate(FunBio = calculate_fungal_biomass(FunL, FunW)) %>%
  ungroup()

# Define bacterial biomass calculation function
# This function calculates bacterial biomass based on bacterial counts and dilution factors.
calculate_bacterial_biomass <- function(counts, dilution) {
  avg_count <- mean(counts, na.rm = TRUE)  # Average bacterial count per field of view
  biomass <- avg_count * dilution * 0.33 / 1e6  # Convert counts to µg/g using bacterial density constant
  biomass_scaled <- biomass * 19  # Apply scaling for 19 drops per mL
  return(biomass_scaled)  # Returns biomass in µg/g
}

# Apply bacterial biomass calculation across each unique ID and date
data <- data %>%
  group_by(ID, Fungal_Date) %>%
  mutate(
    BacBio = mapply(calculate_bacterial_biomass, list(c(Bac1, Bac2, Bac3, Bac4, Bac5)), Bacterial.Dilution)
  ) %>%
  ungroup()

# Protozoa calculations (flagellates and amoebae)
# Calculates total protozoa counts (scaled by dilution and drop factor).
data <- data %>%
  mutate(
    Flag = Flagellates * Main.Dilution * 19,
    Amoe = Amoeba * Main.Dilution * 19,
    Proto = Flag + Amoe  # Total protozoa count
  )

# Nematode calculations
# Calculates scaled counts for different nematode types.
data <- data %>%
  mutate(
    BfNem = Bf.Nem * Main.Dilution * 19,
    FfNem = Ff.Nem * Main.Dilution * 19,
    PNem = Pred.Nem * Main.Dilution * 19,
    RfNem = Rf.Nem * Main.Dilution * 19
  )

# Summarize results by ID and date
# Calculates fungal-to-bacterial biomass ratio (F:B) and other metrics.
summary_results <- data %>%
  group_by(ID, Fungal_Date) %>%
  summarize(
    BacBio = mean(BacBio, na.rm = TRUE),  # Average bacterial biomass
    FunBio = sum(FunBio, na.rm = TRUE),  # Total fungal biomass
    `F:B` = FunBio / BacBio,  # Fungal-to-bacterial biomass ratio
    Proto = sum(Proto, na.rm = TRUE),  # Total protozoa count
    Flag = sum(Flag, na.rm = TRUE),
    Amoe = sum(Amoe, na.rm = TRUE),
    BfNem = sum(BfNem, na.rm = TRUE),  # Total bacterial-feeding nematodes
    FfNem = sum(FfNem, na.rm = TRUE),  # Total fungal-feeding nematodes
    PNem = sum(PNem, na.rm = TRUE),  # Total predatory nematodes
    RfNem = sum(RfNem, na.rm = TRUE),  # Total root-feeding nematodes
    .groups = "drop"
  )

# Save results to a CSV file
# Outputs results to the "output" directory.
write.csv(summary_results, "output/Final_Summary_Results.csv", row.names = FALSE)

# Display the summary results to check
print("Final Summary Results:")
print(summary_results)
