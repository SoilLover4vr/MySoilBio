# MySoilBio: Open-Source Soil Biomass Calculations in R

MySoilBio is an open-source repository providing R scripts for calculating microbial biomass (including bacterial and fungal biomass), protozoa counts, and nematode abundance in soil samples using shadowing microscopy and dilution scaling. This tool is designed to support soil health research by offering an accessible, reproducible alternative to proprietary tools.

---

## Features
- Calculate **bacterial biomass** (µg/g soil).
- Calculate **fungal biomass** (µg/g soil).
- Estimate **protozoa counts** (flagellates and amoebae).
- Enumerate **nematodes** (bacterial-feeding, fungal-feeding, predatory, and root-feeding).
- Includes example datasets for reproducibility and testing.

---

## Repository Structure
```
MySoilBio/
  README.md           # Project overview and instructions (this file)
  LICENSE             # Licensing information
  scripts/            # R scripts for analysis
    calculate_biomass.R
    calculate_protozoa_nematodes.R
    load_metadata.R
    load_fungal_data.R
    run_analysis.R
  data/               # Example datasets
    metadata.csv
    fungal_data.csv
  output/             # Placeholder for output results (empty by default)
```

---

## Getting Started

### Prerequisites
- **R (Version 4.0 or later)**
- R packages:
  - `dplyr`

### Installation
1. Clone this repository:
   ```bash
   git clone https://github.com/SoilLover4vr/MySoilBio.git
   cd MySoilBio
   ```
2. Place your input data files (`metadata.csv` and `fungal_data.csv`) in the `data/` folder.
3. Open the R scripts in your preferred IDE (e.g., RStudio).

---

## Usage

### Example Workflow
1. **Prepare the Input Data Files**:
   - `metadata.csv`: Contains sample-level metadata (dilution factors, bacterial counts, nematode counts, etc.).
   - `fungal_data.csv`: Contains fungal fragment measurements (length and width).

2. **Source the R Scripts**:
   - Run the following commands in R to load the required functions:
     ```R
     source("scripts/load_metadata.R")
     source("scripts/load_fungal_data.R")
     source("scripts/calculate_biomass.R")
     source("scripts/calculate_protozoa_nematodes.R")
     source("scripts/run_analysis.R")
     ```

3. **Execute the Analysis**:
   - Run the `run_analysis()` function to process the input data and generate the results:
     ```R
     results <- run_analysis("data/metadata.csv", "data/fungal_data.csv")
     ```

4. **Export Results**:
   - Save the results to a CSV file:
     ```R
     write.csv(results, "output/Final_Summary_Results.csv", row.names = FALSE)
     ```

5. **Review Results**:
   - The output file will include metrics such as:
     - Bacterial biomass (BacBio)
     - Fungal biomass (FunBio)
     - Fungal-to-bacterial biomass ratio (F:B)
     - Protozoa counts (Proto)
     - Nematode counts (BfNem, FfNem, RfNem, PNem)

---

## Example Data
Example datasets (`metadata.csv` and `fungal_data.csv`) are provided in the `data/` folder for testing. These files contain:
- **metadata.csv**:
  - Sample IDs, dilution factors, bacterial counts, and nematode counts.
- **fungal_data.csv**:
  - Sample IDs, fungal fragment lengths, and widths.

---

## Contributions
Contributions to this project are welcome! To contribute:
1. Fork the repository.
2. Create a new branch for your feature or bug fix.
3. Submit a pull request with your changes.

---

## License
This project is licensed under the MIT License. See the `LICENSE` file for details.

---

## Contact
For questions or suggestions, please open an issue in this repository or reach out to the maintainers.