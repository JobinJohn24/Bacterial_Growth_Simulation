# Data Folder

## Overview
place your experimental bacterial growth data here for model fitting & comparison.

## Expected Format
CSV files with following columns: 
```csv
time, desnity, concentration
```
## Column Definitions
- ** time **: Time points in hours(e.g., 0-48).
- ** desnity **: Bacterial count in CFU/mL.
- ** concentration **: Antibiotic concentration in µg/mL.

## Data Collection Guidelines
- Measure at regulat time intervals (e.g., every 0.5-1 hour).
- Use consistent bacterial species and medium.
- Include control (no antibiotic) samples
- Include multiple antibiotic concentrations.

# Using the data

1. ** Prepare the data: ** - ensure CSV with required columns.
2. ** Rename CSV: ** `growth.csv`
3. ** Place CSV in folder: ** `data/growth.csv`
4. ** Run simulation: ** `Rscript src/main./R`

