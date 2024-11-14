# Bumblebee Latitudinal Distribution and Climate Change Analysis

This project investigates the impact of climate change on bumblebee latitudinal distribution, focusing on the impact of global temperature variations. By integrating biodiversity data from the Barcode of Life Data System (BOLD) with an external temperature dataset, potential shifts in bumblebee populations and their adaptations to climate change was analyzed. 

## Table of Contents
- [Introduction](#introduction)
- [Project Structure](#project-structure)
- [Data Sources](#data-sources)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Analysis Workflow](#analysis-workflow)
- [Results](#results)
- [References](#references)
- [Attribution](#attribution)

## Introduction
Understanding how climate change affects bumblebee populations is critical for biodiversity conservation. By studying latitudinal distribution and temperature data, this project examines how temperature shifts impact the distribution and genetic diversity of bumblebee populations across Canada and USA.

## Project Structure
- `Data/`: Raw and processed data, including BOLD and external temperature data.
- `R Script/`: Script file for data processing, analysis, and visualization.
- `Plots/`:  Significant figures, and visualizations.
- Assignment 1 PDF: Written evaluationa and analysis. 
- `README.md`: Project overview and instructions.

## Data Sources
1. **BOLD Database**: Provides genetic data on bumblebee populations across various locations.
2. **External Temperature Data**: Temperature data from climate sources such as NOAA or the WorldClim database, detailing monthly and annual averages across latitudes.

## Requirements
- R (version 4.0 or higher)
- R packages: `tidyverse`, `dplyr`, `vegan`, `dplyr`, `plotly`, `sf`, `maps`, `mapdata`, `ggplot2`, `cowplot`

## Usage
1. **Load Data**: Import BOLD sequences and external temperature data into R.
2. **Data Preprocessing**: Clean and format data as needed.
3. **Analysis**: Use scripts in `scripts/` to conduct analyses, including:
   - Correlation between latitude, species diversity, and temperature
   - Trends in latitudinal distribution related to temperature changes over time
4. **Visualization**: Generate graphs to illustrate trends and correlations.

## Analysis Workflow
1. **Data Import and Cleaning**:
   - Load BOLD genetic data and temperature data.
   - Filter and preprocess for relevant regions and years.

2. **Analysis**:
   - Perform spatial analysis to examine latitudinal shifts.
   - Correlate species diversity with temperature data.

3. **Visualization**:
   - Map latitudinal distribution.
   - Plot trends over time.

## References
- [BOLD Systems](http://www.boldsystems.org/)
- Climate data source: [NASA Climate Change](https://climate.nasa.gov/vital-signs/global-temperature/?intent=121)

## Attribution
- **Primary Author**: Isha Baxi for developing Bumblebee Latitude Analysis 
- **Secondary Author**: [Vivian Phung](https://github.com/vivianp17) for improving data efficiency for species richness, expanding on filtering steps, creating an additional figure (Bubble plot for species richness and geographically) and other minor changes. 
