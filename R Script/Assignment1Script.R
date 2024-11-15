
# BINF6210 - Software Tools
# Assignment 1 
# Isha Baxi

#### DATA EXPLORATION ####

## ---- Load Packages ----

# Assume all libraries have been installed. 
#If not, install.packages("nameOfpackage")

library(tidyverse)
library(vegan) 
library(dplyr) # Used for data manipulations
library(ggplot2) # Used for visualizations 
library(plotly)  # Used for visualizations (specifically ggplotly())
library(sf)
library(maps)
library(mapdata)
library(cowplot)  # For dual-axis plotting

## ---- Load NASA Temp Data and Clean -----

# Read and Load NASA Average Temperature Data
nasa_data <- read.csv("../Data/nasa_data.csv", header = TRUE, stringsAsFactors = FALSE)

# Convert the data frame into proper columns
nasa_data <- as.data.frame(nasa_data)

nasa_data <- nasa_data %>% 
  filter(!is.na(X)& !is.na(X.1) & !is.na(X.2)) # Remove NA values

# Rename columns
names(nasa_data)[names(nasa_data) == "X"] <- "year"
names(nasa_data)[names(nasa_data) == "X.1"] <- "No_Smoothing"
names(nasa_data)[names(nasa_data) == "X.2"] <- "Lowess(5)"

nasa_data <- nasa_data[-1, ] # Remove first row

## ---- Load Data from BOLD API  ----

# Extracting Bombus data from Canada + USA in TSV data format
# Extracted data in TSV is also saved in /data folder 
BOLDdata <- read_tsv(file = "http://v4.boldsystems.org/index.php/API_Public/specimen?taxon=Bombus&geo=Canada,United%20States&format=tsv") 

# Note: BOLD is currently transitioning from version 4 to version 5. To maintain compatibility with RStudio, the link was adjusted to access the previous version. Future updates to the BOLD system may require further link adjustments to ensure continued access.

BOLDdata_filtered <- BOLDdata %>% 
  filter(country == "Canada" | country == "United States")

## ---- Data Cleaning + Processing  ----

# Creating data frame with necessary variables from BOLDdata
bombus <- BOLDdata_filtered %>% # Extracting BOLDdata for new bombus data frame
  select(processid, species_name, lat, lon, province_state) %>% # Selecting processid, species_name, lat, lon, province_state variables to add into bombus 
  mutate(across(where(is.character), ~ na_if(., " "))) %>% # Changing any empty cells in bombus to NA
  filter(!is.na(processid)& !is.na(species_name) & !is.na(lat) & !is.na(lon) & !is.na(province_state)) # Removing NA cells

# --- Summarize Key Variables and Check for Data Errors ---

# Summarize Numerical Variables 

# Latitude mean, max, min and range
mean(bombus[["lat"]])
max(bombus[["lat"]])
min(bombus[["lat"]])
range(bombus[['lat']])

# Longitude mean, max, min and range
mean(bombus[["lon"]])
max(bombus[["lon"]])
min(bombus[["lon"]])
range(bombus[['lon']])
# The ranges for latitude and longitude are accurate, as referenced by Google Maps. 

# Summarize species count for each specimen 
species_count <- bombus %>%
  count(species_name, sort = TRUE)
summarize(species_count)

sum(species_count$n) # Check to ensure that n = 3432, for each specimen in dataset being classified

# Location
location <- bombus %>%
  count(province_state, sort = TRUE)

sum(location$n) # Check to ensure that n = 3432, for each specimen in dataset being classified into location


# Filter out latitudes and longitudes that are out of range
bombus <- bombus %>%
  filter(lat >= 24 & lat <= 83) %>%  # North American latitudes range from roughly 24 to 83 
  filter(lon >= -172 & lon <= -52)   # North American longitudes range from roughly -172 to -52

# Boxplots for outliers of numerical variables 
boxplot(bombus$lat, main = "Latitude")
boxplot(bombus$lon, main = "Longitude")

# Plot Histograms 
ggplot(bombus, aes(x = species_name)) +
  geom_bar(fill = "skyblue", color = "black") +
  labs(title = "Histogram of Species Types",
       x = "Species",
       y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 0.5))

ggplot(bombus, aes(x = province_state)) +
  geom_bar(fill = "skyblue", color = "black") +
  labs(title = "Histogram of Province/States",
       x = "Province/State",
       y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 0.5))

## ---- Organizing into Data Frames + Matrices  ----

# Classify latitude into 3 groups (Low, Med, High)
latitude_group <- cut(bombus$lat,
                      breaks = c(25, 44, 63, 82), 
                      labels = c("Low", "Mid", "High"),
                      include.lowest = TRUE)
bombus$latitude_group <- latitude_group # Add coloumn to bombus data frame

# Using dplyr's count to improve efficiency in calculating species richness by latitude group
species_richness1 <- bombus %>%
  count(latitude_group, species_name) %>% # Count occurrences of each species in latitude groups
  spread(key = species_name, value = n, fill = 0) # Convert counts into a wide format

# Extract year from last 2 digits of processid (Assuming that the year was when the specimen was found)
year <- str_extract(bombus$processid, "\\d{2}$")
bombus$year <- year
bombus$year <- (paste0("20", bombus$year))  # Ensures that year is numeric when adding to bombus data frame

# Compare species name with frequency count of each species with average latitude
species_count <- bombus %>%
  group_by(species_name) %>%
  summarise(
    species_freq = n(),           # Count the frequency of each species
    avg_lat = mean(lat, na.rm = TRUE)  # Calculate average latitude for each species
  )

# Find the year with the most specimen caught
year <- bombus %>%
  count(year, sort = TRUE)
year

# Counts total species records for each year
species_freq_per_year <- bombus %>%
  group_by(year) %>%
  summarise(species_count = n())  

# Data frame merging average temperature data + species frequency by years in Common
combined_df <- merge(nasa_data, species_freq_per_year, by = "year")
combined_df$No_Smoothing <- as.numeric(as.character(combined_df$No_Smoothing)) # Convert Average to Numeric

# Data cleaning for combined data frame
combined_df <- combined_df %>% # Remove Lowess(5) Column
  select(-"Lowess(5)")

#### ANALYSIS TO ADDRESS QUESTION ####

## ---- Statistical Testing  ----

shannon_diversity <- diversity(species_richness, index = "shannon") # Conduct Shannon Index Test for 3 lat categories
shannon_diversity

# Spearman correlation pre-requisites
str(combined_df) # Check the structure of the combined_df dataframe

# Convert species_richness to numeric
combined_df$species_count <- as.numeric(as.character(combined_df$species_count))
combined_df$No_Smoothing <- as.numeric(as.character(combined_df$No_Smoothing))

# Boxplot to visualize potential outliers
boxplot(combined_df$species_count, main = "Species Richness")
boxplot(combined_df$No_Smoothing, main = "Average Temperature Change")

# Calculate Spearman correlation
spearman_corr <- cor(combined_df$species_count, combined_df$No_Smoothing, method = "spearman")
print(spearman_corr)

## ---- Visualizations ----

### PLOT 1: Species Richness Bubble Plot

# Calculate species richness per latitude, longitude, and species
species_richness_data <- bombus %>%
  group_by(lat,lon,species_name) %>% 
  summarise(SpeciesRichness = n(), .groups = 'drop')

#Getting world map
world <- st_as_sf(map("world", plot = FALSE, fill = TRUE))

# Create the bubble plot 
bubbleplot <- ggplot(data = species_richness_data) +
  geom_sf(data = world, fill = "grey90", color = "white") +
  geom_point(aes(x = lon, y = lat, size = SpeciesRichness, color = species_name), alpha = 0.8) +
  scale_size(range = c(2, 10), name = "Species Richness") +  # Adjust bubble size range
  scale_color_viridis_d(name = "Species") +        # Adjust color scale
  coord_sf(xlim = c(-150, -50), ylim = c(20, 85), expand = FALSE) + 
  labs(title = "Bubble Map of Bee Species Richness by Location",
       subtitle = "Point size and color indicate richness level",
       x = "Longitude",
       y = "Latitude") +
  theme_minimal() +
  theme(legend.position = "none")
new_plot <- ggplotly(bubbleplot, tooltip = c("x", "y", "size","colour")) # converting ggplot to plotly object
print(new_plot) # Display New Plot

### PLOT 2: Species Density Heat Map
species_location_data <- bombus %>%
  filter(!is.na(lat) & !is.na(lon) & !is.na(province_state))  # Filter out missing values

# Plot the map with a heat map
plot_2 <- ggplot() +
  # Add the map of Canada and the USA
  geom_sf(data = world, fill = "grey90", color = "white") +
  
  # Add heat map layer based on species density
  stat_density_2d(
    data = species_location_data,
    aes(x = lon, y = lat, fill = ..level..),  # 'fill' for heat map
    geom = "polygon",  # Use polygons for filled density
    alpha = 0.5       # Adjust transparency
  ) +
  
  # Add contour lines 
  geom_density_2d(
    data = species_location_data, 
    aes(x = lon, y = lat), 
    color = "red"   
  ) +
  
  # Add points for species locations
  geom_point(data = species_location_data, aes(x = lon, y = lat), size = 0.5, alpha = 0.7) +
  
  # Customize the labels and theme
  labs(
    title = "Species Locations with Density Heat Map in Canada and USA",
    # subtitle = "Visualizing Species Concentrations with Contour Lines",
    x = "Longitude (°W)",
    y = "Latitude (°N)",
    fill = "Density",
    color = "Province/State"
  ) +
  
  # Focus on Canada and USA region
  coord_sf(xlim = c(-150, -50), ylim = c(20, 85)) + 
  
  scale_fill_viridis_c(option = "plasma", name = "Density") +
  scale_color_viridis_d(name = "Province/State", option = "plasma") +
  
  # Improving overall aesthetics of map
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust=0.5),
    plot.subtitle = element_text(size = 12, face = "italic"),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    legend.position = "bottom",
    legend.key.width = unit(2, "cm")
  )

print(plot_2)

### PLOT 3: Time-Series Species Richness vs. Avg. Temperature

# Species Richness over time
species_richness_plot <- ggplot(combined_df, aes(x = year, y = species_count)) +
  geom_line(color = "blue", size = 1) +
  geom_point(color = "blue", size = 2) +
  scale_y_continuous(name = "Species Richness") +  # Y-axis for Species Richness
  theme_minimal() +
  theme(axis.title.y = element_text(color = "blue", size = 12))

# Second plot: Temperature Change over time (on the same x-axis)
temperature_change_plot <- ggplot(combined_df, aes(x = year, y = No_Smoothing)) +
  geom_line(color = "red", size = 1) +
  geom_point(color = "red", size = 2) +
  scale_y_continuous(name = "Temperature Change (°C)", 
                     sec.axis = sec_axis(~ ., name = "Temperature Change (°C)")) +  # Secondary axis for temperature
  theme_minimal() +
  theme(axis.title.y.right = element_text(color = "red", size = 12))

# Combine the plots
combined_plot <- plot_grid(
  species_richness_plot + theme(axis.title.x = element_blank()),  # Remove x-axis label from first plot
  temperature_change_plot + theme(axis.title.x = element_blank()),  # Remove x-axis label from second plot
  ncol = 1, align = "v", axis = "l"  # Align both plots vertically
)

# Add title
title <- ggdraw() + draw_label("Species Richness and Temperature Change Over Time", 
                               fontface = 'bold', size = 14)

# Add shared x-axis
x_axis <- ggdraw() + draw_label("Year", size = 12)

# Adding a annotation to Figure 3
annotation_text <- "Note: Trends in species richness may reflect sampling effort changes over time."
annotation<- ggdraw() + 
  draw_label(annotation_text, fontface = 'italic', size = 10, hjust = 0.5, color = "gray40")

# Final plot with annotation
plot_3 <- plot_grid(title, combined_plot, x_axis, annotation, ncol = 1, rel_heights = c(0.1, 1, 0.05, 0.1))
print(plot_3)
