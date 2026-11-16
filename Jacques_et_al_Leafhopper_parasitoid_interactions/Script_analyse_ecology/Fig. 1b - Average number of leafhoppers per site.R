# Clear all objects from the global environment
rm(list = ls())

# Set the file path for the data
setwd("C:/Users/Adminwork/OneDrive - Université Laval/Maîtrise/Analyses dans R/Analyses_finales")

# Load the libraries
#remotes::install_github("teunbrand/ggh4x") # Force the installation of an old package
library(readxl)
library(ggplot2)
library(dplyr)
library(openxlsx)
library(forcats)
library(tidyverse)
library(ggh4x)
library(patchwork)
library(grid)
library(gridExtra)

# Figure 1 - Average number of leafhoppers by year and by month ----

# Load the data from the Excel file (sheet 5)
leafhoppers_raw <- read_excel("Cicadelles_final.xlsx", sheet = 5)

# Create a dataframe with the regions and the number of sites
new_data <- data.frame(
  year = c("2023", "2024"),
  nb_sites = c(18, 16))

# Filter the data to exclude the "FHG" farm
# Group the data by year, week of trap removal (experimental week) and group
# Calculate the average number of leafhopper by year, week and group considering the nb. of site each year
# Here the number of the week is based on 7 days week, starting from May 1st. 
# So May 1st to May 7th is week 1, May 8th to May 14th is week 2, May 15th to May 21th is week 3,
# May 22th to May 28th is week 4, May 29th to June 4th is week 5, ect.

leafhoppers1 <- leafhoppers_raw %>%
  filter(farm != "FHG") %>%  # Exclude the "FHG" farm
  mutate(year = as.character(year)) %>%  # Convert "year" to character
  group_by(year, week, group) %>%  # Group by year, week, and group
  summarise(abundance = sum(number), .groups = 'drop') %>%  # Sum the abundance for each group
  left_join(new_data, by = "year") %>%  # Merge with "nb_site" data
  mutate(average = abundance / nb_sites)  # Calculate the average abundance per site

# Establish the order of the species for the graph
leafhoppers1 <- leafhoppers1 %>%
  mutate(group = fct_relevel(group, "Empoasca fabae", "Macrosteles quadrilineatus", "Local species"))

# Assign each week to its corresponding month
leafhoppers1 <- leafhoppers1 %>%
  mutate(month = case_when(
    week >= 1  & week <= 4  ~ "May",
    week >= 5  & week <= 8  ~ "Jun",
    week >= 9  & week <= 13 ~ "Jul",
    week >= 14 & week <= 17 ~ "Aug",
    week >= 18 & week <= 21 ~ "Sept",
    week >= 22 & week <= 26 ~ "Oct",
    week >= 27 & week <= 30 ~ "Nov",
    TRUE ~ NA_character_  # Assign NA to any week outside the expected range
  ))

# Verify if each week had a month attributed
table(leafhoppers1$month)

# Convert "month" to a factor to ensure correct ordering in the plot
leafhoppers1 <- leafhoppers1 %>%
  mutate(month = factor(month, levels = c("May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov")))


# Plot 1 - Comparison between species
Figure_1_final <- ggplot(leafhoppers1, aes(x = week, y = average, color = group)) +
  geom_point(size = 3) +  # Point size
  geom_line(aes(group = group), linewidth = 1) +  # Line width
  ggh4x::facet_grid2(group ~ year, scales = "free_y", independent = "y", space = "free_x") +  # Allows a scale adapted to each facet
  scale_x_continuous(
    breaks = c(1, 5, 9, 14, 18, 22,27),  
    labels = c("May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov")  
  ) +
  scale_color_manual(
    values =  c("Empoasca fabae" = "#7fbf5b",   # Red
                "Local species" = "#1f74b4",   # Blue
                "Macrosteles quadrilineatus" = "#e27a1c")
  ) +
  labs(
    title = "",
    x = "Month",
    y = "Average number of leafhoppers per site",
    color = "Species"  # Legend with the title "Species"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),  # Facet title size
    axis.text = element_text(size = 10),  # Axis text size
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  # Larger and centered main title
    plot.subtitle = element_text(size = 12, hjust = 0.5),  # Subtitle
    axis.title = element_text(size = 12, face = "bold"),  # More visible axis titles
    legend.position = "none",  # Remove legend
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank())   # Remove minor gridlines


Figure_1_final

ggsave("Abundance_per_site.tiff", width = 220, height = 190, units = "mm",dpi=600, compression = "lzw")

