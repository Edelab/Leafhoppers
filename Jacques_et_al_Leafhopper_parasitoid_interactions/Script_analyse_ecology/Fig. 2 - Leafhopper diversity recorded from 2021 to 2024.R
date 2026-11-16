rm(list = ls())

setwd('C:/Users/JOJAC54/OneDrive - Université Laval/Jacques et al. 2026/Codes GitHub')


# Packages
library(dplyr)
library(tidyr)
library(iNEXT)
library(ggplot2)
library(readxl)
library(tidyverse)
library(readxl)


# Biodiversity among years
data<-read_excel("div_4.xlsx")
data$Year<-as.factor(data$Year)
data$Farm<-as.factor(data$Farm)
str(data)

# Keep only genus data
species_cols <- names(data)[5:ncol(data)]

# Summarize data per year
data_yearly_aggregated <- data %>%
  dplyr::select(Year, all_of(species_cols)) %>% # Select only Year and species columns
  group_by(Year) %>% # Group by year
  summarise(across(all_of(species_cols), sum), .groups = 'drop') # Sum abundances for each species within each year

# Matrix for iNEXT
species_matrix_yearly <- data_yearly_aggregated %>%
  pivot_longer(
    cols = all_of(species_cols),
    names_to = "Species",
    values_to = "Abundance"
  ) %>%
  pivot_wider(
    names_from = Year, # Now the columns will be '2023', '2024'
    values_from = Abundance,
    values_fill = 0 # Ensure any missing species in a year get 0 abundance
  ) %>%
  column_to_rownames("Species") # Make species names the row names

# Check data
print(head(species_matrix_yearly))
print(dim(species_matrix_yearly)) # Should be [number of species] x [number of unique years, e.g., 2]

# Hill numbers (q=0 for richness, q=1 for Shannon, q=2 for Simpson)
iNEXT_yearly_results <- iNEXT(x = species_matrix_yearly,
                              q = c(0, 1, 2), # Calculate richness, Shannon (exp(H')), and Simpson (1/D)
                              datatype = "abundance")

# Basic data information
iNEXT_yearly_results$DataInfo

# See the statistical results
iNEXT_yearly_results$iNextEst

# Asymptotic estimates (the 'total' diversity for each year and order)
iNEXT_yearly_results$AsyEst


# To improve the plot

# Extract the coverage-based estimation results data frame directly
plot_data_coverage_based <- iNEXT_yearly_results$iNextEst$coverage_based


plot_data_coverage_based <- iNEXT_yearly_results$iNextEst$coverage_based



# Enhance the data for better plotting labels
plot_data_coverage_based_enhanced <- plot_data_coverage_based %>%
  mutate(
    # Create a descriptive label for Diversity Order (q)
    Diversity_Order = factor(Order.q,
                             levels = c(0, 1, 2),
                             labels = c("q = 0 (Genus Richness)",
                                        "q = 1 (Shannon Diversity)",
                                        "q = 2 (Simpson Diversity)"))
  )

# Filter the data to include only 'Rarefaction' and 'Extrapolation' methods
plot_data_filtered <- plot_data_coverage_based_enhanced %>%
  filter(Method %in% c("Rarefaction", "Extrapolation"))

year_colors <- c(
  "2021" = "#1b9e77",
  "2022" = "#d95f02",
  "2023" = "dodgerblue",
  "2024" = "#7570b3"
)

figdiv<-ggplot(plot_data_filtered, aes(x = SC, y = qD, color = Assemblage)) +
  geom_line(aes(linetype = Method), linewidth=0.8) +
  geom_ribbon(aes(ymin = qD.LCL, ymax = qD.UCL, fill = Assemblage), alpha = 0.2) +
  facet_wrap(~Diversity_Order, scales = "free_y") +
  labs(x = "Sample coverage", y = "Genus diversity") +
  scale_color_manual(name = "Year", values = year_colors) +
  scale_fill_manual(name = "Year", values = year_colors) +
  scale_linetype_manual(
    name = "Method", # Name of the legend title for line types
    values = c("Rarefaction" = "solid", "Extrapolation" = "dashed")
  ) +
  theme_bw(12) +
  xlim(0.2,1)+
  theme(legend.position = "top") 

figdiv

figdiv2<-ggplot(plot_data_filtered, aes(x = m, y = qD, color = Assemblage)) +
  geom_line(aes(linetype = Method), linewidth=0.8) +
  geom_ribbon(aes(ymin = qD.LCL, ymax = qD.UCL, fill = Assemblage), alpha = 0.2) +
  facet_wrap(~Diversity_Order, scales = "free_y") +
  labs(x = "Number of leafhoppers", y = "Genus diversity") +
  scale_color_manual(name = "Year", values = year_colors) +
  scale_fill_manual(name = "Year", values = year_colors) +
  scale_linetype_manual(
    name = "Method", # Name of the legend title for line types
    values = c("Rarefaction" = "solid", "Extrapolation" = "dashed")
  ) +
  theme_bw(12) +
  theme(legend.position = "top") 

figdiv2



ggsave("figdiv2.tiff", width = 180, height = 100, units = "mm",dpi=600, compression = "lzw")


# The total diversity estimated and the observed

# 1. Prepare data for Asymptotic Estimates
asy_est_data <- iNEXT_yearly_results$AsyEst %>%
  mutate(
    Diversity_Order = factor(Diversity,
                             levels = c("Species richness", "Shannon diversity", "Simpson diversity"),
                             labels = c("q = 0 (Genus Richness)",
                                        "q = 1 (Shannon Diversity)",
                                        "q = 2 (Simpson Diversity)")),
    Type = "Asymptotic"
  ) %>%
  rename(Value = Estimator, LCL_val = LCL, UCL_val = UCL) %>%
  dplyr::select(Assemblage, Diversity_Order, Value, LCL_val, UCL_val, Type)

# 2. Prepare data for Observed Values
obs_data <- iNEXT_yearly_results$iNextEst$size_based %>%
  filter(Method == "Observed") %>%
  mutate(
    Diversity_Order = factor(Order.q,
                             levels = c(0, 1, 2),
                             labels = c("q = 0 (Genus Richness)",
                                        "q = 1 (Shannon Diversity)",
                                        "q = 2 (Simpson Diversity)")),
    Type = "Observed"
  ) %>%
  rename(Value = qD, LCL_val = qD.LCL, UCL_val = qD.UCL) %>%
  dplyr::select(Assemblage, Diversity_Order, Value, LCL_val, UCL_val, Type)

# 3. Combine both data frames
combined_diversity_data <- bind_rows(asy_est_data, obs_data)

estimate_type_colors <- c("Asymptotic" = "grey60", "Observed" = "black") # Choose colors you like
estimate_type_shapes <- c("Asymptotic" = 16, "Observed" = 16)

# 4. Plot the combined data with points and error bars
figest<-ggplot(combined_diversity_data, aes(x = factor(Assemblage), y = Value, color = Type)) +
  geom_point(aes(shape = Type),
             position = position_dodge(width = 0.8), # Dodge points to prevent overlap
             size = 3) + # Adjust point size for visibility
  geom_errorbar(aes(ymin = LCL_val, ymax = UCL_val),
                position = position_dodge(width = 0.8), # Dodge error bars to match points
                width = 0.25, # Width of the horizontal caps on the error bars
                linewidth = 0.8) + # Thickness of the error bar lines
  facet_wrap(~ Diversity_Order, scales = "free_y") + # Separate plots for each diversity order
  labs(x = "Year",
       y = "Genus diversity") +
  scale_color_manual(name = "Estimate type", values = estimate_type_colors) +
  scale_shape_manual(name = "Estimate type", values = estimate_type_shapes) +

theme_classic(base_size = 12) +
  theme(axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 0, hjust = 0.5), # Rotate x-axis labels for readability
        legend.position = "top",
        axis.title = element_text(size = 12, face = "bold"),
        text = element_text(size = 12),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        strip.background = element_rect(fill = "grey80", color = "black", linewidth = 0.3),
        panel.border = element_rect(color = "black", fill = NA))


figest


ggsave("figest.tiff", width = 180, height = 100, units = "mm",dpi=600, compression = "lzw")


#### The five most abundance genus per year
data<-read_excel("div_4.xlsx")
data$Year<-as.factor(data$Year)
data$Farm<-as.factor(data$Farm)
str(data)


genus_cols_names <- names(data)[!names(data) %in% c("Year", "Farm", "Traps", "Days")]

top_genera_data <- data %>%
  pivot_longer(
    cols = all_of(genus_cols_names), 
    names_to = "Genus",              
    values_to = "Count"              
  ) %>%
  group_by(Year, Genus) %>%
  summarise(Total_Count = sum(Count, na.rm = TRUE), .groups = 'drop') %>%
  filter(Total_Count > 0) %>%
  group_by(Year) %>%
  arrange(desc(Total_Count)) %>%  
  slice_head(n = 5) %>%
  ungroup()               

# Custom colors for the genera
custom_colors <- c(
  "Empoasca" = "#1f78b4",
  "Macrosteles" = "#e31a1c",
  "Graminella" = "#33a02c",
  "Hebata" = "#fb9a99",
  "Xestocephalus" = "#ff7f00",
  "Erythroneura" = "#b2df8a",
  "Paraphlepsius" = "#fdbf6f",
  "Draeculacephala" = "#a6cee3"
)

figtop5 <- ggplot(top_genera_data, aes(x = reorder(Genus, Total_Count), y = Total_Count, fill = Genus)) +
  geom_col() + 
  facet_wrap(~ Year, scales = "free_y") +
  coord_flip() + 
  labs(x = "Genus", y = "Total abundance") +
  scale_fill_manual(values = custom_colors) +  # ici !
  theme_bw(base_size = 12) + 
  theme(
    legend.position = "none", 
    axis.text.y = element_text(size = 9, face = "italic"), 
    axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        text = element_text(size = 12))

figtop5

ggsave("figtop5.tiff", width = 160, height = 100, units = "mm",dpi=600, compression = "lzw")


#### Presence-Absence----
library(dplyr)   
library(tidyr)   
library(ggplot2) 

genus_cols_names <- names(data)[!names(data) %in% c("Year", "Farm", "Traps", "Days")]

genus_pa_data <- data %>%
  pivot_longer(
    cols = all_of(genus_cols_names),
    names_to = "Genus",
    values_to = "Count"
  ) %>%
  group_by(Year, Genus) %>%
  summarise(Total_Count = sum(Count, na.rm = TRUE), .groups = 'drop') %>% # Sum counts, handle NAs, drop grouping
  mutate(Presence = ifelse(Total_Count > 0, 1, 0)) %>%
  mutate(Year = factor(Year)) %>%
  mutate(Genus = factor(Genus, levels = sort(unique(Genus))))

figabsprs<-ggplot(genus_pa_data, aes(x = Year, fct_rev(fct_inorder(Genus)), fill = factor(Presence))) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_manual(
    values = c("0" = "lightgrey", "1" = "grey40"),
    labels = c("Absent", "Present"),
    name = "Presence/Absence"
  ) +
  labs(x = "Year", y = "Genus"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
    axis.text.y = element_text(size = 10, face="italic"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_fixed(ratio = 0.3)

figabsprs

ggsave("figabsprs.tiff", width = 120, height = 250, units = "mm",dpi=600, compression = "lzw")

fig_abundance_presence <- ggplot(genus_pa_data, aes(x = Year, y = fct_rev(fct_inorder(Genus)), fill = factor(Presence))) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(Total_Count > 0, Total_Count, "")), color = "black", size = 3) + # Add this line
  scale_fill_manual(
    values = c("0" = "black", "1" = "gray90"),
    labels = c("Absent", "Present"),
    name = "Presence/Absence"
  ) +
  labs(x = "Year", y = "Genus") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 10, face = "italic"),
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_fixed(ratio = 0.3)

fig_abundance_presence


ggsave("fig_abundance_presence.tiff", width = 120, height = 250, units = "mm",dpi=600, compression = "lzw")


# Heatmap

genus_relative_abundance <- genus_pa_data %>%
  group_by(Year) %>%
  mutate(Relative_Abundance = Total_Count / sum(Total_Count, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(Relative_Abundance = ifelse(is.nan(Relative_Abundance), 0, Relative_Abundance))


fheat<-ggplot(genus_relative_abundance, aes(x = Year, y = fct_rev(fct_inorder(Genus)), fill = log10(Relative_Abundance + 1e-6))) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient(low = "white", high = "black",
                      name = "Log10",
                      labels = scales::label_scientific(digits = 1)) +  labs(x = "Year",
                                                                             y = "Genus") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, face = "italic", color = "black"),
        panel.grid = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))
fheat

ggsave("fheat.tiff", width = 100, height = 250, units = "mm",dpi=600, compression = "lzw")

