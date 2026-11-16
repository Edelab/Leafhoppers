rm(list = ls())

setwd("C:/Users/Adminwork/OneDrive - Université Laval/Maîtrise/Analyses dans R/Analyses_finales")

#library
library(car)
library(lme4)
library(DHARMa)
library(readxl)
library(MuMIn)
library(ggplot2)
library(tidyr)
library(dplyr)
library(glmmTMB)
library(splines)
library(performance)
library(spdep)
library(ggeffects)
library(esquisse)
library(lubridate)
library(broom.mixed)
library(mgcv)
library(patchwork)

#The model with species interaction had a high variance inflation factor
#Even changing the polynomial terms the VIF was > 42 (not acceptable)
#So we will model each species 
#The first model will be the abundance overtime and test the abundance between years

# Abundance of each group x year ----
data<-read_excel("data1_clean.xlsx", 2)
data$Year<-as.factor(data$Year)
data$Site<-as.factor(data$Site)
data$Farm<-as.factor(data$Farm)

# Convert wide format to long format
long_data <- data %>%
  pivot_longer(cols = c(Mqa, Efa, Localsp), 
               names_to = "species", 
               values_to = "abundance")

long_data <- long_data %>%
  mutate(NumTraps = 5)

head(long_data)

long_data_wea <- long_data %>%
  group_by(Farm, Year, Day) %>%
  mutate(TotalAbundance = sum(abundance)) %>%  # Compute total abundance across species
  ungroup() %>%
  group_by(Farm, Year, species, Day) %>%
  summarize(
    abundance = sum(abundance),  # Sum species abundance per farm per day
    NumTraps = sum(NumTraps),  # Sum number of traps per farm per day
    duration = mean(duration),  # Average trap duration
    T2M_L0= mean (T2M_L0),
    T2M_L7= mean (T2M_L7),
    T2M_L15= mean (T2M_L15),
    T2M_L30= mean (T2M_L30),
    Prec_L0= mean (Prec_L0),
    Prec_L7= mean (Prec_L7),
    Prec_L15= mean (Prec_L15),
    Prec_L30= mean (Prec_L30),
    WS2M_L0= mean (WS2M_L0),
    WS2M_L7= mean (WS2M_L7),
    WS2M_L15= mean (WS2M_L15),
    WS2M_L30= mean (WS2M_L30),
    lat = mean(lat),  # Mean latitude for the farm
    long = mean(long),  # Mean longitude for the farm
    .groups = "drop"
  )

long_data_wea$species<-as.factor(long_data_wea$species)

# Subset the data for each group of species

defa<-subset(long_data_wea, species=="Efa")
dloc<-subset(long_data_wea, species=="Localsp")
dmqa<-subset(long_data_wea, species=="Mqa")

# Raw abundance----
defa %>%
  group_by(Farm, Year, Day) %>%
  summarise(TotalAbundance = sum(abundance), .groups = "drop") %>%
  ggplot(aes(x = Day, y = TotalAbundance, color = Year, group = Year)) +
  geom_line(size = 0.9) +
  facet_wrap(~Farm, scales = "free_y") +
  labs(title = "Empoasca fabae",
       x = "Julian day", y = "Abundance") +
  theme_minimal(12)

dloc %>%
  group_by(Farm, Year, Day) %>%
  summarise(TotalAbundance = sum(abundance), .groups = "drop") %>%
  ggplot(aes(x = Day, y = TotalAbundance, color = Year, group = Year)) +
  geom_line(size = 0.9) +
  facet_wrap(~Farm, scales = "free_y") +
  labs(title = "Local species",
       x = "Julian day", y = "Abundance") +
  theme_minimal(12)

dmqa %>%
  group_by(Farm, Year, Day) %>%
  summarise(TotalAbundance = sum(abundance), .groups = "drop") %>%
  ggplot(aes(x = Day, y = TotalAbundance, color = Year, group = Year)) +
  geom_line(size = 0.9) +
  facet_wrap(~Farm, scales = "free_y") +
  labs(title = "Macrosteles quadrilineatus",
       x = "Julian day", y = "Abundance") +
  theme_minimal(12)

# Seasonal abundance

# Empoasca----

mod1.efa<-glmmTMB(abundance ~ ns(Day, df=2) + 
                  offset(log(NumTraps * duration)) + (1 | Farm)
                  + (1 | Year),
                family = nbinom2, data = defa)

mod2.efa<-glmmTMB(abundance ~ ns(Day, df=3) + 
                    offset(log(NumTraps * duration)) + (1 | Farm)
                  + (1 | Year),
                  family = nbinom2, data = defa)

mod3.efa<-glmmTMB(abundance ~ ns(Day, df=4) + 
                    offset(log(NumTraps * duration)) + (1 | Farm)
                  + (1 | Year),
                  family = nbinom2, data = defa)

compare_performance(mod1.efa, mod2.efa, mod3.efa, rank=TRUE)

# mod 3 is the best
res_mod3.efa <- simulateResiduals(mod3.efa)
testDispersion(res_mod3.efa)
testZeroInflation(res_mod3.efa)
testUniformity(res_mod3.efa)

# Spatial autocorrelation
location <- paste0(defa$long, "_", defa$lat)
res_grouped <- recalculateResiduals(res_mod3.efa, group = location)
testSpatialAutocorrelation(res_grouped, x = unique(defa$long), y = unique(defa$lat))

# Temporal autocorrelation
res_agg <- recalculateResiduals(res_mod3.efa, group = defa$Day)
testTemporalAutocorrelation(res_agg, time = unique(defa$Day))

# Local----
mod1.loc<-glmmTMB(abundance ~ ns(Day, df=2) + 
                    offset(log(NumTraps * duration)) + (1 | Farm)
                  + (1 | Year),
                  family = nbinom2, data = dloc)

mod2.loc<-glmmTMB(abundance ~ ns(Day, df=3) + 
                    offset(log(NumTraps * duration)) + (1 | Farm)
                  + (1 | Year),
                  family = nbinom2, data = dloc)

mod3.loc<-glmmTMB(abundance ~ ns(Day, df=4) + 
                    offset(log(NumTraps * duration)) + (1 | Farm)
                  + (1 | Year),
                  family = nbinom2, data = dloc)

compare_performance(mod1.loc, mod2.loc, mod3.loc, rank=TRUE)

# mod 3 is the best
res_mod3.loc <- simulateResiduals(mod3.loc)
testDispersion(res_mod3.loc)
testZeroInflation(res_mod3.loc)
testUniformity(res_mod3.loc)

# Spatial autocorrelation
location <- paste0(dloc$long, "_", dloc$lat)
res_grouped <- recalculateResiduals(res_mod3.loc, group = location)
testSpatialAutocorrelation(res_grouped, x = unique(dloc$long), y = unique(dloc$lat))

# Temporal autocorrelation
res_agg <- recalculateResiduals(res_mod3.loc, group = dloc$Day)
testTemporalAutocorrelation(res_agg, time = unique(dloc$Day))

# We have a problem with the residuals

mod3.p.loc <- glmmTMB(abundance ~ poly(Day, 3) + 
                      offset(log(NumTraps * duration)) + 
                      (1 | Farm) + 
                      (1 | Year),
                    family = nbinom2, data = dloc)

res_mod3.p.loc <- simulateResiduals(mod3.p.loc)
testDispersion(res_mod3.p.loc)
testZeroInflation(res_mod3.p.loc)
testUniformity(res_mod3.p.loc)

# Spatial autocorrelation
location <- paste0(dloc$long, "_", dloc$lat)
res_grouped <- recalculateResiduals(res_mod3.p.loc, group = location)
testSpatialAutocorrelation(res_grouped, x = unique(dloc$long), y = unique(dloc$lat))

# Temporal autocorrelation
res_agg <- recalculateResiduals(res_mod3.p.loc, group = dloc$Day)
testTemporalAutocorrelation(res_agg, time = unique(dloc$Day))

summary(mod3.p.loc)

compare_performance(mod1.loc, mod2.loc, mod3.loc, mod3.p.loc, rank=TRUE)

# Keep mod3.loc
res_mod3.p.loc <- simulateResiduals(mod3.p.loc)
testDispersion(res_mod3.p.loc)
testZeroInflation(res_mod3.p.loc)
testUniformity(res_mod3.p.loc)

# Spatial autocorrelation
location <- paste0(dloc$long, "_", dloc$lat)
res_grouped <- recalculateResiduals(res_mod3.p.loc, group = location)
testSpatialAutocorrelation(res_grouped, x = unique(dloc$long), y = unique(dloc$lat))

# Temporal autocorrelation
res_agg <- recalculateResiduals(res_mod3.p.loc, group = dloc$Day)
testTemporalAutocorrelation(res_agg, time = unique(dloc$Day))


# Macrosteles----

mod1.mqa<-glmmTMB(abundance ~ ns(Day, df=2) + 
                    offset(log(NumTraps * duration)) + (1 | Farm)
                  + (1 | Year),
                  family = nbinom2, data = dmqa)

mod2.mqa<-glmmTMB(abundance ~ ns(Day, df=3) + 
                    offset(log(NumTraps * duration)) + (1 | Farm)
                  + (1 | Year),
                  family = nbinom2, data = dmqa)

mod3.mqa<-glmmTMB(abundance ~ ns(Day, df=4) + 
                    offset(log(NumTraps * duration)) + (1 | Farm)
                  + (1 | Year),
                  family = nbinom2, data = dmqa)

compare_performance(mod1.mqa, mod2.mqa, mod3.mqa, rank=TRUE)

# mod 3 is the best
res_mod3.mqa <- simulateResiduals(mod3.mqa)
testDispersion(res_mod3.mqa)
testZeroInflation(res_mod3.mqa)
testUniformity(res_mod3.mqa)

# Spatial autocorrelation
location <- paste0(dloc$long, "_", dmqa$lat)
res_grouped <- recalculateResiduals(res_mod3.mqa, group = location)
testSpatialAutocorrelation(res_grouped, x = unique(dmqa$long), y = unique(dmqa$lat))

# Temporal autocorrelation
res_agg <- recalculateResiduals(res_mod3.mqa, group = dmqa$Day)
testTemporalAutocorrelation(res_agg, time = unique(dmqa$Day))

# Exporting model results for supp. material and manuscript

# Extract fixed coefficients fixos CI
tab_efa <- tidy(mod3.efa, effects = "fixed", conf.int = TRUE) |> mutate(species = "Efa")
tab_loc <- tidy(mod3.p.loc, effects = "fixed", conf.int = TRUE) |> mutate(species = "Localsp")
tab_mqa <- tidy(mod3.mqa, effects = "fixed", conf.int = TRUE) |> mutate(species = "Mqa")

# Combine 
tab_all <- bind_rows(tab_efa, tab_loc, tab_mqa)

# Extract R2 from each model
r2_all <- bind_rows(
  r2(mod3.efa) |> as.data.frame() |> mutate(species = "Efa"),
  r2(mod3.p.loc) |> as.data.frame() |> mutate(species = "Localsp"),
  r2(mod3.mqa) |> as.data.frame() |> mutate(species = "Mqa")
)

# Combine with the first table
tab_final <- left_join(tab_all, r2_all, by = "species")

# Rename columns to be more easier
tab_final <- tab_final |> 
  rename(
    Term = term,
    Estimate = estimate,
    Std_Error = std.error,
    CI_lower = conf.low,
    CI_upper = conf.high,
    P_value = p.value,
    R2_marginal = R2_marginal,
    R2_conditional = R2_conditional
  ) |> 
  select(species, Term, Estimate, Std_Error, CI_lower, CI_upper, P_value, R2_marginal, R2_conditional)

# Export
write.csv(tab_final, "model_results_with_r2.csv", row.names = FALSE)

# Plot----

# Create a grid with fixed offset (e.g.,: 10 traps per 7 days)
newdata <- data.frame(
  Day = seq(min(defa$Day), max(defa$Day), by = 1),
  NumTraps = 10,
  duration = 7,
  Farm = NA,
  Year = NA
)

# Predictions for Empoasca fabae
pred_efa <- predict(mod3.efa, newdata = newdata, type = "response", se.fit = TRUE)
efa_df <- newdata %>%
  mutate(species = "Efa",
         fit = pred_efa$fit,
         se = pred_efa$se.fit,
         lower = pmax(fit - 1.96 * se, 0),
         upper = fit + 1.96 * se)

# Predictions for Local
pred_loc <- predict(mod3.p.loc, newdata = newdata, type = "response", se.fit = TRUE)
loc_df <- newdata %>%
  mutate(species = "Localsp",
         fit = pred_loc$fit,
         se = pred_loc$se.fit,
         lower = fit - 1.96 * se,
         upper = fit + 1.96 * se)

# Predictions for Macrosteles
pred_mqa <- predict(mod3.mqa, newdata = newdata, type = "response", se.fit = TRUE)
mqa_df <- newdata %>%
  mutate(species = "Mqa",
         fit = pred_mqa$fit,
         se = pred_mqa$se.fit,
         lower = fit - 1.96 * se,
         upper = fit + 1.96 * se)

# Combine all predictions
pred_all <- bind_rows(efa_df, loc_df, mqa_df)

# Change to the name of the species
pred_all$species <- dplyr::recode(pred_all$species,
                           "Efa" = "Empoasca fabae",
                           "Mqa" = "Macrosteles quadrilineatus",
                           "Localsp" = "Local species",)

# Define colors
cores_especies <- c("Empoasca fabae" = "#7fbf5b",
                    "Macrosteles quadrilineatus" = "#e27a1c",
                    "Local species" = "#1f74b4")


# All species in the same graphic

fig1 <- ggplot(pred_all, aes(x = Day, y = fit, color = species, fill = species)) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
  scale_color_manual(values = cores_especies, 
                     breaks = c("Empoasca fabae", "Macrosteles quadrilineatus", "Local species")) +
  scale_fill_manual(values = cores_especies, 
                    breaks = c("Empoasca fabae", "Macrosteles quadrilineatus", "Local species")) +
  scale_x_continuous(
    name = "Month",
    breaks = c(121, 152, 182, 213, 244, 274, 305),
    labels = c("May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov")
  ) +
  scale_y_continuous(
    name = "Abundance per farm",
    breaks = seq(0, 300, by = 50),
    limits = c(0, 300)
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.title = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "top",
    text = element_text(size = 12),
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank()   # Remove minor gridlines
  )

fig1

ggsave("Abundance_per_month2.tiff", width = 170, height = 100, units = "mm",dpi=600, compression = "lzw")

ggsave("Abundance_per_month_large2.tiff", width = 280, height = 120, units = "mm",dpi=600, compression = "lzw")



# Abundance x weather----

#### Clean the environment and left defa, deloc dmqa

# Raw plot with the weather at 30 days before
defa$species <- "Efa"
dloc$species <- "Localsp"
dmqa$species <- "Mqa"

# Combine them all
all_data <- rbind(defa, dloc, dmqa)

all_data$norm_abund <- all_data$abundance / (all_data$NumTraps * all_data$duration)

ggplot(all_data, aes(x = T2M_L30, y = norm_abund)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", se = TRUE, color = "#ee7474") +
  facet_wrap(~ species, nrow = 3) +
  xlim(15,23)+
  ylab("Abundance per trap-day") +
  xlab("Temperature (°C, 30-day lag)") +
  theme_bw()

ggplot(all_data, aes(x = Prec_L30, y = norm_abund)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", se = TRUE, color = "#8bc2d1") +
  facet_wrap(~ species, nrow = 3) +
  ylab("Abundance per trap-day") +
  xlab("Precipitation (mm, 30-day lag)") +
  theme_bw()

# Empoasca fabae----
defa$DayF <- factor(defa$Day) 

mod1.efaw <- glmmTMB(abundance ~ T2M_L0 + Prec_L0 +
                           poly(Day, 2)+
                          offset(log(NumTraps * duration)) + 
                       ar1(DayF + 0 | Farm) + (1 | Year),
                        family = nbinom2, data = defa)

mod2.efaw <- glmmTMB(abundance ~ T2M_L7 + Prec_L7 +
                       poly(Day, 2)+
                       offset(log(NumTraps * duration)) + 
                       ar1(DayF + 0 | Farm) + (1 | Year),
                     family = nbinom2, data = defa)

mod3.efaw <- glmmTMB(abundance ~ T2M_L15 + Prec_L15 +
                       poly(Day, 2)+
                       offset(log(NumTraps * duration)) + 
                       ar1(DayF + 0 | Farm) + (1 | Year),
                     family = nbinom2, data = defa)

mod4.efaw <- glmmTMB(abundance ~ T2M_L30 + Prec_L30 +
                       poly(Day, 2)+
                       offset(log(NumTraps * duration)) + 
                       ar1(DayF + 0 | Farm) + (1 | Year),
                     family = nbinom2, data = defa)

AICc(mod1.efaw, mod2.efaw, mod3.efaw, mod4.efaw)

# mod4
res_mod4 <- simulateResiduals(mod4.efaw)
testDispersion(res_mod4)
testZeroInflation(res_mod4)
testUniformity(res_mod4)
check_collinearity(mod4.efaw)

# Spatial autocorrelation
location <- paste0(defa$long, "_", defa$lat)
res_grouped <- recalculateResiduals(res_mod4, group = location)
testSpatialAutocorrelation(res_grouped, x = unique(defa$long), y = unique(defa$lat))

# Temporal autocorrelation
res_agg <- recalculateResiduals(res_mod4 , group = defa$Day)
testTemporalAutocorrelation(res_agg, time = unique(defa$Day))

summary(mod4.efaw)

# Local----
dloc$DayF <- factor(dloc$Day) 

mod1.locw <- glmmTMB(abundance ~ T2M_L0 + Prec_L0 +
                       poly(Day, 2)+
                       offset(log(NumTraps * duration)) + 
                       ar1(DayF + 0 | Farm) + (1 | Year),
                     family = nbinom2, data = dloc)

mod2.locw <- glmmTMB(abundance ~ T2M_L7 + Prec_L7 +
                       poly(Day, 2)+
                       offset(log(NumTraps * duration)) + 
                       ar1(DayF + 0 | Farm) + (1 | Year),
                     family = nbinom2, data = dloc)

mod3.locw <- glmmTMB(abundance ~ T2M_L15 + Prec_L15 +
                       poly(Day, 2)+
                       offset(log(NumTraps * duration)) + 
                       ar1(DayF + 0 | Farm) + (1 | Year),
                     family = nbinom2, data = dloc)

mod4.locw <- glmmTMB(abundance ~ T2M_L30 + Prec_L30 +
                       poly(Day, 2)+
                       offset(log(NumTraps * duration)) + 
                       ar1(DayF + 0 | Farm) + (1 | Year),
                     family = nbinom2, data = dloc)

AICc(mod1.locw, mod2.locw, mod3.locw, mod4.locw)

# mod4
res_mod4 <- simulateResiduals(mod4.locw)
testDispersion(res_mod4)
testZeroInflation(res_mod4)
testUniformity(res_mod4)
check_collinearity(mod4.locw)

# Spatial autocorrelation
location <- paste0(dloc$long, "_", dloc$lat)
res_grouped <- recalculateResiduals(res_mod4, group = location)
testSpatialAutocorrelation(res_grouped, x = unique(dloc$long), y = unique(dloc$lat))

# Temporal autocorrelation
res_agg <- recalculateResiduals(res_mod4 , group = dloc$Day)
testTemporalAutocorrelation(res_agg, time = unique(dloc$Day))

summary(mod4.locw)

# Macrosteles quadrilineatus----
dmqa$DayF <- factor(dmqa$Day) 

mod1.mqaw <- glmmTMB(abundance ~ T2M_L0 + Prec_L0 +
                       poly(Day, 2)+
                       offset(log(NumTraps * duration)) + 
                       ar1(DayF + 0 | Farm) + (1 | Year),
                     family = nbinom2, data = dmqa)

mod2.mqaw <- glmmTMB(abundance ~ T2M_L7 + Prec_L7 +
                       poly(Day, 2)+
                       offset(log(NumTraps * duration)) + 
                       ar1(DayF + 0 | Farm) + (1 | Year),
                     family = nbinom2, data = dmqa)

mod3.mqaw <- glmmTMB(abundance ~ T2M_L15 + Prec_L15 +
                       poly(Day, 2)+
                       offset(log(NumTraps * duration)) + 
                       ar1(DayF + 0 | Farm) + (1 | Year),
                     family = nbinom2, data = dmqa)

mod4.mqaw <- glmmTMB(abundance ~ T2M_L30 + Prec_L30 +
                       poly(Day, 2)+
                       offset(log(NumTraps * duration)) + 
                       ar1(DayF + 0 | Farm) + (1 | Year),
                     family = nbinom2, data = dmqa)

AICc(mod1.mqaw, mod2.mqaw, mod3.mqaw, mod4.mqaw)


# mod4
res_mod4 <- simulateResiduals(mod4.mqaw)
testDispersion(res_mod4)
testZeroInflation(res_mod4)
testUniformity(res_mod4)
check_collinearity(mod4.mqaw)

# Spatial autocorrelation
location <- paste0(dmqa$long, "_", dmqa$lat)
res_grouped <- recalculateResiduals(res_mod4, group = location)
testSpatialAutocorrelation(res_grouped, x = unique(dmqa$long), y = unique(dmqa$lat))

# Temporal autocorrelation
res_agg <- recalculateResiduals(res_mod4 , group = dmqa$Day)
testTemporalAutocorrelation(res_agg, time = unique(dmqa$Day))

# Overall summary of the models
summary(mod4.efaw)
summary(mod4.locw)
summary(mod4.mqaw)

mod4.efaw_simple <- glmmTMB(
  abundance ~ T2M_L30 + Prec_L30 + poly(Day, 2) +
    offset(log(NumTraps * duration)) +
    (1 | Farm) + (1 | Year),
  family = nbinom2, data = defa)

r.squaredGLMM(mod4.efaw_simple)

mod4.locw_simple <- glmmTMB(
  abundance ~ T2M_L30 + Prec_L30 + poly(Day, 2) +
    offset(log(NumTraps * duration)) +
    (1 | Farm) + (1 | Year),
  family = nbinom2, data = dloc)

r.squaredGLMM(mod4.locw_simple)

mod4.mqaw_simple <- glmmTMB(
  abundance ~ T2M_L30 + Prec_L30 + poly(Day, 2) +
    offset(log(NumTraps * duration)) +
    (1 | Farm) + (1 | Year),
  family = nbinom2, data = dmqa)

r.squaredGLMM(mod4.mqaw_simple)


# Extract coefficients from each model
model_results <- list(
  efaw = tidy(mod4.efaw, conf.int = TRUE),
  locw = tidy(mod4.locw, conf.int = TRUE),
  mqaw = tidy(mod4.mqaw, conf.int = TRUE)
)

all_results <- bind_rows(model_results, .id = "species")

species_colors <- c("Empoasca fabae" = "#7fbf5b",   # Red
                    "Local species" = "#1f74b4",   # Blue
                    "Macrosteles quadrilineatus" = "#e27a1c"               
)

species_names <- c(
  "mqaw" = "Macrosteles quadrilineatus",
  "locw" = "Local species",
  "efaw" = "Empoasca fabae"
)

all_results <- bind_rows(model_results, .id = "species_code") %>%
  mutate(species = factor(species_names[species_code],
                          levels = c("Local species",  # Premier
                                     "Macrosteles quadrilineatus",  # Deuxième
                                     "Empoasca fabae")))  # Troisième

# Filter for just the variables of interest
temp_effects <- all_results %>% 
  filter(term %in% c("T2M_L30")) %>% 
  mutate(variable = "Temperature")

precip_effects <- all_results %>% 
  filter(term %in% c("Prec_L30")) %>% 
  mutate(variable = "Precipitation")

temp_plot <- ggplot(temp_effects, aes(x = estimate, y = species, 
                                      color = species)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), 
                 height = 0.2, size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = species_colors) +
  labs(x = "Coefficient (Erreur standard)", 
       y = NULL,
       title = "Température (30 jours avant)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        axis.text = element_text(size = 10),  # Taille des axes
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 12, face = "bold"),  # Titre principal plus grand et centré
        axis.title = element_text(size = 12))  # Titres des axes plus visibles

temp_plot

ggsave("Abundance_per_temperature.tiff", width = 120, height = 90, units = "mm",dpi=600, compression = "lzw")

precip_plot <- ggplot(precip_effects, aes(x = estimate, y = species, 
                                          color = species)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), 
                 height = 0.2, size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = species_colors) +
  labs(x = "Coefficient (Erreur standard)", 
       y = NULL,
       title = "Précipitations (30 jours avant)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        axis.text = element_text(size = 10),  # Taille des axes
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 12, face = "bold"),  # Titre principal plus grand et centré
        axis.title = element_text(size = 12))  # Titres des axes plus visibles

precip_plot

ggsave("Abundance_per_precipitations.tiff", width = 120, height = 90, units = "mm",dpi=600, compression = "lzw")

combined_plot <- temp_plot + precip_plot + 
  plot_layout(nrow = 1) +
  plot_annotation(tag_levels = 'A')

combined_plot

ggsave("Abundance_temp_prec_no_bold_axis_francais.tiff", width = 220, height = 80, units = "mm",dpi=600, compression = "lzw")


# Leafhopper abundance plot by temperature ----

# Generate a temperature grid with fixed covariates
new_temp <- data.frame(
  T2M_L30 = seq(min(defa$T2M_L30, na.rm = TRUE), max(defa$T2M_L30, na.rm = TRUE), length.out = 100),
  Prec_L30 = mean(defa$Prec_L30, na.rm = TRUE),
  Day = median(defa$Day, na.rm = TRUE),
  DayF = as.factor(round(median(defa$Day, na.rm = TRUE))),
  NumTraps = 10,
  duration = 7,
  Farm = NA,
  Year = NA
)

# Add offset
new_temp$logoffset <- log(new_temp$NumTraps * new_temp$duration)

# Predictions for each species
pred_efa_temp <- predict(mod4.efaw, newdata = new_temp, type = "response", se.fit = TRUE)
pred_mqa_temp <- predict(mod4.mqaw, newdata = new_temp, type = "response", se.fit = TRUE)
pred_loc_temp <- predict(mod4.locw, newdata = new_temp, type = "response", se.fit = TRUE)

# Data frames for each species
efa_temp_df <- new_temp %>%
  mutate(species = "Empoasca fabae",
         fit = pred_efa_temp$fit,
         se = pred_efa_temp$se.fit,
         lower = pmax(fit - 1.96 * se, 0),
         upper = fit + 1.96 * se)

mqa_temp_df <- new_temp %>%
  mutate(species = "Macrosteles quadrilineatus",
         fit = pred_mqa_temp$fit,
         se = pred_mqa_temp$se.fit,
         lower = pmax(fit - 1.96 * se, 0),
         upper = fit + 1.96 * se)

loc_temp_df <- new_temp %>%
  mutate(species = "Local species",
         fit = pred_loc_temp$fit,
         se = pred_loc_temp$se.fit,
         lower = pmax(fit - 1.96 * se, 0),
         upper = fit + 1.96 * se)

# Define species order
ordre_especes <- c("Empoasca fabae", "Macrosteles quadrilineatus", "Local species")

# Apply factor levels in the desired order
efa_temp_df$species <- factor(efa_temp_df$species, levels = ordre_especes)
mqa_temp_df$species <- factor(mqa_temp_df$species, levels = ordre_especes)
loc_temp_df$species <- factor(loc_temp_df$species, levels = ordre_especes)

# Combine predictions
pred_temp_all <- bind_rows(efa_temp_df, mqa_temp_df, loc_temp_df)

# Colors
cores_especes <- c("Empoasca fabae" = "#7fbf5b",
                   "Macrosteles quadrilineatus" = "#e27a1c",
                   "Local species" = "#1f74b4")

# Temperature prediction plot
fig_temp <- ggplot(pred_temp_all, aes(x = T2M_L30, y = fit, color = species, fill = species)) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
  scale_color_manual(values = cores_especes) +
  scale_fill_manual(values = cores_especes) +
  labs(x = "Average temperature (°C)", y = "Abundance per farm") +
  theme_bw(base_size = 12) +
  theme(legend.title = element_blank(),
        legend.position = "top",
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        panel.grid.major = element_blank(),  # remove major grid lines
        panel.grid.minor = element_blank())  # remove minor grid lines

fig_temp

ggsave("Abundance_by_temperature2.tiff", width = 160, height = 100, units = "mm", dpi = 600, compression = "lzw")

# Leafhopper abundance plot by precipitation ----

# Generate a precipitation grid with other covariates fixed
new_prec <- data.frame(
  T2M_L30 = mean(defa$T2M_L30, na.rm = TRUE),
  Prec_L30 = seq(min(defa$Prec_L30, na.rm = TRUE), max(defa$Prec_L30, na.rm = TRUE), length.out = 100),
  Day = median(defa$Day, na.rm = TRUE),
  DayF = as.factor(round(median(defa$Day, na.rm = TRUE))),
  NumTraps = 10,
  duration = 7,
  Farm = NA,
  Year = NA
)

# Add offset
new_prec$logoffset <- log(new_prec$NumTraps * new_prec$duration)

# Predictions for each species
pred_efa_prec <- predict(mod4.efaw, newdata = new_prec, type = "response", se.fit = TRUE)
pred_mqa_prec <- predict(mod4.mqaw, newdata = new_prec, type = "response", se.fit = TRUE)
pred_loc_prec <- predict(mod4.locw, newdata = new_prec, type = "response", se.fit = TRUE)

# Data frames for each species
efa_prec_df <- new_prec %>%
  mutate(species = "Empoasca fabae",
         fit = pred_efa_prec$fit,
         se = pred_efa_prec$se.fit,
         lower = pmax(fit - 1.96 * se, 0),
         upper = fit + 1.96 * se)

mqa_prec_df <- new_prec %>%
  mutate(species = "Macrosteles quadrilineatus",
         fit = pred_mqa_prec$fit,
         se = pred_mqa_prec$se.fit,
         lower = pmax(fit - 1.96 * se, 0),
         upper = fit + 1.96 * se)

loc_prec_df <- new_prec %>%
  mutate(species = "Local species",
         fit = pred_loc_prec$fit,
         se = pred_loc_prec$se.fit,
         lower = pmax(fit - 1.96 * se, 0),
         upper = fit + 1.96 * se)

efa_prec_df$species <- factor(efa_prec_df$species, levels = ordre_especes)
mqa_prec_df$species <- factor(mqa_prec_df$species, levels = ordre_especes)
loc_prec_df$species <- factor(loc_prec_df$species, levels = ordre_especes)

# Combine predictions
pred_prec_all <- bind_rows(efa_prec_df, mqa_prec_df, loc_prec_df)

# Colors
cores_especes <- c("Empoasca fabae" = "#7fbf5b",
                   "Macrosteles quadrilineatus" = "#e27a1c",
                   "Local species" = "#1f74b4")

# Precipitation prediction plot
fig_prec <- ggplot(pred_prec_all, aes(x = Prec_L30, y = fit, color = species, fill = species)) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
  scale_color_manual(values = cores_especes) +
  scale_fill_manual(values = cores_especes) +
  labs(x = "Accumulated precipitations (mm)", y = "Abundance per farm") +
  theme_classic(base_size = 12) +
  theme(legend.title = element_blank(),
        legend.position = "top",
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        panel.grid.major = element_blank(),  # remove major grid lines
        panel.grid.minor = element_blank())  # remove minor grid lines

fig_prec

ggsave("Abundance_by_precipitations.tiff", width = 160, height = 100, units = "mm", dpi = 600, compression = "lzw")