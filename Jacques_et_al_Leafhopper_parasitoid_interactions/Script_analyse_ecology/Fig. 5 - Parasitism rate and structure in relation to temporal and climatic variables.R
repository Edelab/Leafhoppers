rm(list = ls())

# library
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
library(multcomp)
library(multcompView)
library(emmeans)

# Parasitism occurrence----

# Empoasca----
dataefa<-read_excel("data_paraefa.xlsx", 2)

# First include the number of traps
dataefa$NumTraps <- 5

# Remove the 0 parasitism data
dataefa_clean <- dataefa %>%  
  filter(  
    !is.na(Efap),          # Remove NA in Efap  
    Efap != 0,              # Remove 0 in Efap  
    !is.na(Efarate),        # Remove NA in Efarate  
    Efarate != "0"          # Remove "0" (string) in Efarate  
  )

# Locals----
datalocal<-read_excel("data_paralocal.xlsx", 2)

# First include the number of traps
datalocal$NumTraps <- 5

# Remove the 0 parasitism data
datalocal_clean <- datalocal %>%  
  filter(  
    !is.na(Localspp),          # Remove NA in Localspp  
    Localspp != 0,              # Remove 0 in Localspp  
    !is.na(Localsprate),        # Remove NA in Localsprate  
    Localsprate != "0"          # Remove "0" (string) in Localsprate  
  )

# Macrosteles----
datamqa<-read_excel("data_paramqa.xlsx", 2)

# First include the number of traps
datamqa$NumTraps <- 5

# Remove the 0 parasitism data
datamqa_clean <- datamqa %>%  
  filter(  
    !is.na(Mqap),          # Remove NA in Mqap  
    Mqap != 0,              # Remove 0 in Mqap  
    !is.na(Mqarate),        # Remove NA in Mqarate  
    Mqarate != "0"          # Remove "0" (string) in Mqarate  
  )

# The rates were too low to analyze the seasonal effect
# so we can model how the probability of parasitism would occur 
# over the season transforming the data in a binary 0 and 1

# Total----
datatotal<-read_excel("data_paratotal.xlsx", 2)

str(datatotal)

# First include the number of traps
datatotal$NumTraps <- 5


datatotal2 <- datatotal %>%
  mutate(
    Parasitism_present = case_when(
      Totalp > 0 | (Totalrate != "0" & !is.na(Totalrate)) ~ 1,  # Parasitism detected
      TRUE ~ 0  # No parasitism (0 or NA)
    )
  )

datatot3 <- datatotal2 %>%
  group_by(Farm, Year, Day) %>%
  summarize(
    # Binary parasitism (1 if ANY trap had parasitism)
    Parasitized = as.integer(any(Totalp > 0 | (Totalrate != "0" & !is.na(Totalrate)))),
    # Total trap effort (NumTraps × duration)
    TotalTraps = sum(NumTraps),  # Sum traps per Farm/Day (5 or 10)
    TotalDuration = mean(duration),  # Mean duration (if variable, else sum)
    # Climate variables (averaged)
    across(c(T2M_L0, T2M_L7, T2M_L15, T2M_L30, 
             Prec_L0, Prec_L7, Prec_L15, Prec_L30,
             WS2M_L0, WS2M_L7, WS2M_L15, WS2M_L30), 
           mean, na.rm = TRUE),
    # Coordinates (constant per Farm)
    lat = first(lat),
    long = first(long),
    .groups = "drop"
  ) %>%
  # Add offset for trap effort (log-transform)
  mutate(
    Offset = log(TotalTraps * TotalDuration)
  )

datatot3$Year<-as.factor(datatot3$Year)
datatot3$Farm<-as.factor(datatot3$Farm)

# Models

m1p <-glmmTMB(
  Parasitized ~ Day + T2M_L0 + Prec_L0 + 
    offset(log(TotalTraps * TotalDuration)) + (1 | Year) + (1 | Farm),  
  data = datatot3,
  family = binomial,
  control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))
)

m2p <- glmmTMB(
  Parasitized ~ Day + T2M_L7 + Prec_L7 + 
    offset(log(TotalTraps * TotalDuration)) + (1 | Year) + (1 | Farm),  
  data = datatot3,
  family = binomial,
  control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))
)

m3p <- glmmTMB(
  Parasitized ~ Day + T2M_L15 + Prec_L15 + 
    offset(log(TotalTraps * TotalDuration)) + (1 | Year) + (1 | Farm),  
  data = datatot3,
  family = binomial,
  control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))
)

m4p <- glmmTMB(
  Parasitized ~ Day + T2M_L30 + Prec_L30 + 
    offset(log(TotalTraps * TotalDuration)) + (1 | Year) + (1 | Farm),  
  data = datatot3,
  family = binomial,
  control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))
)

summary(m4p)

r2_nakagawa(m4p)

AIC(m1p, m2p, m3p, m4p)

AICc(m1p, m2p, m3p, m4p)

# Model 4 is better

# Residuals
sim_resids <- simulateResiduals(m4p, plot = FALSE)
testDispersion(sim_resids)
testZeroInflation(sim_resids)
testUniformity(sim_resids)
check_collinearity(m4p)

# Spatial autocorrelation
location <- paste0(datatot3$long, "_", datatot3$lat)
res_grouped <- recalculateResiduals(sim_resids, group = location)
testSpatialAutocorrelation(res_grouped, x = unique(datatot3$long), y = unique(datatot3$lat))

# Temporal autocorrelation
res_agg <- recalculateResiduals(sim_resids, group = datatot3$Day)
testTemporalAutocorrelation(res_agg, time = unique(datatot3$Day))

summary(m4p)
r2(m4p)

# Day

# Generate predictions for Day (holding T2M_L30 at its mean)

# The model estimates how the probability of at least one trap detecting parasitism 
# changes over the season, taking into account climate and sampling effort.
# The probability to detect parasitism is expressed by julian day and by farm, 
# meaning that if the probability is 0,6 in June (x julian day), it represent that
# you have 60% of chance of having parasitism occurence (at least one parasitized
# leafhoppers) in the 10 traps that stayed in the field 7 days, at the mean temperature.
# The fact that we set 10 traps means that we predict the probability approximately
# by farm.

pred_day <- ggpredict(
  m4p,
  terms = "Day [all]",  # Vary Day across all observed values
  condition = c(
    TotalTraps = 10,     # Fix traps to 10
    TotalDuration = 7,   # Fix duration to 7 days
    T2M_L30 = mean(datatot3$T2M_L30, na.rm = TRUE)  # Hold temperature at mean
  ),
  bias_correction = TRUE  # Account for random effects
)

# Prediction plot for day
figpara <- ggplot(pred_day, aes(x, predicted)) +
  geom_line(linewidth = 2, color = "#7A7A7A") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3, fill = "#7A7A7A") +
  labs(y = "Probability of parasitism" ) +
  scale_y_continuous(
    limits = c(0.2, 0.9),  # Set y-axis range
    breaks = seq(0.2, 0.8, by = 0.2)  # Set intervals
  ) +
  scale_x_continuous(
    name = "Month",
    breaks = c(121, 152, 182, 213, 244, 274, 305),
    labels = c("May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov")) +
  theme_classic(base_size = 12) +
  theme(
    legend.title = element_blank(),
    axis.text = element_text(size = 10),  # Axis text size
    axis.title = element_text(size = 12, face = "bold"),  # More prominent axis titles
    legend.position = "right",
    text = element_text(size = 12),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()   # Remove minor grid lines
  )

figpara

ggsave("Parasitism_month.tiff", width = 160, height = 90, units = "mm",dpi=600, compression = "lzw")

# Generate predictions for T2M_L30
pred_temp <- ggpredict(
  m4p,
  terms = "T2M_L30 [all]",  # Vary temperature
  condition = c(
    TotalTraps = 10,        # Fix traps to 10
    TotalDuration = 7,      # Fix duration to 7 days
    Day = median(datatot3$Day, na.rm = TRUE),  # Hold Day at median
    Prec_L30 = mean(datatot3$Prec_L30, na.rm = TRUE)  # Hold precipitation at mean
  ),
  bias_correction = TRUE
)

# Prediction plot for temperature
figtemp <- ggplot(pred_temp, aes(x, predicted)) +
  geom_line(linewidth = 2, color = "#7A7A7A") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3, fill = "#7A7A7A") +
  scale_y_continuous(
    limits = c(0.0, 0.9),  # Set y-axis range
    breaks = seq(0.2, 0.8, by = 0.2)  # Set intervals
  ) +
  labs(x = "Temperature (°C)", y = "Probability of parasitism") +
  theme_bw(base_size = 12) +
  theme(
    legend.title = element_blank(),
    axis.text = element_text(size = 10),  # Axis text size
    axis.title = element_text(size = 12, face = "bold"),  # More prominent axis titles
    legend.position = "right",
    text = element_text(size = 12),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()   # Remove minor grid lines
  )

figtemp

ggsave("Parasitism_temp2.tiff", width = 160, height = 100, units = "mm",dpi=600, compression = "lzw")
  
# Generate predictions for Prec_L30
pred_prec <- ggpredict(
  m4p,
  terms = "Prec_L30 [all]",  # Vary precipitation
  condition = c(
    TotalTraps = 10,
    TotalDuration = 7,
    Day = median(datatot3$Day, na.rm = TRUE),
    T2M_L30 = mean(datatot3$T2M_L30, na.rm = TRUE)  # Also hold temperature constant
  ),
  bias_correction = TRUE
)

# Prediction plot for precipitation
figprec <- ggplot(pred_prec, aes(x, predicted)) +
  geom_line(linewidth = 2, color = "#7A7A7A") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3, fill = "#7A7A7A") +
  scale_y_continuous(
    limits = c(0.0, 0.9),
    breaks = seq(0.2, 0.8, by = 0.2)
  ) +
  labs(x = "Accumulated precipitations (mm)", y = "Probability of parasitism") +
  theme_bw(base_size = 12) +
  theme(
    legend.title = element_blank(),
    axis.text = element_text(size = 10),  # Axis text size
    axis.title = element_text(size = 12, face = "bold"),  # More prominent axis titles
    legend.position = "right",
    text = element_text(size = 12),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()   # Remove minor grid lines
  )


figprec


ggsave("Parasitism_prec2.tiff", width = 160, height = 100, units = "mm",dpi=600, compression = "lzw")

#Parasitism rate FINAL----

#Model by day including all the zeros rate----
datat<-read_excel("data_paratotal2.xlsx", 3)
datat$Year<-as.factor(datat$Year)
datat$Farm<-as.factor(datat$Farm)
datat$Prate<-as.numeric(datat$Prate)
datat$Totalab<-as.numeric(datat$Totalab) 



#candidate models

mod1.p2 <- glmmTMB(Prate ~ ns(Month, df=2) +
                     offset(log(Traps * Duration)) +
                     (1 | Farm) + (1 | Year),
                   ziformula = ~ 1, 
                   weights = Totalab,
                   family = betabinomial, data = datat)

mod2.p2 <- glmmTMB(Prate ~ ns(Month, df=3)+
                     offset(log(Traps * Duration)) +
                     (1 | Farm) + (1 | Year),
                   ziformula = ~ 1, 
                   weights = Totalab,
                   family = betabinomial, data = datat)

mod3.p2 <- glmmTMB(Prate ~  ns(Month, df=4) +
                     offset(log(Traps * Duration)) +
                     (1 | Farm) + (1 | Year),
                   ziformula = ~ 1, 
                   weights = Totalab,
                   family = betabinomial, data = datat)


compare_performance(mod1.p2, mod2.p2, mod3.p2, rank=TRUE)

summary(mod3.p2)

#residuals
sim_resids1 <- simulateResiduals(mod3.p2, plot = FALSE)
testDispersion(sim_resids1)
testUniformity(sim_resids1)
testZeroInflation(sim_resids1)

#spatial autocorrelation
location <- paste0(datat$long, "_", datat$lat)
res_grouped <- recalculateResiduals(sim_resids1, group = location)
testSpatialAutocorrelation(res_grouped, x = unique(datat$long), y = unique(datat$lat))

#temporal
res_agg <- recalculateResiduals(sim_resids1, group = datat$Day)
testTemporalAutocorrelation(res_agg, time = unique(datat$Day))

#plot
#esquisser()

new_data <- data.frame(Month = seq(5, 10, length.out = 100)) 

if (exists("datat")) { # Check if your original data is loaded
  new_data$Traps <- mean(datat$Traps, na.rm = TRUE)
  new_data$Duration <- mean(datat$Duration, na.rm = TRUE)
  new_data$Totalab <- mean(datat$Totalab, na.rm = TRUE) # Add Totalab here
} else {
  # If 'datat' is not available, you must set these manually based on your data.
  # These are placeholders; adjust them to sensible values from your data.
  new_data$Traps <- 10        # Example mean traps
  new_data$Duration <- 7     # Example mean duration
  new_data$Totalab <- 100    # Example mean Totalab (assuming it's a count/total)
  message("Warning: 'datat' not found. Using example values for Traps, Duration, and Totalab in predictions.")
}

new_data$Farm <- factor("dummy_farm") # Needs to be a factor with a placeholder level
new_data$Year <- factor("dummy_year") 

predictions_with_se <- predict(mod3.p2, newdata = new_data, type = "response", se.fit = TRUE, re.form = NA)
new_data$Predicted_Prate <- predictions_with_se$fit  # Extrait les valeurs numériques
new_data$SE <- predictions_with_se$se.fit

new_data$LowerCI <- new_data$Predicted_Prate - 1.96 * new_data$SE
new_data$UpperCI <- new_data$Predicted_Prate + 1.96 * new_data$SE

new_data$LowerCI <- pmax(0, new_data$LowerCI)
new_data$UpperCI <- pmin(1, new_data$UpperCI)


fignew<-ggplot(data = datat, aes(x = Month, y = Prate)) +
  geom_point(alpha = 1, position = position_jitter(width = 0.1), 
             size = 2, shape = 21, fill = "white",color = "#7A7A7A",stroke = 0.6) +
  geom_ribbon(data = new_data, aes(y = Predicted_Prate, ymin = LowerCI, ymax = UpperCI),
              fill = "#7A7A7A", alpha = 0.3) +
  geom_line(data = new_data, aes(y = Predicted_Prate), color = "#7A7A7A", linewidth = 1) +
  labs(x = "Month", y = "Parasitism rate (%)") +
  scale_x_continuous(breaks = c(5, 6, 7, 8, 9, 10),
                     labels = c("May", "Jun", "Jul", "Aug", "Sept", "Oct")) +
  scale_y_continuous(limits = c(0, 0.40), labels = scales::percent) + # Adjusted max to 0.45 for clarity
  theme_bw(base_size = 12) +
  theme(
    legend.title = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "right",
    text = element_text(size = 12),
    panel.grid.major = element_blank(),  # supprime quadrillage principal
    panel.grid.minor = element_blank())

fignew

ggsave("Parasitism_rate_month_new2.tiff", width = 160, height = 100, units = "mm",dpi=600, compression = "lzw")


# Parasitism proportion by month and group (Empoasca, Macrosteles and Locals)----

# Load data
datap <- read_excel("ptotal.xlsx")

# Display data structure
str(datap)

# Replace month names with their abbreviations
datap$Month <- dplyr :: recode(datap$Month, 
                               "May" = "May", 
                               "June" = "Jun", 
                               "July" = "Jul", 
                               "August" = "Aug", 
                               "September" = "Sept", 
                               "October" = "Oct")

# Restore month levels in the desired order
month_levels <- c("May", "Jun", "Jul", "Aug", "Sept", "Oct")
datap$Month <- factor(datap$Month, levels = month_levels)

# Proportion by group and month (overall)
# For the first plot
figprop <- ggplot(datap, aes(x = Month, y = Prop, fill = factor(Group, levels = c("Efa", "Mqa", "Localsp")))) +
  geom_col(position = "stack", width = 0.5, color='black') +
  geom_text(
    data = datap %>% group_by(Month) %>% summarise(Total = sum(N)),
    aes(x = Month, y = 1.05, label = paste("n =", Total), fill = NULL),
    vjust = 0, color = "black", size = 3  # Adjust label size
  ) +
  scale_fill_manual(
    #name = "Species Group",
    values = c("Efa" = "#7fbf5b", "Localsp" = "#1f74b4", "Mqa" = "#e27a1c"),
    labels = c("Empoasca fabae", "Macrosteles quadrilineatus", "Local species")
  ) +
  labs(x = "Month", y = "Proportion of parasitized leafhoppers") +
  theme_classic(base_size = 12) +
  theme(legend.title = element_blank(),
        axis.text = element_text(size = 10),  # Axis text size
        axis.title = element_text(size = 12, face = "bold"),  # More visible axis titles
        legend.position = "top",
        text = element_text(size = 12),
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank())   # Remove minor grid lines

figprop

figprop + coord_cartesian(ylim = c(0.97, 1))

ggsave("Parasitism_RAW_rate_months.tiff", width = 160, height = 100, units = "mm",dpi=600, compression = "lzw")

ggsave("Parasitism_RAW_rate_month_zoom97.2.tiff", width = 160, height = 100, units = "mm",dpi=600, compression = "lzw")

