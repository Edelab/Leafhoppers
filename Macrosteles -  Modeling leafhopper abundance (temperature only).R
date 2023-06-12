# ********************************************************************
# Macrosteles -  Modeling leafhopper abundance (temperature only) ####
# ********************************************************************
# Co-written by  Jeanne Durivage [1] et Gaétan Daigle [2], March 13, 2023
# [1] : Laboratoire d'Edel Pérez-Lopez, Département de phytologie, Faculté des sciences de l'Agriculture et de l'Alimentation, Université Laval
# [2] : Consultant spécialisé en évaluation ou en statistique, Département de mathématiques et statistique, Faculté de Sciences et de Génie


rm(list=ls())
setwd("/Users/Jeanne/Documents/EdeLab/Donnees_et_graphiques_d'abondance")
load("populations") 
library(dplyr); library(lubridate)

## This imports the "pop" df. It includes weekly catches broken down by species. 
## For modeling purposes, we're interested in leafhoppers in general, so the species doesn't matter.
## The first block of code is used to aggregate the catch data.
## It also removes the numerous columns that are unnecessary for this modeling.


df.Mq <- pop %>% 
  filter(genus == "Macrosteles" & specie == "quadrilineatus") %>% 
  aggregate(
    cbind(nb, weighted) ~ region + site + pose + withdrawal + traps + mean.temp +
      base10 + base8 + base5 + base2 + mean.date,
    FUN = sum) %>% 
  mutate(exposition.length = as.numeric(withdrawal - pose),
         day = yday(mean.date),
         region = as.factor(region),
         site = as.factor(site) )  %>% 
  dplyr::select( -c(pose, withdrawal, mean.date)) 


names(df.Mq)

# The dataset contains the following columns:
##nb: Number of leafhoppers caught during a survey. ** This is the response variable we want to model. **
##traps: Number of sticky traps used for the survey. This variable is used to "weight" nb, which is the number of captures.
##exposure.length: Number of days of exposure of sticky traps for the survey. Also used to "weight" nb.
## site: Site where the survey was carried out. We're not really interested in knowing the extent of its impact for this model.
## region: The 6 sites have been grouped into 4 regions. Again, we're not really interested in the extent of its impact for this model.
## day: number of days since January 1st. We want to try to explain "nb" with this variable.
## temperature, expressed in several ways. We don't know which of these indices is best for modeling leafhopper abundance of leafhoppers. 
#The code will therefore have to select the model with the best BIC to choose which index to keep.
#So we want to know whether temperature influences the abundance of leafhoppers caught, and which index is best.
#Here are the "candidate" indices
#mean.temp: average temperature (°C) at the time of the survey.
#base2: degree-days base 2°C accumulated while traps were exposed.
#base5: degree-days base 5°C...
#base8: degree-days base 8°C...
#base10: degree-days base 10°C...
## weighted: number of captures weighted by number of traps and days of exposure. Used only to illustrate results (last code block).

# *****************************************
# Finding the best model for abundance ####
# *****************************************
# Creation of different models (linear and 1st to 3rd degree polynomial)
library(lme4); library(cAIC4)

# Strategy: Each degree (0 to 3) is tested for each candidate temperature index.

## MEAN TEMPERATURE ##
Mq.mt.0 <- glmer.nb(nb ~ 1 + (1 | site) + offset(log(exposition.length * traps)), 
                 control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE), 
                 data = df.Mq)
Mq.mt.1 <- glmer.nb(nb ~ poly(day, 1) + poly(mean.temp, 1) + (1 | site) +
                   offset(log(exposition.length * traps)), 
                 control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                 data = df.Mq)
Mq.mt.2 <- glmer.nb(nb ~ poly(day, 2) + poly(mean.temp, 2) + (1 | site) +
                   offset(log(exposition.length * traps)), 
                 control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                 data = df.Mq)
Mq.mt.3 <- glmer.nb(nb ~ poly(day, 3) + poly(mean.temp, 3) + (1 | site) +
                   offset(log(exposition.length * traps)), 
                 control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                 data = df.Mq)

## base2 ##
Mq.b2.0 <- glmer.nb(nb ~ 1 + (1 | site) + offset(log(exposition.length * traps)), 
                 control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE), 
                 data = df.Mq)
Mq.b2.1 <- glmer.nb(nb ~ poly(day, 1) + poly(base2, 1) + (1 | site) +
                   offset(log(exposition.length * traps)), 
                 control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                 data = df.Mq)
Mq.b2.2 <- glmer.nb(nb ~ poly(day, 2) + poly(base2, 2) + (1 | site) +
                   offset(log(exposition.length * traps)), 
                 control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                 data = df.Mq)
Mq.b2.3 <- glmer.nb(nb ~ poly(day, 3) + poly(base2, 3) + (1 | site) +
                   offset(log(exposition.length * traps)), 
                 control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                 data = df.Mq)

## base5 ##
Mq.b5.0 <- glmer.nb(nb ~ 1 + (1 | site) + offset(log(exposition.length * traps)), 
                 control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE), 
                 data = df.Mq)
Mq.b5.1 <- glmer.nb(nb ~ poly(day, 1) + poly(base5, 1) + (1 | site) +
                   offset(log(exposition.length * traps)), 
                 control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                 data = df.Mq)
Mq.b5.2 <- glmer.nb(nb ~ poly(day, 2) + poly(base5, 2) + (1 | site) +
                   offset(log(exposition.length * traps)), 
                 control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                 data = df.Mq)
Mq.b5.3 <- glmer.nb(nb ~ poly(day, 3) + poly(base5, 3) + (1 | site) +
                   offset(log(exposition.length * traps)), 
                 control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                 data = df.Mq)

## base8 ##
Mq.b8.0 <- glmer.nb(nb ~ 1 + (1 | site) + offset(log(exposition.length * traps)), 
                 control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE), 
                 data = df.Mq)
Mq.b8.1 <- glmer.nb(nb ~ poly(day, 1) + poly(base8, 1) + (1 | site) +
                   offset(log(exposition.length * traps)), 
                 control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                 data = df.Mq)
Mq.b8.2 <- glmer.nb(nb ~ poly(day, 2) + poly(base8, 2) + (1 | site) +
                   offset(log(exposition.length * traps)), 
                 control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                 data = df.Mq)
Mq.b8.3 <- glmer.nb(nb ~ poly(day, 3) + poly(base8, 3) + (1 | site) +
                   offset(log(exposition.length * traps)), 
                 control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                 data = df.Mq)

## base10 ##
Mq.b10.0 <- glmer.nb(nb ~ 1 + (1 | site) + offset(log(exposition.length * traps)), 
                  control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE), 
                  data = df.Mq)
Mq.b10.1 <- glmer.nb(nb ~ poly(day, 1) + poly(base10, 1) + (1 | site) +
                    offset(log(exposition.length * traps)), 
                  control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                  data = df.Mq)
Mq.b10.2 <- glmer.nb(nb ~ poly(day, 2) + poly(base10, 2) + (1 | site) +
                    offset(log(exposition.length * traps)), 
                  control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                  data = df.Mq)
Mq.b10.3 <- glmer.nb(nb ~ poly(day, 3) + poly(base10, 3) + (1 | site) +
                    offset(log(exposition.length * traps)), 
                  control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                  data = df.Mq)



# Pre-selection of the best model
library(AICcmodavg)
Mq.modeles <- list(Mq.mt.0, Mq.mt.1, Mq.mt.2, Mq.mt.3, Mq.b10.0, Mq.b10.1, Mq.b10.2, Mq.b10.3,
                Mq.b8.0, Mq.b8.1, Mq.b8.2, Mq.b8.3, Mq.b5.0, Mq.b5.1, Mq.b5.2, Mq.b5.3, 
                Mq.b2.0, Mq.b2.1, Mq.b2.2, Mq.b2.3)
Mq.noms <- c("Mq.mt.0", "Mq.mt.1", "Mq.mt.2", "Mq.mt.3", "Mq.b10.0", "Mq.b10.1", "Mq.b10.2", "Mq.b10.3",
          "Mq.b8.0", "Mq.b8.1", "Mq.b8.2", "Mq.b8.3", "Mq.b5.0", "Mq.b5.1", "Mq.b5.2", "Mq.b5.3", 
          "Mq.b2.0", "Mq.b2.1", "Mq.b2.2", "Mq.b2.3")


Mq.choix.BIC <- bictab(cand.set = Mq.modeles, modnames = Mq.noms)
head(Mq.choix.BIC)
# Average temperature, base 2 and base 5 all have the same value. This can be explained by the fact that when the species immigrates to Quebec, the weather is already very warm.
# Confirmation of this theory:
    range(df.Mq$mean.temp)
    range(df.Mq$base10)
    # The most basic temperature is very high indeed. mean.temp, base2 and base5 are equivalent.

# Observation of the best pre-selected model
summary(Mq.mt.1)

# From the observations in the last line of code, it seems that the model is less good for temperature.
# From the chosen index, try combining several degrees to find the best combination.

Mq.mt.12 <- glmer.nb(nb ~ poly(day, 1) + poly(mean.temp, 2) + (1 | site) +
                    offset(log(exposition.length * traps)), 
                  control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                  data = df.Mq)
Mq.mt.13 <- glmer.nb(nb ~ poly(day, 1) + poly(mean.temp, 3) + (1 | site) +
                    offset(log(exposition.length * traps)), 
                  control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                  data = df.Mq)
Mq.mt.21 <- glmer.nb(nb ~ poly(day, 2) + poly(mean.temp, 1) + (1 | site) +
                    offset(log(exposition.length * traps)), 
                  control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                  data = df.Mq)
Mq.mt.23 <- glmer.nb(nb ~ poly(day, 2) + poly(mean.temp, 3) + (1 | site) +
                    offset(log(exposition.length * traps)), 
                  control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                  data = df.Mq)  
Mq.mt.31 <- glmer.nb(nb ~ poly(day, 3) + poly(mean.temp, 1) + (1 | site) +
                    offset(log(exposition.length * traps)), 
                  control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                  data = df.Mq)
Mq.mt.32 <- glmer.nb(nb ~ poly(day, 3) + poly(mean.temp, 2) + (1 | site) +
                    offset(log(exposition.length * traps)), 
                  control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                  data = df.Mq)


Mq.b8.12 <- glmer.nb(nb ~ poly(day, 1) + poly(base8, 2) + (1 | site) +
                    offset(log(exposition.length * traps)), 
                  control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                  data = df.Mq)
Mq.b8.13 <- glmer.nb(nb ~ poly(day, 1) + poly(base8, 3) + (1 | site) +
                    offset(log(exposition.length * traps)), 
                  control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                  data = df.Mq)
Mq.b8.21 <- glmer.nb(nb ~ poly(day, 2) + poly(base8, 1) + (1 | site) +
                    offset(log(exposition.length * traps)), 
                  control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                  data = df.Mq)
Mq.b8.23 <- glmer.nb(nb ~ poly(day, 2) + poly(base8, 3) + (1 | site) +
                    offset(log(exposition.length * traps)), 
                  control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                  data = df.Mq)  
Mq.b8.31 <- glmer.nb(nb ~ poly(day, 3) + poly(base8, 1) + (1 | site) +
                    offset(log(exposition.length * traps)), 
                  control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                  data = df.Mq)
Mq.b8.32 <- glmer.nb(nb ~ poly(day, 3) + poly(base8, 2) + (1 | site) +
                    offset(log(exposition.length * traps)), 
                  control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                  data = df.Mq)


Mq.b10.12 <- glmer.nb(nb ~ poly(day, 1) + poly(base10, 2) + (1 | site) +
                    offset(log(exposition.length * traps)), 
                  control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                  data = df.Mq)
Mq.b10.13 <- glmer.nb(nb ~ poly(day, 1) + poly(base10, 3) + (1 | site) +
                    offset(log(exposition.length * traps)), 
                  control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                  data = df.Mq)
Mq.b10.21 <- glmer.nb(nb ~ poly(day, 2) + poly(base10, 1) + (1 | site) +
                    offset(log(exposition.length * traps)), 
                  control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                  data = df.Mq)
Mq.b10.23 <- glmer.nb(nb ~ poly(day, 2) + poly(base10, 3) + (1 | site) +
                    offset(log(exposition.length * traps)), 
                  control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                  data = df.Mq)  
Mq.b10.31 <- glmer.nb(nb ~ poly(day, 3) + poly(base10, 1) + (1 | site) +
                    offset(log(exposition.length * traps)), 
                  control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                  data = df.Mq)
Mq.b10.32 <- glmer.nb(nb ~ poly(day, 3) + poly(base10, 2) + (1 | site) +
                    offset(log(exposition.length * traps)), 
                  control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                  data = df.Mq)

# Selection of the best model from the combinations created.
Mq.modeles2 <- list(Mq.mt.12, Mq.mt.13, Mq.mt.21, Mq.mt.23, Mq.mt.31, Mq.mt.32,
                    Mq.mt.2, Mq.mt.3, Mq.mt.1,
                    Mq.b10.12, Mq.b10.13, Mq.b10.21, Mq.b10.23, Mq.b10.31, Mq.b10.32,
                    Mq.b10.2, Mq.b10.3, Mq.b10.1, 
                    Mq.b8.12, Mq.b8.13, Mq.b8.21, Mq.b8.23, Mq.b8.31, Mq.b8.32,
                    Mq.b8.2, Mq.b8.3, Mq.b8.1)
Mq.noms2 <- c("Mq.mt.12", "Mq.mt.13", "Mq.mt.21", "Mq.mt.23", "Mq.mt.31", "Mq.mt.32", 
              "Mq.mt.1", "Mq.mt.2", "Mq.mt.3",
              "Mq.b8.12", "Mq.b8.13", "Mq.b8.21", "Mq.b8.23", "Mq.b8.31", "Mq.b8.32", 
               "Mq.b8.1", "Mq.b8.2", "Mq.b8.3",
               "Mq.b10.12", "Mq.b10.13", "Mq.b10.21", "Mq.b10.23", "Mq.b10.31", "Mq.b10.32", 
              "Mq.b10.1", "Mq.b10.2", "Mq.b10.3")

Mq.choix.BIC2 <- bictab(cand.set = Mq.modeles2, modnames = Mq.noms2)
head(Mq.choix.BIC2)
  # Model Mq.b10.12 and model Mq.mt.12 have a similar result.Validation :
    
summary(Mq.mt.12)
summary(Mq.b10.12)

  # base10 better models the "day" factor. Mq.mt.12 better models the "temperature" factor. The difference is slight.
  # I've decided to keep Mq.mt.12, as we're more interested in temperature modeling results in our study.

  
#Selected model
summary(Mq.mt.12)


# Cleannig environnement
rm(list=setdiff(ls(), c("Mq.mt.12", "df.Mq")))  






# ********************************************************
# Illustrate the model, with its confidence intervals ####
# ********************************************************
library(ggplot2)
setwd("/Users/Jeanne/Documents/EdeLab/Modelisation_de_l'abondance/Graphiques_de_modelisation")

# A) Depending on the temperature index selected

# Model prediction on NEWDATA
range(df.Mq$mean.temp)

newdata.temp <- expand.grid(
  mean.temp = seq( 11.5, 24.6, 0.1),
  day = quantile(df.Mq$day, c(0, 0.50, 1)),
  exposition.length = 1, 
  traps = 1)

newdata.temp$pred = predict(Mq.mt.12, newdata=newdata.temp, re.form=NA, type="response")

# Bootstrap CI (confidence intervals) :  
boot = bootMer(Mq.mt.12, 
               function(x) predict(x, newdata = newdata.temp, re.form = NA, type = 'response'), 
               nsim = 200, 
               use.u = FALSE, 
               type = 'parametric', 
               parallel = 'multicore', 
               ncpus = 4)

newdata.temp$lowerCI = apply(boot$t, 2, function(x) quantile(x, 0.025))
newdata.temp$upperCI = apply(boot$t, 2, function(x) quantile(x, 0.975))

# Modeling illustration, with weighted data
unique(newdata.temp$day) 

ggplot() +
  
  geom_line(data=newdata.temp[newdata.temp$day==140,], 
            aes(x=mean.temp, y=pred), color="cadetblue3", linewidth=1.3, lty=1 ) +
  geom_line(data=newdata.temp[newdata.temp$day==140,], 
            aes(x=mean.temp, y=lowerCI), color="cadetblue3", linewidth=1.3, lty=2) +
  geom_line(data=newdata.temp[newdata.temp$day==140,], 
            aes(x=mean.temp, y=upperCI), color="cadetblue3", linewidth=1.3, lty=2 ) +
  
  geom_line(data=newdata.temp[newdata.temp$day==202,], 
            aes(x=mean.temp, y=pred), color="gray30", linewidth=1.3, lty=1 ) +
  geom_line(data=newdata.temp[newdata.temp$day==202,], 
            aes(x=mean.temp, y=lowerCI), color="gray30", linewidth=1.3, lty=2) +
  geom_line(data=newdata.temp[newdata.temp$day==202,], 
            aes(x=mean.temp, y=upperCI), color="gray30", linewidth=1.3, lty=2 ) +
  
  geom_line(data=newdata.temp[newdata.temp$day==247,], 
            aes(x=mean.temp, y=pred), color="red", linewidth=1.3, lty=1 ) +
  geom_line(data=newdata.temp[newdata.temp$day==247,], 
            aes(x=mean.temp, y=lowerCI), color="red", linewidth=1.3, lty=2) +
  geom_line(data=newdata.temp[newdata.temp$day==247,], 
            aes(x=mean.temp, y=upperCI), color="red", linewidth=1.3, lty=2 ) +
  
  geom_point(data=df.Mq, aes(x=mean.temp, y=weighted), size = 1.75) +
  
  theme_bw() +
  theme(axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 18, face="italic"),
        axis.title = element_text(size = 16))    +
  labs(y = "Weighted captures", x = "Temperature (°C)", title = "Macrosteles quadrilineatus") +
  scale_x_continuous(breaks = seq(10,25, 2))
        

#B) Depending on the time of year
range(df.Mq$day)

newdata.day = expand.grid(
  day = seq(140, 247, 1), 
  mean.temp = quantile(df.Mq$mean.temp, c(0, 0.50, 1)),
  exposition.length = 1, 
  traps = 1) 

newdata.day$pred = predict(Mq.mt.12, newdata=newdata.day, re.form=NA, type="response")

# Bootstrap CI:  
boot = bootMer(Mq.mt.12, 
               function(x) predict(x, newdata = newdata.day, re.form = NA, type = 'response'), 
               nsim = 200, 
               use.u = FALSE, 
               type = 'parametric', 
               parallel = 'multicore', 
               ncpus = 4)

newdata.day$lowerCI = apply(boot$t, 2, function(x) quantile(x, 0.025)) # Générer l'intervalle de confiance inférieur
newdata.day$upperCI = apply(boot$t, 2, function(x) quantile(x, 0.975)) # Générer l'intervalle de confiance supérieur


# Modeling illustration, with weighted data
mean.temp = sort(unique(newdata.day$mean.temp)) 
round(mean.temp,2)

ggplot() +
  
  geom_line(data=newdata.day[newdata.day$mean.temp==mean.temp[1],], 
            aes(x=day, y=pred), color="cadetblue3", linewidth=1.3, lty=1 ) +
  geom_line(data=newdata.day[newdata.day$mean.temp==mean.temp[1],], 
            aes(x=day, y=lowerCI), color="cadetblue3", linewidth=1.3, lty=2) +
  geom_line(data=newdata.day[newdata.day$mean.temp==mean.temp[1],], 
            aes(x=day, y=upperCI), color="cadetblue3", linewidth=1.3, lty=2 ) +
  
  geom_line(data=newdata.day[newdata.day$mean.temp==mean.temp[2],], 
            aes(x=day, y=pred), color="gray30", linewidth=1.3, lty=1 ) +
  geom_line(data=newdata.day[newdata.day$mean.temp==mean.temp[2],], 
            aes(x=day, y=lowerCI), color="gray30", linewidth=1.3, lty=2) +
  geom_line(data=newdata.day[newdata.day$mean.temp==mean.temp[2],], 
            aes(x=day, y=upperCI), color="gray30", linewidth=1.3, lty=2 ) +
  
  geom_line(data=newdata.day[newdata.day$mean.temp==mean.temp[3],], 
            aes(x=day, y=pred), color="red", linewidth=1.3, lty=1 ) +
  geom_line(data=newdata.day[newdata.day$mean.temp==mean.temp[3],], 
            aes(x=day, y=lowerCI), color="red", linewidth=1.3, lty=2) +
  geom_line(data=newdata.day[newdata.day$mean.temp==mean.temp[3],], 
            aes(x=day, y=upperCI), color="red", linewidth=1.3, lty=2 ) +
  
  geom_point(data=df.Mq, aes(x= day, y=weighted), size = 1.75) +
  
  theme_bw() +
  theme(axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 18, face="italic"),
        axis.title = element_text(size = 16))    +
  labs(x = " ", y = "Weighted captures", title = "Macrosteles quadrilineatus") +
  scale_x_continuous(breaks=c(123, 153, 183, 214, 245),
                     labels=c("May", "June", "July", "August", "September")) 
  

