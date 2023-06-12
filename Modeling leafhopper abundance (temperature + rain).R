# *******************************************************
# Modeling leafhopper abundance (temperature + rain) ####
# *******************************************************
# Co-written by  Jeanne Durivage [1] et Gaétan Daigle [2], March 13, 2023
# [1] : Laboratoire d'Edel Pérez-Lopez, Département de phytologie, Faculté des sciences de l'Agriculture et de l'Alimentation, Université Laval
# [2] : Consultant spécialisé en évaluation ou en statistique, Département de mathématiques et statistique, Faculté de Sciences et de Génie

rm(list=ls())
setwd("/Users/Jeanne/Documents/EdeLab/Donnees_et_graphiques_d'abondance")
load("populations") 
  

## This imports the "pop" df. It includes weekly catches broken down by species. 
## For modeling purposes, we're interested in leafhoppers in general, so the species doesn't matter.
## The first block of code is used to aggregate the catch data.
## It also removes the numerous columns that are unnecessary for this modeling.


  library(dplyr); library(lubridate)

  df.pluie <- pop %>% 
    aggregate( cbind(nb, weighted) ~ region + site + pose + withdrawal + 
                                     traps + mean.rain + base10 + mean.date,
             FUN = sum) %>% 
     mutate(exposition.length = as.numeric(withdrawal - pose),
            day = yday(mean.date),
            region = as.factor(region),
            site = as.factor(site) )  %>% 
     dplyr::select( -c(pose, withdrawal, mean.date) ) 
  

names(df.pluie)
unique(df.pluie$site)

# Creation of different models (linear and 1st to 3rd degree polynomial)
  library(lme4); library(cAIC4)

# Here's the best model when you don't consider rain:
nbm.base <- glmer.nb(nb ~ poly(day, 2) + poly(base10, 1) + (1 | site) +
                   offset(log(exposition.length * traps)), 
                 control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                 data = df.pluie)

  # Strategy: each degree (0 to 3) is tested for rain + with or without rain
  ## MEAN TEMPERATURE ##  
  nbm.0 <- glmer.nb(nb ~ 1 + (1 | site) + offset(log(exposition.length * traps)), 
                   control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE), 
                  data = df.pluie)
  nbm.1 <- glmer.nb(nb ~ poly(day, 2) + poly(base10, 1) + poly(mean.rain, 1) +
                     (1 | site) +
                    offset(log(exposition.length * traps)), 
                  control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                  data = df.pluie)
  nbm.2 <- glmer.nb(nb ~ poly(day, 2) + poly(base10, 1) + poly(mean.rain, 2) +
                     (1 | site) +
                    offset(log(exposition.length * traps)), 
                  control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                  data = df.pluie)
  nbm.3 <- glmer.nb(nb ~ poly(day, 2) + poly(base10, 1) + poly(mean.rain, 3) +
                     (1 | site) +
                    offset(log(exposition.length * traps)), 
                  control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                  data = df.pluie)

# Best model choice
  library(AICcmodavg)
  modeles <- list(nbm.0, nbm.1, nbm.2, nbm.3, nbm.base)
  noms <- c("nbm.0", "nbm.1", "nbm.2", "nbm.3", "nbm.base")
  bictab(cand.set = modeles, modnames = noms)
    ## The model without precipitation ('base') is better according to the BIC criterion
  
    ## Let's take a look at the best model with precipitation, mt.2
    summary(nbm.1)
    ## Rain is not significant at all, even in the best model.
    ## The BIC criterion evaluates that the best model does not consider rain.

# ***********************************************************
# Visualization of rainfall data with weighted abundance ####
# ***********************************************************
  library(ggplot2)
    ggplot() +
      geom_point(data=df.pluie, aes(x = mean.rain, y = weighted, color = region)) +
      theme_bw() +
      theme(panel.grid = element_blank())
  
# ************************************************
# Exploring the rain-temperature relationship ####
# ************************************************
    ggplot() +
      geom_point(data=df.pluie, aes(x = mean.rain, y = base10, color = region)) +
      theme_bw() +
      theme(panel.grid = element_blank())
    
    
    cor(df.pluie$mean.rain, df.pluie$base10)
        # Low correlation    