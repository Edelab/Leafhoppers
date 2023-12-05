# *****************************************************
# Modeling leafhopper abundance (temperature only) ####
# *****************************************************
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

  df <- pop %>% 
      aggregate(
     cbind(nb, weighted) ~ region + site + pose + withdrawal + traps + mean.temp +
             base10 + base8 + base5 + base2 + mean.date,
             FUN = sum) %>% 
     mutate(exposition.length = as.numeric(withdrawal - pose),
            day = yday(mean.date),
            region = as.factor(region),
            site = as.factor(site) )  %>% 
     dplyr::select( -c(pose, withdrawal, mean.date)) 
  

names(df)
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
  library(lme4)
  library(cAIC4)

# Strategy: Each degree (0 to 3) is tested for each candidate temperature index.

## MEAN TEMPERATURE ##  
  mt.0 <- glmer.nb(nb ~ 1 + (1 | site) + offset(log(exposition.length * traps)), 
                  control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE), 
                  data = df)
  mt.1 <- glmer.nb(nb ~ poly(day, 1) + poly(mean.temp, 1) + (1 | site) +
                    offset(log(exposition.length * traps)), 
                  control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                  data = df)
  mt.2 <- glmer.nb(nb ~ poly(day, 2) + poly(mean.temp, 2) + (1 | site) +
                    offset(log(exposition.length * traps)), 
                  control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                  data = df)
  mt.3 <- glmer.nb(nb ~ poly(day, 3) + poly(mean.temp, 3) + (1 | site) +
                    offset(log(exposition.length * traps)), 
                  control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                  data = df)

  ## base2 ##
  b2.0 <- glmer.nb(nb ~ 1 + (1 | site) + offset(log(exposition.length * traps)), 
                   control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE), 
                   data = df)
  b2.1 <- glmer.nb(nb ~ poly(day, 1) + poly(base2, 1) + (1 | site) +
                     offset(log(exposition.length * traps)), 
                   control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                   data = df)
  b2.2 <- glmer.nb(nb ~ poly(day, 2) + poly(base2, 2) + (1 | site) +
                     offset(log(exposition.length * traps)), 
                   control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                   data = df)
  b2.3 <- glmer.nb(nb ~ poly(day, 3) + poly(base2, 3) + (1 | site) +
                     offset(log(exposition.length * traps)), 
                   control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                   data = df)
  
  ## base5 ##
  b5.0 <- glmer.nb(nb ~ 1 + (1 | site) + offset(log(exposition.length * traps)), 
                   control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE), 
                   data = df)
  b5.1 <- glmer.nb(nb ~ poly(day, 1) + poly(base5, 1) + (1 | site) +
                     offset(log(exposition.length * traps)), 
                   control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                   data = df)
  b5.2 <- glmer.nb(nb ~ poly(day, 2) + poly(base5, 2) + (1 | site) +
                     offset(log(exposition.length * traps)), 
                   control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                   data = df)
  b5.3 <- glmer.nb(nb ~ poly(day, 3) + poly(base5, 3) + (1 | site) +
                     offset(log(exposition.length * traps)), 
                   control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                   data = df)  
  
  ## base8 ##
  b8.0 <- glmer.nb(nb ~ 1 + (1 | site) + offset(log(exposition.length * traps)), 
                  control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE), 
                  data = df)
  b8.1 <- glmer.nb(nb ~ poly(day, 1) + poly(base8, 1) + (1 | site) +
                    offset(log(exposition.length * traps)), 
                  control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                  data = df)
  b8.2 <- glmer.nb(nb ~ poly(day, 2) + poly(base8, 2) + (1 | site) +
                    offset(log(exposition.length * traps)), 
                  control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                  data = df)
  b8.3 <- glmer.nb(nb ~ poly(day, 3) + poly(base8, 3) + (1 | site) +
                    offset(log(exposition.length * traps)), 
                  control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                  data = df)

  ## base10 ##
  b10.0 <- glmer.nb(nb ~ 1 + (1 | site) + offset(log(exposition.length * traps)), 
                    control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE), 
                    data = df)
  b10.1 <- glmer.nb(nb ~ poly(day, 1) + poly(base10, 1) + (1 | site) +
                      offset(log(exposition.length * traps)), 
                    control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                    data = df)
  b10.2 <- glmer.nb(nb ~ poly(day, 2) + poly(base10, 2) + (1 | site) +
                      offset(log(exposition.length * traps)), 
                    control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                    data = df)
  b10.3 <- glmer.nb(nb ~ poly(day, 3) + poly(base10, 3) + (1 | site) +
                      offset(log(exposition.length * traps)), 
                    control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                    data = df)
  
    
  
  # Best model choice
  library(AICcmodavg)
  modeles <- list(mt.0, mt.1, mt.2, mt.3, b10.0, b10.1, b10.2, b10.3,
                  b8.0, b8.1, b8.2, b8.3, b5.0, b5.1, b5.2, b5.3, 
                  b2.0, b2.1, b2.2, b2.3)
  noms <- c("mt.0", "mt.1", "mt.2", "mt.3", "b10.0", "b10.1", "b10.2", "b10.3",
            "b8.0", "b8.1", "b8.2", "b8.3", "b5.0", "b5.1", "b5.2", "b5.3", 
            "b2.0", "b2.1", "b2.2", "b2.3")
  
  choix.BIC <- bictab(cand.set = modeles, modnames = noms)

  head(choix.BIC)

# Observation of the pre-selected model
  summary(b10.2)

  # The 'base10' temperature index is the best, closely followed by 'base8'.
  # From the observations of the last line of code, it seems that the 2nd degree is less good for the temperature part of the model.
  # From the indices chosen, I'll try to combine several degrees to find the best combination. 
  
  b8.12 <- glmer.nb(nb ~ poly(day, 1) + poly(base8, 2) + (1 | site) +
                     offset(log(exposition.length * traps)), 
                   control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                   data = df)
  b8.13 <- glmer.nb(nb ~ poly(day, 1) + poly(base8, 3) + (1 | site) +
                     offset(log(exposition.length * traps)), 
                   control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                   data = df)
  b8.21 <- glmer.nb(nb ~ poly(day, 2) + poly(base8, 1) + (1 | site) +
                     offset(log(exposition.length * traps)), 
                   control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                   data = df)
  b8.23 <- glmer.nb(nb ~ poly(day, 2) + poly(base8, 3) + (1 | site) +
                     offset(log(exposition.length * traps)), 
                   control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                   data = df)  
  b8.31 <- glmer.nb(nb ~ poly(day, 3) + poly(base8, 1) + (1 | site) +
                     offset(log(exposition.length * traps)), 
                   control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                   data = df)
  b8.32 <- glmer.nb(nb ~ poly(day, 3) + poly(base8, 2) + (1 | site) +
                     offset(log(exposition.length * traps)), 
                   control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                   data = df)
  
  
  
  b10.12 <- glmer.nb(nb ~ poly(day, 1) + poly(base10, 2) + (1 | site) +
                      offset(log(exposition.length * traps)), 
                    control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                    data = df)
  b10.13 <- glmer.nb(nb ~ poly(day, 1) + poly(base10, 3) + (1 | site) +
                      offset(log(exposition.length * traps)), 
                    control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                    data = df)
  b10.21 <- glmer.nb(nb ~ poly(day, 2) + poly(base10, 1) + (1 | site) +
                      offset(log(exposition.length * traps)), 
                    control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                    data = df)
  b10.23 <- glmer.nb(nb ~ poly(day, 2) + poly(base10, 3) + (1 | site) +
                      offset(log(exposition.length * traps)), 
                    control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                    data = df)  
  b10.31 <- glmer.nb(nb ~ poly(day, 3) + poly(base10, 1) + (1 | site) +
                      offset(log(exposition.length * traps)), 
                    control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                    data = df)
  b10.32 <- glmer.nb(nb ~ poly(day, 3) + poly(base10, 2) + (1 | site) +
                      offset(log(exposition.length * traps)), 
                    control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                    data = df)
  
  # Selection of the best model from the combinations created.
  modeles2 <- list(b10.12, b10.13, b10.21, b10.23, b10.31, b10.32,
                   b10.2, b10.3, b10.1, 
                   b8.12, b8.13, b8.21, b8.23, b8.31, b8.32,
                   b8.2, b8.3, b8.1)
  noms2 <- c("b8.12", "b8.13", "b8.21", "b8.23", "b8.31", "b8.32", 
            "b8.1", "b8.2", "b8.3",
            "b10.12", "b10.13", "b10.21", "b10.23", "b10.31", "b10.32", 
            "b10.1", "b10.2", "b10.3")
  choix2.BIC <- bictab(cand.set = modeles2, modnames = noms2)
  head(choix2.BIC)
  # The b10.21 model is the best.
  
  # Visualize the final model
  summary(b10.21)

#The following lines were add in 05/12/2023

#R2 conditional and marginal
library(performance)
model_performance(b10.21)

#AGVIF
library(car)
vif(b10.21)

#model suitability DHARMa
library(DHARMa)
modsuit<- simulateResiduals(fittedModel = b10.21)
testResiduals(modsuit)


# Cleannig environnement
rm(list=setdiff(ls(), c("b10.21", "df")))  
  
  
  
# ********************************************************
# Illustrate the model, with its confidence intervals ####
# ********************************************************
    library(ggplot2)
    setwd("/Users/Jeanne/Documents/EdeLab/Modelisation_de_l'abondance/Graphiques_de_modelisation")
  
  # A) Depending on the temperature index selected
    range(df$base10)
    
    # Model prediction on NEWDATA
    newdata.temp <- expand.grid(
      base10 = seq( 0, 14.6, 0.1),
      day = quantile(df$day, c(0, 0.50, 1)),
      exposition.length = 1, 
      traps = 1)
    
    newdata.temp$pred <- predict(b10.21, newdata=newdata.temp, re.form=NA, type="response")
    
    # Bootstrap CI (confidence intervals) :  
    boot <- bootMer(b10.21, 
                   function(x) predict(x, newdata = newdata.temp, re.form = NA, type = 'response'), 
                   nsim = 200, 
                   use.u = FALSE, 
                   type = 'parametric', 
                   parallel = 'multicore', 
                   ncpus = 4)
    
    newdata.temp$lowerCI <- apply(boot$t, 2, function(x) quantile(x, 0.025))
    newdata.temp$upperCI <- apply(boot$t, 2, function(x) quantile(x, 0.975))
    
    # Modeling illustration, with weighted data
    unique(newdata.temp$day)
    
    ggplot() +
      
      geom_line(data=newdata.temp[newdata.temp$day==118,], 
                aes(x=base10, y=pred), color="cadetblue3", linewidth=1.3, lty=1 ) +
      geom_line(data=newdata.temp[newdata.temp$day==118,], 
                aes(x=base10, y=lowerCI), color="cadetblue3", linewidth=1.3, lty=2) +
      geom_line(data=newdata.temp[newdata.temp$day==118,], 
                aes(x=base10, y=upperCI), color="cadetblue3", linewidth=1.3, lty=2 ) +
      
      geom_line(data=newdata.temp[newdata.temp$day==186,], 
                aes(x=base10, y=pred), color="gray30", linewidth=1.3, lty=1 ) +
      geom_line(data=newdata.temp[newdata.temp$day==186,], 
                aes(x=base10, y=lowerCI), color="gray30", linewidth=1.3, lty=2) +
      geom_line(data=newdata.temp[newdata.temp$day==186,], 
                aes(x=base10, y=upperCI), color="gray30", linewidth=1.3, lty=2 ) +
      
      geom_line(data=newdata.temp[newdata.temp$day==247,], 
                aes(x=base10, y=pred), color="red", linewidth=1.3, lty=1 ) +
      geom_line(data=newdata.temp[newdata.temp$day==247,], 
                aes(x=base10, y=lowerCI), color="red", linewidth=1.3, lty=2) +
      geom_line(data=newdata.temp[newdata.temp$day==247,], 
                aes(x=base10, y=upperCI), color="red", linewidth=1.3, lty=2 ) +
      
      geom_point(data=df, aes(x=base10, y=weighted), size = 1.75) +
      
      theme_bw() +
      labs(y = "Weighted captures", x = " ", title = "All leafhopper species") +
      theme(axis.text.y = element_text(size = 14, colour = "black"),
            axis.text.x.top = element_blank(),
            axis.text.x.bottom = element_text(size = 14, colour = "black"),
            axis.ticks.x.top = element_blank(),
            plot.title = element_text(size = 18),
            axis.title = element_text(size = 16))     +
      scale_x_continuous(position = "top",                 
                        sec.axis = sec_axis(~ . + 10,  # J'ai effacé l'axe x principal et mis un 2e axe x pour montrer les températures, plutôt que les degrés-jours base10
                                            name = paste("Temperature (°C)"),
                                            breaks = seq(10,25, 2))) 
      
    
    #B) Depending on time of year ##
    
    newdata.day <- expand.grid(
      day = seq(118, 247, 1), 
      base10 = quantile(df$base10, c(0, 0.50, 1)),
      exposition.length = 1, 
      traps = 1) 
    
    newdata.day$pred <- predict(b10.21, newdata=newdata.day, re.form=NA, type="response")
    
    # Bootstrap CI:  
    boot <- bootMer(b10.21, 
                   function(x) predict(x, newdata = newdata.day, re.form = NA, type = 'response'), 
                   nsim = 200, 
                   use.u = FALSE, 
                   type = 'parametric', 
                   parallel = 'multicore', 
                   ncpus = 4)
    
    newdata.day$lowerCI <- apply(boot$t, 2, function(x) quantile(x, 0.025)) # Générer l'intervalle de confiance inférieur
    newdata.day$upperCI <- apply(boot$t, 2, function(x) quantile(x, 0.975)) # Générer l'intervalle de confiance supérieur
    
    
    # Modeling illustration, with weighted data
    base10 <- sort(unique(newdata.day$base10)) 
    round(base10, 2)
    
    ggplot() +
      
      geom_line(data=newdata.day[newdata.day$base10==base10[1],], 
                aes(x=day, y=pred), color="cadetblue3", linewidth=1.3, lty=1 ) +
      geom_line(data=newdata.day[newdata.day$base10==base10[1],], 
                aes(x=day, y=lowerCI), color="cadetblue3", linewidth=1.3, lty=2) +
      geom_line(data=newdata.day[newdata.day$base10==base10[1],], 
                aes(x=day, y=upperCI), color="cadetblue3", linewidth=1.3, lty=2 ) +
      
      geom_line(data=newdata.day[newdata.day$base10==base10[2],], 
                aes(x=day, y=pred), color="gray30", linewidth=1.3, lty=1 ) +
      geom_line(data=newdata.day[newdata.day$base10==base10[2],], 
                aes(x=day, y=lowerCI), color="gray30", linewidth=1.3, lty=2) +
      geom_line(data=newdata.day[newdata.day$base10==base10[2],], 
                aes(x=day, y=upperCI), color="gray30", linewidth=1.3, lty=2 ) +
      
      geom_line(data=newdata.day[newdata.day$base10==base10[3],], 
                aes(x=day, y=pred), color="red", linewidth=1.3, lty=1 ) +
      geom_line(data=newdata.day[newdata.day$base10==base10[3],], 
                aes(x=day, y=lowerCI), color="red", linewidth=1.3, lty=2) +
      geom_line(data=newdata.day[newdata.day$base10==base10[3],], 
                aes(x=day, y=upperCI), color="red", linewidth=1.3, lty=2 ) +
      
      geom_point(data=df, aes(x= day, y=weighted), size = 1.75) +
      
      theme_bw() +
      labs(color = "Région", x = " ", y = "Weighted captures", title = "All leafhopper species") +
      theme(axis.text = element_text(size = 14, colour = "black"),
            plot.title = element_text(size = 18),
            axis.title = element_text(size = 16))    +
      scale_x_continuous(breaks=c(123, 153, 183, 214, 245),
                    labels=c("May","June", "July", "August", "September")) 

    
   
