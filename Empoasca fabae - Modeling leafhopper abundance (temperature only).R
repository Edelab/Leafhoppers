# ****************************************************************
# Empoasca - Modeling leafhopper abundance (temperature only) ####
# ****************************************************************
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

df.Ef <- pop %>% 
  filter(genus == "Empoasca" & specie == "fabae") %>% 
  aggregate(
    cbind(nb, weighted) ~ region + site + pose + withdrawal + traps + mean.temp +
      base10 + base8 + base5 + base2 + mean.date,
    FUN = sum) %>% 
  mutate(exposition.length = as.numeric(withdrawal - pose),
         day = yday(mean.date),
         region = as.factor(region),
         site = as.factor(site) )  %>% 
  dplyr::select( -c(pose, withdrawal, mean.date)) 


names(df.Ef)

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
Ef.mt.0 <- glmer.nb(nb ~ 1 + (1 | site) + offset(log(exposition.length * traps)), 
                    control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE), 
                    data = df.Ef)
Ef.mt.1 <- glmer.nb(nb ~ poly(day, 1) + poly(mean.temp, 1) + (1 | site) +
                      offset(log(exposition.length * traps)), 
                    control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                    data = df.Ef)
Ef.mt.2 <- glmer.nb(nb ~ poly(day, 2) + poly(mean.temp, 2) + (1 | site) +
                      offset(log(exposition.length * traps)), 
                    control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                    data = df.Ef)
Ef.mt.3 <- glmer.nb(nb ~ poly(day, 3) + poly(mean.temp, 3) + (1 | site) +
                      offset(log(exposition.length * traps)), 
                    control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                    data = df.Ef)

## base2 ##
Ef.b2.0 <- glmer.nb(nb ~ 1 + (1 | site) + offset(log(exposition.length * traps)), 
                    control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE), 
                    data = df.Ef)
Ef.b2.1 <- glmer.nb(nb ~ poly(day, 1) + poly(base2, 1) + (1 | site) +
                      offset(log(exposition.length * traps)), 
                    control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                    data = df.Ef)
Ef.b2.2 <- glmer.nb(nb ~ poly(day, 2) + poly(base2, 2) + (1 | site) +
                      offset(log(exposition.length * traps)), 
                    control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                    data = df.Ef)
Ef.b2.3 <- glmer.nb(nb ~ poly(day, 3) + poly(base2, 3) + (1 | site) +
                      offset(log(exposition.length * traps)), 
                    control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                    data = df.Ef)

## base5 ##
Ef.b5.0 <- glmer.nb(nb ~ 1 + (1 | site) + offset(log(exposition.length * traps)), 
                    control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE), 
                    data = df.Ef)
Ef.b5.1 <- glmer.nb(nb ~ poly(day, 1) + poly(base5, 1) + (1 | site) +
                      offset(log(exposition.length * traps)), 
                    control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                    data = df.Ef)
Ef.b5.2 <- glmer.nb(nb ~ poly(day, 2) + poly(base5, 2) + (1 | site) +
                      offset(log(exposition.length * traps)), 
                    control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                    data = df.Ef)
Ef.b5.3 <- glmer.nb(nb ~ poly(day, 3) + poly(base5, 3) + (1 | site) +
                      offset(log(exposition.length * traps)), 
                    control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                    data = df.Ef)

## base8 ##
Ef.b8.0 <- glmer.nb(nb ~ 1 + (1 | site) + offset(log(exposition.length * traps)), 
                    control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE), 
                    data = df.Ef)
Ef.b8.1 <- glmer.nb(nb ~ poly(day, 1) + poly(base8, 1) + (1 | site) +
                      offset(log(exposition.length * traps)), 
                    control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                    data = df.Ef)
Ef.b8.2 <- glmer.nb(nb ~ poly(day, 2) + poly(base8, 2) + (1 | site) +
                      offset(log(exposition.length * traps)), 
                    control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                    data = df.Ef)
Ef.b8.3 <- glmer.nb(nb ~ poly(day, 3) + poly(base8, 3) + (1 | site) +
                      offset(log(exposition.length * traps)), 
                    control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                    data = df.Ef)

## base10 ##
Ef.b10.0 <- glmer.nb(nb ~ 1 + (1 | site) + offset(log(exposition.length * traps)), 
                     control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE), 
                     data = df.Ef)
Ef.b10.1 <- glmer.nb(nb ~ poly(day, 1) + poly(base10, 1) + (1 | site) +
                       offset(log(exposition.length * traps)), 
                     control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                     data = df.Ef)
Ef.b10.2 <- glmer.nb(nb ~ poly(day, 2) + poly(base10, 2) + (1 | site) +
                       offset(log(exposition.length * traps)), 
                     control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                     data = df.Ef)
Ef.b10.3 <- glmer.nb(nb ~ poly(day, 3) + poly(base10, 3) + (1 | site) +
                       offset(log(exposition.length * traps)), 
                     control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                     data = df.Ef)



# Best model choice
library(AICcmodavg)
Ef.modeles <- list(Ef.mt.0, Ef.mt.1, Ef.mt.2, Ef.mt.3, Ef.b10.0, Ef.b10.1, Ef.b10.2, Ef.b10.3,
                   Ef.b8.0, Ef.b8.1, Ef.b8.2, Ef.b8.3, Ef.b5.0, Ef.b5.1, Ef.b5.2, Ef.b5.3, 
                   Ef.b2.0, Ef.b2.1, Ef.b2.2, Ef.b2.3)
Ef.noms <- c("Ef.mt.0", "Ef.mt.1", "Ef.mt.2", "Ef.mt.3", "Ef.b10.0", "Ef.b10.1", "Ef.b10.2", "Ef.b10.3",
             "Ef.b8.0", "Ef.b8.1", "Ef.b8.2", "Ef.b8.3", "Ef.b5.0", "Ef.b5.1", "Ef.b5.2", "Ef.b5.3", 
             "Ef.b2.0", "Ef.b2.1", "Ef.b2.2", "Ef.b2.3")

Ef.choix.BIC <- bictab(cand.set = Ef.modeles, modnames = Ef.noms)

head(Ef.choix.BIC)
# Average temperature, base 2 and base 5 all have the same value. This can be explained by the fact that when the species immigrates to Quebec, it is already very warm.
# Confirmation of this theory:
    range(df.Ef$mean.temp)
    # The lowest temperature is very high


# Observation of pre-selected models
summary(Ef.b8.3)
summary(Ef.b10.3)
summary(Ef.mt.3)
  # Similar summaries

# Using the selected indices, try combining several degrees to find the best combination.

Ef.mt.12 <- glmer.nb(nb ~ poly(day, 1) + poly(mean.temp, 2) + (1 | site) +
                       offset(log(exposition.length * traps)), 
                     control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                     data = df.Ef)
Ef.mt.13 <- glmer.nb(nb ~ poly(day, 1) + poly(mean.temp, 3) + (1 | site) +
                       offset(log(exposition.length * traps)), 
                     control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                     data = df.Ef)
Ef.mt.21 <- glmer.nb(nb ~ poly(day, 2) + poly(mean.temp, 1) + (1 | site) +
                       offset(log(exposition.length * traps)), 
                     control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                     data = df.Ef)
Ef.mt.23 <- glmer.nb(nb ~ poly(day, 2) + poly(mean.temp, 3) + (1 | site) +
                       offset(log(exposition.length * traps)), 
                     control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                     data = df.Ef)  
Ef.mt.31 <- glmer.nb(nb ~ poly(day, 3) + poly(mean.temp, 1) + (1 | site) +
                       offset(log(exposition.length * traps)), 
                     control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                     data = df.Ef)
Ef.mt.32 <- glmer.nb(nb ~ poly(day, 3) + poly(mean.temp, 2) + (1 | site) +
                       offset(log(exposition.length * traps)), 
                     control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                     data = df.Ef)


Ef.b8.12 <- glmer.nb(nb ~ poly(day, 1) + poly(base8, 2) + (1 | site) +
                       offset(log(exposition.length * traps)), 
                     control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                     data = df.Ef)
Ef.b8.13 <- glmer.nb(nb ~ poly(day, 1) + poly(base8, 3) + (1 | site) +
                       offset(log(exposition.length * traps)), 
                     control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                     data = df.Ef)
Ef.b8.21 <- glmer.nb(nb ~ poly(day, 2) + poly(base8, 1) + (1 | site) +
                       offset(log(exposition.length * traps)), 
                     control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                     data = df.Ef)
Ef.b8.23 <- glmer.nb(nb ~ poly(day, 2) + poly(base8, 3) + (1 | site) +
                       offset(log(exposition.length * traps)), 
                     control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                     data = df.Ef)  
Ef.b8.31 <- glmer.nb(nb ~ poly(day, 3) + poly(base8, 1) + (1 | site) +
                       offset(log(exposition.length * traps)), 
                     control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                     data = df.Ef)
Ef.b8.32 <- glmer.nb(nb ~ poly(day, 3) + poly(base8, 2) + (1 | site) +
                       offset(log(exposition.length * traps)), 
                     control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                     data = df.Ef)


Ef.b10.12 <- glmer.nb(nb ~ poly(day, 1) + poly(base10, 2) + (1 | site) +
                        offset(log(exposition.length * traps)), 
                      control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                      data = df.Ef)
Ef.b10.13 <- glmer.nb(nb ~ poly(day, 1) + poly(base10, 3) + (1 | site) +
                        offset(log(exposition.length * traps)), 
                      control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                      data = df.Ef)
Ef.b10.21 <- glmer.nb(nb ~ poly(day, 2) + poly(base10, 1) + (1 | site) +
                        offset(log(exposition.length * traps)), 
                      control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                      data = df.Ef)
Ef.b10.23 <- glmer.nb(nb ~ poly(day, 2) + poly(base10, 3) + (1 | site) +
                        offset(log(exposition.length * traps)), 
                      control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                      data = df.Ef)  
Ef.b10.31 <- glmer.nb(nb ~ poly(day, 3) + poly(base10, 1) + (1 | site) +
                        offset(log(exposition.length * traps)), 
                      control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                      data = df.Ef)
Ef.b10.32 <- glmer.nb(nb ~ poly(day, 3) + poly(base10, 2) + (1 | site) +
                        offset(log(exposition.length * traps)), 
                      control = glmerControl(optimizer="bobyqa", calc.derivs = FALSE),
                      data = df.Ef)

# Selection of the best model from the combinations created.
Ef.modeles2 <- list(Ef.mt.12, Ef.mt.13, Ef.mt.21, Ef.mt.23, Ef.mt.31, Ef.mt.32,
                    Ef.mt.2, Ef.mt.3, Ef.mt.1,
                    Ef.b10.12, Ef.b10.13, Ef.b10.21, Ef.b10.23, Ef.b10.31, Ef.b10.32,
                    Ef.b10.2, Ef.b10.3, Ef.b10.1, 
                    Ef.b8.12, Ef.b8.13, Ef.b8.21, Ef.b8.23, Ef.b8.31, Ef.b8.32,
                    Ef.b8.2, Ef.b8.3, Ef.b8.1)
Ef.noms2 <- c("Ef.mt.12", "Ef.mt.13", "Ef.mt.21", "Ef.mt.23", "Ef.mt.31", "Ef.mt.32", 
              "Ef.mt.1", "Ef.mt.2", "Ef.mt.3",
              "Ef.b8.12", "Ef.b8.13", "Ef.b8.21", "Ef.b8.23", "Ef.b8.31", "Ef.b8.32", 
              "Ef.b8.1", "Ef.b8.2", "Ef.b8.3",
              "Ef.b10.12", "Ef.b10.13", "Ef.b10.21", "Ef.b10.23", "Ef.b10.31", "Ef.b10.32", 
              "Ef.b10.1", "Ef.b10.2", "Ef.b10.3")

Ef.choix.BIC2 <- bictab(cand.set = Ef.modeles2, modnames = Ef.noms2)
head(Ef.choix.BIC2)
# The mt31 model is the best, but the others, with different temperatures, are similar. This can be explained by the fact that this is a migratory species. 
summary(Ef.mt.31)

#The following lines were add on 05/12/2023

#R2 conditional and marginal
library(performance)
model_performance(Ef.mt.31)

#AGVIF
library(car)
vif(Ef.mt.31)

#model suitability DHARMa
library(DHARMa)
modsuit<- simulateResiduals(fittedModel = Ef.mt.31)
testResiduals(modsuit)


# Cleannig environnement
rm(list=setdiff(ls(), c("Ef.mt.31", "df.Ef")))  





# ********************************************************
# Illustrate the model, with its confidence intervals ####
# ********************************************************
library(ggplot2)
setwd("/Users/Jeanne/Documents/EdeLab/Modelisation_de_l'abondance/Graphiques_de_modelisation")

# A) Depending on the temperature index selected
range(df.Ef$mean.temp)

# Model prediction on NEWDATA
newdata.temp <- expand.grid(
  mean.temp = seq( 10.9, 24.6, 0.1),
  day = quantile(df.Ef$day, c(0, 0.50, 1)),
  exposition.length = 1, 
  traps = 1)

newdata.temp$pred = predict(Ef.mt.31, newdata=newdata.temp, re.form=NA, type="response")

# Bootstrap CI (confidence intervals) :  
boot = bootMer(Ef.mt.31, 
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
  
  geom_line(data=newdata.temp[newdata.temp$day==138,], 
            aes(x=mean.temp, y=pred), color="cadetblue3", linewidth=1.3, lty=1 ) +
  geom_line(data=newdata.temp[newdata.temp$day==138,], 
            aes(x=mean.temp, y=lowerCI), color="cadetblue3", linewidth=1.3, lty=2) +
  geom_line(data=newdata.temp[newdata.temp$day==138,], 
            aes(x=mean.temp, y=upperCI), color="cadetblue3", linewidth=1.3, lty=2 ) +
  
  geom_line(data=newdata.temp[newdata.temp$day==196,], 
            aes(x=mean.temp, y=pred), color="gray30", linewidth=1.3, lty=1 ) +
  geom_line(data=newdata.temp[newdata.temp$day==196,], 
            aes(x=mean.temp, y=lowerCI), color="gray30", linewidth=1.3, lty=2) +
  geom_line(data=newdata.temp[newdata.temp$day==196,], 
            aes(x=mean.temp, y=upperCI), color="gray30", linewidth=1.3, lty=2 ) +
  
  geom_line(data=newdata.temp[newdata.temp$day==247,], 
            aes(x=mean.temp, y=pred), color="red", linewidth=1.3, lty=1 ) +
  geom_line(data=newdata.temp[newdata.temp$day==247,], 
            aes(x=mean.temp, y=lowerCI), color="red", linewidth=1.3, lty=2) +
  geom_line(data=newdata.temp[newdata.temp$day==247,], 
            aes(x=mean.temp, y=upperCI), color="red", linewidth=1.3, lty=2 ) +
  
  geom_point(data=df.Ef, aes(x=mean.temp, y=weighted), color = "black", size = 1.75) +
  
  theme_bw() +
  theme(axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 18, face="italic"),
        axis.title = element_text(size = 16))    +
  labs(x = "Temperature (°C)", y = "Weighted captures", title = "Empoasca fabae") +
  scale_x_continuous(breaks = seq(10,25, 2)) 


#B) Depending on the time of year
range(df.Ef$day)

newdata.day = expand.grid(
  day = seq(138, 247, 1), 
  mean.temp = quantile(df.Ef$mean.temp, c(0, 0.50, 1)),
  exposition.length = 1, 
  traps = 1) 

newdata.day$pred = predict(Ef.mt.31, newdata=newdata.day, re.form=NA, type="response")

# Bootstrap CI:  
boot = bootMer(Ef.mt.31, 
               function(x) predict(x, newdata = newdata.day, re.form = NA, type = 'response'), 
               nsim = 200, 
               use.u = FALSE, 
               type = 'parametric', 
               parallel = 'multicore', 
               ncpus = 4)

newdata.day$lowerCI = apply(boot$t, 2, function(x) quantile(x, 0.025))
newdata.day$upperCI = apply(boot$t, 2, function(x) quantile(x, 0.975))


# Modeling illustration, with weighted data
mean.temp = sort(unique(newdata.day$mean.temp)) 
round(mean.temp, 2)

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
  
  geom_point(data=df.Ef, aes(x= day, y=weighted), size = 1.75, color = "black") +
  
  theme_bw() +
  theme(axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 18, face="italic"),
        axis.title = element_text(size = 16))    +
  labs(x = " ", y = "Weighted captures", title = "Empoasca fabae") +
  scale_x_continuous(breaks=c(153, 183, 214, 245),
                     labels=c("June", "July", "August", "September")) 





