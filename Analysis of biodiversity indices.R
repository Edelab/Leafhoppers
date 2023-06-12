# *********************************************************
# Analysis of biodiversity indices: Leafhopper project ####
# *********************************************************

# Co-written by  Jeanne Durivage [1] et Gaétan Daigle [2], March 13, 2023
# [1] : Laboratoire d'Edel Pérez-Lopez, Département de phytologie, Faculté des sciences de l'Agriculture et de l'Alimentation, Université Laval
# [2] : Consultant spécialisé en évaluation ou en statistique, Département de mathématiques et statistique, Faculté de Sciences et de Génie

# Code structure :
  # 1a) Calculation of alpha diversity indices
  # 1b) Illustration of indices
  # 2) Shannon repeated measures Anova
  # 3) Simpson repeated measures Anova
  # 4) Repeated measures Anova on species richness
  # 5) Simpson's inverse repeated measures Anova
  # 6) Repeated measures Anova on diversity index
  # 7) Repeated measures Anova on Pielou


rm(list=ls())
setwd("/Users/Jeanne/Documents/EdeLab/Indices_de_biodiversite")
load("pour_vegan") 
str(veg)


# ******************************************
# Calculation of alpha diversity indices ###
# ******************************************
  
  library(vegan); library(dplyr); library(tidyr); library(tibble)
  
  # Creation of a dataframe containing diversity indices
  indices <- data.frame( 
    shannon = diversity(veg, index = "shannon"),      
    simpson = diversity(veg, index = "simpson"),
    inverse.simpson = diversity(veg, index = "invsimpson"),
    richesse.specifique = specnumber(veg),
    diversity.index = specnumber(veg)/rowSums(veg),
    pielou = diversity(veg, index = "shannon") / log(specnumber(veg)) ) %>%
    rownames_to_column("reg_an") %>%
  separate(col=reg_an, into=c("site", "annee"), sep="_") %>% 
  mutate(region = case_when(
      site == 'TR' ~ 'Mauricie' ,
      site %in% c('OP1', 'OP2', 'FG') ~ 'Capitale-Nationale',
      site %in% c('VPLF', 'STD') ~ 'Montérégie',
      site == 'DEM' ~ 'Chaudière-Appalache'),
    region_code = case_when(
      region == 'Chaudière-Appalache' ~ 'CA',
      region == 'Capitale-Nationale' ~ 'CN',
      region == 'Mauricie' ~ 'MAUR',
      region == 'Montérégie' ~ 'MONT'
    ))

  # Dataframe formatting for analysis
  indices2 = within(indices[indices$annee != "both" & indices$site !="Province",],{
    region = as.factor(region)
    annee = as.factor(annee)
    an = as.numeric(as.character(annee))
    site = as.factor(site)
  })
  
  str(indices2)

  # Structure in data
  with(indices2, table(paste(region, site, sep="-"), annee))

  # Cleaning the environment
  rm(veg, indices)
  
  
# ****************************
# Illustration of indices ####
# ****************************
  
  library(ggplot2); library(ggsignif)

regions <- indices2 %>% 
  dplyr::select(region_code, shannon, simpson) %>% 
  rename("facteur" = "region_code") %>% 
  mutate(type = "region")

annees <-   indices2 %>% 
  dplyr::select(annee, shannon, simpson)  %>%   
  rename("facteur" = "annee") %>% 
  mutate(type = "annee")

long_graph <- rbind(regions, annees) %>% 
  pivot_longer(cols=c("shannon", "simpson"))
    

col_indices <- c("gray50","gray75", "green4", "indianred2", "darkorange3", "deepskyblue3")
type <- c(annee = "Année", region = "Région")
variable <- c(shannon = "Shannon", simpson = "Simpson")


 bx<-long_graph %>% 
    ggplot( aes(x = facteur, y = value, fill = facteur) ) +
    geom_boxplot (show.legend = FALSE)  +
    facet_grid(name ~ type, scales = "free", space = "free_x",
               labeller=labeller(type=type, variable=variable)) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    scale_fill_manual(values=col_indices) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.3)))
 bx
 
 long_graph %>% 
   filter(name == "shannon") %>% 
   ggplot( aes(x = facteur, y = value, fill = facteur) ) +
   geom_boxplot (show.legend = FALSE)  +
   facet_grid(name ~ type, scales = "free", space = "free_x",
              labeller=labeller(type=type, variable=variable)) +
   theme_bw() +
   theme(axis.title.x = element_blank(),
         axis.title.y = element_blank()) +
   scale_fill_manual(values=col_indices) +
   scale_y_continuous(expand = expansion(mult = c(0.05, 0.3)))
 
 long_graph %>% 
   filter(name == "simpson") %>% 
   ggplot( aes(x = facteur, y = value, fill = facteur) ) +
   geom_boxplot (show.legend = FALSE)  +
   facet_grid(name ~ type, scales = "free", space = "free_x",
              labeller=labeller(type=type, variable=variable)) +
   theme_bw() +
   theme(axis.title.x = element_blank(),
         axis.title.y = element_blank()) +
   scale_fill_manual(values=col_indices) +
   scale_y_continuous(expand = expansion(mult = c(0.05, 0.3)))
 

# ************************************
# Shannon repeated measures Anova ####
# ************************************
  
  library(nlme); library(emmeans); library(multcomp)
  fit1 = lme(shannon ~ region*annee, random=~1 | site, data=indices2)
  fit2 = update(fit1, weights = varIdent(form = ~ 1 | region))
  fit3 = update(fit1, weights = varIdent(form = ~ 1 | annee))
  
  AIC(fit1, fit2, fit3)
  
  joint_tests(fit1)
  ## The region has a significant influence on biodiversity.
  ## The year has a significant influence on biodiversity.
  ## Year and region do not interact on the biodiversity index result.
  
  # Multiple comparaisons between regions
  a1 = emmeans(fit1, ~ region)
  cld(a1, Letters=letters)
  
  # Multiple comparaisons between years;
  a1 = emmeans(fit1, ~ annee)
  cld(a1, Letters=letters)
  
  # Normality
  r = residuals(fit1,type="normalized", level=0)
  hist(r,freq=F)
  xfit<-seq(min(r),max(r),length=40)
  yfit<-dnorm(xfit, mean=mean(r), sd=sd(r)) 
  lines(xfit, yfit,col="red",lwd=2) 
  shapiro.test(r) 
  e1071::kurtosis(r)  
    ## Shapiro et Kurtosis adéquats.
  
  # Variance homogeneity
  plot(fitted(fit1, level=0), r, pch=16, ylab="Normalized residuals", xlab="Predicted values")
  abline(h=0, lty=2)
  
  
# ************************************
# Simpson repeated measures Anova ####
# ************************************
  
  fit1 = lme( simpson ~ region*annee, random=~1 | site, data=indices2)
  fit2 = update(fit1, weights = varIdent(form = ~ 1 | region))
  fit3 = update(fit1, weights = varIdent(form = ~ 1 | annee))
  
  AIC(fit1, fit2, fit3)
  
      # Fit1 Normality
      r = residuals(fit1,type="normalized", level=0)
      hist(r,freq=F)
      xfit<-seq(min(r),max(r),length=40)
      yfit<-dnorm(xfit, mean=mean(r), sd=sd(r))
      lines(xfit, yfit,col="red",lwd=2) 
      shapiro.test(r) 
      e1071::kurtosis(r)  
      
      # Fit1 Variance homogeneity
      plot(fitted(fit1, level=0), r, pch=16, ylab="Normalized residuals", xlab="Predicted values")
      abline(h=0, lty=2)
      
  
      # Fit2 Normality
      r = residuals(fit2,type="normalized", level=0)
      hist(r,freq=F)
      xfit<-seq(min(r),max(r),length=40)
      yfit<-dnorm(xfit, mean=mean(r), sd=sd(r))
      lines(xfit, yfit,col="red",lwd=2) 
      shapiro.test(r) 
      e1071::kurtosis(r)  
      
      # Fit2 Variance homogeneity
      plot(fitted(fit1, level=0), r, pch=16, ylab="Normalized residuals", xlab="Predicted values")
      abline(h=0, lty=2)
      
      ### Since fit1 and fit2 have very similar AIC scores,
      ### but fit2 does not respect normality at all, we use fit1.

  joint_tests(fit1)
  ## The region does not have a significant influence on biodiversity.
  ## The year has a significant influence on biodiversity.
  ## Year and region do not interact on the biodiversity index result.
  
  
# ***********************************************
# Repeated measures Anova on species richness####
# ***********************************************
  
  fit1 = lme( richesse.specifique ~ region*annee, random=~1 | site, data=indices2)
  fit2 = update(fit1, weights = varIdent(form = ~ 1 | region))
  fit3 = update(fit1, weights = varIdent(form = ~ 1 | annee))
  
  AIC(fit1, fit2, fit3)
  
  joint_tests(fit1)
  ## The region does not have a significant influence on biodiversity.
  ## The year has no significant influence on biodiversity.
  ## The year and the region do not interact on the biodiversity index result.
  
  # Normality
  r = residuals(fit1,type="normalized", level=0)
  hist(r,freq=F)
  xfit<-seq(min(r),max(r),length=40)
  yfit<-dnorm(xfit, mean=mean(r), sd=sd(r))
  lines(xfit, yfit,col="red",lwd=2) 
  shapiro.test(r) 
  e1071::kurtosis(r)
  
  # Variance homogeneity
  plot(fitted(fit1, level=0), r, pch=16, ylab="Normalized residuals", xlab="Predicted values")
  abline(h=0, lty=2)
  

# **********************************************
# Simpson's inverse repeated measures Anova ####
# **********************************************
    
  fit1 = lme( inverse.simpson ~ region*annee, random=~1 | site, data=indices2)
  fit2 = update(fit1, weights = varIdent(form = ~ 1 | region))
  fit3 = update(fit1, weights = varIdent(form = ~ 1 | annee))
  
  AIC(fit1, fit2, fit3)

  
  joint_tests(fit1)
  ## Region has no significant influence on biodiversity.
  ## The year has a significant influence.
  ## Year and region do not interact on the biodiversity index result.
  
  # Normality
  r = residuals(fit1,type="normalized", level=0)
  hist(r,freq=F)
  xfit<-seq(min(r),max(r),length=40)
  yfit<-dnorm(xfit, mean=mean(r), sd=sd(r))
  lines(xfit, yfit,col="red",lwd=2) 
  shapiro.test(r) 
  e1071::kurtosis(r)
  
  # Variance homogeneity
  plot(fitted(fit1, level=0), r, pch=16, ylab="Normalized residuals", xlab="Predicted values")
  abline(h=0, lty=2)
  
  
# ***********************************************
# Repeated measures Anova on diversity index ####
# ***********************************************
  
  fit1 = lme( diversity.index ~ region*annee, random = ~1 | site, data=indices2)
  fit2 = update(fit1, weights = varIdent(form = ~ 1 | region)) # variances heterogenes selon la region
  fit3 = update(fit1, weights = varIdent(form = ~ 1 | annee))  # variances heterogenes selon l'annee
  
  AIC(fit1, fit2, fit3)
  
  joint_tests(fit1)
  ## The region has no significant influence on biodiversity.
  ## The year has a significant influence on biodiversity.
  ## Year and region do not interact on the biodiversity index result.
  
  # Normality
  r = residuals(fit1,type="normalized", level=0)
  hist(r,freq=F)
  xfit<-seq(min(r),max(r),length=40)
  yfit<-dnorm(xfit, mean=mean(r), sd=sd(r))
  lines(xfit, yfit,col="red",lwd=2) 
  shapiro.test(r) 
  e1071::kurtosis(r)  

  
  # Variance homogeneity
  plot(fitted(fit1, level=0), r, pch=16, ylab="Normalized residuals", xlab="Predicted values")
  abline(h=0, lty=2)

  
  
# **************************************
# Repeated measures Anova on Pielou ####
# **************************************
  
  
  fit1 = lme( pielou ~ region*annee, random=~1 | site, data=indices2)
  fit2 = update(fit1, weights = varIdent(form = ~ 1 | region)) # variances heterogenes selon la region
  fit3 = update(fit1, weights = varIdent(form = ~ 1 | annee))  # variances heterogenes selon l'annee
  
  AIC(fit1, fit2, fit3)
  
  
  joint_tests(fit1)
  ## The region does not have a significant influence on biodiversity.
  ## The year has no significant influence on biodiversity.
  ## The year and the region do not interact on the biodiversity index result.
  
  # Normalite
  r = residuals(fit1,type="normalized", level=0)
  hist(r,freq=F)
  xfit<-seq(min(r),max(r),length=40)
  yfit<-dnorm(xfit, mean=mean(r), sd=sd(r))
  lines(xfit, yfit,col="red",lwd=2) 
  shapiro.test(r) 
  e1071::kurtosis(r)
  
  # Variance homogeneity
  plot(fitted(fit1, level=0), r, pch=16, ylab="Normalized residuals", xlab="Predicted values")
  abline(h=0, lty=2)
    
  