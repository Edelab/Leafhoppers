# *************************************************************
# Impact of insecticide treatments on leafhopper abundance ####
# *************************************************************
# By Jeanne Durivage, under the supervision of Edel Pérez-Lopez, Phytology department, Laval University.


# Structure of this code :
# 1) ANOVA with all species
# 2) ANOVA with Macrosteles quadrilineatus
# 3) ANOVA with Empoasca fabae

### There will be some data manipulation at the beginning of each section. 

# ******************************
# 1) Anova with all species ####
# ******************************

rm(list=ls())
setwd("/Users/Jeanne/Documents/EdeLab/Insecticides")
load("insecticides") 
  str(insecticides)
  
  # Definition of factor variables
    insecticides2 = within(insecticides,{
      id = as.factor(id)
      site = as.factor(site)
      matiere.active = as.factor(matiere.active)
      traitement = as.factor(traitement)
      temps = traitement })
    str(insecticides2)
 

  # Structure of data 
    with(insecticides2, table(site, matiere.active))
  
  # Investigating a transformation to the response variable
    library(MASS)
    boxcox = boxcox(pondere ~ matiere.active*temps , data=insecticides2,
                    lambda = seq(-2, 2, length = 1000))
    # Log transformation required
    
  # Analyse;
  library(nlme); library(emmeans); library(multcomp)
  fit = lme(log(pondere) ~ matiere.active*temps, random = ~ 1 | site / id, data=insecticides2)
  
  joint_tests(fit)
  
  # Averages by active ingredient:time (all leafhopper species)
  a1 = emmeans(fit, ~ temps | matiere.active, type="response")
  cld(a1)
  
  
  # Normality
  r = residuals(fit,type="normalized", level=0)
  hist(r,freq=F)
  xfit<-seq(min(r),max(r),length=40)
  yfit<-dnorm(xfit, mean=mean(r), sd=sd(r))
  lines(xfit, yfit,col="red",lwd=2) 
  shapiro.test(r) 
  e1071::kurtosis(r)  
  
  
  # Variance homogeneity
  plot(fitted(fit, level=0), r, pch=16, ylab="Normalized residuals", xlab="Predicted values")
  abline(h=0, lty=2)
  
  
  
  
  # *********************************************
  # 2) Anova with Macrosteles quadrilineatus ####
  # *********************************************
  
  # Definition of factor variables
  insecticides2_macro = within(insecticides_macro,{
    id = as.factor(id)
    site = as.factor(site)
    matiere.active = as.factor(matiere.active)
    traitement = as.factor(traitement)
    temps = traitement })
  str(insecticides2_macro)
  
  
  # Investigating a transformation to the response variable
  boxcox = boxcox((pondere + 0.0000001) ~ matiere.active*temps , data=insecticides2_macro,
                  lambda = seq(-2, 2, length = 1000))
  # Square root transformation required
  
  # Analyse;
  fit.Mq = lme(sqrt(pondere) ~ matiere.active*temps, random = ~ 1 | site / id, data=insecticides2_macro)
  
  joint_tests(fit.Mq)
  
  # Averages by active ingredient:time (all leafhopper species)
  a1.Mq = emmeans(fit.Mq, ~ temps | matiere.active, type="response")
  cld(a1.Mq)
  
  
  # Normality
  r.Mq = residuals(fit.Mq,type="normalized", level=0)
  hist(r.Mq,freq=F)
  xfit.Ef<-seq(min(r.Mq),max(r.Mq),length=40)
  yfit.Ef<-dnorm(xfit.Ef, mean=mean(r.Mq), sd=sd(r.Mq))
  lines(xfit.Ef, yfit.Ef,col="red",lwd=2) 
  shapiro.test(r.Mq) 
  e1071::kurtosis(r.Mq)  
  
  
  # Variance homogeneity
  plot(fitted(fit.Mq, level=0), r.Mq, pch=16, ylab="Normalized residuals", xlab="Predicted values")
  abline(h=0, lty=2)
  
  
  
  
  # *********************************
  # 3) Anova with Empoasca fabae ####
  # *********************************

  
  # Definition of factor variables
  insecticides2_empoasca = within(insecticides_empoasca,{
    id = as.factor(id)
    site = as.factor(site)
    matiere.active = as.factor(matiere.active)
    traitement = as.factor(traitement)
    temps = traitement})
  str(insecticides2_empoasca)
  
  # Investigating a transformation to the response variable
  boxcox = boxcox(pondere ~ matiere.active*temps , data=insecticides2_empoasca,
                  lambda = seq(-2, 2, length = 1000))
  # Square root transformation required
  
  # Analyse;
  fit.Ef = lme(sqrt(pondere) ~ matiere.active*temps, random = ~ 1 | site / id, data=insecticides2_empoasca)
  
  joint_tests(fit.Ef)
  
  # Averages by active ingredient:time (all leafhopper species)
  a1.Ef = emmeans(fit.Ef, ~ temps | matiere.active, type="response")
  cld(a1.Ef)
  
  
  # Normality
  r.Ef = residuals(fit.Ef,type="normalized", level=0)
  hist(r.Ef,freq=F)
  xfit.Ef<-seq(min(r.Ef),max(r.Ef),length=40)
  yfit.Ef<-dnorm(xfit.Ef, mean=mean(r.Ef), sd=sd(r.Ef))
  lines(xfit.Ef, yfit.Ef,col="red",lwd=2) 
  shapiro.test(r.Ef) 
  e1071::kurtosis(r.Ef)  
  # Normality is OK.
  
  # Variance homogeneity
  plot(fitted(fit.Ef, level=0), r.Ef, pch=16, ylab="Normalized residuals", xlab="Predicted values")
  abline(h=0, lty=2)
  # OK
  
  
  