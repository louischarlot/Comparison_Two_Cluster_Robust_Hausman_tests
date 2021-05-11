# Dans ce code R, on va essayer de reproduire la 2nde méthode proposé par Miller et Cameron (2015)
# Cette méthode se décompose en plusieurs étapes: 
#                         - utilisation du paired cluster bootstrap pour estimer V_hat dans B subsamples
#                         - Haussman test for RE et FE  


# POSSIBILITÉ 1: reproduction du code stata rhaussman 
# POSSIBILITÉ 2: utilisation du code source R haussman test et modif V_hat  



#############################################################################################################
##############################                                                 ##############################  
##############################  TEST DE REPROCUDCTION DU CODE STATA RHAUSSMAN  ##############################
##############################                                                 ##############################  
#############################################################################################################





#############################################################################################################
##############################                                                 ##############################  
##############################       CODE SUR R A PARTIR DE LA DEFINITION      ##############################
##############################                                                 ##############################  
#############################################################################################################


### 1. Install and load the packages
#install.packages("plm")
#install.packages("rsample")
#install.packages("survey")
#install.packages("sandwich")
#install.packages("nlme")
#install.packages("lme4")


library(plm)
library(rsample)
library(survey)
library(sandwich)
library(nlme)
library(lme4)

### 2. Load the data 
data("Gasoline", package = "plm")
data = Gasoline
form <- lgaspcar ~ lincomep + lrpmg + lcarpcap
#characteristics for model 
x = form
model = c("within", "random")


set.seed(1)
B = 399 #bootstrap coefficient 

### 3. Generate FE and RE models 
fe_mod <- plm(formula = x, data = data, model = "within") #V_hat cannot be computed with plm
fe_mod_man <- lm(lgaspcar ~lincomep + lrpmg + lcarpcap -1 + factor(country), data = data) 
summary(fe_mod)
summary(fe_mod_man)


re_mod <- plm(formula = x, data = data, model = "random")
re_mod_man <- lmer(lgaspcar ~ lincomep + lrpmg + lcarpcap - 1 + (1|country), data = data)
re_mod_man_2 <- lme(lgaspcar ~ lincomep + lrpmg + lcarpcap - 1, data = data, random = ~1|country)

summary(re_mod)
summary(re_mod_man)

### 4. Compute the bootstrapped V_hat
V_hat_FE <- vcovBS(fe_mod_man, cluster = data$country, R = 399, type = "xy")

V_hat_RE <- vcovBS(re_mod_man_2, cluster = data$country, R = 399)

V_hat = V_hat_FE - V_hat_RE 

### 5. Obtaining beta_1_FE and beta_1_RE

### 6. Compute the Haussman test statistic 

T_hauss <- (beta_1_FE - beta_1_RE)V^-1(beta_1_FE - beta_1_RE), 