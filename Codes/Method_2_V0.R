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

#0. Install and load usefull packages

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

#1. Model FE and RE (parametrisation)
    ### a) general parametrisation 
data("Gasoline", package = "plm")
data = Gasoline
form <- lgaspcar ~ lincomep + lrpmg + lcarpcap
#characteristics for model 
x = form
y = data$lgaspcar
model = c("within", "random")

  ### b) bootstrap parametrisation 
set.seed(1)

B = 399 #bootstrap coefficient

#generate the starting values
fe_0 <- plm(formula = x , data = data, model = "within", vcov = vcovHC)
beta_0_fe <- fe_0$coefficients
re_0 <- plm(formula = x , data = data, model = "random", vcov = vcovHC)
beta_0_re <- re_0$coefficients

k_fe = length(beta_0_fe)
k_re = length(beta_0_re) #beta_0_re has 4 coefficient because the first is the intercept

beta_fe_boot = matrix(0,B,k_fe)
beta_re_boot = matrix(0,B,k_re)

#2. Start boostrap code

for (b in 1:B) {
  ### a) resample ? 


  index_b <- sample(length(y),length(y),replace=TRUE)
  x_b <- lgaspcar[index_b] ~ lincomep[index_b] + lrpmg[index_b] + lcarpcap[index_b]
  
  ### b) FE model   
  fe_mod <- plm(formula = x_b , data = data, model = "within", vcov = vcovHC)
  beta_fe_boot[b,1:k_fe] <- fe_mod$coefficients
  
  ### c) RE model
  re_mod <- plm(formula = x_b , data = data, model = "random", vcov = vcovHC)
  beta_re_boot[b,1:k_re] <- re_mod$coefficients
  

}

# re coefficients include time varying intercept while fe does not 
# to have the same vector size, drop the intercept
beta_0_re <- beta_0_re[2:4]
beta_re_boot<-beta_re_boot[1:399,2:4]

# 3. Create the haussman statistic 

  ### a) Generate a vector of differences in coefficients
diff_beta0_hat <- beta_0_fe - beta_0_re
  ### b) Generate bootstrapped differences in coefficients
diff_betaboot_hat <- beta_fe_boot - beta_re_boot
  ### c) generate covariance matrix of bootstrapped differences
# we can use that V_hat(FE-RE) = V_hat_FE - V_hat_RE 
# use the funciton vcov?
# we can also use the formula in Cameron Miller

V_beta_boot_hat = 1/(399-1)*sum((diff_betaboot_hat)^2)

  ### d) generate the Hausman test statistic
H= diff_betaboot_hat*(V_beta_boot_hat)^(-1)*diff_betaboot_hat

# PROPOSITION RACHEL FOR 3.
# 3. Create the haussman statistic 

# we can use that V_hat(FE-RE) = V_hat_FE - V_hat_RE so we first calculate the two V_hat separately
betahat_bar_FE = 1/B*colSums(beta_fe_boot)
varhat_betahat_FE <- 1/(B-1)*colSums(t(beta_fe_boot-betahat_bar_FE)%*%(beta_fe_boot-betahat_bar_FE))

betahat_bar_RE = 1/B*colSums(beta_re_boot)
varhat_betahat_RE <- 1/(B-1)*colSums(t(beta_re_boot-betahat_bar_RE)%*%(beta_re_boot-betahat_bar_RE))

# The pairs cluster bootstrap variance is:
V_hat_FE_RE = varhat_betahat_FE - varhat_betahat_RE

# We calculate the difference between the estimated coefficients:
diff_beta0_hat = beta_0_fe - beta_0_re

### d) generate the Hausman test statistic for each coefficient
H1 = diff_beta0_hat[1]*(V_hat_FE_RE[1]^(-1))*diff_beta0_hat[1]
H2 = diff_beta0_hat[2]*(V_hat_FE_RE[2]^(-1))*diff_beta0_hat[2]
H3 = diff_beta0_hat[3]*(V_hat_FE_RE[3]^(-1))*diff_beta0_hat[3]

# 4. Test the statistic 

#############################################################################################################
##############################                                                 ##############################  
##############################       CODE SUR R A PARTIR DE LA DEFINITION      ##############################
##############################                                                 ##############################  
#############################################################################################################


# 0. Install and load the packages
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

# 1. Load the data 
data("Gasoline", package = "plm")
data = Gasoline
form <- lgaspcar ~ lincomep + lrpmg + lcarpcap
#characteristics for model 
x = form
model = c("within", "random")


set.seed(1)
B = 399 #bootstrap coefficient 

# 2. Generate FE and RE models 
fe_mod <- plm(formula = x, data = data, model = "within") #V_hat cannot be computed with plm
fe_mod_man <- lm(lgaspcar ~lincomep + lrpmg + lcarpcap -1 + factor(country), data = data) 
summary(fe_mod)
summary(fe_mod_man)

re_mod <- plm(formula = x, data = data, model = "random")
re_mod_man <- lmer(lgaspcar ~ lincomep + lrpmg + lcarpcap - 1 + (1|country), data = data)
re_mod_man_2 <- lme(lgaspcar ~ lincomep + lrpmg + lcarpcap - 1, data = data, random = ~1|country)

summary(re_mod)
summary(re_mod_man)

# 3. Compute the bootstrapped V_hat
V_hat_FE <- vcovBS(fe_mod_man, cluster = data$country, R = 399, type = "xy")

V_hat_RE <- vcovBS(re_mod_man_2, cluster = data$country, R = 399)

V_hat = V_hat_FE - V_hat_RE 

# 4. Obtaining beta_1_FE and beta_1_RE

# 5. Compute the Haussman test statistic 

T_hauss <- (beta_1_FE - beta_1_RE)V^-1(beta_1_FE - beta_1_RE)
