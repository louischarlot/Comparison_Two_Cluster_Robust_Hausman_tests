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
  ### a) resample ? => PROPOSITION RACHEL

  # get a vector with all clusters
  c <- sort(unique(data$country))
  
  # group the data points per cluster
  clust.group <- function(c) {
    data[data$country==c,]
  }
  
  clust.list <- lapply(c,clust.group)
  
  # resample clusters with replacement
  c.sample <- sample(c, replace=T)  
  
  clust.sample <- clust.list[c.sample]
  
  clust.size <- 19
  
  # combine the cluster list back to a single data matrix
  clust.bind <- function(c) {
    matrix(unlist(c),nrow=clust.size)
  }
  
  c.boot <- do.call(rbind,lapply(clust.sample,clust.bind)) # c.boot = the new data set (single bootstrap replicate)
  
  # Just to maintain columns name
  colnames(c.boot) <- names(data)
  
### Here the problem is that in the new sample dataset c.boot, some 'clusters' have been sampled more than once so we need to rename them otherwise plm does not work
  c.boot <- as.data.frame(c.boot)
  c.boot <- c.boot[order(c.boot$country),]
  c.boot$country <- rep(1:18, each=19) # no more duplicates !
    
###______________________end proposition Rachel for cluster resampling    

  index_b <- sample(length(y),length(y),replace=TRUE)
  x_b <- lgaspcar[index_b] ~ lincomep[index_b] + lrpmg[index_b] + lcarpcap[index_b]
# x_b <- lgaspcar ~ lincomep + lrpmg + lcarpcap     ## if using above cluster resampling method  
  
  ### b) FE model   
  fe_mod <- plm(formula = x_b , data = data, model = "within", vcov = vcovHC)
# fe_mod <- plm(formula = x_b , data = c.boot, model = "within", vcov = vcovHC)  ## if using above cluster resampling method  
  beta_fe_boot[b,1:k_fe] <- fe_mod$coefficients
  
  ### c) RE model
  re_mod <- plm(formula = x_b , data = data, model = "random", vcov = vcovHC)
# re_mod <- plm(formula = x_b , data = c.boot, model = "random", vcov = vcovHC)  ## if using above cluster resampling method  
  beta_re_boot[b,1:k_re] <- re_mod$coefficients
}

# re coefficients include time varying intercept while fe does not 
# to have the same vector size, drop the intercept
beta_0_re <- beta_0_re[2:4]
beta_re_boot<-beta_re_boot[1:399,2:4]


# 3. Create the haussman statistic 
# we can use that V_hat(FE-RE) = V_hat_FE - V_hat_RE so we first calculate the two V_hat separately

# Calculate V_hat_FE using the formula in Cameron and Miller (2015)
betahat_bar_FE = 1/B*colSums(beta_fe_boot)
beta_fe_boot_demeaned <- sweep(beta_fe_boot, 2, betahat_bar_FE, "-")

list0_FE <- lapply(1:399, matrix, data=c(0,0,0,0,0,0,0,0,0), nrow=3, ncol=3)
for (i in 1:399){
  list0_FE[[i]] <- beta_fe_boot_demeaned[i,]%*%t(beta_fe_boot_demeaned[i,])
}
sum_mat_fe <- Reduce("+", list0_FE)

varhat_betahat_FE <- 1/(B-1)*sum_mat_fe

# Calculate V_hat_RE similarly 
betahat_bar_RE = 1/B*colSums(beta_re_boot)
beta_re_boot_demeaned <- sweep(beta_re_boot, 2, betahat_bar_RE, "-")

list0_RE <- lapply(1:399, matrix, data=c(0,0,0,0,0,0,0,0,0), nrow=3, ncol=3)
for (i in 1:399){
  list0_RE[[i]] <- beta_re_boot_demeaned[i,]%*%t(beta_re_boot_demeaned[i,])
}
sum_mat_re <- Reduce("+", list0_RE)

varhat_betahat_RE <- 1/(B-1)*sum_mat_re

# The pairs cluster bootstrap variance is:
V_hat_FE_RE = varhat_betahat_FE - varhat_betahat_RE

# We calculate the difference between the estimated coefficients:
diff_beta0_hat = beta_0_fe - beta_0_re

### d) generate the Hausman test statistic => follows a chi-square distribution with 3 degrees of freedom (because 3 coefs estimated)
H1 = t(diff_beta0_hat)%*%(V_hat_FE_RE^(-1))%*%diff_beta0_hat


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
