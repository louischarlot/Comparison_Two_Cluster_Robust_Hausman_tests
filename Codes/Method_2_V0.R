# Dans ce code R, on va essayer de reproduire la 2nde méthode proposé par Miller et Cameron (2015)
# Cette méthode se décompose en plusieurs étapes: 
#                         - utilisation du paired cluster bootstrap pour estimer V_hat dans B subsamples
#                         - Haussman test for RE et FE  


#############################################################################################################
#############################                                                   #############################  
#############################  CODE METHODE 2: PAIRS CLUSTER BOOTSTRAP VARIANCE #############################
#############################                                                   #############################  
#############################################################################################################

#0. Install and load usefull packages

install.packages("plm")
install.packages("rsample")
install.packages("survey")
install.packages("sandwich")
install.packages("BBmisc")
install.packages("purrr")

library(plm)
library(rsample)
library(survey)
library(sandwich)
library(BBmisc)
library(purrr)

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

# compute the estimated coefficients 
fe_0 <- plm(formula = x , data = data, model = "within")
beta_0_fe <- fe_0$coefficients
re_0 <- plm(formula = x , data = data, model = "random")
beta_0_re <- re_0$coefficients

# Initialization for the bootstrap
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
     
  x_b <- lgaspcar ~ lincomep + lrpmg + lcarpcap     
  
  ### b) FE model   
  fe_mod <- plm(formula = x_b , data = c.boot, model = "within")  
  beta_fe_boot[b,1:k_fe] <- fe_mod$coefficients
  
  ### c) RE model
  re_mod <- plm(formula = x_b , data = c.boot, model = "random") 
  beta_re_boot[b,1:k_re] <- re_mod$coefficients
}

# re coefficients include time varying intercept while fe does not 
# to have the same vector size, drop the intercept
beta_0_re <- beta_0_re[2:4]
beta_re_boot<-beta_re_boot[1:B,2:4]


# 3. Create the haussman statistic 
# we can use that V_hat(FE-RE) = V_hat_FE - V_hat_RE so we first calculate the two V_hat separately

# Calculate V_hat_FE using the formula in Cameron and Miller (2015)
betahat_bar_FE = 1/B*colSums(beta_fe_boot)
beta_fe_boot_demeaned <- sweep(beta_fe_boot, 2, betahat_bar_FE, "-")

list0_FE <- lapply(1:B, matrix, data=c(0,0,0,0,0,0,0,0,0), nrow=3, ncol=3)
for (i in 1:B){
  list0_FE[[i]] <- beta_fe_boot_demeaned[i,]%*%t(beta_fe_boot_demeaned[i,])
}
sum_mat_fe <- Reduce("+", list0_FE)

varhat_betahat_FE <- 1/(B-1)*sum_mat_fe

# Calculate V_hat_RE similarly 
betahat_bar_RE = 1/B*colSums(beta_re_boot)
beta_re_boot_demeaned <- sweep(beta_re_boot, 2, betahat_bar_RE, "-")

list0_RE <- lapply(1:B, matrix, data=c(0,0,0,0,0,0,0,0,0), nrow=3, ncol=3)
for (i in 1:B){
  list0_RE[[i]] <- beta_re_boot_demeaned[i,]%*%t(beta_re_boot_demeaned[i,])
}
sum_mat_re <- Reduce("+", list0_RE)

varhat_betahat_RE <- 1/(B-1)*sum_mat_re

# The pairs cluster bootstrap variance is:
V_hat_FE_RE = varhat_betahat_FE - varhat_betahat_RE

# We calculate the difference between the estimated coefficients:
diff_beta0_hat = beta_0_fe - beta_0_re

### d) generate the Hausman test statistic => follows a chi-square distribution with 3 degrees of freedom (because 3 coefs estimated)
H = t(diff_beta0_hat)%*%(V_hat_FE_RE^(-1))%*%diff_beta0_hat


# 4. Test the statistic 
# we want to compare the H statistic to a Chi-square distribution with 3 degrees of freedom
#generate the p-value of H: 
p_value_H <- pchisq(H, df=3, lower.tail=FALSE)
p_value_H

