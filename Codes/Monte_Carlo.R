

# Title : Monte-Carlo comparison of the 2 methods (Wooldridge (2010) and Miller et Cameron (2015)) of Cluster-Robust 
#         Hausman test for RE vs FE.

# We realize these 2 tests with CORRELATION = "NO"  or "YES", which means NO correlation Corr(c_i,x_it) != 0 (Respect of Assumption RE.1.b of Wooldridge (2010)) 
# or correlation Corr(c_i,x_it) = 0 (Failure of RE.1.b).

# => More details available in the pdf.





#########################################################################################################
#########################################################################################################
# SETTINGS :         
#########################################################################################################
#########################################################################################################
rm(list=ls()) 		# Clear workspace

#install.packages("plm")
#install.packages("rsample")
#install.packages("survey")
#install.packages("sandwich")
#install.packages("BBmisc")
#install.packages("purrr")

library(plm)
library(rsample)
library(survey)
library(sandwich)
library(BBmisc)
library(purrr)


# Set working directory 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load the 2 METHODS
source("Code_methode_1.R")
source("Code_methode_2.R")


set.seed(1) 		    # Set seed for random number generator



#########################################################################################################
# CHOOSE CORRELATION OR NOT  :         
#########################################################################################################

#CORRELATION <- "NO"             # Case 1 : NO correlation Corr(c_i,x_it) != 0 (Respect of Assumption RE.1.b of Wooldridge (2010)) 
CORRELATION <- "YES"           # Case 2: correlation Corr(c_i,x_it) = 0 (Failure of RE.1.b)





#########################################################################################################
#########################################################################################################
# SETTINGS FOR  DGP :         
#########################################################################################################
#########################################################################################################

number_clusters = 10        # Set the number of clusters
number_times = 10
n = number_clusters*number_times      # The sample size


B = 399 # Set number of bootstrap iterations

num = 50			# Number of Monte Carlo iterations

delta = 5

# Number of w_it variables:
NUMBER_MEAN_VARIABLES <- 1

# Variance for Robust Wald statistic:
vcov_chosen <- vcovHC # Heteroscedasticity-consistent estimation of the covariance matrix of the coefficient estimates in regression models.







#########################################################################################################
#########################################################################################################
# MONTE-CARLO:         
#########################################################################################################
#########################################################################################################




#Initialise !!!!!!!!!!!!!!!!
p_value_method_1_sum = 0
p_value_method_2_sum = 0
Hausman_stat_method_1_sum = 0
Hausman_stat_method_2_sum = 0
number_times_p_value_method_1_inf_005 = 0
number_times_p_value_method_2_inf_005 = 0
number_times_p_value_method_1_inf_001 = 0
number_times_p_value_method_2_inf_001 = 0


for (iter in 1:num) {
  
  
  #########################################################################################################
  # DGP :         
  #########################################################################################################
  x_it = rnorm(n,1,1) 
  u_it = rnorm(n,0,1)  
  
  # Set the clusters i:
  cluster <- rep(1:number_clusters, times=1, each=n/number_clusters)
  #cluster <- as.character(rep(1:number_clusters, times=1, each=n/number_clusters))
  #cluster <- paste(cluster, "cluster", sep=".")
  
  # Set the times t:
  time <- as.character(rep(1:number_times, times=n/number_times, each=1))
  time <- paste(time, "time", sep=".")
  
  # Cluster names:
  cluster_names <- as.character(rep(1:number_clusters, times=1))
  cluster_names <- paste(cluster_names, "cluster", sep=".")
  
  # w_i is the meab of x_it in each cluster i: 
  data_to_determine_wi <- data.frame(x_it, cluster)
  data_to_determine_wi$w_i <- ave(data_to_determine_wi$x_it, data_to_determine_wi$cluster)
  w_i <- data_to_determine_wi$w_i
  
  
  
  
  if(CORRELATION == "NO"){
    ###########################################################################################################
    # Case 1 : NO correlation Corr(c_i,x_it) != 0 (Respect of RE.1.b) #########################################
    ###########################################################################################################
    # Same c_i (taken randomly) within a cluster:
    c_i = rep(rnorm(number_clusters,0,1), times=1, each=n/number_clusters)
    ###########################################################################################################
  }
  
  if(CORRELATION == "YES"){
    ###########################################################################################################
    # Case 2: correlation Corr(c_i,x_it) = 0 (Failure of RE.1.b)  #############################################
    ###########################################################################################################
    # There is correlation with the averages w_i (of x_it): 
    c_i = 0.01 * rep(rnorm(number_clusters,0,1), times=1, each=n/number_clusters)  +  0.99* w_i
    ###########################################################################################################
  }
  
  
  y_it = x_it * delta + c_i + u_it
  
  # Create the dataframe of observations:
  data0 <- data.frame(y_it = y_it, x_it = x_it, cluster = cluster, time = time)
  
  # We put the data into a panel-dataframe:
  data <- pdata.frame(data0, index = c("cluster","time"))
  
  # w_i is the meab of x_it in each cluster i:
  data$w_i <- ave(data$x_it, data$cluster)
  
  
  
  
  
  #########################################################################################################
  # Method 1:         
  #########################################################################################################
  M1 <- Method_1(y_it,x_it,w_i,data, vcov_chosen, NUMBER_MEAN_VARIABLES)
  
  p_value_method_1_sum = p_value_method_1_sum + M1$p.value
  Hausman_stat_method_1_sum = Hausman_stat_method_1_sum + M1$statistic
  
  if(M1$p.value < 0.05){
    number_times_p_value_method_1_inf_005 = number_times_p_value_method_1_inf_005 + 1
  }
  if(M1$p.value < 0.01){
    number_times_p_value_method_1_inf_001 = number_times_p_value_method_1_inf_001 + 1
  }
  
  
  
  #########################################################################################################
  # Method 2:         
  #########################################################################################################
  M2 <- Method_2(y_it,x_it,data, B)
  
  p_value_method_2_sum = p_value_method_2_sum + M2[1]
  Hausman_stat_method_2_sum = Hausman_stat_method_2_sum + M2[2]
  
  if(M2[1] < 0.05){
    number_times_p_value_method_2_inf_005 = number_times_p_value_method_2_inf_005 + 1
  }
  
  if(M2[1] < 0.01){
    number_times_p_value_method_2_inf_001 = number_times_p_value_method_2_inf_001 + 1
  }
  
  
}

# Results for method 1:

p_value_method_1_sum/num

Hausman_stat_method_1_sum/num

number_times_p_value_method_1_inf_005/num

number_times_p_value_method_1_inf_001/num


# Results for method 2:

p_value_method_2_sum/num

Hausman_stat_method_2_sum/num

number_times_p_value_method_2_inf_005/num

number_times_p_value_method_2_inf_001/num




#############
# n=1000  (WITH number_clusters = 100 and number_times = 10), B=399, num = 50

### 
# Case 1 : NO correlation Corr(c_i,x_it) = 0 (Respect of RE.1.b)  
# p_value_method_1_sum/num                      0.4664159
# Hausman_stat_method_1_sum/num                 1.453387 
# number_times_p_value_method_1_inf_005/num     0.14
# number_times_p_value_method_1_inf_001/num     0.04

# p_value_method_2_sum/num                      0.4930527
# Hausman_stat_method_2_sum/num                 1.086238
# number_times_p_value_method_2_inf_005/num     0.08
# number_times_p_value_method_2_inf_001/num     0.02

###
# Case 2 :  correlation Corr(c_i,x_it) != 0 (Failure of RE.1.b) 
# p_value_method_1_sum/num                      1.732134e-10
# Hausman_stat_method_1_sum/num                 91.66104 
# number_times_p_value_method_1_inf_005/num     1
# number_times_p_value_method_1_inf_001/num     1

# p_value_method_2_sum/num                      1.640855e-05
# Hausman_stat_method_2_sum/num                 33.11566
# number_times_p_value_method_2_inf_005/num     1
# number_times_p_value_method_2_inf_001/num     1





#############
# n=100  (WITH number_clusters = 10 and number_times = 10), B=399, num = 50
### 
# Case 1 : NO correlation Corr(c_i,x_it) = 0 (Respect of RE.1.b)  
# p_value_method_1_sum/num                      0.3955996
# Hausman_stat_method_1_sum/num                 2.753435  
# number_times_p_value_method_1_inf_005/num     0.18
# number_times_p_value_method_1_inf_001/num     0.12

# p_value_method_2_sum/num                      0.7279315
# Hausman_stat_method_2_sum/num                 0.2181224
# number_times_p_value_method_2_inf_005/num     0
# number_times_p_value_method_2_inf_001/num     0

###
# Case 2 :  correlation Corr(c_i,x_it) != 0 (Failure of RE.1.b) 
# p_value_method_1_sum/num                      0.04652651
# Hausman_stat_method_1_sum/num                 20.36881 
# number_times_p_value_method_1_inf_005/num     0.84
# number_times_p_value_method_1_inf_001/num     0.76

# p_value_method_2_sum/num                      0.1065866
# Hausman_stat_method_2_sum/num                 5.08248
# number_times_p_value_method_2_inf_005/num     0.58
# number_times_p_value_method_2_inf_001/num     0.32











