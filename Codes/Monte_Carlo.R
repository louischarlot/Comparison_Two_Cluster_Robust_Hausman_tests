

# Title : Monte-Carlo comparison of 2 methods





#########################################################################################################
#########################################################################################################
# SETTINGS :         
#########################################################################################################
#########################################################################################################


#install.packages("plm")
library (plm)


# Set working directory 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))







set.seed(1) 		    # Set seed for random number generator








#########################################################################################################
#########################################################################################################
# DGP :         
#########################################################################################################
#########################################################################################################

n = 100      			  # Set the sample size
number_clusters = 10        # Set the number of clusters
number_times = 10

delta = 5
x_it = rnorm(n,1,1) 
u_it = rnorm(n,0,1)  




# Set the clusters i:
cluster <- as.character(rep(1:number_clusters, times=1, each=n/number_clusters))
cluster <- paste(cluster, "cluster", sep=".")

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





###########################################################################################################
# Case 1 : NO correlation Corr(c_i,x_it) != 0 (Respect of RE.1.b) #########################################
###########################################################################################################

# Same c_i (taken randomly) within a cluster:
# c_i = rep(rnorm(number_clusters,0,1), times=1, each=n/number_clusters)

###########################################################################################################
# Case 2: correlation Corr(c_i,x_it) = 0 (Failure of RE.1.b)  #############################################
###########################################################################################################

# There is correlation with the averages w_i (of x_it): 
c_i = 0.01 * rep(rnorm(number_clusters,0,1), times=1, each=n/number_clusters)  +  0.99* w_i


  
  
  
  
y_it = x_it * delta + c_i + u_it





# Create the dataframe of observations:
data0 <- data.frame(y_it = y_it, x_it = x_it, cluster = cluster, time = time)

# We put the data into a panel-dataframe:
data <- pdata.frame(data0, index = c("cluster","time"))



#########################################################################################################
#########################################################################################################
# MONTE-CARLO METHOD 1 :         
#########################################################################################################
#########################################################################################################

# Variance for Robust Wald statistic:
vcov_chosen <- vcovHC # CHECK BETTER !!!!!!!!!!!!

# w_i is the meab of x_it in each cluster i:
data$w_i <- ave(data$x_it, data$cluster)
NUMBER_MEAN_VARIABLES <- 1




Method_1(y_it,x_it,w_i,data, vcov_chosen, NUMBER_MEAN_VARIABLES)



# chisq = 16.913, df = 1, p-value = 3.914e-05






for (it in 1:num) {
  # Data generating process
  x = rnorm(n,1,1) 
  sigma_epsilon = sqrt(x^2*sigma2_beta + sigma2_u)
  epsilon = rnorm(n,0,1)*sigma_epsilon
  y = x * beta + epsilon
  y[y<0] <- 0
  
  result <- optim(par = theta_start, RC_l, y = y, x = x, method = c("L-BFGS-B"), lower = c(-Inf,0,0.1), upper = c(Inf,Inf,Inf), hessian=TRUE)
  theta_hat = result$par
  
  
  # For the HISTOGRAM ##############################################################
  
  
}









# Case 1 : NO correlation Corr(c_i,x_it) != 0 (Respect of RE.1.b)  : chisq = 1.0029, df = 1, p-value = 0.3166


# Case 2 :  correlation Corr(c_i,x_it) = 0 (Failure of RE.1.b) 
















