

# Title : Monte-Carlo comparison of 2 methods





#########################################################################################################
#########################################################################################################
# SETTINGS :         
#########################################################################################################
#########################################################################################################

# Set working directory 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

rm(list=ls()) 		  # Clear workspace

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
c_i = rep(rnorm(number_clusters,0,1), times=1, each=n/number_clusters)  + 0.5 * w_i

  
  
  
  
  
y_it = x_it * delta + c_i + u_it





# Create the dataframe of observations:
data0 <- data.frame(y_it = y_it, x_it = x_it, cluster = cluster, time = time)

# We put the data into a panel-dataframe:
data <- pdata.frame(data0, index = c("cluster","time"))



#########################################################################################################
#########################################################################################################
# METHOD 1 :         
#########################################################################################################
#########################################################################################################




# w_i is the meab of x_it in each cluster i:
data$w_i <- ave(data$x_it, data$cluster)
NUMBER_MEAN_VARIABLES <- 1



# Variance for Robust Wald statistic:
vcov_chosen <- vcovHC # CHECK BETTER !!!!!!!!!!!!





# We run the pooled regression proposed by (Wooldridge, 2010): y_it = beta*x_it + chi*w_i + (a_i + u_it):
auxfm <- y_it ~ x_it + w_i
auxmod <- plm(formula = auxfm, data = data, model = "pooling") 

# Number of "mean" variables (in our example: lincomep.mean, lrpmg.mean and lcarpcap.mean)
nvars <- NUMBER_MEAN_VARIABLES
# Identity matrix of dimension = (nvars x nvars)
Id <- diag(1, nvars)
# Vector of zeros of dimension = nvars
Zeros <- rep(0, nvars) # here just for clarity of illustration

# Covariance matrix for the auxiliary regression for the "mean" variables (in our example: lincomep.mean, lrpmg.mean and lcarpcap.mean):
# => IT CAN (AND IN OUR CASE SHOULD) BE ROBUSTIFIED !!!
Covariance_mean <- vcov_chosen(auxmod)[(nvars+2):(nvars*2+1),
                                       (nvars+2):(nvars*2+1)]

# Operation that finally gives the coefficients of the "mean" variables in the auxiliary regression in the right form
# (in our example: lincomep.mean, lrpmg.mean and lcarpcap.mean) 
Estimates_mean <- Id %*% coef(auxmod)[(nvars+2):(nvars*2+1)] - Zeros


# We calculate our Cluster-Robust Wald statistic => SO HERE IT IS: t(Estimates_mean) %*% (Covariance_mean)^(-1) %*% Estimates_mean
# "crossprod(A,B)" is the cross-product of matrices A and B gives t(A) %*% B.
# "solve(a,b)" will solve the equation a %*% x = b for x, where b can be either a vector or a matrix.
Wald_stat <- as.numeric(crossprod(Estimates_mean, solve(Covariance_mean, Estimates_mean)))

# We calculate the p-value of our Cluster-Robust Hausman test:
# "pchisq" gives the probability that a chi2(df) > Haussman_stat (with df number of degrees of freedom) => CHECK ONCE AGAIN !!!!
pWald <- pchisq(Wald_stat, df = nvars, lower.tail = FALSE)

# We name "df" the degrees of freedom and "chisq" the calculated chi-squared:
df <- nvars
names(df) <- "df"
names(Wald_stat) <- "chisq"

# If "vcov" function is not the default one, we display which one it is in our final result: => NE MARCHE PAS MAIS PAS TRÈS GRAVE !!!!!!!!!!!!!!
if (!is.null(vcov)) {
  vcov_chosen <- paste(", vcov: ",
                       paste(deparse(substitute(vcov_chosen))),
                       sep="")
}

# We display the results of our Cluster-Robust Hausman test:
haus_robust <- list(statistic   = Wald_stat,
                    p.value     = pWald,
                    parameter   = df,
                    method      = paste("Regression-based Hausman test", vcov_chosen, sep=""),
                    alternative = "one model is inconsistent",
                    data.name   = paste(deparse(substitute(x))))
class(haus_robust) <- "htest"
haus_robust
haus_robust$statistic






# Case 1 : NO correlation Corr(c_i,x_it) != 0 (Respect of RE.1.b)  : chisq = 1.0029, df = 1, p-value = 0.3166


# Case 2 :  correlation Corr(c_i,x_it) = 0 (Failure of RE.1.b) 












