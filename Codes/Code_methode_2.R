
# In this R code, we use the 2nde méthod to calculate tehe Cluster-Robst Hausman test, proposed by Miller et Cameron (2015)
# This method is divided un many steps: 
#                         - utilization of paired cluster bootstrap to estimate V_hat in B subsamples
#                         - Haussman test for RE vs FE  


#############################################################################################################
#############################                                                   #############################  
#############################  CODE METHODE 2: PAIRS CLUSTER BOOTSTRAP VARIANCE #############################
#############################                                                   #############################  
#############################################################################################################



Method_2 <- function (y_it,x_it,data, B) {
  
  
  # compute the estimated coefficients 
  reg <- y_it ~ x_it
  fe_0 <- plm(formula = reg, data = data, model = "within")
  beta_0_fe <- fe_0$coefficients
  re_0 <- plm(formula = reg , data = data, model = "random")
  beta_0_re <- re_0$coefficients
  
  # Initialization for the bootstrap
  k_fe = length(beta_0_fe)
  k_re = length(beta_0_re) 
  
  beta_fe_boot = matrix(0,B,k_fe)
  beta_re_boot = matrix(0,B,k_re)
  
  #2. Start boostrap code
  
  for (b in 1:B) {
    ### a) resample ? => PROPOSITION RACHEL
    
    # get a vector with all clusters
    c <- sort(unique(data0$cluster))
    # c <- as.character(c)
    
    # group the data points per cluster
    clust.group <- function(c) {
      data0[data0$cluster==c,]             #### là ça ne fonctionne plus
    }
    
    clust.list <- lapply(c,clust.group)  #### et du coup là non plus
    
    # resample clusters with replacement
    c.sample <- sample(c, replace=T)  
    
    clust.sample <- clust.list[c.sample]
    
    clust.size <- number_times
    
    # combine the cluster list back to a single data matrix
    clust.bind <- function(c) {
      matrix(unlist(c),nrow=clust.size)
    }
    
    c.boot <- do.call(rbind,lapply(clust.sample,clust.bind)) # c.boot = the new data set (single bootstrap replicate)
    
    # Just to maintain columns name
    colnames(c.boot) <- names(data0)
    
    ### Here the problem is that in the new sample dataset c.boot, some 'clusters' have been sampled more than once so we need to rename them otherwise plm does not work
    c.boot <- as.data.frame(c.boot)
    c.boot <- c.boot[order(c.boot$cluster),]
    c.boot$cluster <- rep(1:number_clusters, each=number_times) # no more duplicates !
    
    ### Reorganize order of columns: "cluster" then "time" then "y_it" then "x_it" and transform y_it and x_it into numbers (because they became characters!)
    c.boot.ordered <- data.frame(cluster = c.boot$cluster, time = c.boot$time, y_it = as.numeric(c.boot$y_it), x_it = as.numeric(c.boot$x_it))
    
    ### b) FE model   
    fe_mod <- plm(formula = y_it ~ x_it, data = c.boot.ordered, model = "within")  
    beta_fe_boot[b,1:k_fe] <- fe_mod$coefficients
    
    ### c) RE model
    re_mod <- plm(formula = reg, data = c.boot.ordered, model = "random") 
    beta_re_boot[b,1:k_re] <- re_mod$coefficients
  }
  
  # re coefficients include time varying intercept while fe does not 
  # to have the same vector size, drop the intercept
  beta_0_re <- beta_0_re[2:k_re]
  beta_re_boot<-beta_re_boot[1:B,2:k_re]
  
  # We have to transpose it as x_it is one Unique variable (and thus it will lead to a line instead of a column as needed)
  beta_re_boot <- as.matrix(beta_re_boot)
  
  
  
  # beta_fe_boot- beta_re_boot :
  beta_fe_minus_re <- beta_fe_boot - beta_re_boot
  
  
  
  
  # 3. Create the haussman statistic 
  # Calculate V_hat using the formula in Cameron and Miller (2015)
  
  betahat_bar = 1/B*colSums(beta_fe_minus_re)
  # We substract to each element of b the mean over B that we just calculated:
  beta_fe_minus_re_demeaned <- sweep(beta_fe_minus_re, 2, betahat_bar, "-")
  
  list0 <- lapply(1:B, matrix, data=NA, nrow=k_fe, ncol=k_fe)
  # we realize the product b*t(b):
  for (i in 1:B){
    list0[[i]] <- beta_fe_minus_re_demeaned[i,]%*%t(beta_fe_minus_re_demeaned[i,])
  }
  # we sum the elements:
  sum_mat <- Reduce("+", list0)
  
  # We finally have the V_hat_FE using the formula in Cameron and Miller (2015)
  varhat_betahat <- 1/(B-1)*sum_mat
  
  
  # We calculate the difference between the estimated coefficients:
  diff_beta0_hat = beta_0_fe - beta_0_re
  
  
  
  
  
  # d) generate the Hausman test statistic => follows a chi-square distribution with 3 degrees of freedom (because 3 coefs estimated)
  H = t(diff_beta0_hat)%*%(varhat_betahat^(-1))%*%diff_beta0_hat
  
  
  
  # 4. Test the statistic
  # we want to compare the H statistic to a Chi-square distribution with 3 degrees of freedom
  # generate the p-value of H:
  p_value_H <- pchisq(H, df=k_fe, lower.tail=FALSE)
  
  
  return(c(p_value_H, H))
  
  
  
}















