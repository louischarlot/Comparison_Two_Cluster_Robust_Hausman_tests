# In this code, we write the function associated to the 1st method to calculate the Cluster-Robst Hausman test, proposed by Wooldridge (2010)

# => More details available in the pdf.


######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
# TITLE : First Method for Cluster Robust Hausman Test RE vs FE (Wooldridge, 2010)
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################

# Some formalisms are inspired by the way functions are coded in plm Package of R. 


Method_1 <- function (y_it,x_it,w_i,data, vcov_chosen, NUMBER_MEAN_VARIABLES) {
  
  
  
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
  
  return(haus_robust)

  
  
}























