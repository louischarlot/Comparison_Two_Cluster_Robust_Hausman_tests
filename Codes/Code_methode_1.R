# Dans ce code R, on va essayer de reproduire la 1ère méthode, expliquée pour Stata
# dans 



#####################################################################################################################################
# ABNDONNED: INSPIRATIONS SUR STATA ? ###############################################################################################
#####################################################################################################################################

# ABNDONNED:
# POSSIBILITÉ 1: "Microeconometrics Using Stata" (2009) à la page 261 ###################################################
# POSSIBILITÉ 2: CODE DE Zachariah Rutledge: voir fichier "Code_methode_1_Zachariah_Rutledge.do" ######################





#####################################################################################################################################
# FONCTION SUR  R EXISTANTE POUR COMMENCER ##########################################################################################
#####################################################################################################################################
# Les deux méthodes implémentées sur Stata sont tirées de (Wooldridge, 2002). Étant donné que notre article cite la version de (Wooldridge, 2010),
# plus récente (donc possiblement améliorée/corrigée par rapport à celle de 2002), il serait peut-être mieux d'essayer d'implémenter 
# celle ci.

# There is a function of the package "plm" that has been developed to calculate the "cluster-robust Hausman test" of (Wooldridge, 2010):


#install.packages("plm")
library (plm)

# Load the data:
data("Gasoline", package = "plm")
form <- lgaspcar ~ lincomep + lrpmg + lcarpcap

# Fixed effects:
fe <- plm(form, data = Gasoline, model = "within") 
# Random effects:
re <- plm(form, data = Gasoline, model = "random") 

# Details of the function:
help(phtest)

# Classical Hausman test:
phtest(fe, re)
phtest(form, data = Gasoline)


# Robust Hausman Tests: comprendre les détails !!!!!!!!!!!!

phtest(x=form, data = Gasoline, model = c("within", "random"), method = "aux")
# robust Hausman test (regression-based)
phtest(x=form, data = Gasoline, model = c("within", "random"), method = "aux", vcov = vcovHC)
# robust Hausman test with vcov supplied as a
# function and additional parameters
phtest(x=form, data = Gasoline, model = c("within", "random"), method = "aux",
       vcov = function(x) vcovHC(x, method="white2", type="HC3"))

# To see the source code of the function:
getAnywhere(phtest)
methods(phtest)
getAnywhere(phtest.formula)



gigi(x=form, data = Gasoline, model = c("within", "random"))

x=form
data = Gasoline
model = c("within", "random")




#####################################################################################################################################
# 1ere TENTATIVE DE CODAGE SUR  R (INSPIRÉ DE "phtest") #############################################################################
#####################################################################################################################################

x = form
data = Gasoline
model = c("within", "random")
# "effect":	the effects introduced in the model, one of "individual", "time", "twoways", or "nested":
# effect = "twoways" or NULL ??? # PAS SÛR => ou RIEN ????? => DANS CE CAS FCT LUI DONNE VALEUR NULL !!!!!!!!!!!!!!!!!!!!!!!!!!!
vcov_chosen = vcovHC # PAS SÛR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#vcov_chosen = NULL
index = NULL



# If the data is not pdata.frame, transform it into a pdata.frame:
# an object of class 'pdata.frame' is a data.frame with an index attribute that describes its individual and time dimensions.
if (!inherits(data, "pdata.frame")) data <- pdata.frame(data, index = index) #, ...)

# Calculate FE and RE model: 
fe_mod <- plm(formula = x, data = data, model = "within", effect = effect)
re_mod <- plm(formula = x, data = data, model = "random", effect = effect)

# Construct (transformed) response of the RE model:
# pmodel.response() gives you the response variable (in this case that would be y) but 
# with the specified transformation (random effects transformation ("quasi demeaning")) applied to it:
reY <- pmodel.response(re_mod)

# Construct the  transformed X for RE and FR => A REVERIFIER !!!!!
# When we have a mode written as "Y=Xβ+ε", the design matrix (or model matrix) is the matrix X.(I THINK => CHECK BETTER !!!)
# The model.matrix methods builds a model matrix with transformations performed as specified by the model and effect arguments 
# (and theta if model = "random" is requested), in this case the supplied data argument should be a model frame created by plm's model.frame method.
reX <- model.matrix(re_mod, cstcovar.rm = "intercept") # "cstcovar.rm": remove the constant columns, one of "none", "intercept", "covariates", "all")
feX <- model.matrix(fe_mod, cstcovar.rm = "all")

# We add ".tilde" at the dimension names of the feX: => VOIR MIEUX POURQUOI ON RENOMME AINSI !!!
dimnames(feX)[[2]] <- paste(dimnames(feX)[[2]], "tilde", sep=".")

# We match the observations of the 2 models (RE and FE), given that the possible dropping of NAs:
# estimated models could have fewer obs (due dropping of NAs) compared to the original data
# => match original data and observations used in estimated models
# (routine adapted from lmtest::bptest)
commonrownames <- intersect(intersect(intersect(row.names(data), names(reY)), row.names(reX)), row.names(feX)) # Rows that are kept for all 3 reY, reX, feX 
# If some rows have been dropped, we take only the rows that are present in all 3 reY, reX, feX: 
if (!(all(c(row.names(data) %in% commonrownames, commonrownames %in% row.names(data))))) { 
  data <- data[commonrownames, ]
  reY  <- reY[commonrownames]
  reX  <- reX[commonrownames, ]
  feX  <- feX[commonrownames, ]
}
# Tests of correct matching of observations we just have done (just for safety ...):
if (!all.equal(length(reY), nrow(data), nrow(reX), nrow(feX)))
  stop("number of cases/observations do not match, most likely due to NAs in \"data\"")
if (any(c(is.na(names(reY)), is.na(row.names(data)), is.na(row.names(reX)), is.na(row.names(feX)))))
  stop("one (or more) rowname(s) is (are) NA")
if (!all.equal(names(reY), row.names(data), row.names(reX), row.names(feX)))
  stop("row.names of cases/observations do not match, most likely due to NAs in \"data\"")


### We construct now data set and formula for auxiliary regression:

# An object of class 'pdata.frame' is a data.frame with an index attribute that describes its individual and time dimensions. => VOIR MIEUX !!!
data_2 <- pdata.frame(cbind(index(data), reY, reX, feX))

# We write the equation for the auxiliary regression, following WOOLDRIDGE (2010)
auxfm <- as.formula(paste("reY~",
                          paste(dimnames(reX)[[2]],
                                collapse="+"), "+",
                          paste(dimnames(feX)[[2]],
                                collapse="+"), sep=""))

# We then run the corresponding pooled regression:
auxmod <- plm(formula = auxfm, data = data_2, model = "pooling") # => ERROR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# Number of "tilde" variables (in our example: lincomep.tilde, lrpmg.tilde and lcarpcap.tilde)
nvars <- dim(feX)[[2]]
# Identity matrix of dimension = (nvars x nvars)
Id <- diag(1, nvars)
# Vector of zeros of dimension = nvars
Zeros <- rep(0, nvars) # here just for clarity of illustration

# Covariance matrix for the auxiliary regression for the "tilde" variables (in our example: lincomep.tilde, lrpmg.tilde and lcarpcap.tilde):
# => IT CAN (AND IN OUR CASE SHOULD) BE ROBUSTIFIED !!!
Covariance_tilde <- vcov_chosen(auxmod)[(nvars+2):(nvars*2+1),
                       (nvars+2):(nvars*2+1)]

# Operation that finally gives the coefficients of the "tilde" variables in the auxiliary regression 
# (in our example: lincomep.tilde, lrpmg.tilde and lcarpcap.tilde) => UNDERSTAND BETTER WHY WE ARE DOING THIS WAY !!!!!!!!!!!!!!!!!!!!!!!!!
Estimates_tilde <- Id %*% coef(auxmod)[(nvars+2):(nvars*2+1)] - Zeros


# We calculate our Cluster-Robust Wald statistic => SO HERE IT IS: t(Estimates_tilde) %*% (Covariance_tilde)^(-1) %*% Estimates_tilde
# "crossprod(A,B)" is the cross-product of matrices A and B gives t(A) %*% B.
# "solve(a,b)" will solve the equation a %*% x = b for x, where b can be either a vector or a matrix.
Wald_stat <- as.numeric(crossprod(Estimates_tilde, solve(Covariance_tilde, Estimates_tilde)))

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

# Wr display the results of our Cluster-Robust Hausman test:
haus2 <- list(statistic   = Wald_stat,
              p.value     = pWald,
              parameter   = df,
              method      = paste("Regression-based Hausman test", vcov_chosen, sep=""),
              alternative = "one model is inconsistent",
              data.name   = paste(deparse(substitute(x))))
class(haus2) <- "htest"
haus2









#####################################################################################################################################
# 2eme TENTATIVE DE CODAGE SUR  R (Wooldridge, 2010) => NON, RESTER SUR LA PREMIERE !!! #############################################
#####################################################################################################################################

# We use:
# y_it : lgaspcar
# x_it : lincomep + lrpmg + lcarpcap + country_i
# cluster : country
# w_it : lincomep + lrpmg + lcarpcap
# w_i : mean(lincomep) + mean(lrpmg) + mean(lcarpcap)


# Fixed-effects manually:
##fixed_manual <-lm(lgaspcar ~lincomep + lrpmg + lcarpcap + factor(country) - 1, data = Gasoline)
##summary(fixed_manual)
# Fixed-effects with plm:
##fe <- plm(form, data = Gasoline, model = "within") 
##summary(fe)



# We add the columns "w_i": mean values by cluster for the different w_it:
##data <- Gasoline
##data$lincomep_mean <- ave(data$lincomep, data$country)
##data$lrpmg_mean <- ave(data$lrpmg, data$country)
##$lcarpcap_mean <- ave(data$lcarpcap, data$country)



# Regression proposed by (Wooldridge, 2010): y_it = beta*x_it + chi*w_i + (a_i + u_it)

##reg_Wooldridge <- lm(lgaspcar ~lincomep + lrpmg + lcarpcap + lincomep_mean + lrpmg_mean + lcarpcap_mean + factor(country) - 1, data)
##summary(reg_Wooldridge)


















