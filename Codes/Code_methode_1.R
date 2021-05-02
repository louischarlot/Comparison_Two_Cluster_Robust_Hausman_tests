# Dans ce code R, on va essayer de reproduire la 1ère méthode, expliquée pour Stata
# dans 



#####################################################################################################################################
# INSPIRATIONS SUR STATA ? ##########################################################################################################
#####################################################################################################################################

# POSSIBILITÉ 1: "Microeconometrics Using Stata" (2009) à la page 261 ###################################################


# POSSIBILITÉ 2: CODE DE Zachariah Rutledge: voir fichier "Code_methode_1_Zachariah_Rutledge.do" ######################





#####################################################################################################################################
# FONCTION SUR  R ###################################################################################################################
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
wi <- plm(form, data = Gasoline, model = "within") 
# Random effects:
re <- plm(form, data = Gasoline, model = "random") 

# Details of the function:
help(phtest)

# Classical Hausman test:
phtest(wi, re)
phtest(form, data = Gasoline)


# Robust Hausman Tests: comprendre les détails !!!!!!!!!!!!

phtest(form, data = Gasoline, method = "aux")
# robust Hausman test (regression-based)
phtest(form, data = Gasoline, method = "aux", vcov = vcovHC)
# robust Hausman test with vcov supplied as a
# function and additional parameters
phtest(form, data = Gasoline, method = "aux",
       vcov = function(x) vcovHC(x, method="white2", type="HC3"))

# To see the source code of the function:
getAnywhere(phtest)
getAnywhere(UseMethod)





# Simplified function of package "plm" is:  ######################################################################

phtest.formula <- function(x, data, model = c("within", "random"),
                           method = c("chisq", "aux"),
                           index = NULL, vcov = NULL, ...) {
  # NB: No argument 'effect' here, maybe introduce?
  #     it gets evaluated due to the eval() call for method="chisq"
  #     and since rev. 305 due to extraction from dots (...) in method="aux" as a quick fix
  #    If introduced as argument, change doc accordingly (currently, effect arg is mentioned in ...)
  
  
  switch(match.arg(method),
         chisq={
           cl <- match.call(expand.dots = TRUE)
           cl$model <- model[1]
           names(cl)[2] <- "formula"
           m <- match(plm.arg, names(cl), 0)
           cl <- cl[c(1,m)]
           cl[[1]] <- as.name("plm")
           plm.model.1 <- eval(cl, parent.frame())
           plm.model.2 <- update(plm.model.1, model = model[2])
           return(phtest(plm.model.1, plm.model.2))
         },
         ################################################## AUX = OPTION ISSUE DE WOODRIDGE (2010) !!!!!!!!
         aux={
           ## some interface checks here
           if (model[1] != "within") {
             stop("Please supply 'within' as first model type")
           }
           
           if (!is.null(vcov) && !is.function(vcov)) stop("argument 'vcov' needs to be a function")
           
           ## set pdata.frame
           if (!inherits(data, "pdata.frame")) data <- pdata.frame(data, index = index) #, ...)
           
           row.names(data) <- NULL # reset rownames of original data set (->numbers rownames in clean sequence) to make rownames
           # comparable for later comparison to obs used in estimation of models (get rid of NA values)
           # [needed because pmodel.response() and model.matrix() do not retain fancy rownames, but rownames]
           
           # rev. 305: quick and dirty fix for missing effect argument in function 
           # signature for formula interface/test="aux": see if effect is in dots and extract
           dots <- list(...)
           
           if (!is.null(dots$effect)) effect <- dots$effect else effect <- NULL
           
           # calculate FE and RE model ###################################################################################################
           fe_mod <- plm(formula = x, data = data, model = "within", effect = effect)
           re_mod <- plm(formula = x, data = data, model = "random", effect = effect)
           
           reY <- pmodel.response(re_mod)
           reX <- model.matrix(re_mod, cstcovar.rm = "intercept")
           feX <- model.matrix(fe_mod, cstcovar.rm = "all")
           
           
           
           
           
           
           
           dimnames(feX)[[2]] <- paste(dimnames(feX)[[2]], "tilde", sep=".")
           ## estimated models could have fewer obs (due dropping of NAs) compared to the original data
           ## => match original data and observations used in estimated models
           ## routine adapted from lmtest::bptest
           commonrownames <- intersect(intersect(intersect(row.names(data), names(reY)), row.names(reX)), row.names(feX))
           if (!(all(c(row.names(data) %in% commonrownames, commonrownames %in% row.names(data))))) {
             data <- data[commonrownames, ]
             reY  <- reY[commonrownames]
             reX  <- reX[commonrownames, ]
             feX  <- feX[commonrownames, ]
           }
           
           # Tests of correct matching of obs (just for safety ...)
           if (!all.equal(length(reY), nrow(data), nrow(reX), nrow(feX)))
             stop("number of cases/observations do not match, most likely due to NAs in \"data\"")
           if (any(c(is.na(names(reY)), is.na(row.names(data)), is.na(row.names(reX)), is.na(row.names(feX)))))
             stop("one (or more) rowname(s) is (are) NA")
           if (!all.equal(names(reY), row.names(data), row.names(reX), row.names(feX)))
             stop("row.names of cases/observations do not match, most likely due to NAs in \"data\"")
           
           ## fetch indices here, check pdata
           ## construct data set and formula for auxiliary regression
           data <- pdata.frame(cbind(index(data), reY, reX, feX))
           auxfm <- as.formula(paste("reY~",
                                     paste(dimnames(reX)[[2]],
                                           collapse="+"), "+",
                                     paste(dimnames(feX)[[2]],
                                           collapse="+"), sep=""))
           auxmod <- plm(formula = auxfm, data = data, model = "pooling")
           nvars <- dim(feX)[[2]]
           R <- diag(1, nvars)
           r <- rep(0, nvars) # here just for clarity of illustration
           omega0 <- vcov(auxmod)[(nvars+2):(nvars*2+1),
                                  (nvars+2):(nvars*2+1)]
           Rbr <- R %*% coef(auxmod)[(nvars+2):(nvars*2+1)] - r
           
           h2t <- as.numeric(crossprod(Rbr, solve(omega0, Rbr)))
           ph2t <- pchisq(h2t, df = nvars, lower.tail = FALSE)
           
           df <- nvars
           names(df) <- "df"
           names(h2t) <- "chisq"
           
           if (!is.null(vcov)) {
             vcov <- paste(", vcov: ",
                           paste(deparse(substitute(vcov))),
                           sep="")
           }
           
           haus2 <- list(statistic   = h2t,
                         p.value     = ph2t,
                         parameter   = df,
                         method      = paste("Regression-based Hausman test", vcov, sep=""),
                         alternative = "one model is inconsistent",
                         data.name   = paste(deparse(substitute(x))))
           class(haus2) <- "htest"
           return(haus2)
         })
}








#####################################################################################################################################
# TENTATIVE DE CODAGE SUR  R ########################################################################################################
#####################################################################################################################################

x = form
data = Gasoline
model = c("within", "random")
# "effect":	the effects introduced in the model, one of "individual", "time", "twoways", or "nested":
effect = "twoways" # PAS SÛR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



# Calculate FE and RE model: 
fe_mod <- plm(formula = x, data = data, model = "within", effect = effect)
re_mod <- plm(formula = x, data = data, model = "random", effect = effect)

# Construct (transformed) response of the RE model:
# pmodel.response() gives you the response variable (in this case that would be y) but 
# with the specified transformation (random effects transformation ("quasi demeaning")) applied to it:
reY <- pmodel.response(re_mod)

# Construct the  design matrix (or model matrix):
# When we have a mode written as "Y=Xβ+ε", the design matrix (or model matrix) is the matrix X.(I THINK => CHECK BETTER !!!)
# => see https://genomicsclass.github.io/book/pages/expressing_design_formula.html
# "cstcovar.rm": remove the constant columns, one of "none", "intercept", "covariates", "all"),
reX <- model.matrix(re_mod, cstcovar.rm = "intercept")
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
# auxmod <- plm(formula = auxfm, data = data_1, model = "pooling") => ERROR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# => BECAUSE OF THIS ERROR, WE RUN THE EQUIVALENT OLS REGRESSION => CHECK THAT IT IS EQUIALENT !!!!!!!!!!!!!
# plm(formula = x, data = data, model = "pooling") <=> lm(formula = x, data = data) => IN OTHER WORDS, CHECK IF THIS IS ALWAYS TRUE !!!!!!!
auxmod <- lm(formula = auxfm, data = data_1)

# Number of variables (in our example: lincomep.tilde, lrpmg.tilde and lcarpcap.tilde)
nvars <- dim(feX)[[2]]
# Identity matrix of dimension = (nvars x nvars)
R <- diag(1, nvars)
# Vector of zeros of dimension = nvars
r <- rep(0, nvars) # here just for clarity of illustration








omega0 <- vcov(auxmod)[(nvars+2):(nvars*2+1),
                       (nvars+2):(nvars*2+1)]
Rbr <- R %*% coef(auxmod)[(nvars+2):(nvars*2+1)] - r

h2t <- as.numeric(crossprod(Rbr, solve(omega0, Rbr)))
ph2t <- pchisq(h2t, df = nvars, lower.tail = FALSE)

df <- nvars
names(df) <- "df"
names(h2t) <- "chisq"

if (!is.null(vcov)) {
  vcov <- paste(", vcov: ",
                paste(deparse(substitute(vcov))),
                sep="")
}

haus2 <- list(statistic   = h2t,
              p.value     = ph2t,
              parameter   = df,
              method      = paste("Regression-based Hausman test", vcov, sep=""),
              alternative = "one model is inconsistent",
              data.name   = paste(deparse(substitute(x))))
class(haus2) <- "htest"
return(haus2)


















