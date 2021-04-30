# Dans ce code R, on va essayer de reproduire la 1ère méthode, expliquée pour Stata
# dans 



#####################################################################################################################################
# INSPIRATIONS SUR STATA ? ##########################################################################################################
#####################################################################################################################################

# POSSIBILITÉ 1: "Microeconometrics Using Stata" (2009) à la page 261 ###################################################


# POSSIBILITÉ 2: CODE DE Zachariah Rutledge: voir fichier "Code_methode_1_Zachariah_Rutledge.do" ######################





#####################################################################################################################################
# TENTATIVE DE CODE R ###############################################################################################################
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



























