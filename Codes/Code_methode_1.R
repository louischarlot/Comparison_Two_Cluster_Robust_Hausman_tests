# Dans ce code R, on va essayer de reproduire la 1ère méthode, expliquée pour Stata
# dans 



#####################################################################################################################################
# 2 INSPIRATIONS SUR STATA ##########################################################################################################
#####################################################################################################################################

# CODE STATA: POSSIBILITÉ 1: "Microeconometrics Using Stata" (2009) à la page 261 ###################################################
# Robust Hausman test using method of Wooldridge (2002)
# quietly xtreg lwage $xlist, re
# scalar theta = e(theta)
# global yandxforhausman lwage exp exp2 wks
# sort id
# foreach x of varlist $yandxforhausman {
#     by id : egen mean "x" = mean("x")
#     generate md "x" = "x" - mean "x"
#     generate red "x" = "x" - theta * mean "x"
#     }
# quietly regress redlwage redexp redexp2 redwks mdexp mdexp2 mdwks
# test mdexp mdexp2 mdwks
# ( 1) mdexp = 0
# ( 2) mdexp2 = 0
# ( 3) mdwks = 0
# F( 3, 4158) 848.39
# Prob > F = 0 . 0000


# CODE STATA: POSSIBILITÉ 2: CODE DE Zachariah Rutledge: voir fichier "Code_methode_1_Zachariah_Rutledge.do" ######################





#####################################################################################################################################
# TENTATIVE DE CODE R ###############################################################################################################
#####################################################################################################################################
# Les deux méthodes implémentées sur Stata sont tirées de (Wooldridge, 2002). Étant donné que notre article cite la version de (Wooldridge, 2010),
# plus récente (donc possiblement améliorée/corrigée par rapport à celle de 2002), il serait peut-être mieux d'essayer d'implémenter 
# celle ci.






























