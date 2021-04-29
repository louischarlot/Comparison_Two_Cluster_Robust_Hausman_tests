# Dans ce code R, on va essayer de reproduire la 1ère méthode, expliquée pour Stata
# dans "Microeconometrics Using Stata" (2009) à la page 261.



# CODE STATA : Robust Hausman test using method of Wooldridge (2002)
# 
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



