clear
clear matrix
clear mata
set maxvar 100000
set matsize 11000
estimates drop _all

use "C:\Users\Zach\Dropbox\Labor Scarcity Study\Hausman Test Files\Robust Hausman for Unbalanced Panels\Robust Hausman Test Sample Data.dta"


***HERE I MAKE THE "T_i" VARIABLE WHICH IDENTIFIES THE NUMBER OF YEARS THAT EACH 
***CROSS SECTIONAL OBSERVATION APPEARS IN THE DATA.  THIS IS CRITICAL AS THIS IS WHAT
***DIFFERS BETWEEN THE BALANCED PANEL SETUP AND THE UNBALANCED PANEL SETUP
***FOR EXAMPLE IF YOU ARE USING AN UNBALANCED PANEL OF STATES BETWEEN 1990 AND 2010
***YOU WOULD SWITCH THE "panel_id" VARIABLE IN MY CODE BELOW FOR THE STATE VARIABLE 
duplicates tag panel_id, gen(T_i)
replace T_i = T_i + 1


qui xtset panel_id year
***RUN THE STANDARD RE MODEL TO RETRIEVE THE SPECIFIC INFORMATION NEEDED TO QUASI-DEMEAN
***THE VARIABLES
qui xtreg ln_p ln_x i.commoditycode i.countycode i.countycode#c.year, re vce(cluster panel_id) theta
di _b[ln_x]

***PULL THE NECESSARY SCALARS OUT OF STATA'S e-MATIX TO CALCULATE "lambda_i"
qui scalar sigma_u_sq = [e(sigma_u)]^2
qui scalar sigma_e_sq = [e(sigma_e)]^2

***GEN THE "lamdda_i" FOR THE UNBALANCED PANEL. NOTE THAT THE CRITICAL DIFFERENCE
***HERE IS THE USE OF "T_i" INSTEAD OF "T".
qui gen lambda_i = 1 - sqrt((sigma_e_sq/(T_i*sigma_u_sq + sigma_e_sq)))


***HERE I CREATE A LIST CALLED "commodity" WHICH IDENTIFIES ALL THE CODES
***FOR MY VARIABLE COMMODITY, WHICH I NEED TO CONSTRUCT DUMMY VARIABLES FOR AND QUASI-DEMEANED
***VARIABLE FOR
global commodity  201119	///
201519	 ///
201999	 ///
202999	 ///
203999	 ///
204999	 ///
205999	 ///
206999	 ///
207999	 ///
209999	 ///
211999	 ///
212199	 ///
212399	 ///
212999	 ///
213199	 ///
214899	 ///
214999	 ///
215199	 ///
215399	 ///
216199	 ///
217999	 ///
218199	 ///
218299	 ///
218399	 ///
218499	 ///
221999	 ///
224999	 ///
226999	 ///
229999	 ///
230639	 ///
234799	 ///
236199	 ///
237199	 ///
237299	 ///
238199	 ///
239999	 ///
301999	 ///
302199	 ///
304199	 ///
304399	 ///
306999	 ///
307189	 ///
307199	 ///
307299	 ///
308999	 ///
309999	 ///
310999	 ///
314189	 ///
314199	 ///
314299	 ///
316189	 ///
316199	 ///
316299	 ///
318999	 ///
323999	 ///
325999	 ///
330999	 ///
331999	 ///
332999	 ///
333999	 ///
337999	 ///
339196	 ///
339999	 ///
340999	 ///
341999	 ///
342999	 ///
343999	 ///
344999	 ///
346999	 ///
348999	 ///
354299	 ///
354999	 ///
355999	 ///
357999	 ///
359999	 ///
363999	 ///
364999	 ///
366999	 ///
372999	 ///
375999	 ///
376999	 ///
378199	 ///
381999	 ///
394199	 ///
394999	 ///
398199	 ///
398399	 ///
398499	 ///
398999	 


***HERE I GENERATE A LIST CALLED "county" WHICH CONTAINS ALL THE CODES FOR THE
***10 COUNTIES I NEED TO MAKE DUMMY VARIABLES AND THE QUASI-DEMEANED DUMMY VARIABLES
global county 	1	 ///
         19 ///
         25 ///
         29 ///
         53 ///
         65 ///
         73 ///
         77 ///
         83 ///
        107 ///
        111 	
	

***HERE I MAKE A DUMMY VARIABLE FOR EACH COUNTY, A TIME-DEMEANED DUMMY VARIABLE
***FOR EACH COUNTY, AND A QUASI-DEMEANED VARIABLE FOR EACH COUNTY
foreach i in $county {	
gen county_`i'_dum = countycode==`i'		
egen mean_county_`i'_dum = mean(county_`i'_dum), by(panel_id)	
gen aug_county_dum_`i' = county_`i'_dum - mean_county_`i'_dum
gen aug2_county_dum_`i' = county_`i'_dum - lambda_i*mean_county_`i'_dum

***HERE I MAKE A COUNTY-SPECIFIC TREND VARIABLE FOR EACH COUNTY, A TIME-DEMEANED
****COUNTY-SPECIFIC TREND VARIABLE FOR EACH COUNTY, AND A QUASI-DEMEANED 
***COUNTY TREND VARIABLE FOR EACH COUNTY
gen county_yr_dum_`i' = county_`i'_dum*year
egen mean_county_yr_`i' = mean(county_yr_dum_`i'), by(panel_id)	
gen aug_county_yr_`i' = county_yr_dum_`i' - mean_county_yr_`i'	
gen aug2_county_yr_`i' = county_yr_dum_`i' - lambda_i*mean_county_yr_`i'			
}


***HERE I MAKE A DUMMY VARIABLE FOR EACH COUNTY, A TIME-DEMEANED DUMMY VARIABLE
***FOR EACH COMMODITY, AND A QUASI-DEMEANED VARIABLE FOR EACH COMMODITY
qui foreach e in $commodity	{	
gen commodity_dum_`e' = commoditycode==`e'	
egen mean_commodity_dum_`e' = mean(commodity_dum_`e'), by(panel_id)	 	
gen aug_commodity_dum_`e' = commodity_dum_`e' - mean_commodity_dum_`e'
gen aug2_commodity_dum_`e' = commodity_dum_`e' - lambda_i*mean_commodity_dum_`e'
}


***HERE I MAKE THE TIME-DEMEANED AND QUASI-DEMEANED DEPENDENT VARIABLES AND 
***MAIN REGRESSOR OF INTEREST
qui egen mean_ln_x = mean(ln_x), by(panel_id)
qui gen aug_ln_x = ln_x-mean_ln_x
qui gen aug2_ln_x = ln_x-lambda_i*mean_ln_x

qui egen mean_ln_p = mean(ln_p), by(panel_id)
qui gen aug_ln_p = ln_p-mean_ln_p
qui gen aug2_ln_p = ln_p - lambda_i*mean_ln_p

***HERE I CONSTRUCT THE QUASI DEMEANED CONSTANT
gen constant = (1-lambda_i)

***THIS IS THE STANDARD RE MODEL
qui xtreg ln_p ln_x i.commoditycode i.countycode i.countycode#c.year, re vce(cluster panel_id) theta 
di _b[ln_x]
di _se[ln_x]

***THIS IS THE MANUALLY QUASI-DEMEANDED RE MODEL...I JUST RUN THIS TO MAKE SURE THAT I 
***CONSTRUCTED THE QUASI DEMEANED VARIABLES CORRECTLY
***NOTE 1: THE COEFFICIENTS FROM THIS MODEL AND THE STANDARD RE MODEL ARE SLIGHTLY DIFFERENT DUE
***TO ROUNDING OF THE MANUALLY QUASI-DEMEANED VARIABLES
***NOTE 2: TO INCLUDE ALL OF THE QUASI-DEMEANED DUMMY VARIABLES AND TRENDS, I AM USING THE "*" OPTION
***WHICH INLCLUDES ALL VARIABLES THAT HAVE THE SAME LETTERS BEFORE THE "*"
qui reg aug2_ln_p aug2_ln_x aug2_commodity_* aug2_county_dum_* aug2_county_yr_* constant , nocons vce(cluster panel_id)
di _b[aug2_ln_x]
di _se[aug2_ln_x]

***THIS IS THE MODEL THAT WOOLDRIDGE (2002) SAYS TO TEST FOR ROBUST RE VS. FE MODEL
***NOTE 1: MY VARIABLE OF INTEREST IS THE "ln_x" VARIABLE, WHICH IS WHY I AM TESTING
***WHETHER "aug_ln_x" (i.e., THE TIME DEMEANED ln_x) IS EQUAL TO ZERO 
qui reg aug2_ln_p aug2_ln_x aug2_commodity_* aug2_county_dum_* constant aug_ln_x, nocons vce(cluster panel_id)
test aug_ln_x=0


***NOTE MY P-VALUE IS .71 SO I WOULD FAIL TO REJECT THE NULL HYPOTHESIS.
***THEREFORE, I DO NOT NEED TO USE AN FE MODEL IN MY CASE.


