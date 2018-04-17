Data structure that we want:
pairid	twinid	mz	dz1	dz2	sex	   age	sbp	  lnsbp   	bmi
12	         1	Yes	No	No	Female	20	109	  4.691348	18.59
12	         2	Yes	No	No	Female	20	115	  4.744932	17.71
23	         1	No	Yes	No	Female	18	103	  4.634729	24.77
23	         2	No	No	Yes	Female	18	94	  4.543295	30.3
36	         1	Yes	No	No	Female	23	103.5	4.639572	24.84
36	         2	Yes	No	No	Female	23	102.5	4.629863	26.35
38	         1	No	Yes	No	Male	  21	125	  4.828314	23.41
38	         2	No	No	Yes	Female	21	116.5	4.757891	21.45

Noting that dz1 and dz2 are redundant (is just twinid in case of dz zyg),
this can be reframed as 
  
pairid	twinid	mz  sex	 age	sbp	  lnsbp     bmi
12	         1	1	  1    20	  109	  4.691348	18.59
12	         2	1	  1    20	  115	  4.744932	17.71
23	         1	0	  1    18	  103	  4.634729	24.77
23	         2	0   1    18	  94	  4.543295	30.3
36	         1	1	  1    23	  103.5	4.639572	24.84
36	         2	1	  1    23	  102.5	4.629863	26.35
38	         1	0	  0	   21	  125	  4.828314	23.41
38	         2	0	  1    21	  116.5	4.757891	21.45


Since the Mz and Dz pops are independent probably it makes most sense to simulate these independently with parameterisation as appropriate.  The samples can later be appended vertically if required, with dummy variable to indicator group type.

So, we have:

pairid	twinid	mz  sex	 age	sbp	  lnsbp     bmi
12	         1	1	  1    20	  109	  4.691348	18.59
12	         2	1	  1    20	  115	  4.744932	17.71
36	         1	1	  1    23	  103.5	4.639572	24.84
36	         2	1	  1    23	  102.5	4.629863	26.35

pairid	twinid	mz  sex	 age	sbp	  lnsbp     bmi
23	         1	0	  1    18	  103	  4.634729	24.77
23	         2	0   1    18	  94	  4.543295	30.3
38	         1	0	  0	   21	  125	  4.828314	23.41
38	         2	0	  1    21	  116.5	4.757891	21.45

Note that pairid will be only unique in each sample, unless an offset  is applied (e.g. of n1/2).


Using the following respective mixed effects model regression results to model a plausible twin dataset including SBP and BMI

# ***************
# . mixed bmi b1.sex, || pairid: || pairid: mz dz1 dz2, covariance(identity) nocons
# ***************
# Mixed-effects ML regression                     Number of obs      =       298
# Group variable: pairid                          Number of groups   =       149
# 
#                                                 Obs per group: min =         2
#                                                                avg =       2.0
#                                                                max =         2
# 
# 
#                                                 Wald chi2(1)       =     30.23
# Log likelihood = -717.67982                     Prob > chi2        =    0.0000
# 
# ------------------------------------------------------------------------------
#          bmi |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
# -------------+----------------------------------------------------------------
#          sex |
#      Female  |  -1.968666   .3580699    -5.50   0.000     -2.67047   -1.266862
#        _cons |   23.74276   .3066701    77.42   0.000      23.1417    24.34382
# ------------------------------------------------------------------------------
# 
# ------------------------------------------------------------------------------
#   Random-effects Parameters  |   Estimate   Std. Err.     [95% Conf. Interval]
# -----------------------------+------------------------------------------------
# pairid: Identity             |
#                   var(_cons) |   4.664747   1.007701      3.054551    7.123752
# -----------------------------+------------------------------------------------
# pairid: Identity             |
#              var(mz dz1 dz2) |    2.31785   .8230305       1.15568    4.648716
# -----------------------------+------------------------------------------------
#                var(Residual) |   2.387788   .4105283      1.704708    3.344578
# ------------------------------------------------------------------------------
# LR test vs. linear regression:       chi2(2) =    77.89   Prob > chi2 = 0.0000
# 
# *************************
# mixed sbp bmi b1.sex, || pairid: || pairid: mz dz1 dz2, covariance(identity) nocons
# ************************* 
# Mixed-effects ML regression                     Number of obs      =       298
# Group variable: pairid                          Number of groups   =       149
# 
#                                                 Obs per group: min =         2
#                                                                avg =       2.0
#                                                                max =         2
# 
# 
#                                                 Wald chi2(2)       =     74.33
# Log likelihood = -1071.5334                     Prob > chi2        =    0.0000
# 
# ------------------------------------------------------------------------------
#          sbp |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
# -------------+----------------------------------------------------------------
#          bmi |   .6584735   .1829424     3.60   0.000     .2999131    1.017034
#      Female  |  -7.784706   1.210266    -6.43   0.000    -10.15678   -5.412628
#        _cons |   105.6484   4.437888    23.81   0.000     96.95034    114.3465
# ------------------------------------------------------------------------------
# 
# ------------------------------------------------------------------------------
#   Random-effects Parameters  |   Estimate   Std. Err.     [95% Conf. Interval]
# -----------------------------+------------------------------------------------
# pairid: Identity             |
#                   var(_cons) |   13.04321   9.196803      3.274917    51.94802
# -----------------------------+------------------------------------------------
# pairid: Identity             |
#              var(mz dz1 dz2) |   30.42725   11.80591      14.22305    65.09275
# -----------------------------+------------------------------------------------
#                var(Residual) |   40.43253   6.868871      28.98176    56.40754
# ------------------------------------------------------------------------------
# LR test vs. linear regression:       chi2(2) =    22.31   Prob > chi2 = 0.0000
# 



# So, in terms of simstudy def this might be: 



mz.pair <- defData(varname = "age", dist = "uniform", formula = "18; 31", id = "id.pair")
mz.pair <- defData(mz.pair,varname = "nTwins", formula = 2, variance = 0, id = "id.pair")
mz.pair <- defData(mz.pair,varname = "female", dist = "binary", formula = 0.6, id = "id.pair")
mz.pair <- defData(mz.pair,varname = "p0_sbp", dist = "normal", formula = 0, variance =  43, id = "id.pair")
mz.pair <- defData(mz.pair,varname = "p0_bmi", dist = "normal", formula = 0, variance = 6.8, id = "id.pair")
mz.pair

#    varname formula variance    dist     link
# 1:     age  18; 31      0.0 uniform identity
# 2:  female     0.6      0.0  binary identity
# 3:  p0_sbp       0     43.0  normal identity
# 4:  p0_bmi       0      6.8  normal identity

dt.mz <- genData(8, mz.pair)
dt.mz
## So these are twin pair level random effects (for Mz twins)
#    idpair      age female    p0_sbp     p0_bmi
# 1:      1 18.27667      1  5.293162 -3.2680182
# 2:      2 19.55972      1  3.831335  8.0497054
# 3:      3 30.94106      1  2.535445 -0.4937271
# 4:      4 24.19007      1  2.450477  0.4954653
# 5:      5 26.79905      0  2.162944 -2.1045156
# 6:      6 18.40369      0 11.015737 -0.9356570
# 7:      7 22.71510      0  8.044301  3.8664231
# 8:      8 22.00455      0 -7.038198  4.2270851

# using summary data from HeartExample_Stata.dta for parameterisation, by running for Mz twins...
# mixed bmi b1.sex if mz==1, || pairid: || pairid: mz
# dz.ind <- defDataAdd(varname = "female", dist = "binary", formula = 0.6)
mz.ind <- defDataAdd(varname = "sbp", dist = "normal",  formula = "120.6 - 6.16*female + p0_sbp", variance = 2.4)
mz.ind <- defDataAdd(mz.ind, varname = "bmi", dist = "normal",  formula = "22.7 -.815*female + p0_bmi", variance = 42)
dt.mz.ind <- genCluster(dt.mz, cLevelVar = "id.pair", numIndsVar = "nTwins", level1ID = "id.twin")
dt.mz.ind  <- addColumns(mz.ind, dt.mz.ind)
dt.mz.ind
dtMz <- addPeriods(dt.mz.ind)
#    varname                      formula variance    dist     link
# 1:  female                          0.5      0.0  binary identity
# 2:     age                    9.5; 10.5      0.0 uniform identity
# 3:     sbp 120.6 - 6.16*female + p0_sbp      2.4  normal identity
# 4:     bmi   22.7 -.815*female + p0_bmi     42.0  normal identity