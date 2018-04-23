# Drawing on demos from
# https://cran.r-project.org/web/packages/simstudy/vignettes/simstudy.html
# with view to apply this to clustered twin data


# install.packages("simstudy")
require(simstudy)

# define variables
def <- defData(     varname = "nr",   dist = "nonrandom", formula = 7, id = "idnum")
def <- defData(def, varname = "x1",   dist = "uniform", formula = "10;20")
def <- defData(def, varname = "y1",   formula = "nr + x1 * 2", variance = 8)
def <- defData(def, varname = "y2",   dist = "poisson", formula = "nr - 0.2 * x1",     link = "log")
def <- defData(def, varname = "xCat", formula = "0.3;0.2;0.5", dist = "categorical")
def <- defData(def, varname = "g1",   dist = "gamma", formula = "5+xCat", variance = 1, link = "log")
def <- defData(def, varname = "a1",   dist = "binary", formula = "-3 + xCat",     link = "logit")

def
# varname       formula variance        dist     link
# 1:      nr             7        0   nonrandom identity
# 2:      x1         10;20        0     uniform identity
# 3:      y1   nr + x1 * 2        8      normal identity
# 4:      y2 nr - 0.2 * x1        0     poisson      log
# 5:    xCat   0.3;0.2;0.5        0 categorical identity
# 6:      g1        5+xCat        1       gamma      log
# 7:      a1     -3 + xCat        0      binary    logit

# generate data
dt <- genData(1000, def)
dt
# idnum nr       x1       y1  y2 xCat         g1 a1
# 1:     1  7 11.71823 30.61837 106    3 2454.58244  1
# 2:     2  7 19.72982 47.37233  16    3 3495.34206  0
# 3:     3  7 11.73166 28.97878 105    3  257.83624  1
# 4:     4  7 11.99212 29.98492 132    2 1196.95586  0
# 5:     5  7 17.72467 44.05648  40    1 1277.82901  0
# ---                                                  
#   996:   996  7 13.96225 33.14931  57    1  455.42894  0
# 997:   997  7 16.52625 37.24160  42    1   89.56957  0
# 998:   998  7 17.83418 42.89763  33    3  408.33425  1
# 999:   999  7 17.23273 44.51284  41    1  185.98869  0
# 1000:  1000  7 17.02450 39.10828  24    3 9837.95798  1


# clustered approach
gen.school <- defData(varname = "s0", dist = "normal", formula = 0, variance = 3, id = "idSchool")
gen.school <- defData(gen.school, varname = "nClasses", dist = "noZeroPoisson", formula = 3)
gen.school
#    idSchool         s0 nClasses
# varname formula variance          dist     link
# 1:       s0       0        3        normal identity
# 2: nClasses       3        0 noZeroPoisson identity

dtSchool <- genData(8, gen.school)
dtSchool <- trtAssign(dtSchool, n = 2)

dtSchool
#    idSchool trtGrp         s0 nClasses
# 1:        1      1 -0.8413091        6
# 2:        2      1  1.1981539        3
# 3:        3      0  1.1729038        4
# 4:        4      0  2.3362275        6
# 5:        5      1  2.1647581        3
# 6:        6      1 -1.2419999        5
# 7:        7      0 -2.4106126        1
# 8:        8      0  0.9470429        1

gen.class <- defDataAdd(varname = "c0", dist = "normal", formula = 0, variance = 2)
gen.class <- defDataAdd(gen.class, varname = "nStudents", dist = "noZeroPoisson", formula = 20)
#      varname formula variance          dist     link
# 1:        c0       0        2        normal identity
# 2: nStudents      20        0 noZeroPoisson identity


dtClass <- genCluster(dtSchool, "idSchool", numIndsVar = "nClasses", level1ID = "idClass")
dtClass <- addColumns(gen.class, dtClass)

head(dtClass, 10)
#     idSchool trtGrp         s0 nClasses idClass          c0 nStudents
#  1:        1      1 -0.8413091        6       1 -3.30814001        12
#  2:        1      1 -0.8413091        6       2 -1.84226332        16
#  3:        1      1 -0.8413091        6       3 -2.05428148        25
#  4:        1      1 -0.8413091        6       4 -0.85623318        17
#  5:        1      1 -0.8413091        6       5 -0.92143165        25
#  6:        1      1 -0.8413091        6       6 -0.37487937        20
#  7:        2      1  1.1981539        3       7 -1.42395635        14
#  8:        2      1  1.1981539        3       8 -0.55502125        22
#  9:        2      1  1.1981539        3       9 -0.83292393        19
# 10:        3      0  1.1729038        4      10 -0.09945386        22

gen.student <- defDataAdd(varname = "Male", dist = "binary", formula = 0.5)
gen.student <- defDataAdd(gen.student, varname = "age", dist = "uniform", formula = "9.5; 10.5")
gen.student <- defDataAdd(gen.student, varname = "test", dist = "normal",  formula = "50 - 5*Male + s0 + c0 + 8 * trtGrp", variance = 2)

gen.student
#    varname                            formula variance    dist     link
# 1:    Male                                0.5        0  binary identity
# 2:     age                          9.5; 10.5        0 uniform identity
# 3:    test 50 - 5*Male + s0 + c0 + 8 * trtGrp        2  normal identity

dtStudent <- genCluster(dtClass, cLevelVar = "idClass", numIndsVar = "nStudents",  level1ID = "idChild")
dtStudent <- addColumns(gen.student, dtStudent)
dtStudent
#      idSchool trtGrp         s0 nClasses idClass        c0 nStudents idChild Male       age     test
#   1:        1      1 -0.8413091        6       1 -3.308140        12       1    1 10.409240 48.39706
#   2:        1      1 -0.8413091        6       1 -3.308140        12       2    1  9.582417 51.51788
#   3:        1      1 -0.8413091        6       1 -3.308140        12       3    1 10.465350 48.92440
#   4:        1      1 -0.8413091        6       1 -3.308140        12       4    1  9.511273 47.39084
#   5:        1      1 -0.8413091        6       1 -3.308140        12       5    1 10.040281 50.15491
# ---                                                                                                
# 574:        8      0  0.9470429        1      29  1.071603        16     574    1  9.873724 47.06487
# 575:        8      0  0.9470429        1      29  1.071603        16     575    0 10.361916 51.06674
# 576:        8      0  0.9470429        1      29  1.071603        16     576    0  9.977288 51.99369
# 577:        8      0  0.9470429        1      29  1.071603        16     577    0  9.945570 52.75790
# 578:        8      0  0.9470429        1      29  1.071603        16     578    1  9.898581 47.99445

# multivariate correlated data (normal, binomial, poisson, gamma, or uniform?)
# mvnormal

# "genCorData requires the user to specify 
#       a mean vector mu, 
#       a single standard deviation or a vector of standard deviations sigma, 
#       and either 
#            a correlation matrix corMatrix or 
#            a correlation coefficient rho and a correlation structure corsrt."

# specifying a specific correlation matrix C
C <- matrix(c(1, 0.7, 0.2, 0.7, 1, 0.8, 0.2, 0.8, 1), nrow = 3)
C
#      [,1] [,2] [,3]
# [1,]  1.0  0.7  0.2
# [2,]  0.7  1.0  0.8
# [3,]  0.2  0.8  1.0

# generate 3 correlated variables with different location and scale for each
# field
dt <- genCorData(1000, mu = c(4, 12, 3), sigma = c(1, 2, 3), corMatrix = C)
dt
#         id       V1        V2         V3
#    1:    1 3.590053 10.467607  1.5365015
#    2:    2 4.390306 13.016778  4.9034960
#    3:    3 5.262963 13.838306  4.6115221
#    4:    4 3.753833  9.309315 -1.8590332
#    5:    5 5.230627 12.562042  0.7605537
#  ---                                   
#  996:  996 2.252507  7.557695 -2.6444042
#  997:  997 4.860714 11.382270 -0.5769707
#  998:  998 4.740323 15.177604  9.5828643
#  999:  999 4.476665 12.935593  4.6317326
# 1000: 1000 4.730124 10.373298 -1.0099788

# estimate correlation matrix
dt[, round(cor(cbind(V1, V2, V3)), 1)]
#     V1  V2  V3
# V1 1.0 0.7 0.2
# V2 0.7 1.0 0.8
# V3 0.2 0.8 1.0
# estimate standard deviation
dt[, round(sqrt(diag(var(cbind(V1, V2, V3)))), 1)]
# V1  V2  V3 
# 1.0 2.1 3.1

# generate 3 correlated variables with different location but same standard
# deviation and compound symmetry (cs) correlation matrix with correlation
# coefficient = 0.4.  Other correlation matrix structures are 'independent'
# ('ind') and 'auto-regressive' ('ar1').

dt <- genCorData(1000, mu = c(4, 12, 3), sigma = 3, rho = 0.4, corstr = "cs", 
                 cnames = c("x0", "x1", "x2"))
dt

# The new data generated by genCorData can be merged with an existing data set. 
# Alternatively, addCorData will do this directly:
# define and generate the original data set
def <- defData(varname = "x", dist = "normal", formula = 0, variance = 1, id = "cid")
dt <- genData(1000, def)

# add new correlate fields a0 and a1 to 'dt'
dt <- addCorData(dt, idname = "cid", mu = c(0, 0), sigma = c(2, 0.2), rho = -0.2, 
                 corstr = "cs", cnames = c("a0", "a1"))
dt

#        cid           x         a0          a1
#    1:    1 -0.23242361  0.3479707 -0.24173786
#    2:    2 -0.08577881  0.6513559 -0.11194365
#    3:    3  1.31791916  1.6415128 -0.06047269
#    4:    4 -0.31450927 -1.5848162  0.19856533
#    5:    5 -1.39704950  0.8907335 -0.04767962
#  ---                                        
#  996:  996  1.17912888  0.4937818 -0.13771564
#  997:  997  1.90045830  1.4829844  0.23149315
#  998:  998 -1.24643005  0.9578953  0.27172403
#  999:  999 -2.34069731 -2.6433781  0.31940490
# 1000: 1000 -0.42068363 -0.5471829  0.13101156

# estimate correlation matrix
dt[, round(cor(cbind(a0, a1)), 1)]
##      a0   a1
## a0  1.0 -0.2
## a1 -0.2  1.0

# estimate standard deviation
dt[, round(sqrt(diag(var(cbind(a0, a1)))), 1)]
##  a0  a1 
## 2.0 0.2


# for non-normal data, can use genCorGen and addCorGen

# Poisson
l <- c(8, 10, 12)  # lambda for each new variable
dx <- genCorGen(1000, nvars = 3, params1 = l, dist = "poisson", rho = 0.3, corstr = "cs", 
                wide = TRUE)
dx
##         id V1 V2 V3
##    1:    1  8 13  8
##    2:    2 12 13  9
##    3:    3 12 17 12
##    4:    4 11 15  9
##    5:    5  7 13 17
##   ---              
##  996:  996  9 11 12
##  997:  997  8 10 10
##  998:  998  6  6 12
##  999:  999  7  8 13
## 1000: 1000  9 10 13
round(cor(as.matrix(dx[, .(V1, V2, V3)])), 2)
##      V1   V2   V3
## V1 1.00 0.32 0.32
## V2 0.32 1.00 0.29
## V3 0.32 0.29 1.00

# Binary data
genCorGen(1000, nvars = 3, params1 = c(0.3, 0.5, 0.7), dist = "binary", rho = 0.8, corstr = "cs", wide = TRUE)
##         id V1 V2 V3
##    1:    1  1  1  1
##    2:    2  0  1  1
##    3:    3  0  1  1
##    4:    4  0  0  0
##    5:    5  0  0  1
##   ---              
##  996:  996  0  1  1
##  997:  997  1  1  1
##  998:  998  1  1  1
##  999:  999  0  1  1
## 1000: 1000  0  0  0

# The gamma distribution requires two parameters - the mean and dispersion. 
# (These are converted into shape and rate parameters more commonly used.)

dx <- genCorGen(1000, nvars = 3, params1 = l, params2 = c(1, 1, 1), dist = "gamma", 
                rho = 0.7, corstr = "cs", wide = TRUE, cnames = "a, b, c")
dx
##         id          a          b          c
##    1:    1  2.2147914  5.7828970 12.4494947
##    2:    2  1.3062344  6.1613218  3.1724069
##    3:    3  4.8616340 29.6643309  6.4130084
##    4:    4  0.5403350  0.6539857  0.2711737
##    5:    5 12.4790275  3.3807058 14.2851746
##   ---                                      
##  996:  996  9.0188944  9.5958862  6.9818731
##  997:  997  1.0142999  2.5976399 13.4244831
##  998:  998  0.3335885  5.6766551  3.0571633
##  999:  999  3.4507202  2.6912672 16.6274046
## 1000: 1000 12.6529929  8.2262360 18.7148794

round(cor(as.matrix(dx[, .(a, b, c)])), 2)
##      a    b    c
## a 1.00 0.67 0.67
## b 0.67 1.00 0.67
## c 0.67 0.67 1.00

# These data sets can be generated in either wide or long form. 
# So far, we have generated wide form data, where there is one row per unique id. 
# Now, we will generate data using the long form, where the correlated data are on different rows,
# so that there are repeated measurements for each id. 
# An id will have multiple records (i.e. one id will appear on multiple rows):

dx <- genCorGen(1000, nvars = 3, params1 = l, params2 = c(1, 1, 1), dist = "gamma", 
                rho = 0.7, corstr = "cs", wide = FALSE, cnames = "NewCol")
dx

##         id period   NewCol
##    1:    1      0 4.618194
##    2:    1      1 6.364941
##    3:    1      2 4.994391
##    4:    2      0 1.616983
##    5:    2      1 2.628199
##   ---                     
## 2996:  999      1 5.019716
## 2997:  999      2 7.262344
## 2998: 1000      0 9.557948
## 2999: 1000      1 5.309435
## 3000: 1000      2 2.403630


round(cor(as.matrix(dx[, .(V1, V2, V3)])), 2)
##      V1   V2   V3
## V1 1.00 0.32 0.32
## V2 0.32 1.00 0.29
## V3 0.32 0.29 1.00
We can also generate correlated binary data by specifying the probabilities:
  
  genCorGen(1000, nvars = 3, params1 = c(0.3, 0.5, 0.7), dist = "binary", rho = 0.8, 
            corstr = "cs", wide = TRUE)
##         id V1 V2 V3
##    1:    1  1  1  1
##    2:    2  0  1  1
##    3:    3  0  1  1
##    4:    4  0  0  0
##    5:    5  0  0  1
##   ---              
##  996:  996  0  1  1
##  997:  997  1  1  1
##  998:  998  1  1  1
##  999:  999  0  1  1
## 1000: 1000  0  0  0
The gamma distribution requires two parameters - the mean and dispersion. (These are converted into shape and rate parameters more commonly used.)

dx <- genCorGen(1000, nvars = 3, params1 = l, params2 = c(1, 1, 1), dist = "gamma", 
                rho = 0.7, corstr = "cs", wide = TRUE, cnames = "a, b, c")
dx
##         id          a          b          c
##    1:    1  2.2147914  5.7828970 12.4494947
##    2:    2  1.3062344  6.1613218  3.1724069
##    3:    3  4.8616340 29.6643309  6.4130084
##    4:    4  0.5403350  0.6539857  0.2711737
##    5:    5 12.4790275  3.3807058 14.2851746
##   ---                                      
##  996:  996  9.0188944  9.5958862  6.9818731
##  997:  997  1.0142999  2.5976399 13.4244831
##  998:  998  0.3335885  5.6766551  3.0571633
##  999:  999  3.4507202  2.6912672 16.6274046
## 1000: 1000 12.6529929  8.2262360 18.7148794
round(cor(as.matrix(dx[, .(a, b, c)])), 2)
##      a    b    c
## a 1.00 0.67 0.67
## b 0.67 1.00 0.67
## c 0.67 0.67 1.00
These data sets can be generated in either wide or long form. So far, we have generated wide form data, where there is one row per unique id. Now, we will generate data using the long form, where the correlated data are on different rows, so that there are repeated measurements for each id. An id will have multiple records (i.e. one id will appear on multiple rows):
  
  dx <- genCorGen(1000, nvars = 3, params1 = l, params2 = c(1, 1, 1), dist = "gamma", 
                  rho = 0.7, corstr = "cs", wide = FALSE, cnames = "NewCol")
dx
##         id period   NewCol
##    1:    1      0 4.618194
##    2:    1      1 6.364941
##    3:    1      2 4.994391
##    4:    2      0 1.616983
##    5:    2      1 2.628199
##   ---                     
## 2996:  999      1 5.019716
## 2997:  999      2 7.262344
## 2998: 1000      0 9.557948
## 2999: 1000      1 5.309435
## 3000: 1000      2 2.403630


# addCorGen allows us to create correlated data from an existing data set, 
# as one can already do using addCorData. In the case of addCorGen, the parameter(s) 
# used to define the distribution are created as a field (or fields) in the dataset. 
# The correlated data are added to the existing data set. In the example below,
# we are going to generate three sets (poisson, binary, and gamma) of correlated data with means 
# that are a function of the variable xbase, which varies by id.

#First we define the data and generate a data set:
  
def <- defData(varname = "xbase", formula = 5, variance = 0.2, dist = "gamma", id = "cid")
def <- defData(def, varname = "lambda", formula = ".5 + .1*xbase", dist = "nonrandom", link = "log")
def <- defData(def, varname = "p", formula = "-2 + .3*xbase", dist = "nonrandom", link = "logit")
def <- defData(def, varname = "gammaMu", formula = ".5 + .2*xbase", dist = "nonrandom",link = "log")
def <- defData(def, varname = "gammaDis", formula = 1, dist = "nonrandom")

dt <- genData(10000, def)
dt
##          cid    xbase   lambda         p  gammaMu gammaDis
##     1:     1 5.698630 2.914980 0.4279032 5.153757        1
##     2:     2 3.564564 2.354802 0.2827968 3.363267        1
##     3:     3 2.664886 2.152196 0.2313802 2.809418        1
##     4:     4 5.678577 2.909141 0.4264312 5.133129        1
##     5:     5 5.369445 2.820585 0.4039179 4.825377        1
##    ---                                                    
##  9996:  9996 8.325959 3.790871 0.6219393 8.716273        1
##  9997:  9997 2.559089 2.129546 0.2257838 2.750596        1
##  9998:  9998 3.525356 2.345587 0.2804172 3.336997        1
##  9999:  9999 2.417845 2.099679 0.2184629 2.673982        1
## 10000: 10000 3.060428 2.239030 0.2531520 3.040694        1

# The Poisson distribution has a single parameter, lambda:
  
dtX1 <- addCorGen(dtOld = dt, idvar = "cid", nvars = 3, rho = 0.1, corstr = "cs", 
                    dist = "poisson", param1 = "lambda", cnames = "a, b, c")
dtX1
##          cid    xbase   lambda         p  gammaMu gammaDis a b c
##     1:     1 5.698630 2.914980 0.4279032 5.153757        1 6 1 2
##     2:     2 3.564564 2.354802 0.2827968 3.363267        1 3 3 3
##     3:     3 2.664886 2.152196 0.2313802 2.809418        1 2 1 3
##     4:     4 5.678577 2.909141 0.4264312 5.133129        1 3 1 4
##     5:     5 5.369445 2.820585 0.4039179 4.825377        1 2 2 3
##    ---                                                          
##  9996:  9996 8.325959 3.790871 0.6219393 8.716273        1 4 4 6
##  9997:  9997 2.559089 2.129546 0.2257838 2.750596        1 1 3 3
##  9998:  9998 3.525356 2.345587 0.2804172 3.336997        1 0 2 1
##  9999:  9999 2.417845 2.099679 0.2184629 2.673982        1 1 2 2
## 10000: 10000 3.060428 2.239030 0.2531520 3.040694        1 3 0 2

# The Bernoulli (binary) distribution has a single parameter, p:
  
dtX2 <- addCorGen(dtOld = dt, idvar = "cid", nvars = 4, rho = 0.4, corstr = "ar1", 
                    dist = "binary", param1 = "p")
dtX2
##          cid    xbase   lambda         p  gammaMu gammaDis V1 V2 V3 V4
##     1:     1 5.698630 2.914980 0.4279032 5.153757        1  0  0  1  0
##     2:     2 3.564564 2.354802 0.2827968 3.363267        1  1  1  0  1
##     3:     3 2.664886 2.152196 0.2313802 2.809418        1  0  0  0  0
##     4:     4 5.678577 2.909141 0.4264312 5.133129        1  1  0  0  0
##     5:     5 5.369445 2.820585 0.4039179 4.825377        1  1  1  1  1
##    ---                                                                
##  9996:  9996 8.325959 3.790871 0.6219393 8.716273        1  1  0  0  1
##  9997:  9997 2.559089 2.129546 0.2257838 2.750596        1  0  0  1  1
##  9998:  9998 3.525356 2.345587 0.2804172 3.336997        1  0  0  1  1
##  9999:  9999 2.417845 2.099679 0.2184629 2.673982        1  1  0  0  1
## 10000: 10000 3.060428 2.239030 0.2531520 3.040694        1  0  0  0  0

# The Gamma distribution has two parameters - in simstudy the mean and dispersion are specified:
  
dtX3 <- addCorGen(dtOld = dt, idvar = "cid", nvars = 4, rho = 0.4, corstr = "cs", 
                    dist = "gamma", param1 = "gammaMu", param2 = "gammaDis")
dtX3
##          cid    xbase   lambda         p  gammaMu gammaDis        V1
##     1:     1 5.698630 2.914980 0.4279032 5.153757        1 2.1924081
##     2:     2 3.564564 2.354802 0.2827968 3.363267        1 8.7328677
##     3:     3 2.664886 2.152196 0.2313802 2.809418        1 2.7422032
##     4:     4 5.678577 2.909141 0.4264312 5.133129        1 4.3125713
##     5:     5 5.369445 2.820585 0.4039179 4.825377        1 4.6422097
##    ---                                                              
##  9996:  9996 8.325959 3.790871 0.6219393 8.716273        1 9.9600065
##  9997:  9997 2.559089 2.129546 0.2257838 2.750596        1 2.3437735
##  9998:  9998 3.525356 2.345587 0.2804172 3.336997        1 0.6124997
##  9999:  9999 2.417845 2.099679 0.2184629 2.673982        1 5.5975817
## 10000: 10000 3.060428 2.239030 0.2531520 3.040694        1 1.5432414
##               V2         V3        V4
##     1: 3.1822804 11.8429426 3.5030663
##     2: 0.8551908  8.4568344 3.7027966
##     3: 1.4664891  9.1577855 3.2776414
##     4: 2.5736640  1.9366756 1.6897873
##     5: 2.9348884  1.7149849 2.0216163
##    ---                               
##  9996: 1.7654044 13.4982386 4.2711241
##  9997: 4.6065786  4.6201099 1.7846534
##  9998: 1.3123053  0.7748302 1.9384933
##  9999: 4.4785144  1.3585290 5.6899324
## 10000: 3.8009607  0.1165655 0.5052408


# If we have data in long form (e.g. longitudinal data), the function will recognize the structure:
  
def <- defData(varname = "xbase", formula = 5, variance = 0.4, dist = "gamma", id = "cid")
def <- defData(def, "nperiods", formula = 3, dist = "noZeroPoisson")
def2 <- defDataAdd(varname = "lambda", formula = ".5+.5*period + .1*xbase", dist = "nonrandom", link = "log")

dt <- genData(1000, def)

dtLong <- addPeriods(dt, idvars = "cid", nPeriods = 3)
dtLong <- addColumns(def2, dtLong)

dtLong
##        cid period    xbase nperiods timeID   lambda
##    1:    1      0 4.419360        1      1 2.564942
##    2:    1      1 4.419360        1      2 4.228875
##    3:    1      2 4.419360        1      3 6.972236
##    4:    2      0 5.941155        2      4 2.986540
##    5:    2      1 5.941155        2      5 4.923972
##   ---                                              
## 2996:  999      1 3.365393        3   2996 3.805850
## 2997:  999      2 3.365393        3   2997 6.274785
## 2998: 1000      0 4.490924        5   2998 2.583364
## 2999: 1000      1 4.490924        5   2999 4.259247
## 3000: 1000      2 4.490924        5   3000 7.022311
### Generate the data

dtX3 <- addCorGen(dtOld = dtLong, idvar = "cid", nvars = 3, rho = 0.6, corstr = "cs", 
                  dist = "poisson", param1 = "lambda", cnames = "NewPois")
dtX3
##        cid period    xbase nperiods timeID   lambda NewPois
##    1:    1      0 4.419360        1      1 2.564942       0
##    2:    1      1 4.419360        1      2 4.228875       3
##    3:    1      2 4.419360        1      3 6.972236       4
##    4:    2      0 5.941155        2      4 2.986540       3
##    5:    2      1 5.941155        2      5 4.923972       4
##   ---                                                      
## 2996:  999      1 3.365393        3   2996 3.805850       7
## 2997:  999      2 3.365393        3   2997 6.274785       9
## 2998: 1000      0 4.490924        5   2998 2.583364       3
## 2999: 1000      1 4.490924        5   2999 4.259247       6
## 3000: 1000      2 4.490924        5   3000 7.022311      10


## APPLICATION FOR TWIN STUDY CONTEXT
# twins are same age, Mz twins are assumed to be same sex; Dz may vary (but modelled seperately
# csv copy of HeartExampleStata.data
heart <- read.csv("C:/Users/Carl/OneDrive/Research/2 - BCA/Research project/bca_rp2/scripts/mbiostatprojectbackgroundreading/HeartExample.csv")
# create compound indicators of mz (2), dz1 (0) or dz2 (1)  -- ie. to pick up on individ' env. re
heart["mz_ind"]<- heart$mz+1-heart$dz1
# run linear mixed effects model
m <- list()
m[["bmi"]] <- lmer(bmi ~ sex + (1|pairid) + (1|pairid:mz_ind), data = heart)
m[["lnsbp"]] <- lmer(lnsbp ~ sex + bmi + (1|pairid) + (1|pairid:mz_ind), data = heart)
summary(m[["lnsbp"]]) # Rough equivalent of stata model: mixed lnsbp sex bmi, || pairid: || pairid: mz dz1 dz2, covariance(identity) nocons
# Linear mixed model fit by REML ['lmerMod']
# Formula: lnsbp ~ sex + bmi + (1 | pairid) + (1 | pairid:mz_ind)
# Data: heart
# 
# REML criterion at convergence: -666.1
# 
# Scaled residuals: 
#   Min       1Q   Median       3Q      Max 
# -2.03202 -0.50601 -0.01025  0.49651  2.08075 
# 
# Random effects:
#   Groups        Name        Variance Std.Dev.
# pairid:mz_ind (Intercept) 0.002165 0.04653 
# pairid        (Intercept) 0.001132 0.03364 
# Residual                  0.002967 0.05447 
# Number of obs: 298, groups:  pairid:mz_ind, 232; pairid, 149
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept) 4.592482   0.035140  130.69
# sex         0.066597   0.010434    6.38
# bmi         0.005725   0.001582    3.62
# 
# Correlation of Fixed Effects:
#   (Intr) sex   
# sex  0.176       
# bmi -0.982 -0.294
summary <- list()
for (var in c("bmi","lnsbp")){
  summary[[var]][["summary"]]         <- summary(m[[var]])
  summary[[var]][["fe"]]              <- summary[[var]][["summary"]][["coefficients"]]
  summary[[var]][["re"]]              <- as.data.frame(VarCorr(m[[var]]))
  summary[[var]][["A"]]               <- summary[[var]][["re"]][1,4] # additional Mz covariance
  summary[[var]][["C"]]               <- summary[[var]][["re"]][2,4] # Dz covariance
  summary[[var]][["E"]]               <- summary[[var]][["re"]][3,4] # residual variance
  summary[[var]][["total_var"]]       <- sum(unlist(summary[[var]][c("A","C","E")]))
  summary[[var]][["cov"]][["mz"]]     <- summary[[var]]$C+summary[[var]]$A
  summary[[var]][["cov"]][["dz"]]     <- summary[[var]]$C
  summary[[var]][["rho"]][["mz"]]     <- summary[[var]][["cov"]]["mz"]/summary[[var]]$total_var
  summary[[var]][["rho"]][["dz"]]     <- summary[[var]][["cov"]]["dz"]/summary[[var]]$total_var
}
cat("Intraclass Correlations, ")
matrix(c(summary$bmi$rho,summary$lnsbp$rho),
       ncol=2,
       dimnames=list(c("Mz","Dz"),c("BMI","lnsbp")))


n.pairs = 149

# Simulate data, drawing on mixed effects models as guide for data generating mechanism
tw.pair <- defData(varname = "n.inds", formula = 2, variance = 0, id = "id.pair") 
tw.pair <- defData(tw.pair,varname = "Mz", dist = "binary", formula = 0.3, id = "id.pair")
tw.pair <- defData(tw.pair,varname = "z0_mz_bmi",   dist = "normal", formula = 0, variance = summary$bmi$cov["mz"],   id = "id.pair")
tw.pair <- defData(tw.pair,varname = "z0_mz_lnsbp", dist = "normal", formula = 0, variance = summary$lnsbp$cov["mz"], id = "id.pair")
tw.pair <- defData(tw.pair,varname = "z0_dz_bmi",   dist = "normal", formula = 0, variance = summary$bmi$cov["dz"],   id = "id.pair")
tw.pair <- defData(tw.pair,varname = "z0_dz_lnsbp", dist = "normal", formula = 0, variance = summary$lnsbp$cov["dz"], id = "id.pair")
tw.pair <- defData(tw.pair,varname = "age", dist = "uniform", formula = "18; 31", id = "id.pair")
tw.pair <- defData(tw.pair,varname = "mz_fem", dist = "binary", formula = 0.5, id = "id.pair")
tw.pair
# varname formula    variance    dist     link
# 1:      n.inds       2 0.000000000  normal identity
# 2:          Mz     0.3 0.000000000  binary identity
# 3:   z0_mz_bmi       0 7.064280080  normal identity
# 4: z0_mz_lnsbp       0 0.003296820  normal identity
# 5:   z0_dz_bmi       0 4.714199627  normal identity
# 6: z0_dz_lnsbp       0 0.001131618  normal identity
# 7:         age  18; 31 0.000000000 uniform identity
# 8:      mz_fem     0.5 0.000000000  binary identity
dt.pair <- genData(n.pairs, tw.pair)
dt  .pair
#      id.pair n.inds Mz  z0_mz_bmi  z0_mz_lnsbp  z0_dz_bmi  z0_dz_lnsbp      age mz_fem
#   1:       1      2  0 -0.6574000  0.064558889 -3.6999346  0.014983375 28.77778      1
#   2:       2      2  0 -1.8892245  0.016784464 -3.0159563 -0.011469922 19.82292      0
#   3:       3      2  1 -1.9766545  0.018124157  0.5023240  0.020114604 28.39733      1
#   4:       4      2  0  0.2970619  0.046058288 -2.4671660  0.040322515 28.64035      0
#   5:       5      2  1  4.5812000 -0.032586016  0.3818789  0.027171521 22.80382      0
# ---                                                                                  
# 145:     145      2  0 -1.8933768  0.021241787 -0.7393719 -0.018024633 23.25968      0
# 146:     146      2  1 -3.0230712  0.070729052 -3.1220987 -0.002867279 26.87104      1
# 147:     147      2  0 -0.3361081  0.102640413  0.8698970  0.016268089 28.48500      0
# 148:     148      2  0  0.7751230  0.009452269 -1.2566928  0.037654063 28.00665      1
# 149:     149      2  0 -6.6343096  0.056883598 -2.0525133  0.011597051 21.54519      1

tw.ind <- defDataAdd(varname = "dz_fem", dist = "binary", formula = 0.5) 
tw.ind <- defDataAdd(tw.ind,varname = "bmi", dist = "normal",  
                     formula = paste0(summary$bmi$fe["(Intercept)","Estimate"]," + ",
                                      summary$bmi$fe["sex","Estimate"],"*((mz_fem*Mz)+(dz_fem*(abs(1-Mz)))) + 
                                      z0_mz_bmi * Mz + z0_dz_bmi * (abs(1-Mz))"), 
                     variance = summary$bmi$E)
tw.ind <- defDataAdd(tw.ind,varname = "lnsbp", dist = "normal",  
                     formula = paste0(summary$lnsbp$fe["(Intercept)","Estimate"]," + ",
                                      summary$lnsbp$fe["sex","Estimate"],"*((mz_fem*Mz)+(dz_fem*(abs(1-Mz)))) +",
                                      summary$lnsbp$fe["bmi","Estimate"],
                                      " * bmi + z0_mz_lnsbp*Mz + z0_dz_lnsbp*(abs(1-Mz))"), 
                     variance = summary$lnsbp$E)

dt.tw.ind <- genCluster(dt.pair, cLevelVar = "id.pair", numIndsVar = "n.inds", level1ID = "id.ind")
dt.tw.ind  <- addColumns(tw.ind, dt.tw.ind)

# create within pair index
# Should be more simple way of doing this; in lieu of knowing a better way, 
# it is calculated as sum of 1 plus 1 if previous value is identical, else 1 + 0
dt.tw.ind[ ,id.twin:= 1+replace((shift(id.pair)==id.pair), is.na(shift(id.pair)), 0)]
# alternately, assuming correctly sorted data: dt.tw.ind["id.twin"] <- 1 + (dt.tw.ind[ ,"id.ind"]%% 2 == 0) 

dt.tw.ind[(dt.tw.ind$Mz==1),"female"] <- dt.tw.ind[(dt.tw.ind$Mz==1),"mz_fem"]
dt.tw.ind[(dt.tw.ind$Mz==0),"female"] <- dt.tw.ind[(dt.tw.ind$Mz==0),"dz_fem"]
dt.tw.ind <- dt.tw.ind[,c("id.pair","id.twin","Mz","age","female","bmi","lnsbp")]
head(dt.tw.ind,10)


# Pearson correlation
cor(dt.tw.ind[dt.tw.ind$Mz==1&dt.tw.ind$id.twin==1,"bmi"],dt.tw.ind[dt.tw.ind$Mz==1&dt.tw.ind$id.twin==2,"bmi"])
# bmi
# bmi 0.7359233
cor(dt.tw.ind[dt.tw.ind$Mz==0&dt.tw.ind$id.twin==1,"bmi"],dt.tw.ind[dt.tw.ind$Mz==0&dt.tw.ind$id.twin==2,"bmi"])
# bmi
# bmi 0.6420553


# Correlated data example
nsim = 1000
# a vector containing the mean(s) of the distributions (for binary, should be prob)
mu <- c(0.8, 0.1)  
# given dist, a parameterisation (e.g. rate [poisson], dispersion [gamma], var [normal] or max [uniform])
# For binary data, set to NULL
param <- c(NULL, NULL) 
distribution = "binary"  # "binary", "poisson" or "gamma", "normal", or "uniform".
rho = 0.78
genCorGen(nsim, nvars = 2, params1 = mu, params2 = param,  dist = distribution, rho = rho,
          corstr = "cs", wide = TRUE)
