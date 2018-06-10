# R script to simulate power for difference in correlations (Pearson, Spearman, and later ... ICC)
# Carl Higgs 2017


# Function to simulate two bivariate normal distributions based on respective population correlations.
# Reports:
#    - the significance of the sample correlation's difference (output 'z_p') 
#    - the power to detect a difference (output 'z_power') given alpha, beta and sidedness parameters.
# Optionally
#    - calculate the power or significance of difference of Pearson, Spearman or Kendall correlations
#    - calculate the significance and power without sampling using the two input correlations (simulation = FALSE).
#    - output a full log of simulation paramaters with results
#    - output a single statistic (z_p, or z_power).
#    - distribution can be one of "poisson", "binary", "gamma", "uniform", "negbinom", "normal"
#    - given dist, a parameterisation (e.g. rate [poisson], dispersion [gamma], var [normal] or max [uniform])
#    - allow parameterisation to explore change in source pop mean and sd (unequal variance, etc)
# To develop
#    - accounting for clustering (ie. icc in twin studies)

# install.packages("Rcpp")
# install.packages("simstudy")
# install.packages("ggplot2")
# install.packages("data.table")
# install.packages('RPostgres')

# may have to set working directory
setwd("C:/Users/Carl/OneDrive/Research/2 - BCA/Research project/bca_rp2/scripts")
require("simstudy") # simstudy is used to generate bivariate data
require(compiler)   # for byte code compilation
require(Rcpp)       # for alternate compilation approach, using C++
sourceCpp('test.cpp') # a C++ script, with more efficient number generators
require(data.table) # Results are stored using data.table
require(parallel) # For parallel processing
require(DBI)  # used to connext to Postgresql using RPostgres
require(config) # used to access config.yml file with Postgres connection parameters

# to deploy R power app (code elsewhere

# to deploy (in R < v3.5)
# library(rsconnect)
# deployApp(appDir = "C:/Users/Carl/OneDrive/Research/2 - BCA/Research project/bca_rp2/scripts/corr_power_app")


## Define tests
# Fishers Z test - no sim
fz_nosim <- function(r1,r2,n1,n2,
                     alpha = 0.05, sidedness=2,method = "pearson",
                     power = TRUE) {
  # Step 1: calculate Fisher's Z
  z1     <- atanh(r1)
  z2     <- atanh(r2)
  # Step 2: take difference
  zdiff  <- z1-z2
  # Step 3: calculate standard error and test statistic
  z_se   <- sqrt(1/(n1-3) + 1/(n2-3))
  z_test <- zdiff/z_se
  # optionally return p-value for observing diff at least this large under H0
  z_p    <- sidedness*pnorm(-abs(z_test))
  if (power == FALSE) return("p" = z_p)
  z_ref   <- qnorm(1-alpha/sidedness)
  z_power <- 1-pnorm(z_ref - abs(z_test))
  return(z_power)
}


# Fishers Z test
fz <- function(a,b,sidedness=2,method = "pearson") {
  # Two samples
  n1 <- nrow(a)
  n2 <- nrow(b)
   
  # Step 1: Compute sample correlation coefficients
  z1     <- atanh(cor(a,method = method)[2,1])
  z2     <- atanh(cor(b,method = method)[2,1])
  zdiff  <- z1-z2
  z_se   <- sqrt(1/(n1-3) + 1/(n2-3))
  z_test <- zdiff/z_se
  z_p    <- sidedness*pnorm(-abs(z_test))
  # return(c(z_p,z_test))
  return(z_p)
}

## Using rccp - but this doesn't play nice with parallelisation
gtv <- function(a,b,M=1e4,method = "pearson") {
  # Two samples
  n1 <- nrow(a)
  n2 <- nrow(b)

  # Step 1: Compute sample correlation coefficients
  r1 <- cor(a,method = method)[2,1]
  r2 <- cor(b,method = method)[2,1]
  r  <- c(r1,r2)

  # Step 2: Generate random numbers
  V2     <- matrix(data=0, nrow = M, ncol = 2)
  V2[,1] <- rchisq_cpp(M, df = n1-1, ncp = 1)
  V2[,2] <- rchisq_cpp(M, df = n2-1, ncp = 1)

  W2     <- matrix(data=0, nrow = M, ncol = 2)
  W2[,1] <- rchisq_cpp(M, df = n1-2, ncp = 1)
  W2[,2] <- rchisq_cpp(M, df = n2-2, ncp = 1)

  Z <-matrix(data = rnorm_cpp(2*M,0,1), nrow=M, ncol = 2)

  # Compute test statistic
  rstar <- r/sqrt(1-r^2)
  top   <- c(sqrt(W2[,1])*rstar[1],sqrt(W2[,2])*rstar[2]) - Z
  G     <- top / sqrt( top^2 + V2 )

  # Compute p value
  Grho <- G[,1] - G[,2];
  p    <- 2*min( mean(Grho<0), mean(Grho>0) );
  return(p)
}

# GTV test statistic, redefined not using rccp
gtv_r <- function(a,b,M=1e4,method = "pearson") {
  # Two samples
  n1 <- nrow(a)
  n2 <- nrow(b)
  
  # Step 1: Compute sample correlation coefficients
  r1 <- cor(a,method = method)[2,1]
  r2 <- cor(b,method = method)[2,1]
  r  <- c(r1,r2)
  
  # Step 2: Generate random numbers
  V2     <- matrix(data=0, nrow = M, ncol = 2)
  V2[,1] <- rchisq(M, df = n1-1, ncp = 1)
  V2[,2] <- rchisq(M, df = n2-1, ncp = 1)
  
  W2     <- matrix(data=0, nrow = M, ncol = 2)
  W2[,1] <- rchisq(M, df = n1-2, ncp = 1)
  W2[,2] <- rchisq(M, df = n2-2, ncp = 1)
  
  Z <-matrix(data = rnorm(2*M), nrow=M, ncol = 2)
  
  # Compute test statistic
  rstar <- r/sqrt(1-r^2)
  top   <- c(sqrt(W2[,1])*rstar[1],sqrt(W2[,2])*rstar[2]) - Z
  G     <- top / sqrt( top^2 + V2 )
  
  # Compute p value
  Grho <- G[,1] - G[,2];
  p    <- 2*min( mean(Grho<0), mean(Grho>0) ); 
  return(p)
}

# SIgned Log-likelihood Ratio test (an 'unmodified' version of test described in Kazemi and Jafari / DiCiccio)
slr <- function(a,b,M=1e4,sidedness=2,method = "pearson") {
  # Two samples
  n  <- c(nrow(a),nrow(b))
  
  # Step 1: Compute sample correlation coefficients
  r  <- c(cor(a,method = method)[2,1], cor(b,method = method)[2,1])
  z  <- atanh(r)
  rf <- tanh(mean(z))
  
  slr <-sign(r[1]-r[2])*sqrt(sum(n*log(((1-rf*r)^2)/((1-r^2)*(1-rf^2)))))
  p    <- 2 * (1 - pnorm(abs(slr))); 
  return(p)
}


pt <- function(a,b,M=1e4,sidedness=2,method = "pearson") {
# Based on Efron and Tibshirani
n  <- c(nrow(a),nrow(b))
r  <- c(cor(a,method = method)[2,1], cor(b,method = method)[2,1])
z  <- atanh(r)
v <- cbind(rank(rbind(a[,1],b[,1]),ties.method = "random"),rank(rbind(a[,2],b[,2]),ties.method = "random"))
rownames(v) <- c(rep("A",n[1]),rep("B",n[2]))
rtest <- numeric(0)
for (i in 1:M){
  permute <- cbind(v,rbinom(sum(n),1,0.5))
  rstar   <- c(cor(permute[permute[,3]==0,c(1,2)],method = method)[2,1],
               cor(permute[permute[,3]==1,c(1,2)],method = method)[2,1])
  zstar   <- atanh(rstar)
  rtest   <- c(rtest,
               abs(zstar[1]-zstar[2]) > abs(z[1]-z[2]))
  } 
p <- mean(rtest)
return(p)
}

zou <- function(a,b,alpha = 0.05,sidedness=2,method = "pearson") {
  # From Zou (2007) and used in Cocor (note typo for U in paper; should be '+')
  #  However, really, this is equivalent to fz test for hypothesis testing purposes
  n  <- c(nrow(a),nrow(b))
  r  <- c(cor(a,method = method)[2,1], cor(b,method = method)[2,1])
  z  <- atanh(r)
  zdiff  <- z[1]-z[2]
  # Step 3: calculate standard error and test statistic
  z_ref  <- qnorm(1-alpha/sidedness)
  z_se   <- sqrt(1/(n-3))
  ci_mat <- matrix(c(-1,-1,1,1),nrow = 2, ncol = 2, dimnames =list(c("Mz","Dz"),c("l","u")))
  z_ci   <- z + ci_mat * z_se * z_ref
  r_ci   <- tanh(z_ci)
  L      <- r[1]-r[2] - sqrt((r[1]      - r_ci[1,1])^2 + (r_ci[2,2] - r[2]     )^2)
  U      <- r[1]-r[2] + sqrt((r_ci[1,2] - r[1]     )^2 + (r[2]      - r_ci[2,1])^2)
  r_diff_ci <- c(L,U)
  ci_test <- (L < 0) && (0 < U)
  return(c(ci_test,r_diff_ci))
}

# compile tests, noting that GTV function is already compiled using rccp (gtv_r is straight r)
fz_ns_compiled <- cmpfun(fz_nosim)
fz_compiled <- cmpfun(fz)
gtv_compiled<- cmpfun(gtv_r)
slr_compiled <- cmpfun(slr)
pt_compiled <- cmpfun(pt)
zou_compiled <- cmpfun(zou)

corr_diff_test <- function(rho = c(.2,.5), n = c(30,90), distr = "normal",
                    param1a = c(0,0), param1b = c(0,0),param2a = c(1,1), param2b = c(1,1),
                    alpha = 0.05, sidedness = 2, test = c("fz","gtv","pt","slr","zou"),
                    method ="pearson", lower.tri = FALSE) {
  # cat("Parameters: ",rho[1],rho[2], n[1],n[2], param1a, param1b, param2a, param2b, distr, alpha, sidedness, method,"\n")
  if(lower.tri==TRUE){
    # only calculate lower matrix half when comparing across all correlation combinations
    if(rho[1] < rho[2]) { 
      return(NA)
    }
  }
  results <- list()
  if ("fz_nosim" %in% test) {
    results[["fz_nosim"]] <- fz_ns_compiled(rho[1],rho[2],n[1],n[2], 
                                      alpha = 0.05, sidedness = 2, method = method, power = FALSE)
    if(length(test)==1) return(results)
  }
  require("simstudy")
  a <- genCorGen(n[1], nvars = 2, params1 = param1a, params2 = param2a,  
                dist = distr, corMatrix = matrix(c(1, rho[1], rho[1], 1), ncol = 2), 
                wide = TRUE)[,2:3]
  b <- genCorGen(n[2], nvars = 2, params1 = param1b, params2 = param2b,  
                dist = distr, corMatrix = matrix(c(1, rho[2], rho[2], 1), ncol = 2), 
                wide = TRUE)[,2:3]
  if ("fz"       %in% test) results[["fz"]]       <- fz_compiled(a,b)
  if ("gtv"      %in% test) results[["gtv"]]      <- gtv(a,b) # uses rccp ; so elsewise compiled
  if ("gtvr"     %in% test) results[["gtvr"]]      <- gtv_compiled(a,b) 
  if ("pt"       %in% test) results[["pt"]]       <- pt_compiled(a,b)
  if ("slr"      %in% test) results[["slr"]]      <- slr_compiled(a,b)
  if ("zou"      %in% test) results[["zou"]]      <- zou_compiled(a,b)[1]
  return(rbind(results[test]))
}

corr_diff_test_compiled <- cmpfun(corr_diff_test)

# Corr power simulation
corr_power <- function(rho = c(.2,.5), n = c(30,90),distr = "normal",
                       param1a = c(0,0), param1b = c(0,0),param2a = c(1,1), param2b = c(1,1),
                       test = c("fz","gtv","slr","zou"),
                       alpha = 0.05, sidedness=2,method="pearson",  
                       nsims = 100,lower.tri = FALSE, power_only = FALSE){
  sim <- list()
  sim[["params"]] <- c("method" = method, "rho_1" = rho[1], "rho_2" = rho[2],"n1" = n[1], "n2" = n[2],
                       "alpha" = alpha, "sidedness" = sidedness, "nsims" = nsims, "distr" = distr)
  sim[["additional"]] <- c("param1a" = param1a, "param1b" = param1b,"param2a" = param2a, "param2b" = param2b)
  if (length(test)==1){
    if ("fz_nosim" %in% test) {
      sim[["analytical"]] <- c("fz_nosim" = fz_nosim(rho[1],rho[2],n[1],n[2], 
                                      alpha = alpha, sidedness = sidedness, method = method, power = TRUE))
      if(power_only==TRUE) return(sim[["analytical"]])
    }
    if (!"fz_nosim" %in% test) {      
      sim_tests <- test
      sim[["power"]] <- mean(replicate(nsims,corr_diff_test(rho = rho, n = n,distr =distr,
                                                              param1a = param1a, param1b = param1b,param2a = param2a, param2b = param2b,
                                                              alpha = alpha, sidedness=sidedness,method=method, 
                                                              test = sim_tests)[,])<alpha)
      if(power_only==TRUE) return(sim[["power"]])
    }
  }
  if (length(test)>1){
    sim[["analytical"]] <- c("fz_nosim" = fz_nosim(rho[1],rho[2],n[1],n[2], 
                                                     alpha = alpha, sidedness = sidedness, method = method, power = TRUE))

    sim_tests <- test[!test %in% "fz_nosim"]
    if (length(sim_tests)==1) {      
      sim[["power"]] <- mean(replicate(nsims,corr_diff_test(rho = rho, n = n,distr =distr,
                                                            param1a = param1a, param1b = param1b,param2a = param2a, param2b = param2b,
                                                            alpha = alpha, sidedness=sidedness,method=method, 
                                                            test = sim_tests)[,])<alpha)
    }
    if (length(sim_tests)>1) {      
    sim[["power"]] <- rowMeans(replicate(nsims,corr_diff_test(rho = rho, n = n,distr =distr,
                                                                param1a = param1a, param1b = param1b,param2a = param2a, param2b = param2b,
                                                                alpha = alpha, sidedness=sidedness,method=method, 
                                                                test = sim_tests)[,])<alpha)
    }
    if ("fz_nosim" %in% test) sim[["power"]] <- c(sim[["analytical"]],sim[["power"]])[test]
    if(power_only==TRUE) return(sim[["power"]])
  }

  cat("\r                                                                \r",sim[["params"]],sim$power,sep="\t")
  flush.console()
  return(sim)
}

corr_power_compiled <- cmpfun(corr_power)


# corrs <- round(seq(-0.90,0.90,0.1),2)
# sim_n <- c(15,30,60,120,240,480,960)
# sim_groups <- expand.grid(sim_n,sim_n) 
# sim_groups["n"] <- rowSums(sim_groups)
# nrow(sim_groups)
# # 49 population size combinations
# sim_corrs  <- expand.grid(corrs,corrs)
# nrow(sim_corrs)
# # 1521 correlation combinations
# 
# group_params <- cbind((sim_groups[,1]/sim_groups[,2]),sim_groups[,3])
# ratios       <- sort(unique(group_params[,1]))
# numbers      <- sort(unique(group_params[,2]))

# Extract results into long dataframe
# > 2*3*5*49*1520
# [1] 2234400  --- actually that's not so bad

# but, we can have the tests in line, and we can amend the resolution (coarser corr seq)
# > 2*3*49*361
# [1] 106134

corrs <- round(seq(-0.90,0.90,0.1),2)
sim_n <- c(15,30,60,120,240,480,960) 
methods <- c("pearson","spearman")
# distributions <- c("normal_0_1","gamm_1.5_0.09","gamma_1_5")
distributions <- c("normal","gamma","gamma")
param1 <- c(0,1.5,1)
param2 <- c(1,0.09,5)
tests <- c( "fz_nosim","fz","gtv","slr","zou")

# construct simulation data table from parameter combinations
dt <- as.data.table(expand.grid(rho2 = corrs,
                    rho1 = corrs,
                    n2   = sim_n,
                    n1   = sim_n,
                    sdist = c(1:3),
                    method = methods))
dt[,("dist"):= distributions[sdist], by = 1:nrow(dt)]
dt[,("p1"):= param1[sdist], by = 1:nrow(dt)]
dt[,("p2"):= param2[sdist], by = 1:nrow(dt)]

setcolorder(dt,c("method","sdist","dist","p1","p2","n1","n2","rho1","rho2"))
dt <- dt[,c("method","dist","p1","p2","n1","n2","rho1","rho2")]


## Initial approach, calculating results within data table chunks using 100 simulations

system.time(dt[, c( "fz_nosim","fz","gtv","slr","zou"):={
  cat("\r",.GRP,"\r")
  as.list(corr_power_compiled(rho = c(rho1,rho2), 
                              n = c(n1,n2),
                              distr = dist,
                              param1a = c(p1,p1), 
                              param1b = c(p1,p1),
                              param2a = c(p2,p2), 
                              param2b = c(p2,p2),
                              test    = tests,
                              alpha   = 0.05, 
                              sidedness=2,
                              method=as.character(method),  
                              nsims = 100,
                              power_only = TRUE))
},by = 1:nrow(dt)] )

# add in extra summary vars
dt[,("ratio"):=n1/n2,by=1:nrow(dt)]
dt[,("n"):=n1+n2,by=1:nrow(dt)]
dt[,("diff"):=abs(rho1-rho2),by=1:nrow(dt)]


### Due to typo issue
## specifying:
#### tests <-  c( "fz_nosim","fz","gtv","slr","zou")
#### dt[,tests :={
## instead of
#### dt[, c( "fz_nosim","fz","gtv","slr","zou"):={
# The 1000 simulation only returned the first test result (no sim), despite processing all results
### Given time constraints I opted to run specific scenarios to process at higher levels of simulation

## How does ratio impact with Pearsons/Spearmans?
# Scenario set A1
#  n1 == 60 and n2 == 120
# nrow(dtk[(n1==60)&(n2==120),)
# 2166 rows
# A1 - 1,000 simulations
dt_sA1 <- dt[(n1==60)&(n2==120),c("method","dist","p1","p2","n1","n2","rho1","rho2")]
system.time(dt_sA1[, c( "fz_nosim","fz","gtv","slr","zou"):={
    cat("\r",.GRP,"\r")
    as.list(corr_power_compiled(rho = c(rho1,rho2),
                                n = c(n1,n2),
                                distr = dist,
                                param1a = c(p1,p1),
                                param1b = c(p1,p1),
                                param2a = c(p2,p2),
                                param2b = c(p2,p2),
                                test    = tests,
                                alpha   = 0.05,
                                sidedness=2,
                                method=as.character(method),
                                nsims = 1000,
                                power_only = TRUE))
  },by = 1:nrow(dt_sA1)] )
# Run on work computer (8 cores)
# user   system  elapsed 
# 48278.95   673.58 48780.27 


# A1 - 10,000 simulations
dt_sA1_10k <- dt[(n1==60)&(n2==120),c("method","dist","p1","p2","n1","n2","rho1","rho2")]
system.time(dt_sA1_10k[, c( "fz_nosim","fz","gtv","slr","zou"):={
  cat("\r",.GRP,"\r")
  as.list(corr_power_compiled(rho = c(rho1,rho2),
                              n = c(n1,n2),
                              distr = dist,
                              param1a = c(p1,p1),
                              param1b = c(p1,p1),
                              param2a = c(p2,p2),
                              param2b = c(p2,p2),
                              test    = tests,
                              alpha   = 0.05,
                              sidedness=2,
                              method=as.character(method),
                              nsims = 10000,
                              power_only = TRUE))
},by = 1:nrow(dt_sA1_10k)] )

# Scenario set A2 
# n1 == 90 and n2 == 90
# Note - this doesn't exist in the dt definition; will manually redefine n1 and n2

# A2 1,000 simulations
dt_sA2 <-  dt[(n1==60)&(n2==120),c("method","dist","p1","p2","n1","n2","rho1","rho2")]
dt_sA2$n1 <- 90
dt_sA2$n2 <- 90
system.time(dt_sA2[, c( "fz_nosim","fz","gtv","slr","zou"):={
  cat("\r",.GRP,"\r")
  as.list(corr_power_compiled(rho = c(rho1,rho2),
                              n = c(n1,n2),
                              distr = dist,
                              param1a = c(p1,p1),
                              param1b = c(p1,p1),
                              param2a = c(p2,p2),
                              param2b = c(p2,p2),
                              test    = tests,
                              alpha   = 0.05,
                              sidedness=2,
                              method=as.character(method),
                              nsims = 1000,
                              power_only = TRUE))
},by = 1:nrow(dt_sA2)] )
# on work computer
# user   system  elapsed 
# 48350.98   736.72 48922.11 

# A2 10,000 simulations
dt_sA2_10k <-  dt[(n1==60)&(n2==120),c("method","dist","p1","p2","n1","n2","rho1","rho2")]
dt_sA2_10k$n1 <- 90
dt_sA2_10k$n2 <- 90
system.time(dt_sA2_10k[, c( "fz_nosim","fz","gtv","slr","zou"):={
  cat("\r",.GRP,"\r")
  as.list(corr_power_compiled(rho = c(rho1,rho2),
                              n = c(n1,n2),
                              distr = dist,
                              param1a = c(p1,p1),
                              param1b = c(p1,p1),
                              param2a = c(p2,p2),
                              param2b = c(p2,p2),
                              test    = tests,
                              alpha   = 0.05,
                              sidedness=2,
                              method=as.character(method),
                              nsims = 10000,
                              power_only = TRUE))
},by = 1:nrow(dt_sA2_10k)] )


## How does sample size impact power estimates?
# Scenario set B1
# rho1 == 0.2 && rho2 == 0.5
# nrow(dtk[(rho1==0.2)&(rho2==0.5),])
# 294 rows

# B 1,000 simulations
dt_sA3 <-  dt[(rho1==0.2)&(rho2==0.5),c("method","dist","p1","p2","n1","n2","rho1","rho2")]

system.time(dt_sA3[, c( "fz_nosim","fz","gtv","slr","zou"):={
  cat("\r",.GRP,"\r")
  as.list(corr_power_compiled(rho = c(rho1,rho2),
                              n = c(n1,n2),
                              distr = dist,
                              param1a = c(p1,p1),
                              param1b = c(p1,p1),
                              param2a = c(p2,p2),
                              param2b = c(p2,p2),
                              test    = tests,
                              alpha   = 0.05,
                              sidedness=2,
                              method=as.character(method),
                              nsims = 1000,
                              power_only = TRUE))
},by = 1:nrow(dt_sA3)] )
# at work
# user  system elapsed 
# 6834.87   92.21 6921.39 

# B 10,000 simulations
dt_sA3_10k <-  dt[(rho1==0.2)&(rho2==0.5),c("method","dist","p1","p2","n1","n2","rho1","rho2")]
system.time(dt_sA3_10k[, c( "fz_nosim","fz","gtv","slr","zou"):={
  cat("\r",.GRP,"\r")
  as.list(corr_power_compiled(rho = c(rho1,rho2),
                              n = c(n1,n2),
                              distr = dist,
                              param1a = c(p1,p1),
                              param1b = c(p1,p1),
                              param2a = c(p2,p2),
                              param2b = c(p2,p2),
                              test    = tests,
                              alpha   = 0.05,
                              sidedness=2,
                              method=as.character(method),
                              nsims = 10000,
                              power_only = TRUE))
},by = 1:nrow(dt_sA3_10k)] )


# additional summary statistics, for each specific scenario
dt_sA1[,("ratio"):=n1/n2,by=1:nrow(dt_sA1)]
dt_sA1[,("n"):=n1+n2,by=1:nrow(dt_sA1)]
dt_sA1[,("diff"):=abs(rho1-rho2),by=1:nrow(dt_sA1)]

dt_sA2[,("ratio"):=n1/n2,by=1:nrow(dt_sA2)]
dt_sA2[,("n"):=n1+n2,by=1:nrow(dt_sA2)]
dt_sA2[,("diff"):=abs(rho1-rho2),by=1:nrow(dt_sA2)]

dt_sA3[,("ratio"):=n1/n2,by=1:nrow(dt_sA3)]
dt_sA3[,("n"):=n1+n2,by=1:nrow(dt_sA3)]
dt_sA3[,("diff"):=abs(rho1-rho2),by=1:nrow(dt_sA3)]

dt_sA3_10k[,("ratio"):=n1/n2,by=1:nrow(dt_sA3)]
dt_sA3_10k[,("n"):=n1+n2,by=1:nrow(dt_sA3)]
dt_sA3_10k[,("diff"):=abs(rho1-rho2),by=1:nrow(dt_sA3)]

# convert to long
# 100 sims
dt.long <- melt(dt, measure.vars = tests, variable.name = "test", value.name = "power")
# 1000 sims
dtk.long <- melt(dtk, measure.vars = tests, variable.name = "test", value.name = "power")
dt_s1.long <- melt(dt_sA1, measure.vars = tests, variable.name = "test", value.name = "power")
dt_s2.long <- melt(dt_sA2, measure.vars = tests, variable.name = "test", value.name = "power")
dt_s3.long <- melt(dt_sA3, measure.vars = tests, variable.name = "test", value.name = "power")

# 10000 sims
dt_s3_10k.long <- melt(dt_sA3_10k, measure.vars = tests, variable.name = "test", value.name = "power")

# overall average power - NOTE - the below are only for the 100 simulation run
## Intent was to produce full sets using 1,000 simulation; as noted in part above two issues occured
# 1. a syntactical error rendered processed results not being retained (described above, full code not shown)
# 2. when reprocessing using an alternate parallelised approach (not shown in this script), 
#    the escape button was accidentally pressed when checking progress.
#  The above two errors led to a further restructing so these real world issues could be mitigated in future
#  1. using parallel processing to improve speed
#  2. using SQL to store results robustly, and multiple processors can write to database simultaneously
######  The use of SQL is important, as results are stored securely external to R: 
########    Solution 1. progress can be checked to ensure results are processing as intended
########    Solution 2. if a power outage occurs, or some crash occurs, results are not lost

#  SO - marginal results using 100 simulations (ie. LIMITED VALIDITY IN THE BELOW; proof of concept)
dt.long[(method=="pearson")&(dist=="normal"),round(mean(power),2),by=test]
#        test   V1
# 1: fz_nosim 0.74
# 2:       fz 0.75
# 3:      gtv 0.75
# 4:      slr 0.83
# 5:      zou 0.74
dt.long[(method=="pearson")&(dist=="gamma")&(p1=="1.5"),round(mean(power),2),by=test]
#        test   V1
# 1: fz_nosim 0.74
# 2:       fz 0.74
# 3:      gtv 0.74
# 4:      slr 0.82
# 5:      zou 0.74
dt.long[(method=="pearson")&(dist=="gamma")&(p1=="1"),round(mean(power),2),by=test]
#        test   V1
# 1: fz_nosim 0.74
# 2:       fz 0.50
# 3:      gtv 0.51
# 4:      slr 0.62
# 5:      zou 0.50


dt.long[(method=="pearson")&
          (dist=="normal")&
          (rho1==0)&
          (rho2==0),round(mean(power),2),by=test]
# test   V1
# 1: fz_nosim 0.03
# 2:       fz 0.05
# 3:      gtv 0.06
# 4:      slr 0.22
# 5:      zou 0.05

dt.long[(method=="pearson")&
          (dist=="normal")&
          (rho1==-0.2)&
          (rho2==-0.2),round(mean(power),2),by=test]
# test   V1
# 1: fz_nosim 0.03
# 2:       fz 0.05
# 3:      gtv 0.06
# 4:      slr 0.20
# 5:      zou 0.05
dt.long[(method=="pearson")&
          (dist=="normal")&
          (rho1==0.9)&
          (rho2==0.9),round(mean(power),2),by=test]
# test   V1
# 1: fz_nosim 0.03
# 2:       fz 0.05
# 3:      gtv 0.05
# 4:      slr 0.20
# 5:      zou 0.05

dt[(method=="pearson")&
          (dist=="normal"),list("mean_n1" = round(mean(n1),0),
                                "mean_n2" = round(mean(n2),0),
                              "fz_nosim" = round(mean(fz_nosim),2),
                               "fz" = round(mean(fz),2),
                               "zou" = round(mean(zou),2),
                               "gtv" = round(mean(gtv),2),
                               "slr" = round(mean(slr),2))
                               ,by=list(ratio)]
#         ratio mean_n1 mean_n2 fz_nosim   fz  zou  gtv  slr
#  1:  1.000000     272     272     0.77 0.77 0.77 0.77 0.78
#  2:  0.500000     158     315     0.77 0.77 0.77 0.78 0.79
#  3:  0.250000      93     372     0.76 0.77 0.77 0.77 0.82
#  4:  0.125000      56     450     0.74 0.75 0.75 0.75 0.85
#  5:  0.062500      35     560     0.71 0.71 0.71 0.72 0.87
#  6:  0.031250      22     720     0.67 0.67 0.67 0.68 0.89
#  7:  0.015625      15     960     0.61 0.62 0.62 0.63 0.91
#  8:  2.000000     315     158     0.77 0.77 0.77 0.78 0.80
#  9:  4.000000     372      93     0.76 0.77 0.77 0.77 0.82
# 10:  8.000000     450      56     0.74 0.74 0.74 0.75 0.84
# 11: 16.000000     560      35     0.71 0.71 0.71 0.72 0.87
# 12: 32.000000     720      22     0.67 0.67 0.67 0.68 0.89
# 13: 64.000000     960      15     0.61 0.62 0.62 0.63 0.91

dt.long[(method=="pearson")&(dist=="normal"),list(p50=round(quantile(power, .50, na.rm=TRUE),2),
                                                  p25=round(quantile(power, .025, na.rm=TRUE),2),
                                                  p75=round(quantile(power, .975, na.rm=TRUE),2)),by=test]
# test  p50  p25 p75
# 1: fz_nosim 0.98 0.03   1
# 2:       fz 0.99 0.04   1
# 3:      gtv 0.99 0.05   1
# 4:      slr 1.00 0.10   1
# 5:      zou 0.99 0.04   1
dt.long[(method=="pearson")&(dist=="gamma")&(p1=="1.5"),list(p50=round(quantile(power, .50, na.rm=TRUE),2),
                                                             p25=round(quantile(power, .25, na.rm=TRUE),2),
                                                             p75=round(quantile(power, .75, na.rm=TRUE),2)),by=test]
# test  p50  p25 p75
# 1: fz_nosim 0.98 0.47   1
# 2:       fz 0.98 0.45   1
# 3:      gtv 0.99 0.46   1
# 4:      slr 1.00 0.72   1
# 5:      zou 0.98 0.45   1
dt.long[(method=="pearson")&(dist=="gamma")&(p1=="1"),list(p50=round(quantile(power, .50, na.rm=TRUE),2),
                                                           p25=round(quantile(power, .025, na.rm=TRUE),2),
                                                           p75=round(quantile(power, .975, na.rm=TRUE),2)),by=test]
# test  p50  p25 p75
# 1: fz_nosim 0.98 0.03   1
# 2:       fz 0.46 0.00   1
# 3:      gtv 0.47 0.00   1
# 4:      slr 0.72 0.00   1
# 5:      zou 0.46 0.00   1


