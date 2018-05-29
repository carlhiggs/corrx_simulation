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

# may have to set working directory
setwd("C:/Users/Carl/OneDrive/Research/2 - BCA/Research project/bca_rp2/scripts")
require(compiler)
require(Rcpp)
require("simstudy")
sourceCpp('test.cpp')
require(data.table)
require(ggplot2)

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
### NOTE!! YOu don't need to do the simulation within this test --- duh
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

gtv_alt <- function(a,b,M=1,method = "pearson") {
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
### WAIT! Don't do this as a f-ing vector--- just single, and return test
# Then treat as per other tests--- do similar refactoring of GTV
# Based on Efron and Tibshirani
n  <- c(nrow(a),nrow(b))
r  <- c(cor(a,method = method)[2,1], cor(b,method = method)[2,1])
z  <- atanh(r)
# g <- c(rep("A",n[1]),rep("B",n[2]))
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

dt_normal <- dt[(dist=="normal"),]
dt_gamma1 <- dt[(dist=="gamma")&(p1==1.5),]
dt_gamma2 <- dt[(dist=="gamma")&(p1==1),]

# run main simulation study
system.time(dt_normal[,  c("fz_nosim","fz","gtv","slr","zou"):={
  cat("\r",.GRP,"\r")
  as.list(corr_power_compiled(rho = c(rho1,rho2), 
                      n = c(n1,n2),
                      distr = dist,
                      param1a = c(p1,p1), 
                      param1b = c(p1,p1),
                      param2a = c(p2,p2), 
                      param2b = c(p2,p2),
                      test    = c( "fz_nosim","fz","gtv","slr","zou"),
                      alpha   = 0.05, 
                      sidedness=2,
                      method=as.character(method),  
                      nsims = 1000,
                      power_only = TRUE))
},by = 1:nrow(dt_normal)] )


system.time(dt_gamma1[,  c("fz_nosim","fz","gtv","slr","zou"):={
  cat("\r",.GRP,"\r")
  as.list(corr_power_compiled(rho = c(rho1,rho2), 
                              n = c(n1,n2),
                              distr = dist,
                              param1a = c(p1,p1), 
                              param1b = c(p1,p1),
                              param2a = c(p2,p2), 
                              param2b = c(p2,p2),
                              test    = c( "fz_nosim","fz","gtv","slr","zou"),
                              alpha   = 0.05, 
                              sidedness=2,
                              method=as.character(method),  
                              nsims = 1000,
                              power_only = TRUE))
},by = 1:nrow(dt_gamma1)] )


system.time(dt_gamma2[,  c("fz_nosim","fz","gtv","slr","zou"):={
  cat("\r",.GRP,"\r")
  as.list(corr_power_compiled(rho = c(rho1,rho2), 
                              n = c(n1,n2),
                              distr = dist,
                              param1a = c(p1,p1), 
                              param1b = c(p1,p1),
                              param2a = c(p2,p2), 
                              param2b = c(p2,p2),
                              test    = c( "fz_nosim","fz","gtv","slr","zou"),
                              alpha   = 0.05, 
                              sidedness=2,
                              method=as.character(method),  
                              nsims = 1000,
                              power_only = TRUE))
},by = 1:nrow(dt_gamma2)] )


## More examples - test runs for now to implement once dt is calculated
# add in extra summary vars
dt[,("ratio"):=n1/n2,by=1:nrow(dt)]
dt[,("n"):=n1+n2,by=1:nrow(dt)]
dt[,("diff"):=abs(rho1-rho2),by=1:nrow(dt)]


### Due to typo issue with simulation (used macro instead of vector to store results in, 
### resulting in main results only stored for first test (the non-simulation one)
### I have now chosen specific scenarios to process at higher levels of simulation

## How does ratio impact with Pearsons/Spearmans?
# Scenario A1
#  n1 == 60 and n2 == 120
# nrow(dtk[(n1==60)&(n2==120),)
# 2166 rows
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
# user   system  elapsed 
# 48278.95   673.58 48780.27 
save.image("C:/Users/Carl/OneDrive/Research/2 - BCA/Research project/bca_rp2/scripts/r_power_work_dt_sA1_20180527.RData")


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

save.image("C:/Users/Carl/OneDrive/Research/2 - BCA/Research project/bca_rp2/scripts/r_power_work_dt_sA1_10k_20180527.RData")
# Scenario A2 
# n1 == 90 and n2 == 90
# Note - this doesn't exist; I'll have to fabricate it
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

# user   system  elapsed 
# 48350.98   736.72 48922.11 
save.image("C:/Users/Carl/OneDrive/Research/2 - BCA/Research project/bca_rp2/scripts/r_power_work_dt_sA2_20180527.RData")

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
save.image("C:/Users/Carl/OneDrive/Research/2 - BCA/Research project/bca_rp2/scripts/r_power_work_dt_sA2_10k_20180527.RData")
# Scenario B1
# rho1 == 0.2 && rho2 == 0.5
# nrow(dtk[(rho1==0.2)&(rho2==0.5),])
# 294 rows
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
# user  system elapsed 
# 6834.87   92.21 6921.39 
save.image("C:/Users/Carl/OneDrive/Research/2 - BCA/Research project/bca_rp2/scripts/r_power_work_dt_sA3_20180527.RData")

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
save.image("C:/Users/Carl/OneDrive/Research/2 - BCA/Research project/bca_rp2/scripts/r_power_work_dt_sA3_10k_20180527.RData")

# Hero attempt
dt_hero <-  dt[,c("method","dist","p1","p2","n1","n2","rho1","rho2")]
system.time(dt_hero[, c( "fz_nosim","fz","gtv","slr","zou"):={
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
save.image("C:/Users/Carl/OneDrive/Research/2 - BCA/Research project/bca_rp2/scripts/r_power_work_dt_sA3_20180527.RData")


### scratch code for getting results processed from environment on work computer
# tmp.env <- new.env() # create a temporary environment
# load("r_power_work_partial_dt.RData", envir=tmp.env) # load workspace into temporary environment
# x <- get("dt", pos=tmp.env) # get the objects you need into your globalenv()
# #x <- tmp.env$x # equivalent to previous line
# rm(tmp.env) # remove the temporary environment to free up memory
# dt <- rbind(other[is.na(fz_nosim)==FALSE,],dt_test)

tmp.env <- new.env() # create a temporary environment
load("r_power_work_dt_normal_20180526.RData", envir=tmp.env) # load workspace into temporary environment
dtk_normal <- get("dt_normal", pos=tmp.env) # get the objects you need into your globalenv()
load("r_power_work_dt_gamma1_20180526.RData", envir=tmp.env) # load workspace into temporary environment
dtk_gamma1 <- get("dt_gamma1", pos=tmp.env) # get the objects you need into your globalenv()
load("r_power_work_dt_gamma2_20180526.RData", envir=tmp.env) # load workspace into temporary environment
dtk_gamma2 <- get("dt_gamma2", pos=tmp.env) # get the objects you need into your globalenv()
rm(tmp.env) # remove the temporary environment to free up memory
#x <- tmp.env$x # equivalent to previous line
dtk <- rbind(dtk_normal,dtk_gamma1,dtk_gamma2)


tmp.env <- new.env() # create a temporary environment
load("r_power_work_dt_sA1_20180527.RData", envir=tmp.env) # load workspace into temporary environment
dt_sA1 <- get("dt_sA1", pos=tmp.env) # get the objects you need into your globalenv()
load("r_power_work_dt_sA2_20180527.RData", envir=tmp.env) # load workspace into temporary environment
dt_sA2 <- get("dt_sA2", pos=tmp.env) # get the objects you need into your globalenv()
load("r_power_work_dt_sA3_20180527.RData", envir=tmp.env) # load workspace into temporary environment
dt_sA3 <- get("dt_sA3", pos=tmp.env) # get the objects you need into your globalenv()
load("r_power_work_dt_sA3_10k_20180527.RData", envir=tmp.env) # load workspace into temporary environment
dt_sA3_10k <- get("dt_sA3_10k", pos=tmp.env) # get the objects you need into your globalenv()
rm(tmp.env) # remove the temporary environment to free up memory
#x <- tmp.env$x # equivalent to previous line

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

# overall average power
dt.long[(method=="pearson")&(dist=="normal"),round(mean(power),2),by=test]
dt.long[(method=="pearson")&(dist=="gamma")&(p1=="1.5"),round(mean(power),2),by=test]
dt.long[(method=="pearson")&(dist=="gamma")&(p1=="1"),round(mean(power),2),by=test]


# concerning that SLR has 22% power to detect difference when no difference....
# Then again, this is beta - probability of not rejecting when false --- 20%..
# hmmm.. no its not, it should be alpha which is 5% - hence other tests
# ie. probability of rejecting when H0 is true
dt.long[(method=="pearson")&
          (dist=="normal")&
          (rho1==0)&
          (rho2==0),round(mean(power),2),by=test]


dt.long[(method=="pearson")&
          (dist=="normal")&
          (rho1==-0.2)&
          (rho2==-0.2),round(mean(power),2),by=test]

dt.long[(method=="pearson")&
          (dist=="normal")&
          (rho1==0.9)&
          (rho2==0.9),round(mean(power),2),by=test]


dt[(method=="pearson")&
          (dist=="normal"),list("mean_n1" = round(mean(n1),0),
                                "mean_n2" = round(mean(n2),0),
                              "fz_nosim" = round(mean(fz_nosim),2),
                               "fz" = round(mean(fz),2),
                               "zou" = round(mean(zou),2),
                               "gtv" = round(mean(gtv),2),
                               "slr" = round(mean(slr),2))
                               ,by=list(ratio)]


dt.long[(method=="pearson")&(dist=="normal"),list(p50=round(quantile(power, .50, na.rm=TRUE),2),
                                                  p25=round(quantile(power, .025, na.rm=TRUE),2),
                                                  p75=round(quantile(power, .975, na.rm=TRUE),2)),by=test]
dt.long[(method=="pearson")&(dist=="gamma")&(p1=="1.5"),list(p50=round(quantile(power, .50, na.rm=TRUE),2),
                                                             p25=round(quantile(power, .25, na.rm=TRUE),2),
                                                             p75=round(quantile(power, .75, na.rm=TRUE),2)),by=test]
dt.long[(method=="pearson")&(dist=="gamma")&(p1=="1"),list(p50=round(quantile(power, .50, na.rm=TRUE),2),
                                                           p25=round(quantile(power, .025, na.rm=TRUE),2),
                                                           p75=round(quantile(power, .975, na.rm=TRUE),2)),by=test]



# Function to produce corrx plots from long data
corrxplot <- function(data.long,method = "pearson",dist = "normal",test="",n1 = 60,n2 = 120,
                      param1a = 0,param1b = 0, param2a = 1,param2b = 1,
                      rho1 = 0.2, rho2 = 0.5, ratio = 0.5,
                      nsims = 100, names = c("Mz twins","Dz twins"), type  = "contour",
                      alpha = 0.05, threshold = 0.8,
                      graph_out = "corrxplot",gwidth = 7,gheight=6.5) {
  
  graph_out <- if(graph_out=="corrxplot") paste("corrx",type,method,dist,
                                                paste0(param1a,'-',param1a),
                                                paste0(param2a,'-',param2b),
                                                paste0('n',n1,'-',n2),
                                                rho1,rho2,".pdf",sep="_") else graph_out
  
  # rename input vars so not ambiguous when referencing data.table fields
  lmethod  <- method
  ldist    <- dist
  ltest    <- test
  lparam1a <- param1a
  lparam1b <- param1b
  lparam2a <- param2a
  lparam2b <- param2b
  ln1      <- n1
  ln2      <- n2
  lrho1    <- rho1
  lrho2    <- rho2
  lratio   <- ratio
  lnsims   <- nsims
  lnames   <- names

  # produce contour plot
  if(type=="contour"){
      # make correlation-power matrix
      corrx_tmp <- matrix(unlist(data.long[(method==lmethod)&
                                          (dist==ldist)&
                                          (p1==lparam1a)&
                                          (p2==lparam2a)&
                                          (n1==ln1)&
                                          (n2==ln2)&
                                          (test==ltest),]$power),
                      nrow=19,ncol=19, byrow = TRUE)
    # lookup table for test names
    test_lookup <- cbind(data.table("test" = c("fz_nosim","fz","gtv","slr","zou")),
                         data.table("name" = c("Fisher's z (no sim)","Fisher's Z","GTV","SLR","Zou's CI")))
      
    # initialise output plot
    pdf(graph_out,width=gwidth,height=gheight)
    # make contour plot
    filled.contour(x = corrs,y = corrs,z = corrx_tmp, nlevels = 100,
                 xlim = c(-.9,.9), ylim = c(-.9,.9), zlim = c(0,1),
                 plot.axes = {contour(x = corrs, y = corrs, z = corrx_tmp,
                                      levels = threshold, at = seq(-.9, .9, 0.2), drawlabels = FALSE, axes = FALSE,
                                      add = TRUE, lwd = 3, col = "steelblue3");
                   abline(v = seq(-9, 0.9, 0.1), lwd = .5, col = "lightgray", lty = 2)
                   abline(h = seq(-9, 0.9, 0.1), lwd = .5, col = "lightgray", lty = 2)
                   axis(1, seq(-.9,.9,0.2))
                   axis(2, seq(-.9,.9,0.2))},
                 key.title = {title(main="Power")},
                 plot.title = title(main = paste0(test_lookup[test==ltest,name]," test ~ ",
                                                  ldist,"((",lparam1a,",",lparam2a,"),(",
                                                  lparam1b,",",lparam2b,"))","\n",
                                                  "Mz = ",ln1,
                                                  ", Dz = ",ln2,
                                                  ", alpha: ",alpha, ", sims: ",lnsims),
                                    xlab = paste0("Correlation in ",lnames[1]),
                                    ylab = paste0("Correlation in ",lnames[2]), adj = 0),
                 color.palette =  colorRampPalette(c("#f7fcf0","#525252")));
    # add arrow to legend to indicate target power level (0.8)
    arrows(0.6, 0.54, 0.77, 0.54, length = 0.14, lwd = 3, col = "steelblue3")  
    # finalise plot
    dev.off()
  }
  
  if(type=="npower"){
    # get data subset
    um <- data.long[(dist==ldist)&
                    (method==lmethod)&
                    (p1==lparam1a)&
                    (p2==lparam2a)&
                    (rho1==lrho1)&
                    (rho2==lrho2)&
                    (ratio==lratio),
                  power,by=list(n,test)]
    
    # interpolate power given log(n)
    um[,"log_n":=log2(n),by=1:nrow(um)]
    fit <- data.table()
    for(x in levels(um$test)) {
      fit <- rbind(fit,
                   cbind("test" = x,
                         "n"    = min(um$n):max(um$n),
                         "power" = splinefun(um[test %in% x,log_n], 
                                             um[test %in% x,power],
                                             method = "monoH.FC")(log2(min(um$n):max(um$n)))))
    }
    
    fit$n <- as.integer(fit$n)
    fit$power <- as.double(fit$power)
    
    # find x intercept give y of 0.8
    cross <- data.table()
    for(x in levels(um$test)) {
      cross <- rbind(cross, 
                     fit[(test %in% x)][which.min(abs(threshold-fit[(test %in% x),power]))])
    }
    
    # associate x intercept (required sample size for 80% power) with test name
    cross$label <- c("FZ (no sim.)",
                     "FZ",
                     "GTV",
                     "SLR",
                     "Zou's CI")
    cross[,"label":= paste0(label," (",n,")"),by=1:nrow(cross)]
    
    # order tests by smallest sample size estimate required to achieve power threshold
    cross <- cross[order(+rank(n))]
    
    # define plot title
    title <- paste0("Power to detect difference in ",method," correlations, by sample size\n",
                    dist,"((",param1a,",",param1b,"),(",param2a,",",param2b,"))","\n",
                    "rho: (",rho1,",",rho2,")",
                    "; Mz to Dz ratio: ",ratio, "; sims: ",nsims)
    # initialise plot
    pdf(graph_out,width=gwidth,height=gheight)
    p <- ggplot(NULL, aes(x = n, y = power, colour = test, group = test))+ 
      scale_x_continuous(trans='log2',bquote(N~(log[2]~scale)), 
                         breaks = unique(um$n),
                         limits = c(30,1920)) +
      scale_y_continuous(bquote(Power~(1-beta)), 
                         breaks = seq(0,1,0.1),
                         limits = c(0,1)) +
      geom_point(data = um)  +
      geom_line(data = fit, lwd = 1) +
      geom_hline(yintercept = 0.8) +
      geom_vline(aes(xintercept = cross$n, colour = cross$test)) +
      scale_colour_discrete(name="Tests (req. n)",
                            breaks=cross$test,
                            labels=cross$label)  +
      theme(panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black")) +
      ggtitle(title)
    
    print(p)
    # finalise plot
    dev.off()
  }
  if(type %in% c("diffpower","diffpowerabs")){
    ln = n1+n2
    # plot power curve by difference given parameters
    um <- data.long[(dist=="normal")&
                    (method=="pearson")&
                    (p1==lparam1a)&
                    (p2==lparam2a)&
                    (n ==ln)&
                    (ratio==lratio),
                  power,by=list(test,rho1,rho2,diff)]
    
    um[,"z1":=atanh(rho1), by = 1:nrow(um)]
    um[,"z2":=atanh(rho2), by = 1:nrow(um)]
    
    if(type=="diffpower") {  
      um[,"diff":=tanh(z1-z2), by = 1:nrow(um)]
      xformula <- bquote("d* ="~tanh(atanh(rho[1])-atanh(rho[2])))
    } else {
      um[,"diff":=abs(tanh(z1-z2)), by = 1:nrow(um)]
      xformula <- bquote("d* = |"~tanh(atanh(rho[1])-atanh(rho[2]))~"|")
    }
    
    # fit by averaging over differences arising from differen rho combinations
    fit <- data.table()
    for(x in levels(um$test)) {
      spline <- with(um[test %in% x,], smooth.spline(diff, 
                                                     power,
                                                     all.knots=seq(0,1,0.01)))
      fit <- rbind(fit,
                   cbind("test" = x,
                         "diff"    = seq(0,1,0.01),
                         "power" = splinefun(spline$x, 
                                             spline$y,
                                             method = "monoH.FC")(seq(0,1,0.01))))
    }
    
    fit$diff <- as.double(fit$diff)
    fit$power <- as.double(fit$power)
    
    cross <- data.table()
    for(x in levels(um$test)) {
      cross <- rbind(cross, 
                     fit[(test %in% x)][which.min(abs(0.8-fit[(test %in% x),power]))])
    }
    
    cross$label <- c("FZ (no sim.)",
                     "FZ",
                     "GTV",
                     "SLR",
                     "Zou's CI")
    cross[,"label":= paste0(label," (",round(diff,2),")"),by=1:nrow(cross)]
    cross <- cross[order(+rank(diff))]
    
    title <- paste0("Power to detect difference in ",method," correlations, by difference*\n",
                    dist,"((",param1a,",",param1b,"),(",param2a,",",param2b,"))","\n",
                    "N: ",ln,
                    "; Mz to Dz ratio: ",ratio, "; sims: ",nsims)
    # initialise plot
    pdf(graph_out,width=gwidth,height=gheight)
    p <- ggplot(NULL, aes(x = diff, y = power, colour = test, group = test))+ 
      scale_x_continuous(xformula, 
                         breaks = seq(0,1,0.1),
                         limits = c(0,1)) +
      scale_y_continuous(bquote(Power~(1-beta)), 
                         breaks = seq(0,1,0.1),
                         limits = c(0,1)) +
      geom_point(data = um)  +
      geom_line(data = fit, lwd = 1) +
      geom_hline(yintercept = 0.8) +
      geom_vline(aes(xintercept = cross$diff, colour = cross$test)) +
      scale_colour_discrete(name="d* for 80% Power",
                            breaks=cross$test,
                            labels=cross$label)  +
      theme(panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black")) +
      ggtitle(title)
    
    print(p)
    # finalise plot
    dev.off()
    
  }
}


# Scenario 1pre - 100 sims 60:120
corrxplot(data.long  =  dt.long,
          method     = "pearson",
          dist       =  "normal",
          test       =  "slr",
          n1         =  60,
          n2         =  120,
          param1a    =  0,
          param1b    =  0,
          param2a    =  1,
          param2b    =  1,
          nsims      =  100,
          names      =  c("Mz twins","Dz twins"),
          type       = "contour",
          graph_out  =  "../figs/corrx_contour_sA1pre_n180_mzdz.5_slr_s100.pdf",
          alpha      =  0.05,
          threshold  =  0.8) 



# Scenario 1  - 1000 sims 60:120
corrxplot(data.long  =  dt_s1.long,
          method     = "pearson",
          dist       =  "normal",
          test       =  "slr",
          n1         =  60,
          n2         =  120,
          param1a    =  0,
          param1b    =  0,
          param2a    =  1,
          param2b    =  1,
          nsims      =  1000,
          names      =  c("Mz twins","Dz twins"),
          type       = "contour",
          graph_out  =  "../figs/corrx_contour_sA1_n180_mzdz.5_slr_s1000.pdf",
          alpha      =  0.05,
          threshold  =  0.8) 


# Scenario 2  - 1000 sims 90:90
corrxplot(data.long  =  dt_s2.long,
          method     = "pearson",
          dist       =  "normal",
          test       =  "slr",
          n1         =  90,
          n2         =  90,
          param1a    =  0,
          param1b    =  0,
          param2a    =  1,
          param2b    =  1,
          nsims      =  1000,
          names      =  c("Mz twins","Dz twins"),
          type       = "contour",
          graph_out  =  "../figs/corrx_contour_sA2_n180_mzdz1_slr_s1000.pdf",
          alpha      =  0.05,
          threshold  =  0.8) 

# plot power curve by N given parameters  - Normal 100 sim
corrxplot(data.long  =  dt.long,
          method     = "pearson",
          dist       =  "normal",
          test       =  "slr",
          rho1       =  0.2,
          rho2       =  0.5,
          ratio      =  0.5,
          param1a    =  0,
          param2a    =  1,
          nsims      =  100,
          type       = "npower",
          graph_out  =  "../figs/corrx_npower_norm_r.2_r.5_mzdz.5_s100.pdf",
          alpha      =  0.05,
          threshold  =  0.8) 



# plot power curve by N given parameters  - Normal - Scenario 3 - 1000sim
corrxplot(data.long  =  dt_s3.long,
          method     = "pearson",
          dist       =  "normal",
          test       =  "slr",
          rho1       =  0.2,
          rho2       =  0.5,
          ratio      =  0.5,
          param1a    =  0,
          param2a    =  1,
          nsims      =  1000,
          type       = "npower",
          graph_out  =  "../figs/corrx_npower_norm_r.2_r.5_mzdz.5_s1000.pdf",
          alpha      =  0.05,
          threshold  =  0.8) 


# plot power curve by N given parameters  - Normal - Scenario 3 - 10000sim
corrxplot(data.long  =  dt_s3_10k.long,
          method     = "pearson",
          dist       =  "normal",
          test       =  "slr",
          rho1       =  0.2,
          rho2       =  0.5,
          ratio      =  0.5,
          param1a    =  0,
          param2a    =  1,
          nsims      =  10000,
          type       = "npower",
          graph_out  =  "../figs/corrx_npower_norm_r.2_r.5_mzdz.5_s10000.pdf",
          alpha      =  0.05,
          threshold  =  0.8) 

# plot power curve by N given parameters  - gamma mild 100 sim
corrxplot(data.long  =  dt.long,
          method     = "pearson",
          dist       =  "gamma",
          test       =  "slr",
          rho1       =  0.2,
          rho2       =  0.5,
          ratio      =  0.5,
          param1a    =  1.5,
          param2a    =  0.09,
          nsims      =  100,
          type       = "npower",
          graph_out  =  "../figs/corrx_npower_mildskew_r.2_r.5_mzdz.5_s100.pdf",
          alpha      =  0.05,
          threshold  =  0.8) 

# plot power curve by N given parameters  - gamma mild - Scenario 3 - 1000sim
corrxplot(data.long  =  dt_s3.long,
          method     = "pearson",
          dist       =  "gamma",
          test       =  "slr",
          rho1       =  0.2,
          rho2       =  0.5,
          ratio      =  0.5,
          param1a    =  1.5,
          param2a    =  0.09,
          nsims      =  1000,
          type       = "npower",
          graph_out  =  "../figs/corrx_npower_mildskew_r.2_r.5_mzdz.5_s1000.pdf",
          alpha      =  0.05,
          threshold  =  0.8) 

# plot power curve by N given parameters  - gamma mild - Scenario 3 - 10000sim
corrxplot(data.long  =  dt_s3_10k.long,
          method     = "pearson",
          dist       =  "gamma",
          test       =  "slr",
          rho1       =  0.2,
          rho2       =  0.5,
          ratio      =  0.5,
          param1a    =  1.5,
          param2a    =  0.09,
          nsims      =  10000,
          type       = "npower",
          graph_out  =  "../figs/corrx_npower_mildskew_r.2_r.5_mzdz.5_s10000.pdf",
          alpha      =  0.05,
          threshold  =  0.8) 



# plot power curve by N given parameters  - gamma extr 100 sim
corrxplot(data.long  =  dt.long,
          method     = "pearson",
          dist       =  "gamma",
          test       =  "slr",
          rho1       =  0.2,
          rho2       =  0.5,
          ratio      =  0.5,
          param1a    =  1,
          param2a    =  5,
          nsims      =  100,
          type       = "npower",
          graph_out  =  "../figs/corrx_npower_extrskew_r.2_r.5_mzdz.5_s100.pdf",
          alpha      =  0.05,
          threshold  =  0.8) 



# plot power curve by N given parameters  - gamma extr - Scenario 3 - 1000sim
corrxplot(data.long  =  dt_s3.long,
          method     = "pearson",
          dist       =  "gamma",
          test       =  "slr",
          rho1       =  0.2,
          rho2       =  0.5,
          ratio      =  0.5,
          param1a    =  1,
          param2a    =  5,
          nsims      =  1000,
          type       = "npower",
          graph_out  =  "../figs/corrx_npower_extrskew_r.2_r.5_mzdz.5_s1000.pdf",
          alpha      =  0.05,
          threshold  =  0.8) 


# plot power curve by N given parameters  - gamma extr - Scenario 3 - 10000sim
corrxplot(data.long  =  dt_s3_10k.long,
          method     = "pearson",
          dist       =  "gamma",
          test       =  "slr",
          rho1       =  0.2,
          rho2       =  0.5,
          ratio      =  0.5,
          param1a    =  1,
          param2a    =  5,
          nsims      =  10000,
          type       = "npower",
          graph_out  =  "../figs/corrx_npower_extrskew_r.2_r.5_mzdz.5_s10000.pdf",
          alpha      =  0.05,
          threshold  =  0.8) 


# difference plot (note the imprecision though - not a perfect approach)
corrxplot(data.long  =  dt.long,
          method     = "pearson",
          dist       =  "normal",
          param1a    =  0,
          param2a    =  1,
          ratio      =  .5,
          n1         =  60,
          n2         =  120,
          nsims      =  100,
          type       = "diffpowerabs",
          graph_out  =  "../figs/corrx_diffpower_normal_60_120_mzdz1_s100.pdf",
          alpha      =  0.05,
          threshold  =  0.8) 


# difference plot (note the imprecision though - not a perfect approach)
corrxplot(data.long  =  dt_s1.long,
          method     = "pearson",
          dist       =  "normal",
          param1a    =  0,
          param2a    =  1,
          ratio      =  .5,
          n1         =  60,
          n2         =  120,
          nsims      =  1000,
          type       = "diffpowerabs",
          graph_out  =  "../figs/corrx_diffpower_normal_60_120_mzdz1_s1000.pdf",
          alpha      =  0.05,
          threshold  =  0.8) 


# difference plot (note the imprecision though - not a perfect approach)
corrxplot(data.long  =  dt_s2.long,
          method     = "pearson",
          dist       =  "normal",
          param1a    =  0,
          param2a    =  1,
          ratio      =  1,
          n1         =  90,
          n2         =  90,
          nsims      =  1000,
          type       = "diffpowerabs",
          graph_out  =  "../figs/corrx_diffpower_normal_90_90_mzdz1_s1000.pdf",
          alpha      =  0.05,
          threshold  =  0.8) 


# difference plot (note the imprecision though - not a perfect approach)
corrxplot(data.long  =  dt_s2.long,
          method     = "pearson",
          dist       = "gamma",
          param1a    =  1,
          param2a    =  5,
          ratio      =  1,
          n1         =  90,
          n2         =  90,
          nsims      =  1000,
          type       = "diffpowerabs",
          graph_out  =  "../figs/corrx_diffpower_extrskew_90_90_mzdz1_s1000.pdf",
          alpha      =  0.05,
          threshold  =  0.8) 


# # plot power curve by difference given parameters - Gamma, extreme 
um <- dt.long[(dist=="gamma")&
                (method=="pearson")&
                (p1==1)&
                (p2==5)&
                (ratio==1)&
                (n==240),
              power,by=list(test,rho1,rho2,diff)]

um[,"z1":=atanh(rho1), by = 1:nrow(um)]
um[,"z2":=atanh(rho2), by = 1:nrow(um)]
um[,"diff1":=tanh(z1-z2), by = 1:nrow(um)]
um[,"diff2":=abs(tanh(z1-z2)), by = 1:nrow(um)]

# fit by averaging over differences arising from differen rho combinations
fit <- data.table()
for(x in levels(um$test)) {
  spline <- with(um[test %in% x,], smooth.spline(diff2, 
                                                 power,
                                                 all.knots=seq(0,1,0.01)))
  fit <- rbind(fit,
               cbind("test" = x,
                     "diff2"    = seq(0,1,0.01),
                     "power" = splinefun(spline$x, 
                                         spline$y,
                                         method = "monoH.FC")(seq(0,1,0.01))))
}

fit$diff2 <- as.double(fit$diff2)
fit$power <- as.double(fit$power)

cross <- data.table()
for(x in levels(um$test)) {
  cross <- rbind(cross, 
                 fit[(test %in% x)][which.min(abs(0.8-fit[(test %in% x),power]))])
}

cross$label <- c("FZ (no sim.)",
                 "FZ",
                 "GTV",
                 "SLR",
                 "Zou's CI")
cross[,"label":= paste0(label," (",round(diff2,2),")"),by=1:nrow(cross)]
cross <- cross[order(+rank(diff2))]

method <- "pearson"
dist <- "gamma"
p <- c(1,5)
rho1 <- 0.2
rho2 <- 0.5
ratio <- 1
n <- 240
title <- paste0("Power to detect difference in ",method," correlations, by difference*\n",
                dist,"((",p[1],",",p[1],"),(",p[2],",",p[2],"))","\n",
                "N: ",n,
                "; Mz to Dz ratio: ",ratio, "; sims: ",100)

ggplot(NULL, aes(x = diff2, y = power, colour = test, group = test))+ 
  scale_x_continuous(bquote("d* = |"~tanh(atanh(rho[1])-atanh(rho[2]))~"|"), 
                     breaks = seq(0,1,0.1),
                     limits = c(0,1)) +
  scale_y_continuous(bquote(Power~(1-beta)), 
                     breaks = seq(0,1,0.1),
                     limits = c(0,1)) +
  geom_point(data = um)  +
  geom_line(data = fit, lwd = 1) +
  geom_hline(yintercept = 0.8) +
  geom_vline(aes(xintercept = cross$diff2, colour = cross$test)) +
  scale_colour_discrete(name="d* for 80% Power",
                        breaks=cross$test,
                        labels=cross$label)  +
  theme(panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  ggtitle(title)

# Non-normal distribution simulation
# normal reference distribution
bivariate_distribution_normal_m0_s1_n50 <- ggExtra::ggMarginal(data = as.data.frame(a), x = "V1", y = "V2")
plot(bivariate_distribution_normal_m0_s1_n50)
# gamma distribution example (mean/shape, and rate/dispersion)
# some positive skew (but distinctly non-normal)
gamma <- genCorGen(1000, nvars = 2, params1 = c(1.5,1.5), params2 = c(0.09,0.09),dist = "gamma",
                   corMatrix = matrix(c(1, 0.5, 0.5, 1), ncol = 2), wide = TRUE)[,2:3]
bivariate_distribution_gamma_m1.5_d0.09_n50 <- ggExtra::ggMarginal(data = as.data.frame(gamma), x = "V1", y = "V2")
plot(bivariate_distribution_gamma_m1.5_d0.09_n50)

# more gamma exploration, using our otherwise defaults
gamma_a <- genCorGen(20, nvars = 2, params1 = c(1,1), params2 = c(0.09,0.09),dist = "gamma",
                     corMatrix = matrix(c(1, 0.5, 0.5, 1), ncol = 2), wide = TRUE)[,2:3]
gamma_b <- genCorGen(90, nvars = 2, params1 = c(1,1), params2 = c(0.09,0.09),dist = "gamma",
                     corMatrix = matrix(c(1, 0.5, 0.5, 1), ncol = 2), wide = TRUE)[,2:3]

ggExtra::ggMarginal(data = as.data.frame(gamma_a), x = "V1", y = "V2")
ggExtra::ggMarginal(data = as.data.frame(gamma_b), x = "V1", y = "V2")



## Alternate fitted contour representation --- a lot bloody easier!
ggplot(dt_s1.long[(n1==60)&
                    (n2==120)&
                    (method=="pearson")&
                    (dist=="gamma")&(p1==1),],aes(x = rho1,y = rho2, z = power, colour = test, group = test))  +
  geom_contour(breaks = 0.8)


ggplot(dt_s2.long[(n1==90)&
                    (n2==90)&
                    (method=="pearson")&
                    (dist=="gamma")&(p1==1),],aes(x = rho1,y = rho2, z = power, colour = test, group = test))  +
  geom_contour(breaks = 0.8)