# R script to simulate power for difference in correlations (Pearson, Spearman, and later ... ICC)
# Carl Higgs 2017
#
# Adapted from https://janhove.github.io/design/2015/04/14/power-simulations-for-comparing-independent-correlations
# Jan Vanhove 14 April 2015

# The approach by Vanhove uses r.test in the psych package to calculate significance of 

# OWN TEST FORMULATION (benefits: not reliant on package; can swap out method for other; can report results)
# Confirm calculation of r.test
r.test(n = 20, r12 = 0.5, n2 = 50, r34 = 0.2)$p
## [1] 0.2207423
z_diff = atanh(0.5) - atanh(0.2)
z_se   = sqrt(1/(20-3) + 1/(50-3))
z_test = (z_diff/z_se)
z_p    = 2*pnorm(-abs(z_test))
z_p
## [1] 0.2207423

# Function to simulate two bivariate normal distributions based on respective population correlations.
# Reports:
#    - the significance of the sample correlation's difference (output 'z_p') 
#    - the power to detect a difference (output 'z_power') given alpha, beta and sidedness parameters.
# Optionally
#    - calculate the power or significance of difference of Pearson, Spearman or Kendall correlations
#    - calculate the significance and power without sampling using the two input correlations (simulation = FALSE).
#    - output a full log of simulation paramaters with results
#    - output a single statistic (z_p, or z_power).
# To develop
#    - accounting for clustering (ie. icc in twin studies)
#    - allow parameterisation to explore change in source pop mean and sd (unequal variance, etc)


# New plan:
#  Incorporate simstudy simulation approach - requires a supplied schema definition 'def' to be set up as below

# Correlated data example

## a vector containing the mean(s) of the distributions (for binary, should be prob)
#mu <- c(0.8, 0.1)  
## given dist, a parameterisation (e.g. rate [poisson], dispersion [gamma], var [normal] or max [uniform])
## For binary data, set to NULL
#param <- c(NULL, NULL) 
#distribution can be one of "poisson", "binary", "gamma", "uniform", "negbinom", "normal"
#rho1 = 0.78

# Testing out compilation of genCorGen according to directions at 
# https://www.r-statistics.com/2012/04/speed-up-your-r-code-using-a-just-in-time-jit-compiler/
# on advice from Koen Simons
require(compiler)
genCorGen_compiled <- cmpfun(genCorGen)

fo <- function() for (i in 1:1000) genCorGen(50, nvars = 2, params1 = 0, params2 = 1, dist = "normal", 
                                             corMatrix = matrix(c(1, 0.5, 0.5, 1), ncol = 2), wide = TRUE)

fo_c <- function() for (i in 1:1000) genCorGen_compiled(50, nvars = 2, params1 = 0, params2 = 1, dist = "normal", 
                                                        corMatrix = matrix(c(1, 0.5, 0.5, 1), ncol = 2), wide = TRUE)

system.time(fo())
system.time(fo_c())
enableJIT(3)
system.time(fo())
system.time(fo_c())
enableJIT(0)
## There doesn't appear to be any clear benefits, at least in basic case with normal dist and small sample size
# > system.time(fo())
# user  system elapsed 
# 14.97    0.42   18.71 
# > system.time(fo_c())
# user  system elapsed 
# 13.76    0.34   17.53 
# > enableJIT(3)
# [1] 0
# > system.time(fo())
# user  system elapsed 
# 14.44    0.36   17.45 
# > system.time(fo_c())
# user  system elapsed 
# 15.25    0.44   19.59 
# > enableJIT(0)
# [1] 3

fo <- function() for (i in 1:1000) genCorGen(500, nvars = 2, params1 = 0, params2 = 1, dist = "gamma", 
                                             corMatrix = matrix(c(1, 0.3, 0.3, 1), ncol = 2), wide = TRUE)

fo_c <- function() for (i in 1:1000) genCorGen_compiled(500, nvars = 2, params1 = 0, params2 = 1, dist = "gamma", 
                                                        corMatrix = matrix(c(1, 0.3, 0.3, 1), ncol = 2), wide = TRUE)

system.time(fo())
system.time(fo_c())
enableJIT(3)
system.time(fo())
system.time(fo_c())
enableJIT(0)
## Likewise, not clear benefits with larger sample size, different correlation and gamma distribution
# > system.time(fo())
# user  system elapsed 
# 13.84    0.14   17.19 
# > system.time(fo_c())
# user  system elapsed 
# 15.09    0.67   20.66 
# > enableJIT(3)
# [1] 0
# > system.time(fo())
# user  system elapsed 
# 13.58    0.41   18.19 
# > system.time(fo_c())
# user  system elapsed 
# 13.88    0.39   17.58 
# > enableJIT(0)
# [1] 3

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

# Fishers Z test - no sim
fz_nosim <- function(r1,r2,n1,n2,sidedness=2,method = "pearson",power = TRUE) {
  # Step 1: Compute sample correlation coefficients
  z1     <- atanh(r1)
  z2     <- atanh(r2)
  zdiff  <- z1-z2
  z_se   <- sqrt(1/(n1-3) + 1/(n2-3))
  z_test <- zdiff/z_se
  z_p    <- sidedness*pnorm(-abs(z_test))
  return(c(z_p,z_test))
  # return(z_p)
}

# GTV test statistic, based on code from Enes
gtv <- function(a,b,M=1e5,method = "pearson") {
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

a <- genCorGen(50, nvars = 2, params1 = 0, params2 = 1, dist = "normal", 
          corMatrix = matrix(c(1, 0.5, 0.5, 1), ncol = 2), wide = TRUE)[,2:3]
b <- genCorGen(50, nvars = 2, params1 = 0, params2 = 1, dist = "normal", 
             corMatrix = matrix(c(1, 0.2, 0.2, 1), ncol = 2), wide = TRUE)[,2:3]
result <- list()
system.time(result[["fz"]]       <- fz_test(a,b))
system.time(result[["fz_nosim"]] <- fz_test_nosim(0.5,0.5,50,50))
system.time(result[["gtv"]]      <- gtv_test(a,b))


fz_compiled <- cmpfun(fz)
fz_ns_compiled <- cmpfun(fz_nosim)
gtv_compiled <- cmpfun(gtv)
fc_fz     <- function() for (i in 1:1000)      fz_test(a,b)
fc_fz_ns  <- function() for (i in 1:1000)  fz_test_nosim(0.5,0.5,50,50)
fc_gtv    <- function() for (i in 1:1000)     gtv_test(a,b)
cfc_fz    <- function() for (i in 1:1000)  fz_compiled(a,b)
cfc_fz_ns <- function() for (i in 1:1000) fz_ns_compiled(0.5,0.5,50,50)
cfc_gtv   <- function() for (i in 1:1000) gtv_compiled(a,b)
system.time(fc_fz())
# > system.time(fc_fz())
# user  system elapsed 
# 1.42    0.02    1.48 
system.time(fc_fz_ns())
# > system.time(fc_fz_ns())
# user  system elapsed 
# 0.00    0.00    0.08 
system.time(fc_gtv())
# > system.time(fc_gtv())
# user  system elapsed 
# 248.23   15.93  286.11 
system.time(cfc_fz())
# > system.time(cfc_fz())
# user  system elapsed 
# 1.07    0.03    1.36 
system.time(cfc_fz_ns())
# > system.time(cfc_fz_ns())
# user  system elapsed 
# 0.01    0.00    0.02 
system.time(cfc_gtv())
# > system.time(cfc_gtv())
# user  system elapsed 
# 208.49   12.06  242.33 

# The self-compiled functions are notable quicker for the simulation tests
# (~15% to 30% reduced run time, for these samples); worth using
# except for analytical test; not worth compiling.

corr_diff_test <- function(rho = c(.2,.5), n = c(30,90), distr = "normal",
                    mu1 = c(0,0), mu2 = c(0,0),param1 = c(1,1), param2 = c(1,1),
                    alpha = 0.05, sidedness = 2, test = c("fz","gtv"),
                    method ="pearson", lower.tri = FALSE) {
  # cat("Parameters: ",rho[1],rho[2], n[1],n[2], mu1, mu2, param1, param2, distr, alpha, sidedness, method,"\n")
  if(lower.tri==TRUE){
    # only calculate lower matrix half when comparing across all correlation combinations
    if(rho[1] < rho[2]) { 
      return(NA)
    }
  }
  results <- list()
  if ("fz_nosim" %in% test) {
    results[["fz_nosim"]] <- fz_nosim(rho[1],rho[2],n[1],n[2], method = method)
    if(length(test)==1) return(results)
  }
  require("simstudy")
  a <- genCorGen(n[1], nvars = 2, params1 = mu1, params2 = param1,  
                dist = distr, corMatrix = matrix(c(1, rho[1], rho[1], 1), ncol = 2), 
                wide = TRUE)[,2:3]
  b <- genCorGen(n[2], nvars = 2, params1 = mu2, params2 = param2,  
                dist = distr, corMatrix = matrix(c(1, rho[2], rho[2], 1), ncol = 2), 
                wide = TRUE)[,2:3]
  if ("fz"       %in% test) results[["fz"]]       <- fz_compiled(a,b)
  if ("gtv"      %in% test) results[["gtv"]]      <- gtv_compiled(a,b)
  return(rbind(results[test]))
}
corr_diff_test(rho = c(.2,.4), n = c(100,300), distr = "normal",test = "fz_nosim") 
corr_diff_test(rho = c(.2,.4), n = c(100,300), distr = "normal",test = "fz") 
corr_diff_test(rho = c(.2,.4), n = c(100,300), distr = "normal",test = "gtv") 
corr_diff_test(rho = c(.2,.4), n = c(100,300), distr = "normal") 

# Corr power simulation
corr_power <- function(rho = c(.2,.5), n = c(30,90),distr = "normal",
                       mu1 = c(0,0),mu2 = c(0,0),param1 = c(1,1), param2 = c(1,1),
                       test = c("fz","gtv"),
                       alpha = 0.05, sidedness=2,method="pearson",  
                       nsims = 100,lower.tri = FALSE, power_only = FALSE){
  sim <- list()
  sim[["params"]]             <- c("method" = method, "rho_1" = rho[1], "rho_2" = rho[2],"n1" = n[1], "n2" = n[2],
                                 "alpha" = alpha, "sidedness" = sidedness, "nsims" = nsims, "distr" = distr)
  sim[["z_ref"]]              <- qnorm(1-alpha/sidedness)
  sim[["z_analytical"]]     <- fz_nosim(rho[1],rho[2],n[1],n[2], method = method)
  sim[["z_analytical_power"]] <- 1-pnorm(sim$z_ref - abs(sim$z_analytical[2]))
  sim[["power"]] <- rowMeans(replicate(nsims,corr_diff_test(rho = rho, n = n,distr =distr,test = test)[,])<alpha)

  #cat('\r',results$params,results$power)
  cat("\r",sim[["params"]],"\t",sim$power)
  flush.console()
  if(power_only==FALSE) return(sim)
  else return(sim[["power"]])
}

# example usage (defaults: 100 sims, rho1 0.2, rho2 0.5, n1 20, n2 50, alpha 0.05, 2 sided, Pearson)
simulation1 <- corr_power()
simulation1
# $params
# method     rho_1     rho_2        n1        n2       sim 
# "pearson"     "0.2"     "0.5"      "20"      "50"    "TRUE" 
# 
# $log
# z_1           z_2          r_diff       z_diff       z_se      z_test       z_ref    z_power    z_p        
# [1,] -0.3049277    0.1421052    -0.436972    -0.4470329   0.2830197 -1.579512    1.959964 0.3518049  0.1142187  
# [2,] 0.1985209     0.08794177   0.1082376    0.1105791    0.2830197 0.3907117    1.959964 0.05829459 0.6960103  
# [3,] 0.101958      0.0836204    0.0181801    0.01833758   0.2830197 0.06479259   1.959964 0.02903485 0.9483391  
# [4,] 0.3320248     0.09623546   0.2243996    0.2357894    0.2830197 0.8331202    1.959964 0.1299043  0.4047769  
# [5,] -0.07972956   -0.1182992   0.03818939   0.03856966   0.2830197 0.1362791    1.959964 0.03409986 0.8916007  
# ....

compute.power(c(0.2,0.5),c(20,50))$power
compute.power(c(0.2,0.5),c(20,50),power_only=TRUE)

corr_power <- function(n.sims = 100, 
                        n = 700,
                        mzdz_ratio = 0.4,
                        distr = "normal",
                        alpha = 0.05, 
                        sidedness = 2,
                        method  = "pearson",
                        log = TRUE,
                        testsim = TRUE,
                        beta = 0.2,
                        res_min = -0.99,
                        res_max = 0.99,
                        res_inc = 0.05,
                        names = c("Population A","Population B"),
                        lower   = FALSE,
                        simulation = TRUE) {
  # Note: for now, method may be 'peasron,'spearman' or 'kendall' ; plan to include icc
  # 'lower' option is trial of only processing lower matrix
  #  -- takes < 0.5x processing time; result could be mirrored, or plotted against delta
  
  start_time <- Sys.time()
  
  # define target power (to mark on plot)
  target <- 1 - beta

  # Sequance of correlations to compare
  corrs <- seq(res_min, res_max, res_inc)
  results <- list()
  results[["params"]]<-c("method" = method, "distr" = distr,
                         "n1" = ceiling(n*(mzdz_ratio)), "n2" = ceiling(n*(1-mzdz_ratio)), "sim" = simulation, "nsim" = n.sims,"\n")
  cat("Simulation for power to detect a difference in bivariate correlations","\n")
  cat("Correlation method: ",results[["params"]][["method"]],"\n")
  cat("Distribution:","bivariate",results[["params"]][["distr"]],"\n")
  cat("Simulation length: ",results[["params"]][["nsim"]],"\n")
  cat("rho1\trho2\tmz\tdz\tpower\n")
  if (simulation==TRUE) {
    # Evaluate pairwise comparisons across nsims for power estimate
    # mapply(compute.power,rho1 = r, rho2 = c,n1 = n1, n2 = n2, threshold=alpha, nsims = n.sims, lower.tri = lower,   power_only=TRUE))
    results <- outer(corrs, corrs, FUN = function(r, c) mapply(compute.power,rho1 = r, 
                                                                             rho2 = c, 
                                                                             n=n,
                                                                             mzdz_ratio=mzdz_ratio,
                                                                             distr=distr,
                                                                             nsims = n.sims, alpha=alpha, sidedness=sidedness, method=method,
                                                                             power_only=TRUE))
  }
  if (simulation==FALSE){
    ## Actually its quite interesting to plot what the probabilities would look like without sampling; 
    results <- outer(corrs, corrs, FUN = function(r, c) mapply(p_zdiff,rho1 = r, 
                                                                       rho2 = c, 
                                                                       n=n,
                                                                       mzdz_ratio=mzdz_ratio,
                                                                       distr=distr,
                                                                       nsims = n.sims, alpha=alpha, sidedness=sidedness, method=method,
                                                                       power_only=TRUE,simulation = testsim,
                                                                       SIMPLIFY = FALSE))
  }
  # format method to proper case for plot
  corr_type <- stringr::str_to_title(method)
  # plot power simulation results
  fig<-filled.contour(x = corrs, y = corrs, z = as.matrix(results), nlevels = 10, 
                 xlim = c(-1,1), ylim = c(-1,1), zlim = c(0,1),
                 plot.axes = {contour(x = corrs, y = corrs, z = as.matrix(results), 
                                      levels = target, at = seq(-1, 1, 0.2), drawlabels = FALSE, axes = FALSE,
                                      add = TRUE, lwd = 3, col = "steelblue3");
                   abline(v = seq(-1, 1, 0.1), lwd = .5, col = "lightgray", lty = 2)
                   abline(h = seq(-1, 1, 0.1), lwd = .5, col = "lightgray", lty = 2)
                   axis(1, seq(-1,1,0.2))
                   axis(2, seq(-1,1,0.2)) },
                 plot.title = title(main = paste0("Power for difference in ",corr_type," correlations","\n",
                                    "Mz = ",results[["params"]][["n1"]],", Dz = ",results[["params"]][["n2"]],", alpha: ",alpha, ", sims: ",n.sims),
                                    xlab = paste0("Correlation in ",names[1]),
                                    ylab = paste0("Correlation in ",names[2]), adj = 0),
                 color.palette =  colorRampPalette(c("#f7fcf0","#525252")));
        arrows(0.63, 0.6, 0.845, 0.6, length = 0.14, lwd = 3, col = "steelblue3")
  

  end_time <- Sys.time()
  
  # display running time
  cat(paste0("Power for difference in ",corr_type," correlations","\n",
             "n1 = ",n[1],", n2 = ",n[2],", alpha: ",alpha, ", sims: ",n.sims,"\n",
             "Processing time: ",end_time - start_time,"\n\n"))
  return(fig)
}


power_4to6 <- corr_power(n.sims    = 100,
                     n       =  700,
                     mzdz_ratio = 0.4,
                     distr   = "normal",
                     alpha   =   0.05, 
                     beta    =   0.2,
                     names   =  c("Mz twins","Dz twins"),
                     method  =  "pearson",
                     lower   =   TRUE)

# an example, iterating over different kinds of correlation
#  -- Although I am not certain implementation is correct for latter two yet.
p <- list()
for(corr in c("pearson","spearman")){
  p[[corr]]  <- corr_power(n.sims  = 100,
                           n1      =  134, 
                           n2      =  134, 
                           alpha   =   0.05, 
                           beta    =   0.2,
                           n1_name = "Mz twins",
                           n2_name = "Dz twins",
                           method  =  corr,
                           lower   =   FALSE) 
}               


# plot combinations of result processing parameters (correlations, sim size, and alpha)
correlations <- c("pearson","spearman")
results <- list()
for(corr in correlations){
  for(n1 in c(30,150)){
    for(n2 in c(30,150)){
results[[corr]][[paste0("n1_",n1)]][[paste0("n2_",n2)]] <- corr_power(n.sims  = 100,
                                                                      n1      =  n1, 
                                                                      n2      =  n2, 
                                                                      alpha   =  0.05, 
                                                                      beta    =  0.2,
                                                                      n1_name =  "Mz twins",
                                                                      n2_name =  "Dz twins",
                                                                      method  =  corr,
                                                                      lower   =  FALSE) 
    }
  }
}


# trick to remove warnings while debugging
assign("last.warning", NULL, envir = baseenv())


coloursets <- list()
coloursets[["pgrn"]] <- c('#276419','#4d9221','#7fbc41','#b8e186','#e6f5d0','#fde0ef','#f1b6da','#de77ae','#c51b7d','#8e0152','#f7f7f7')
coloursets[["BrBG"]] <- c('#003c30','#01665e','#35978f','#80cdc1','#c7eae5','#f5f5f5','#f6e8c3','#dfc27d','#bf812d','#8c510a','#543005')
coloursets[["RdYlBu"]] <- c('#313695','#4575b4','#74add1','#abd9e9','#e0f3f8','#ffffbf','#fee090','#fdae61','#f46d43','#d73027','#a50026')


# Simulation data

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
