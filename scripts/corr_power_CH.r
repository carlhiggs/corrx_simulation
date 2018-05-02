# R script to simulate power for difference in correlations (Pearson, Spearman, and later ... ICC)
# Carl Higgs 2017
#
# Initial inspiration for corr_power_plot from https://janhove.github.io/design/2015/04/14/power-simulations-for-comparing-independent-correlations
# Jan Vanhove 14 April 2015


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
# To develop
#    - accounting for clustering (ie. icc in twin studies)
#    - allow parameterisation to explore change in source pop mean and sd (unequal variance, etc)


# Testing out compilation of genCorGen according to directions at 
# https://www.r-statistics.com/2012/04/speed-up-your-r-code-using-a-just-in-time-jit-compiler/
# on advice from Koen Simons
require(compiler)
install.packages("simstudy")
require("simstudy")
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

# Modified SIgned Log-likelihood Ratio test (Kazemi & Jafari, 2015; maybe based on Barndorff Nielsen 1991
# note limitation that Rf (here, rstar.f)
mslr <- function(a,b,M=1e4,sidedness=2,method = "pearson") {
  # Two samples
  n1 <- nrow(a)
  n2 <- nrow(b)
  
  # Step 1: Compute sample correlation coefficients
  r1      <- cor(a,method = method)[2,1]
  r2      <- cor(b,method = method)[2,1]
  z1     <- atanh(r1)
  z2     <- atanh(r2)
  rf     <- tanh(mean(z1,z2))
  
  # Step 2: Generate random numbers
  V2     <- matrix(data=0, nrow = M, ncol = 2)
  V2[,1] <- rchisq(M, df = n1-1, ncp = 1)
  V2[,2] <- rchisq(M, df = n2-1, ncp = 1)
  
  W2     <- matrix(data=0, nrow = M, ncol = 2)
  W2[,1] <- rchisq(M, df = n1-2, ncp = 1)
  W2[,2] <- rchisq(M, df = n2-2, ncp = 1)
  
  N <-matrix(data = rnorm(2*M), nrow=M, ncol = 2)
  
  
  # Compute test statistic
  rstar.f  <- rf/sqrt(1-rf^2)
  rstar.top <- rstar.f * c(V2[,1],V2[,2]) + N
  rstar     <- rstar.top / sqrt( rstar.top^2 + W2^2 )

  # non looping approach (but otherwise as per formula) which does not result in vector of slr
  #  but a vector is required to calculate slr mean and variance
  slr <- sign(rstar[,1]-rstar[,2])*sqrt(sum(c(n1*log(((1-rf*r1)^2)/((1-r1^2)*(1-rf^2))),
                                              n2*log(((1-rf*r2)^2)/((1-r2^2)*(1-rf^2))))))
  
  ## non-working attempt using rstar to recalculate rf as vector
  # rf     <- tanh(rowMeans(atanh(rstar)))
  # slr <- numeric(length = 0)
  # for (i in 1:length(rf))   slr <- c(slr,sign(r1-r2)*sqrt(sum(c(n1*log(((1-rf[i]*r1)^2)/((1-r1^2)*(1-rf[i]^2))),
  #                                                              n1*log(((1-rf[i]*r2)^2)/((1-r1^2)*(1-rf[i]^2)))))))
  
  ## non-working attempt, using rstar as vector
  # slr <- sign(rstar[,1]-rstar[,2])*sqrt(sum(c(n1*log(((1-rf*rstar[,1])^2)/((1-rstar[,1]^2)*(1-rf^2))),
  #                                            n2*log(((1-rf*rstar[,2])^2)/((1-rstar[,2]^2)*(1-rf^2))))))

  # 6. Repeat steps 3-5 for a large number of times (say M = 10,000).
  ###### Shouldn't that be 2-5, so the RVs are redrawn?? Or is it just each i of length M?
  # 7. Compute the sample mean and sample variance of SLR and compute   the MSLR in (9).
  mslr      <- (slr - mean(slr))/sqrt(var(slr))
  # 8. Determine the p-value 2 as:  2 * pnorm(abs(mslr))
  # Compute p value
  p    <- 2 * (1 - pnorm(abs(mslr))); 
  return(p)
}
  


a <- genCorGen(50, nvars = 2, params1 = 0, params2 = 1, dist = "normal", 
          corMatrix = matrix(c(1, 0.5, 0.5, 1), ncol = 2), wide = TRUE)[,2:3]
b <- genCorGen(50, nvars = 2, params1 = 0, params2 = 1, dist = "normal", 
             corMatrix = matrix(c(1, 0.2, 0.2, 1), ncol = 2), wide = TRUE)[,2:3]
result <- list()
system.time(result[["fz"]]       <- fz(a,b))
system.time(result[["fz_nosim"]] <- fz_nosim(0.5,0.5,50,50))
system.time(result[["gtv"]]      <- gtv(a,b))
mslr(a,b, M=10000)

fz_compiled <- cmpfun(fz)
fz_ns_compiled <- cmpfun(fz_nosim)
gtv_compiled <- cmpfun(gtv)
fc_fz     <- function() for (i in 1:1000)      fz(a,b)
fc_fz_ns  <- function() for (i in 1:1000)  fz_nosim(0.5,0.5,50,50)
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
                    param1a = c(0,0), param1b = c(0,0),param2a = c(1,1), param2b = c(1,1),
                    alpha = 0.05, sidedness = 2, test = c("fz","gtv"),
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
    results[["fz_nosim"]] <- fz_nosim(rho[1],rho[2],n[1],n[2], 
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
  if ("gtv"      %in% test) results[["gtv"]]      <- gtv_compiled(a,b)
  return(rbind(results[test]))
}
corr_diff_test(rho = c(.2,.4), n = c(100,300), distr = "normal",test = "fz_nosim") 
corr_diff_test(rho = c(.2,.4), n = c(100,300), distr = "normal",test = "fz") 
corr_diff_test(rho = c(.2,.4), n = c(100,300), distr = "normal",test = "gtv") 
corr_diff_test(rho = c(.2,.4), n = c(100,300), distr = "normal") 

corr_diff_test_compiled <- cmpfun(corr_diff_test)
cdt  <- function() for (i in 1:1000)  corr_diff_test(rho = c(.2,.4), n = c(100,300), distr = "normal",test = "fz") 
cdt_c  <- function() for (i in 1:1000)  corr_diff_test_compiled(rho = c(.2,.4), n = c(100,300), distr = "normal",test = "fz") 
system.time(cdt())
system.time(cdt_c())
## No apparent speed benefits to compilation at the individual test stage
# > system.time(cdt())
# user  system elapsed 
# 34.62    1.12   42.71 
# > system.time(cdt_c())
# user  system elapsed 
# 35.02    1.22   44.27 

# Corr power simulation
corr_power <- function(rho = c(.2,.5), n = c(30,90),distr = "normal",
                       param1a = c(0,0), param1b = c(0,0),param2a = c(1,1), param2b = c(1,1),
                       test = c("fz","gtv"),
                       alpha = 0.05, sidedness=2,method="pearson",  
                       nsims = 100,lower.tri = FALSE, power_only = FALSE){
  sim <- list()
  sim[["params"]] <- c("method" = method, "rho_1" = rho[1], "rho_2" = rho[2],"n1" = n[1], "n2" = n[2],
                       "alpha" = alpha, "sidedness" = sidedness, "nsims" = nsims, "distr" = distr)
  
  sim[["additional"]] <- c("param1a" = param1a, "param1b" = param1b,"param2a" = param2a, "param2b" = param2b)
  sim[["analytical"]] <- c("fz_nosim" = fz_nosim(rho[1],rho[2],n[1],n[2], 
                                    alpha = alpha, sidedness = sidedness, method = method, power = TRUE))
  sim_tests <- test[!test %in% "fz_nosim"]
  sim[["power"]] <- rowMeans(replicate(nsims,corr_diff_test(rho = rho, n = n,distr =distr,
                                                            param1a = param1a, param1b = param1b,param2a = param2a, param2b = param2b,
                                                            alpha = alpha, sidedness=sidedness,method=method, 
                                                            test = sim_tests)[,])<alpha)
  if ("fz_nosim" %in% test) sim[["power"]] <- c(sim[["analytical"]],sim[["power"]])[test]

  #cat('\r',results$params,results$power)
  cat("\r                                                                \r",sim[["params"]],sim$power,sep="\t")
  flush.console()
  if(power_only==FALSE) return(sim)
  else return(sim[["power"]])
}

corr_power_compiled <- cmpfun(corr_power)
cp  <- function() for (i in 1:10)  corr_power(power_only=TRUE)
cp_c  <- function() for (i in 1:10)  corr_power_compiled(power_only=TRUE)
# system.time(cp())
# system.time(cp_c())
# > system.time(cp())
# pearson 0.2 0.5 30 90 0.05 2 100 normal 	 0.38 0.4   
# user  system elapsed 
# 266.84   12.50  313.37 
# > system.time(cp_c())
# pearson 0.2 0.5 30 90 0.05 2 100 normal 	 0.35 0.37
# user  system elapsed 
# 263.76   12.39  307.74 
## Pretty similar timing; still, i'll use just in case.

corr_power_compiled(power_only=TRUE, nsim = 10, test=c("fz","fz_nosim","gtv"))
# fz  fz_nosim       gtv                                     
# 0.1000000 0.3494663 0.1000000 
## NOTE - issue if only one simulation is processed - the rowMeans subcommand bugs out
# corr_power_compiled(power_only=TRUE, nsim = 10, test=c("fz","fz_nosim"))
# Error in rowMeans(replicate(nsims, corr_diff_test(rho = rho, n = n, distr = distr,  : 
# 'x' must be an array of at least two dimensions


corr_power_plot <- function(nsims = 100, res_min = -0.95, res_max = 0.95, res_inc = 0.05, 
                            n = c(30,90),distr = "normal",
                            param1a = c(0,0), param1b = c(0,0),param2a = c(1,1), param2b = c(1,1),
                            tests = c("fz_nosim","fz","gtv"),
                            alpha = 0.05, beta = 0.2, sidedness=2,method="pearson",  
                            names = c("Population A","Population B"),
                            lower.tri = FALSE){
  cat("\n","Correlation power plot simulation commenced at",as.character(Sys.time()),"\n")
  results <- list()
  results[["params"]] <- c("method" = method,"n1" = n[1], "n2" = n[2],
                           "param1a" = param1a, "param1b" = param1b,"param2a" = param2a, "param2b" = param2b,
                           "alpha" = alpha, "sidedness" = sidedness, 
                           "nsims" = nsims, "distr" = distr)
  results[["tests"]]  <- tests
  results[["z_ref"]]  <- qnorm(1-alpha/sidedness)
  corrs <- round(seq(res_min, res_max, res_inc),2)
  
  # create result holder matrices per test
  for (test in results[["tests"]]) {
    results[[test]] <- matrix(data = NA, nrow = length(corrs), ncol = length(corrs))
    colnames(results[[test]]) <- rownames(results[[test]]) <- format(corrs, trim=TRUE)
  }
  # print header for corr_power output
  cat("\tmethod","rho_1","rho_2","n1","n2","alpha","sides","nsims","distr","PowerXtests","\n",sep="\t")
  cat("\t",rep("-",45),"\n")
  # calculate pairwise results
  for (r in corrs){
    for (c in corrs){
      temp <- corr_power_compiled(rho = c(r,c), n = c(n[1],n[2]),distr = distr,
                                  param1a = param1a, param1b = param1b,param2a = param2a, param2b = param2b,
                                  test = tests,
                                  alpha = alpha, sidedness=sidedness,method=method,  
                                  nsims = nsims,lower.tri = lower.tri, power_only = TRUE)  
      for (test in results[["tests"]]) {
        results[[test]][r == corrs,c == corrs] <- temp[test]
      }
    }         
  }
  
  # format method to proper case for plot
  corr_type <- stringr::str_to_title(method)
  # define target power (to mark on plot)
  target <- 1 - beta

  # plot power simulation results
   for (test in results[["tests"]]) {
    results[["fig"]][[test]]<-filled.contour(x = corrs, y = corrs, z = as.matrix(results[[test]]), nlevels = 10,
                                xlim = c(-1,1), ylim = c(-1,1), zlim = c(0,1),
                                plot.axes = {contour(x = corrs, y = corrs, z = as.matrix(results[[test]]),
                                                     levels = target, at = seq(-1, 1, 0.2), drawlabels = FALSE, axes = FALSE,
                                                     add = TRUE, lwd = 3, col = "steelblue3");
                                  abline(v = seq(-1, 1, 0.1), lwd = .5, col = "lightgray", lty = 2)
                                  abline(h = seq(-1, 1, 0.1), lwd = .5, col = "lightgray", lty = 2)
                                  axis(1, seq(-1,1,0.2))
                                  axis(2, seq(-1,1,0.2))},
                                plot.title = title(main = paste0(test," test","\n",
                                                                 "Mz = ",results[["params"]][["n1"]],
                                                                 ", Dz = ",results[["params"]][["n2"]],
                                                                 ", alpha: ",alpha, ", sims: ",nsims),
                                                   xlab = paste0("Correlation in ",names[1]),
                                                   ylab = paste0("Correlation in ",names[2]), adj = 0),
                                color.palette =  colorRampPalette(c("#f7fcf0","#525252")));
   # arrows(0.63, 0.6, 0.845, 0.6, length = 0.14, lwd = 3, col = "steelblue3")                          
   }
  cat("\n","Completed at ",as.character(Sys.time()),"\n")
  return(results)
}

corr_pplot_compiled <- cmpfun(corr_power_plot)

system.time(results<- corr_pplot_compiled(nsims = 10, res_min = -.3, res_max = 0.3, res_inc = 0.1, n = c(30,90)))
# Correlation power plot simulation commenced at 2018-05-01 21:57:57 
# method	rho_1	rho_2	n1	n2	alpha	sides	nsims	distr	PowerXtests	
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#   pearson	0.3	0.3	30	90	0.05	2	10	normal	0.025	0	0
# Completed at  2018-05-01 21:59:00 
# user  system elapsed 
# 58.14    5.07   63.17

system.time(results<- corr_pplot_compiled())
## Note - this time test was run on more powerful work computer
# Correlation power plot simulation commenced at 2018-05-01 22:05:06 
# method	rho_1	rho_2	n1	n2	alpha	sides	nsims	distr	PowerXtests	
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#   pearson	0.95	0.95	30	90	0.05	2	100	normal	0.025	0.09	0.09
# Completed at  2018-05-02 02:06:36 
# user   system  elapsed 
# 14425.09    94.50 14489.76
#
# There was some error viewing all plots some time later - not sure why; however, results should be saved
corrs<- seq(-0.95,0.95,0.05)
nsims = 100 
n = c(30,90)
distr = "normal"
tests = c("fz_nosim","fz","gtv")
alpha = 0.05
beta = 0.2
sidedness=2
target = 1-beta
method="pearson"  
names = c("Population A","Population B")
for (test in results[["tests"]]) {
results[["fig"]][[test]]<-filled.contour(x = corrs, y = corrs, z = as.matrix(results[[test]]), nlevels = 10,
                                         xlim = c(-1,1), ylim = c(-1,1), zlim = c(0,1),
                                         plot.axes = {contour(x = corrs, y = corrs, z = as.matrix(results[[test]]),
                                                              levels = target, at = seq(-1, 1, 0.2), drawlabels = FALSE, axes = FALSE,
                                                              add = TRUE, lwd = 3, col = "steelblue3");
                                           abline(v = seq(-1, 1, 0.1), lwd = .5, col = "lightgray", lty = 2)
                                           abline(h = seq(-1, 1, 0.1), lwd = .5, col = "lightgray", lty = 2)
                                           axis(1, seq(-1,1,0.2))
                                           axis(2, seq(-1,1,0.2))},
                                         plot.title = title(main = paste0(test," test","\n",
                                                                          "Mz = ",results[["params"]][["n1"]],
                                                                          ", Dz = ",results[["params"]][["n2"]],
                                                                          ", alpha: ",alpha, ", sims: ",nsims),
                                                            xlab = paste0("Correlation in ",names[1]),
                                                            ylab = paste0("Correlation in ",names[2]), adj = 0),
                                         color.palette =  colorRampPalette(c("#f7fcf0","#525252")));
arrows(0.63, 0.6, 0.845, 0.6, length = 0.14, lwd = 3, col = "steelblue3")    
}



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
fz(gamma_a,gamma_b)
gtv(gamma_a,gamma_b)   
corr_power(dist = "gamma",param1a = c(1.5,1.5), param1b = c(1.5,1.5), param2a = c(.09,.09), param2b = c(.09,.09))
# pearson	0.2	0.5	30	90	0.05	2	100	gamma	0.37	0.37$params        
# method     rho_1     rho_2        n1        n2     alpha sidedness     nsims     distr 
# "pearson"     "0.2"     "0.5"      "30"      "90"    "0.05"       "2"     "100"   "gamma" 
# 
# $additional
# param1a1 param1a2 param1b1 param1b2 param2a1 param2a2 param2b1 param2b2 
# 1.50     1.50     1.50     1.50     0.09     0.09     0.09     0.09 
# 
# $analytical
# fz_nosim 
# 0.3494663 
# 
# $power
# fz  gtv 
# 0.37 0.37 

# simulate using non-normal, positively skewed distribution
system.time(res_gamma<- corr_pplot_compiled(dist = "gamma",param1a = c(1.5,1.5), param1b = c(1.5,1.5),param2a = c(.09,.09),param2b = c(.09,.09)))






#other mslrt exploration
# param1 are the probabilities
binary <- genCorGen(50, nvars = 2, params1 = c(.3, .5), dist = "binary", 
                    corMatrix = matrix(c(1, 0.8, 0.8, 1), ncol = 2), wide = TRUE)


