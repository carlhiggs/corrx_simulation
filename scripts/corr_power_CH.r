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




p_zdiff <- function(rho = c(.2,.5), n = 700, mzdz_ratio = 0.4, distr = "normal",
                    alpha = 0.05, sidedness = 2,
                    method = "pearson", log = TRUE, output  = "z_p",simulation = TRUE,
                    mu1 = c(0,0), mu2 = c(0,0),param1 = c(1,1), param2 = c(1,1)) {
  # print(paste(rho, n, mu1, mu2, param1, param2, distr, alpha, sidedness, method,log, simulation))  
  require("simstudy")
  sim <- list()
  sim[["z_method"]] <- method
  sim[["rho_1"]]  <- rho[1]
  sim[["rho_2"]]  <- rho[2]   
  sim[["n1"]]     <- n*mzdz_ratio
  sim[["n2"]]     <- n*(1-mzdz_ratio)
  # print(matrix(c(1, sim$rho_1, sim$rho_1, 1)))
  # print(matrix(c(1, sim$rho_2, sim$rho_2, 1)))
  
  if(simulation==TRUE){
    sim[["sim"]] <- TRUE
    sim[["z_1"]]    <- atanh(cor(genCorGen(sim$n1, nvars = 2, params1 = mu1, params2 = param1,  
                                           dist = distr, corMatrix = matrix(c(1, sim$rho_1, sim$rho_1, 1), ncol = 2), wide = TRUE)[,2:3], 
                                           method = method)[1,2])
    sim[["z_2"]]    <- atanh(cor(genCorGen(sim$n2, nvars = 2, params1 = mu2, params2 = param2,  
                                           dist = distr, corMatrix = matrix(c(1, sim$rho_2, sim$rho_2, 1), ncol = 2), wide = TRUE)[,2:3], 
                                           method = method)[1,2])
    sim[["r_diff"]] <- tanh(sim[["z_1"]]) - tanh(sim[["z_2"]])
  } else {
    sim[["sim"]] <- FALSE
    sim[["z_1"]]    <- atanh(sim$rho_1)
    sim[["z_2"]]    <- atanh(sim$rho_2)
    sim[["r_diff"]] <- sim$rho_1 - sim$rho_2
  }
  sim[["z_diff"]] <- sim[["z_1"]] - sim[["z_2"]]
  sim[["z_se"]]   <- sqrt(1/(sim$n1-3) + 1/(sim$n2-3))
  sim[["z_test"]] <- sim[["z_diff"]]/sim[["z_se"]]
  sim[["z_ref"]]  <- qnorm(1-alpha/sidedness)
  sim[["z_power"]] <- 1-pnorm(sim[["z_ref"]]-abs(sim[["z_test"]]))
  sim[["z_p"]]    <- sidedness*pnorm(-abs(sim[["z_test"]]))
  # print(paste(rho, n, mu1, mu2, param1, param2, distr, alpha, sidedness, method,log, simulation, sim$z_power,"/n"))
  if (log==FALSE){
    return(sim[[output]]) 
  }
  return(sim)
}

rbind(p_zdiff(rho=c(0.2,0.5),n = c(20,50  ), simulation = FALSE),
      p_zdiff(rho=c(0.2,0.5),n = c(200,500), method="pearson", simulation = FALSE),
      p_zdiff(rho=c(0.2,0.5),n = c(20,50  ), method="pearson" ),
      p_zdiff(rho=c(0.2,0.5),n = c(200,500), method="pearson" ),
      p_zdiff(rho=c(0.2,0.5),n = c(20,50  ), method="spearman", simulation = FALSE),
      p_zdiff(rho=c(0.2,0.5),n = c(200,500), method="spearman", simulation = FALSE),
      p_zdiff(rho=c(0.2,0.5),n = c(20,50  ), method="spearman"),
      p_zdiff(rho=c(0.2,0.5),n = c(200,500), method="spearman"),
      p_zdiff(rho=c(0.2,0.5),n = c(20,50  ), method="kendall", simulation = FALSE),
      p_zdiff(rho=c(0.2,0.5),n = c(200,500), method="kendall", simulation = FALSE),
      p_zdiff(rho=c(0.2,0.5),n = c(20,5   ), method="kendall" ),
      p_zdiff(rho=c(0.2,0.5),n = c(200,500), method="kendall" ))
#       z_method   rho_1 rho_2 n1  n2  sim   z_1        z_2       r_diff     z_diff     z_se       z_test     z_ref    z_power    z_p         
#  [1,] "pearson"  0.2   0.5   20  50  FALSE 0.2027326  0.5493061 -0.3       -0.3465736 0.2830197  -1.224557  1.959964 0.2310457  0.2207423   
#  [2,] "pearson"  0.2   0.5   200 500 FALSE 0.2027326  0.5493061 -0.3       -0.3465736 0.08419154 -4.11649   1.959964 0.9844787  3.846864e-05
#  [3,] "pearson"  0.2   0.5   20  50  TRUE  0.1717027  0.5955662 -0.363852  -0.4238635 0.2830197  -1.497647  1.959964 0.3219269  0.134225    
#  [4,] "pearson"  0.2   0.5   200 500 TRUE  0.2824285  0.4998857 -0.1868761 -0.2174572 0.08419154 -2.582887  1.959964 0.7333324  0.009797744 
#  [5,] "spearman" 0.2   0.5   20  50  FALSE 0.2027326  0.5493061 -0.3       -0.3465736 0.2830197  -1.224557  1.959964 0.2310457  0.2207423   
#  [6,] "spearman" 0.2   0.5   200 500 FALSE 0.2027326  0.5493061 -0.3       -0.3465736 0.08419154 -4.11649   1.959964 0.9844787  3.846864e-05
#  [7,] "spearman" 0.2   0.5   20  50  TRUE  0.02105574 0.60847   -0.5219966 -0.5874143 0.2830197  -2.075525  1.959964 0.5459997  0.03793793  
#  [8,] "spearman" 0.2   0.5   200 500 TRUE  0.3152196  0.5729482 -0.212343  -0.2577286 0.08419154 -3.061217  1.959964 0.8646068  0.002204391 
#  [9,] "kendall"  0.2   0.5   20  50  FALSE 0.2027326  0.5493061 -0.3       -0.3465736 0.2830197  -1.224557  1.959964 0.2310457  0.2207423   
# [10,] "kendall"  0.2   0.5   200 500 FALSE 0.2027326  0.5493061 -0.3       -0.3465736 0.08419154 -4.11649   1.959964 0.9844787  3.846864e-05
# [11,] "kendall"  0.2   0.5   20  5   TRUE  0.137706   0.4236489 -0.2631579 -0.2859429 0.747545   -0.3825093 1.959964 0.05734547 0.7020836   
# [12,] "kendall"  0.2   0.5   200 500 TRUE  0.1685118  0.3645807 -0.1823078 -0.196069  0.08419154 -2.328844  1.959964 0.6438914  0.01986733  


# Compute power 
compute.power <- function(rho = c(.2,.5), n = 700, mzdz_ratio = 0.4,distr = "normal",
                          alpha = 0.05, sidedness=2,method="pearson", rho1 = NULL, rho2 = NULL, n1 = NULL, n2 = NULL, 
                          nsims = 1000,lower.tri = FALSE,log=TRUE, simulation=TRUE,power_only = FALSE,
                          mu1 = c(0,0),mu2 = c(0,0),param1 = c(1,1), param2 = c(1,1)) {
  if((is.null(rho1)==FALSE)&&(is.null(rho2)==FALSE)) rho <- c(rho1,rho2)
  # if((is.null(mu1a)==FALSE)&&(is.null(mu1b)==FALSE)) mu1 <- c(mu1a,mu1b)
  # if((is.null(mu2a)==FALSE)&&(is.null(mu2b)==FALSE)) mu2 <- c(mu2a,mu2b)
  # if((is.null(param1a)==FALSE)&&(is.null(param1b)==FALSE)) param1 <- c(param1a,param1b)
  # if((is.null(param2a)==FALSE)&&(is.null(param2b)==FALSE)) param2 <- c(param2a,param2b)
  # if((is.null(distr1)==FALSE)&&(is.null(distr2)==FALSE)) distr <- c(distr1,distr2)
  if(lower.tri==TRUE){
    # only calculate lower matrix half when comparing across all correlation combinations
    if(rho[1] < rho[2]) { 
      return(NA)
    }
  }
  results <- list()
  results[["params"]]<-c("method" = method, "rho_1" = rho[1], "rho_2" = rho[2], 
                         "n1" = ceiling(n*(mzdz_ratio)), "n2" = ceiling(n*(1-mzdz_ratio)), "sim" = simulation)
  results[["log"]] <- t(replicate(nsims, p_zdiff(rho=rho, n=n,  distr=distr, alpha=alpha, sidedness=sidedness, method=method,log=log, simulation=simulation,
                                                 mu1=mu1, mu2=mu2, param1=param1, param2=param2))[7:15,])
  results[["power"]]<- mean(unlist(results[["log"]][,"z_p"]) < alpha)
  #cat('\r',results$params,results$power)
  cat("\r",rho[1],"\t",rho[2],"\t",results[["params"]][["n1"]],"\t",results[["params"]][["n2"]],"\t",results$power)
  flush.console()
  if(power_only==FALSE) return(results)
  else return(results[["power"]])
}

# example usage (defaults: 1000 sims, rho1 0.5, rho2 0.2, n1 20, n2 50, alpha 0.05, 2 sided, Pearson)
simulation1 <- compute.power()
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
