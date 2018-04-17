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

p_zdiff <- function(rho1, rho2, n1, n2, alpha = 0.05, sidedness = 2,
                    method = "pearson", log = TRUE, output  = "z_p",simulation = TRUE) {
  require("MASS")
  sim <- list()
  sim[["z_method"]] <- method
  sim[["rho_1"]]  <- rho1
  sim[["rho_2"]]  <- rho2    
  sim[["n1"]]     <- n1
  sim[["n2"]]     <- n2
  if(simulation==TRUE){
    sim[["sim"]] <- TRUE
    sim[["z_1"]]    <- atanh(cor(mvrnorm(n1, mu = c(0, 0),Sigma = matrix(c(1, rho1, rho1, 1), ncol = 2)), method = method)[1,2])
    sim[["z_2"]]    <- atanh(cor(mvrnorm(n2, mu = c(0, 0),Sigma = matrix(c(1, rho2, rho2, 1), ncol = 2)), method = method)[1,2])
    sim[["r_diff"]] <- tanh(sim[["z_1"]]) - tanh(sim[["z_2"]])
  } else {
    sim[["sim"]] <- FALSE
    sim[["z_1"]]    <- atanh(rho1)
    sim[["z_2"]]    <- atanh(rho2)
    sim[["r_diff"]] <- rho1 - rho2
  }
  sim[["z_diff"]] <- sim[["z_1"]] - sim[["z_2"]]
  sim[["z_se"]]   <- sqrt(1/(n1-3) + 1/(n2-3))
  sim[["z_test"]] <- sim[["z_diff"]]/sim[["z_se"]]
  sim[["z_ref"]]  <- qnorm(1-alpha/sidedness)
  sim[["z_power"]] <- 1-pnorm(sim[["z_ref"]]-abs(sim[["z_test"]]))
  sim[["z_p"]]    <- sidedness*pnorm(-abs(sim[["z_test"]]))
  if (log==FALSE){
    return(sim[[output]]) 
  }
  return(sim)
}

rbind(p_zdiff(0.2,0.5,20,50, simulation = FALSE),
      p_zdiff(0.2,0.5,200,500, method="pearson", simulation = FALSE),
      p_zdiff(0.2,0.5,20,50,   method="pearson" ),
      p_zdiff(0.2,0.5,200,500, method="pearson" ),
      p_zdiff(0.2,0.5,20,50,   method="spearman", simulation = FALSE),
      p_zdiff(0.2,0.5,200,500, method="spearman", simulation = FALSE),
      p_zdiff(0.2,0.5,20,50,   method="spearman"),
      p_zdiff(0.2,0.5,200,500, method="spearman"),
      p_zdiff(0.2,0.5,20,50,   method="kendall", simulation = FALSE),
      p_zdiff(0.2,0.5,200,500, method="kendall", simulation = FALSE),
      p_zdiff(0.2,0.5,20,50,   method="kendall" ),
      p_zdiff(0.2,0.5,200,500, method="kendall" ))
#       z_method   rho_1 rho_2 n1  n2  sim   z_1        z_2       r_diff     z_diff     z_se       z_test    z_ref    z_power   z_p         
#  [1,] "pearson"  0.2   0.5   20  50  FALSE 0.2027326  0.5493061 -0.3       -0.3465736 0.2830197  -1.224557 1.959964 0.2310457 0.2207423   
#  [2,] "pearson"  0.2   0.5   200 500 FALSE 0.2027326  0.5493061 -0.3       -0.3465736 0.08419154 -4.11649  1.959964 0.9844787 3.846864e-05
#  [3,] "pearson"  0.2   0.5   20  50  TRUE  0.08106579 0.5163649 -0.3940009 -0.4352991 0.2830197  -1.538053 1.959964 0.3365449 0.1240357   
#  [4,] "pearson"  0.2   0.5   200 500 TRUE  0.1877089  0.6038379 -0.35424   -0.4161291 0.08419154 -4.942647 1.959964 0.9985713 7.706881e-07
#  [5,] "spearman" 0.2   0.5   20  50  FALSE 0.2027326  0.5493061 -0.3       -0.3465736 0.2830197  -1.224557 1.959964 0.2310457 0.2207423   
#  [6,] "spearman" 0.2   0.5   200 500 FALSE 0.2027326  0.5493061 -0.3       -0.3465736 0.08419154 -4.11649  1.959964 0.9844787 3.846864e-05
#  [7,] "spearman" 0.2   0.5   20  50  TRUE  0.1147872  0.4721117 -0.3256182 -0.3573245 0.2830197  -1.262543 1.959964 0.2427697 0.2067535   
#  [8,] "spearman" 0.2   0.5   200 500 TRUE  0.1908226  0.529357  -0.2963497 -0.3385343 0.08419154 -4.021002 1.959964 0.9803503 5.79511e-05 
#  [9,] "kendall"  0.2   0.5   20  50  FALSE 0.2027326  0.5493061 -0.3       -0.3465736 0.2830197  -1.224557 1.959964 0.2310457 0.2207423   
# [10,] "kendall"  0.2   0.5   200 500 FALSE 0.2027326  0.5493061 -0.3       -0.3465736 0.08419154 -4.11649  1.959964 0.9844787 3.846864e-05
# [11,] "kendall"  0.2   0.5   20  50  TRUE  0.1163111  0.444228  -0.3013534 -0.3279168 0.2830197  -1.158636 1.959964 0.211471  0.2466045   
# [12,] "kendall"  0.2   0.5   200 500 TRUE  0.05281287 0.3782312 -0.3084065 -0.3254183 0.08419154 -3.865214 1.959964 0.9716262 0.0001109919

# Compute power 
compute.power <- function(rho1 = 0.5, rho2 = 0.2, n1 = 20, n2 = 50, 
                          threshold = 0.05, sidedness=2,method="pearson",
                          nsims = 1000,lower.tri = FALSE,log=TRUE, simulation=TRUE,power_only = FALSE) {
  if(lower.tri==TRUE){
    # only calculate lower matrix half when comparing across all correlation combinations
    if(rho1 < rho2) { 
      return(NA)
    }
  }
  results <- list()
  results[["params"]]<-c("method" = method, "rho_1" = rho1, "rho_2" = rho2, "n1" = n1, "n2" = n2, "sim" = simulation)
  results[["log"]] <- t(replicate(nsims, p_zdiff(rho1 = rho1, rho2 = rho2, n1 = n1, n2 = n2,
                                               alpha = threshold,sidedness = sidedness,
                                               method = method,log = log, simulation = simulation))[7:15,])
  results[["power"]]<- mean(unlist(results[["log"]][,"z_p"]) < threshold)
  if(power_only==FALSE) return(results)
  else return(results[["power"]])
}

# example usage (defaults: 1000 sims, rho1 0.5, rho2 0.2, n1 20, n2 50, alpha 0.05, 2 sided, Pearson)
simulation1 <- compute.power()
#       z_method   rho_1 rho_2 n1  n2  sim   z_1        z_2       r_diff     z_diff     z_se       z_test    z_ref    z_power   z_p         
#  [1,] "pearson"  0.2   0.5   20  50  FALSE 0.2027326  0.5493061 -0.3       -0.3465736 0.2830197  -1.224557 1.959964 0.2310457 0.2207423   
#  [2,] "pearson"  0.2   0.5   200 500 FALSE 0.2027326  0.5493061 -0.3       -0.3465736 0.08419154 -4.11649  1.959964 0.9844787 3.846864e-05
#  [3,] "pearson"  0.2   0.5   20  50  TRUE  0.08106579 0.5163649 -0.3940009 -0.4352991 0.2830197  -1.538053 1.959964 0.3365449 0.1240357   
#  [4,] "pearson"  0.2   0.5   200 500 TRUE  0.1877089  0.6038379 -0.35424   -0.4161291 0.08419154 -4.942647 1.959964 0.9985713 7.706881e-07
#  [5,] "spearman" 0.2   0.5   20  50  FALSE 0.2027326  0.5493061 -0.3       -0.3465736 0.2830197  -1.224557 1.959964 0.2310457 0.2207423   
#  [6,] "spearman" 0.2   0.5   200 500 FALSE 0.2027326  0.5493061 -0.3       -0.3465736 0.08419154 -4.11649  1.959964 0.9844787 3.846864e-05
#  [7,] "spearman" 0.2   0.5   20  50  TRUE  0.1147872  0.4721117 -0.3256182 -0.3573245 0.2830197  -1.262543 1.959964 0.2427697 0.2067535   
#  [8,] "spearman" 0.2   0.5   200 500 TRUE  0.1908226  0.529357  -0.2963497 -0.3385343 0.08419154 -4.021002 1.959964 0.9803503 5.79511e-05 
#  [9,] "kendall"  0.2   0.5   20  50  FALSE 0.2027326  0.5493061 -0.3       -0.3465736 0.2830197  -1.224557 1.959964 0.2310457 0.2207423   
# [10,] "kendall"  0.2   0.5   200 500 FALSE 0.2027326  0.5493061 -0.3       -0.3465736 0.08419154 -4.11649  1.959964 0.9844787 3.846864e-05
# [11,] "kendall"  0.2   0.5   20  50  TRUE  0.1163111  0.444228  -0.3013534 -0.3279168 0.2830197  -1.158636 1.959964 0.211471  0.2466045   
# [12,] "kendall"  0.2   0.5   200 500 TRUE  0.05281287 0.3782312 -0.3084065 -0.3254183 0.08419154 -3.865214 1.959964 0.9716262 0.0001109919
compute.power(0.1,0.5,50,50)$power
compute.power(0.1,0.5,50,50,power_only=TRUE)

corr_power <- function(n.sims = 100, 
                        n1=20, 
                        n2=20, 
                        alpha = 0.05, 
                        beta = 0.2,
                        sidedness = 2,
                        res_min = -1,
                        res_max = 1,
                        res_inc = 0.05,
                        n1_name = "Population A",
                        n2_name = "Population B",
                        method  = "pearson",
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
  if (simulation==TRUE) {
    # Evaluate pairwise comparisons across nsims for power estimate
    results <- outer(corrs, corrs, FUN = function(r, c) mapply(compute.power,rho1 = r, rho2 = c,n1 = n1, n2 = n2,
                                                               threshold=alpha, nsims = n.sims, lower.tri = lower,
                                                               power_only=TRUE))
  }
  if (simulation==FALSE){
    ## Actually its quite interesting to plot what the probabilities would look like without sampling; 
    results <- outer(corrs, corrs, FUN = function(r, c) mapply(p_zdiff,rho1 = r, rho2 = c,n1, n2,log=FALSE))
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
                                    "n1 = ",n1,", n2 = ",n2,", alpha: ",alpha, ", sims: ",n.sims),
                                    xlab = paste0("Correlation in ",n1_name),
                                    ylab = paste0("Correlation in ",n2_name), adj = 0),
                 color.palette =  colorRampPalette(c("#f7fcf0","#525252")));
        arrows(0.63, 0.6, 0.845, 0.6, length = 0.14, lwd = 3, col = "steelblue3")
  

  end_time <- Sys.time()
  
  # display running time
  cat(paste0("Power for difference in ",corr_type," correlations","\n",
             "n1 = ",n1,", n2 = ",n2,", alpha: ",alpha, ", sims: ",n.sims,"\n",
             "Processing time: ",round((end_time - start_time)/60,1)," minutes\n\n"))
  return(fig)
}

corr_power(n.sims  = 100,
           n1      =  134, 
           n2      =  134, 
           alpha   =   0.05, 
           beta    =   0.2,
           n1_name = "Mz twins",
           n2_name = "Dz twins",
           method  =  "pearson",
           lower   =   FALSE) 

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
