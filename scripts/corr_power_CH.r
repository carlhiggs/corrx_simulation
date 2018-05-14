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

# install.packages("Rcpp")
# install.packages("simstudy")
require(compiler)
require(Rcpp)
require("simstudy")
sourceCpp('test.cpp')


# # deploy to shinyapps.io
# library(rsconnect)
# rsconnect::deployApp('C:/Users/Carl/OneDrive/Research/2 - BCA/Research project/bca_rp2/scripts/corr_power_app')
# 

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


  # cocor_test <- cocor(~V1+V2|V1+V2,list(a,b))
  # zouci <- unlist(cocor_test@zou2007)
  # zou <- ((zouci[1]) < 0 && (0 < zouci[2])) 


# # quick comparison of tests - fz is de facto gold standard
# n1 <- 90
# n2 <- 90
# r1 <- -0.2
# r2 <- 0.5
# p1 <- 0
# p2 <- 1
# dist <- "normal"
# method <- "pearson"
# a <- genCorGen(n1, nvars = 2, params1 = p1, params2 = p2, dist = dist, 
#                corMatrix = matrix(c(1, r1, r1, 1), ncol = 2), wide = TRUE)[,2:3]
# b <- genCorGen(n2, nvars = 2, params1 = p1, params2 = p2, dist = dist, 
#                corMatrix = matrix(c(1, r2, r2, 1), ncol = 2), wide = TRUE)[,2:3]
# hist(unlist(a[,1]))
# corr <- c(cor(a)[1,2],cor(b)[1,2])
# cat("corr: ",corr,"diff: ", corr[1] - corr[2])
# 
# cat("test","\t","p"      ,"\n",
#     "fz_a","\t", fz_nosim(r1,r2,n1,n2,method = method,power=FALSE),"\n",
#     "fz"  ,"\t", fz(a,b,method = method),"\n",
#     "gtv" ,"\t", gtv(a,b,method = method),"\n",
#     "pt"  ,"\t", pt(a,b,method = method),"\n",
#     "slr" ,"\t", slr(a,b,method = method),"\n",
#     "zou" ,"\t", zou(a,b,method = method))



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
# corr_diff_test(rho = c(.2,.4), n = c(100,300), distr = "normal",test = "fz_nosim") 
# corr_diff_test(rho = c(.2,.4), n = c(100,300), distr = "normal",test = "fz") 
# corr_diff_test(rho = c(.2,.4), n = c(100,300), distr = "normal",test = "gtv") 
# corr_diff_test(rho = c(.2,.4), n = c(100,300), distr = "normal",test = "pt") 
# corr_diff_test(rho = c(.2,.4), n = c(100,300), distr = "normal",test = "slr") 
# corr_diff_test(rho = c(.2,.4), n = c(100,300), distr = "normal",test = "zou") 
# corr_diff_test(rho = c(.2,.4), n = c(100,300), distr = "normal") 

corr_diff_test_compiled <- cmpfun(corr_diff_test)
# cdt  <- function() for (i in 1:1000)  corr_diff_test(rho = c(.2,.4), n = c(100,300), distr = "normal",test = "fz") 
# cdt_c  <- function() for (i in 1:1000)  corr_diff_test_compiled(rho = c(.2,.4), n = c(100,300), distr = "normal",test = "fz") 
# system.time(cdt())
# system.time(cdt_c())
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
# cp  <- function() for (i in 1:10)  corr_power(power_only=TRUE)
# cp_c  <- function() for (i in 1:10)  corr_power_compiled(power_only=TRUE)
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

# corr_power_compiled(power_only=TRUE, nsim = 10, test=c("fz","fz_nosim","gtv"))
# corr_power_compiled(power_only=TRUE, nsim = 10, test=c("gtv","gtvr"))
# # fz  fz_nosim       gtv                                     
# # 0.1000000 0.3494663 0.1000000 
# ## NOTE - issue if only one simulation is processed - the rowMeans subcommand bugs out
# # corr_power_compiled(power_only=TRUE, nsim = 10, test=c("fz","fz_nosim"))
# # Error in rowMeans(replicate(nsims, corr_diff_test(rho = rho, n = n, distr = distr,  : 
# # 'x' must be an array of at least two dimensions

corr_power_plot_plot <- function(data, nsims = 100, res_min = -0.95, res_max = 0.95, res_inc = 0.05, 
                            n = c(30,90),distr = "normal",
                            param1a = c(0,0), param1b = c(0,0),param2a = c(1,1), param2b = c(1,1),
                            test = c("fz","gtv","pt","slr","zou"),
                            alpha = 0.05, beta = 0.2, sidedness=2,method="pearson",  
                            names = c("Population A","Population B"),
                            lower.tri = FALSE){
  results <- plot
  print(results)
  # format method to proper case for plot
  corr_type <- stringr::str_to_title(method)
  # define target power (to mark on plot)
  target <- 1 - beta
  corrs <- round(seq(res_min, res_max, res_inc),2)
  # plot power simulation results
  out<-filled.contour(x = corrs, y = corrs, z = as.matrix(results), nlevels = 10,
                                             xlim = c(-1,1), ylim = c(-1,1), zlim = c(0,1),
                                             plot.axes = {contour(x = corrs, y = corrs, z = as.matrix(results),
                                                                  levels = target, at = seq(-1, 1, 0.2), drawlabels = FALSE, axes = FALSE,
                                                                  add = TRUE, lwd = 3, col = "steelblue3");
                                               abline(v = seq(-1, 1, 0.1), lwd = .5, col = "lightgray", lty = 2)
                                               abline(h = seq(-1, 1, 0.1), lwd = .5, col = "lightgray", lty = 2)
                                               axis(1, seq(-1,1,0.2))
                                               axis(2, seq(-1,1,0.2))},
                                             plot.title = title(main = paste0(test," test ~ ",
                                                                              distr,"((",param1a,",",param2a,"),(",
                                                                              param1b,",",param2b,"))","\n",
                                                                              "Mz = ",n[1],
                                                                              ", Dz = ",n[2],
                                                                              ", alpha: ",alpha, ", sims: ",nsims),
                                                                xlab = paste0("Correlation in ",names[1]),
                                                                ylab = paste0("Correlation in ",names[2]), adj = 0),
                                             color.palette =  colorRampPalette(c("#f7fcf0","#525252")));
    # arrows(0.63, 0.6, 0.845, 0.6, length = 0.14, lwd = 3, col = "steelblue3")                          
  return(out)
}

cp_plotplot <- cmpfun(corr_power_plot_plot)

corr_power_plot <- function(nsims = 100, res_min = -0.95, res_max = 0.95, res_inc = 0.05, 
                            n = c(30,90),distr = "normal",
                            param1a = c(0,0), param1b = c(0,0),param2a = c(1,1), param2b = c(1,1),
                            tests = c("fz_nosim","fz","gtv","slr","zou"),
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
                                plot.title = title(main = paste0(test," test ~ ",
                                                                 distr,"((",param1a,",",param2a,"),(",
                                                                            param1b,",",param2b,"))","\n",
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

# run simulations
sim_res <- list()

# first run complete - Pearson, Normal(0,1), all tests except PT
sim_n <- c(15,30,60,120,240,480,960)
dist <- "normal"
method <-"pearson"
sim_tests <- c("fz_nosim","fz","gtv","slr","zou")
param1a <- param1b <- c(0,0)
param2a <- param2b <- c(1,1)
for (n1 in sim_n) {
  for (n2 in sim_n) {
      system.time(sim_res[[method]][[dist]][[paste0(param1a,"_",param1b,"_",param2a,"_",param2b)[1]]][[paste0(n1,"_",n2)]] <- 
                    corr_pplot_compiled(nsims = 100, res_min = -.95, res_max = 0.95, res_inc = 0.05, 
                                        n = c(n1,n2),
                                        param1a = param1a,
                                        param1b = param1b,
                                        param2a = param2a,
                                        param2b = param2b,
                                        dist    = dist,
                                        tests   = sim_tests,
                                        method  = method))
  }
}


# 
# # test using parallel processing 20180512
# # install.packages("foreach")
# # install.packages("doParallel")
# library(foreach)
# library(doParallel)
# 
# #setup parallel backend to use 8 processors
# cl<-makeCluster(8)
# registerDoParallel(cl)
# 
# # test a list to test this out
# sim_corr_res2 <- list()
# 
# 
# 
# 
# nsims <- 100
# res <- c(-0.95,0.95,0.05)
# sim_n <- c(15,30,60,120,240,480,960)
# sim_groups <- expand.grid(sim_n,sim_n) 
# # rccp compiled commands need to be treated in special way as per
# # https://stackoverflow.com/questions/25062383/cant-run-rcpp-function-in-foreach-null-value-passed-as-symbol-address
# # so using gtvr (i couldn't get the foreach .noexport option for rccp functions to work)
# sim_tests <- c("fz_nosim","fz","gtvr","slr","zou")
# 
# # second run - Spearman, Normal(0,1), all tests except PT
# method <-"spearman"
# dist <- "normal"
# param1a <- param1b <- c(0,0)
# param2a <- param2b <- c(1,1)
# paramstr <- paste0(param1a,"_",param1b,"_",param2a,"_",param2b)[1]
# # foreach(i =1:nrow(sim_groups)) %do%
# sim_corr_res2[[method]][[dist]][[paramstr]][[]] <- foreach(i =1:nrow(sim_groups))  %dopar% {
#     corr_pplot_compiled(nsims = nsims, res_min = res[1], res_max = res[2], res_inc = res[3], 
#                                       n = as.numeric(sim_groups[i,]),
#                                       param1a = param1a,
#                                       param1b = param1b,
#                                       param2a = param2a,
#                                       param2b = param2b,
#                                       dist    = dist,
#                                       tests   = sim_tests,
#                                       method  = method)
# }
# 
# 
# method <-"pearson"
# dist <- "gamma"
# param1a <- param1b <- c(1.5,1.5)
# param2a <- param2b <- c(0.09,0.09)
# paramstr <- paste0(param1a,"_",param1b,"_",param2a,"_",param2b)[1]
# # foreach(i =1:nrow(sim_groups)) %do%
# sim_corr_res2[[method]][[dist]][[paramstr]] <- foreach(i =1:nrow(sim_groups))  %dopar% {
#   # rccp compiled commands need to be treated in special way as per
#   # https://stackoverflow.com/questions/25062383/cant-run-rcpp-function-in-foreach-null-value-passed-as-symbol-address
#   
#   corr_pplot_compiled(nsims = nsims, res_min = res[1], res_max = res[2], res_inc = res[3], 
#                       n = as.numeric(sim_groups[i,]),
#                       param1a = param1a,
#                       param1b = param1b,
#                       param2a = param2a,
#                       param2b = param2b,
#                       dist    = dist,
#                       tests   = sim_tests,
#                       method  = method)
# }
# 
# method <-"spearman"
# dist <- "gamma"
# param1a <- param1b <- c(1.5,1.5)
# param2a <- param2b <- c(0.09,0.09)
# paramstr <- paste0(param1a,"_",param1b,"_",param2a,"_",param2b)[1]
# # foreach(i =1:nrow(sim_groups)) %do%
# sim_corr_res2[[method]][[dist]][[paramstr]] <- foreach(i =1:nrow(sim_groups))  %dopar% {
#   # rccp compiled commands need to be treated in special way as per
#   # https://stackoverflow.com/questions/25062383/cant-run-rcpp-function-in-foreach-null-value-passed-as-symbol-address
#   
#   corr_pplot_compiled(nsims = nsims, res_min = res[1], res_max = res[2], res_inc = res[3], 
#                       n = as.numeric(sim_groups[i,]),
#                       param1a = param1a,
#                       param1b = param1b,
#                       param2a = param2a,
#                       param2b = param2b,
#                       dist    = dist,
#                       tests   = sim_tests,
#                       method  = method)
# }
# 
# method <-"pearson"
# dist <- "gamma"         
# param1a <- param1b <- c(1,1)
# param2a <- param2b <- c(5,5)
# paramstr <- paste0(param1a,"_",param1b,"_",param2a,"_",param2b)[1]
# # foreach(i =1:nrow(sim_groups)) %do%
# sim_corr_res2[[method]][[dist]][[paramstr]] <- foreach(i =1:nrow(sim_groups))  %dopar% {
#   # rccp compiled commands need to be treated in special way as per
#   # https://stackoverflow.com/questions/25062383/cant-run-rcpp-function-in-foreach-null-value-passed-as-symbol-address
#   
#   corr_pplot_compiled(nsims = nsims, res_min = res[1], res_max = res[2], res_inc = res[3], 
#                       n = as.numeric(sim_groups[i,]),
#                       param1a = param1a,
#                       param1b = param1b,
#                       param2a = param2a,
#                       param2b = param2b,
#                       dist    = dist,
#                       tests   = sim_tests,
#                       method  = method)
# }
# method <-"spearman"
# dist <- "gamma"         
# param1a <- param1b <- c(1,1)
# param2a <- param2b <- c(5,5)
# paramstr <- paste0(param1a,"_",param1b,"_",param2a,"_",param2b)[1]
# # foreach(i =1:nrow(sim_groups)) %do%
# sim_corr_res2[[method]][[dist]][[paramstr]] <- foreach(i =1:nrow(sim_groups))  %dopar% {
#   # rccp compiled commands need to be treated in special way as per
#   # https://stackoverflow.com/questions/25062383/cant-run-rcpp-function-in-foreach-null-value-passed-as-symbol-address
#   
#   corr_pplot_compiled(nsims = nsims, res_min = res[1], res_max = res[2], res_inc = res[3], 
#                       n = as.numeric(sim_groups[i,]),
#                       param1a = param1a,
#                       param1b = param1b,
#                       param2a = param2a,
#                       param2b = param2b,
#                       dist    = dist,
#                       tests   = sim_tests,
#                       method  = method)
# }
# stopCluster(cl)

# 
# 
# sim_n <- c(15,30,60,120,240,480,960)
# sim_dist <- c("normal","gamma")
# sim_method <- c("pearson","spearman")
# sim_tests <- c("fz_nosim","fz","gtv","slr","zou")
# sim_res <- list()
# for (method in sim_method) {
#   for (dist in sim_dist) {
#     for (n1 in sim_n) {
#       for (n2 in sim_n) {
#         if (dist=="normal") {
#           param1a <- param1b <- c(0,0)
#           param2a <- param2b <- c(1,1)
#           system.time(sim_res[[method]][[dist]][[paste0(param1a,"_",param1b,"_",param2a,"_",param2b)[1]]][[paste0(n1,"_",n2)]] <- 
#             corr_pplot_compiled(nsims = 100, res_min = -.95, res_max = 0.95, res_inc = 0.05, 
#                                 n = c(n1,n2),
#                                 param1a = param1a,
#                                 param1b = param1b,
#                                 param2a = param2a,
#                                 param2b = param2b,
#                                 dist    = dist,
#                                 tests   = sim_tests,
#                                 method  = method))
#         }
#         if (dist=="gamma") {
#           param1a <- param1b <- c(1.5,1.5)
#           param2a <- param2b <- c(0.09,0.09)
#           system.time(sim_res[[method]][[dist]][[paste0(param1a,"_",param1b,"_",param2a,"_",param2b)]][[paste0(n1,"_",n2)]] <- 
#                         corr_pplot_compiled(nsims = 100, res_min = -.95, res_max = 0.95, res_inc = 0.05, 
#                                             n = c(n1,n2),
#                                             param1a = param1a,
#                                             param1b = param1b,
#                                             param2a = param2a,
#                                             param2b = param2b,
#                                             dist    = dist,
#                                             method  = method))
#           param1a <- param1b <- c(1,1)
#           param2a <- param2b <- c(5,5)
#           system.time(sim_res[[method]][[dist]][[paste0(param1a,"_",param1b,"_",param2a,"_",param2b)]][[paste0(n1,"_",n2)]] <- 
#                         corr_pplot_compiled(nsims = 100, res_min = -.95, res_max = 0.95, res_inc = 0.05, 
#                                             n = c(n1,n2),
#                                             param1a = param1a,
#                                             param1b = param1b,
#                                             param2a = param2a,
#                                             param2b = param2b,
#                                             dist    = dist,
#                                             method  = method))
#         }
#       }
#     }
#   }
# }



corrs <- round(seq(-0.90,0.90,0.1),2)
sim_n <- c(15,30,60,120,240,480,960)
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
# test3 <- dt[1:5,]
system.time(dt[, tests:={
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



# 
# corr_power_compiled(rho = test2[1, c(rho1,rho2)], 
#                     n = test2[1,c(n1,n2)],
#                     distr = test2[1,dist],
#                     param1a = test2[1,c(p1,p1)], 
#                     param1b = test2[1,c(p1,p1)],
#                     param2a = test2[1,c(p2,p2)], 
#                     param2b = test2[1,c(p2,p2)],
#                     test    = test2[1,test],
#                     alpha   = 0.05, 
#                     sidedness=2,
#                     method =method,  
#                     nsims = 100,
#                     power_only = TRUE)
# 
# corr_power_compiled(rho = c(-0.95,-0.95), 
#                     n = c(30,90),
#                     distr = "normal",
#                     param1a = c(0,0), 
#                     param1b = c(0,0),
#                     param2a = c(1,1), 
#                     param2b = c(1,1),
#                     test    = "fz_nosim",
#                     alpha   = 0.05, 
#                     sidedness=2,
#                     method ="pearson",  
#                     nsims = 100,
#                     power_only = TRUE)

# for (c = 1:nrow(sim_groups)) {
#   c.a <- as.character(format(sim_groups[c,1],nsmall=2))
#   c.b <- as.character(format(sim_groups[c,2],nsmall=2))
#   # ratio_vec <- vector(length = 0)
#   for(k = 1:nrow(sim_groups)) {
#     i <- sim_groups[k,1]
#     j <- sim_groups[k,2]
#     n_str <- paste0(i,'_',j)
#     # ratio_vec <- c(ratio_vec,paste0('(',i,',',j,')'))
#     test_rec <- sim_res$pearson$normal$`0_0_1_1`[[n_str]]$tests
#     for(test in test_rec){
#       # cat(paste(n_str,test,c.a,c.b,'\n'))
#       power_vec <- rbind(power_vec,
#                          cbind(test = test,
#                                method = method,
#                                dist = dist,
#                                n1 = i,
#                                n2 = j,
#                                n = i+j,
#                                ratio = format(i/j,nsmall=2)))
# nstr = paste0('(',i,',',j,')'),
# power = sim_res$pearson$normal$`0_0_1_1`[[n_str]][[test]][c.a,c.b]))
#   }
# }
# outtab <- matrix(power_vec,nrow=length(test_rec),ncol=length(sim_n))
# rownames(outtab) <- test_rec
# colnames(outtab) <- ratio_vec
# 
# library(ggplot2)
# p <- ggplot(power_vec)  
# p2 <- p + geom_point(aes(x = n2, y = power, colour = test), size = 3) 
# p3 <- p2 + geom_smooth(aes(x = n2, y = power))
# p3
# 
# system.time(results<- corr_pplot_compiled(nsims = 10, res_min = -.3, res_max = 0.3, res_inc = 0.1, n = c(30,90)))
# # Correlation power plot simulation commenced at 2018-05-01 21:57:57 
# # method	rho_1	rho_2	n1	n2	alpha	sides	nsims	distr	PowerXtests	
# # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# #   pearson	0.3	0.3	30	90	0.05	2	10	normal	0.025	0	0
# # Completed at  2018-05-01 21:59:00 
# # user  system elapsed 
# # 58.14    5.07   63.17
# 
# system.time(results<- corr_pplot_compiled())
# ## Note - this time test was run on more powerful work computer
# # Correlation power plot simulation commenced at 2018-05-01 22:05:06 
# # method	rho_1	rho_2	n1	n2	alpha	sides	nsims	distr	PowerXtests	
# # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# #   pearson	0.95	0.95	30	90	0.05	2	100	normal	0.025	0.09	0.09
# # Completed at  2018-05-02 02:06:36 
# # user   system  elapsed 
# # 14425.09    94.50 14489.76
# #
# # There was some error viewing all plots some time later - not sure why; however, results should be saved
# corrs<- seq(-0.95,0.95,0.05)
# nsims = 100 
# n = c(30,90)
# distr = "normal"
# tests = c("fz_nosim","fz","gtv")
# alpha = 0.05
# beta = 0.2
# sidedness=2
# target = 1-beta
# method="pearson"  
# names = c("Population A","Population B")
# for (test in results[["tests"]]) {
# results[["fig"]][[test]]<-filled.contour(x = corrs, y = corrs, z = as.matrix(results[[test]]), nlevels = 10,
#                                          xlim = c(-1,1), ylim = c(-1,1), zlim = c(0,1),
#                                          plot.axes = {contour(x = corrs, y = corrs, z = as.matrix(results[[test]]),
#                                                               levels = target, at = seq(-1, 1, 0.2), drawlabels = FALSE, axes = FALSE,
#                                                               add = TRUE, lwd = 3, col = "steelblue3");
#                                            abline(v = seq(-1, 1, 0.1), lwd = .5, col = "lightgray", lty = 2)
#                                            abline(h = seq(-1, 1, 0.1), lwd = .5, col = "lightgray", lty = 2)
#                                            axis(1, seq(-1,1,0.2))
#                                            axis(2, seq(-1,1,0.2))},
#                                          plot.title = title(main = paste0(test," test ~ ",
#                                                                  distr,"((",param1a,",",param2a,"),(",
#                                                                                                     param1b,",",param2b,"))","\n",
#                                                                           "Mz = ",results[["params"]][["n1"]],
#                                                                           ", Dz = ",results[["params"]][["n2"]],
#                                                                           ", alpha: ",alpha, ", sims: ",nsims),
#                                                             xlab = paste0("Correlation in ",names[1]),
#                                                             ylab = paste0("Correlation in ",names[2]), adj = 0),
#                                          color.palette =  colorRampPalette(c("#f7fcf0","#525252")));
# arrows(0.63, 0.6, 0.845, 0.6, length = 0.14, lwd = 3, col = "steelblue3")    
# }
# 
# 
# 
# # Non-normal distribution simulation
# # normal reference distribution
# bivariate_distribution_normal_m0_s1_n50 <- ggExtra::ggMarginal(data = as.data.frame(a), x = "V1", y = "V2") 
# plot(bivariate_distribution_normal_m0_s1_n50)
# # gamma distribution example (mean/shape, and rate/dispersion)
# # some positive skew (but distinctly non-normal)
# gamma <- genCorGen(1000, nvars = 2, params1 = c(1.5,1.5), params2 = c(0.09,0.09),dist = "gamma", 
#                    corMatrix = matrix(c(1, 0.5, 0.5, 1), ncol = 2), wide = TRUE)[,2:3]
# bivariate_distribution_gamma_m1.5_d0.09_n50 <- ggExtra::ggMarginal(data = as.data.frame(gamma), x = "V1", y = "V2") 
# plot(bivariate_distribution_gamma_m1.5_d0.09_n50)
# 
# # more gamma exploration, using our otherwise defaults
# gamma_a <- genCorGen(20, nvars = 2, params1 = c(1,1), params2 = c(0.09,0.09),dist = "gamma", 
#                      corMatrix = matrix(c(1, 0.5, 0.5, 1), ncol = 2), wide = TRUE)[,2:3]
# gamma_b <- genCorGen(90, nvars = 2, params1 = c(1,1), params2 = c(0.09,0.09),dist = "gamma", 
#                      corMatrix = matrix(c(1, 0.5, 0.5, 1), ncol = 2), wide = TRUE)[,2:3]
# 
# ggExtra::ggMarginal(data = as.data.frame(gamma_a), x = "V1", y = "V2") 
# ggExtra::ggMarginal(data = as.data.frame(gamma_b), x = "V1", y = "V2") 
# fz(gamma_a,gamma_b)
# gtv(gamma_a,gamma_b)   
# corr_power(dist = "gamma",param1a = c(1.5,1.5), param1b = c(1.5,1.5), param2a = c(.09,.09), param2b = c(.09,.09))
# # pearson	0.2	0.5	30	90	0.05	2	100	gamma	0.37	0.37$params        
# # method     rho_1     rho_2        n1        n2     alpha sidedness     nsims     distr 
# # "pearson"     "0.2"     "0.5"      "30"      "90"    "0.05"       "2"     "100"   "gamma" 
# # 
# # $additional
# # param1a1 param1a2 param1b1 param1b2 param2a1 param2a2 param2b1 param2b2 
# # 1.50     1.50     1.50     1.50     0.09     0.09     0.09     0.09 
# # 
# # $analytical
# # fz_nosim 
# # 0.3494663 
# # 
# # $power
# # fz  gtv 
# # 0.37 0.37 
# 
# # simulate using non-normal, Extreme positively skewed distribution
# system.time(res_gamma<- corr_pplot_compiled(dist = "gamma",param1a = c(1.5,1.5), param1b = c(1.5,1.5),param2a = c(.09,.09),param2b = c(.09,.09)))
# ## Note - run on work computer, so time is optimistic
# # Correlation power plot simulation commenced at 2018-05-02 22:25:35 
# # method	rho_1	rho_2	n1	n2	alpha	sides	nsims	distr	PowerXtests	
# # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# #   pearson	0.95	0.95	30	90	0.05	2	100	gamma	0.025	0.07	0.07
# # Completed at  2018-05-03 02:22:07 
# # user   system  elapsed 
# # 14169.60    64.91 14192.36 
# ## about 3 hours and 54 minutes
# # cp_plotplot(data = results[["fz"]],dist = "gamma",param1a = c(1.5,1.5), param1b = c(1.5,1.5),param2a = c(.09,.09),param2b = c(.09,.09))
# corrs<- seq(-0.95,0.95,0.05)
# nsims = 100 
# n = c(30,90)
# dist = "gamma"
# param1a = c(1.5,1.5)
# param1b = c(1.5,1.5)
# param2a = c(.09,.09)
# param2b = c(.09,.09)
# tests = c("fz_nosim","fz","gtv")
# alpha = 0.05
# beta = 0.2
# sidedness=2
# target = 1-beta
# method="pearson"  
# names = c("Population A","Population B")
# for (test in res_gamma[["tests"]]) {
#   res_gamma[["fig"]][[test]]<-filled.contour(x = corrs, y = corrs, z = as.matrix(res_gamma[[test]]), nlevels = 10,
#                                            xlim = c(-1,1), ylim = c(-1,1), zlim = c(0,1),
#                                            plot.axes = {contour(x = corrs, y = corrs, z = as.matrix(res_gamma[[test]]),
#                                                                 levels = target, at = seq(-1, 1, 0.2), drawlabels = FALSE, axes = FALSE,
#                                                                 add = TRUE, lwd = 3, col = "steelblue3");
#                                              abline(v = seq(-1, 1, 0.1), lwd = .5, col = "lightgray", lty = 2)
#                                              abline(h = seq(-1, 1, 0.1), lwd = .5, col = "lightgray", lty = 2)
#                                              axis(1, seq(-1,1,0.2))
#                                              axis(2, seq(-1,1,0.2))},
#                                            plot.title = title(main = paste0(test," test ~ ",
#                                                                             dist,"((",param1a,",",param2a,"),(",
#                                                                             param1b,",",param2b,"))","\n",
#                                                                             "Mz = ",n[1],
#                                                                             ", Dz = ",n[2],
#                                                                             ", alpha: ",alpha, ", sims: ",nsims),
#                                                               xlab = paste0("Correlation in ",names[1]),
#                                                               ylab = paste0("Correlation in ",names[2]), adj = 0),
#                                            color.palette =  colorRampPalette(c("#f7fcf0","#525252")));
#   arrows(0.63, 0.6, 0.845, 0.6, length = 0.14, lwd = 3, col = "steelblue3")    
# }
# 
# gamma_c <- genCorGen(90, nvars = 2, params1 = c(1,1), params2 = c(5,5),dist = "gamma", 
#                      corMatrix = matrix(c(1, 0.2, 0.2, 1), ncol = 2), wide = TRUE)[,2:3]
# gamma_d <- genCorGen(90, nvars = 2, params1 = c(1,1), params2 = c(5,5),dist = "gamma", 
#                      corMatrix = matrix(c(1, 0.5, 0.5, 1), ncol = 2), wide = TRUE)[,2:3]
# 
# ggExtra::ggMarginal(data = as.data.frame(gamma_c), x = "V1", y = "V2") 
# ggExtra::ggMarginal(data = as.data.frame(gamma_d), x = "V1", y = "V2") 
# 
# system.time(res_gamma<- corr_pplot_compiled(dist = "gamma",param1a = c(1,1), param1b = c(1,1),param2a = c(5,5),param2b = c(5,5)))
# 
# 
# 
# 
# #other mslrt exploration
# # param1 are the probabilities
# binary <- genCorGen(50, nvars = 2, params1 = c(.3, .5), dist = "binary", 
#                     corMatrix = matrix(c(1, 0.8, 0.8, 1), ncol = 2), wide = TRUE)


