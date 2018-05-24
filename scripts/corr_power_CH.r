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

# example of calculating power for all simulation combinations
test3 <- dt[1:100,]
system.time(test3[,  c( "fz_nosim","fz","gtv","slr","zou"):={
  cat("\r",.GRP,"\r")
  as.list(corr_power_compiled(rho = c(rho1,rho2), 
                              n = c(n1,n2),
                              distr = dist,
                              param1a = c(p1,p1), 
                              param1b = c(p1,p1),
                              param2a = c(p2,p2), 
                              param2b = c(p2,p2),
                              test    =  c( "fz_nosim","fz","gtv","slr","zou"),
                              alpha   = 0.05, 
                              sidedness=2,
                              method=as.character(method),  
                              nsims = 100,
                              power_only = TRUE))
},by = 1:nrow(test3)] )

#     method   dist p1 p2 n1 n2 rho1 rho2  fz_nosim   fz  gtv  slr  zou
# 1: pearson normal  0  1 15 15 -0.9 -0.9 0.0250000 0.12 0.12 0.14 0.12
# 2: pearson normal  0  1 15 15 -0.9 -0.8 0.1480538 0.16 0.17 0.19 0.16
# 3: pearson normal  0  1 15 15 -0.9 -0.7 0.3162464 0.31 0.33 0.35 0.31
# 4: pearson normal  0  1 15 15 -0.9 -0.6 0.4794100 0.54 0.57 0.63 0.54
# 5: pearson normal  0  1 15 15 -0.9 -0.5 0.6181794 0.65 0.66 0.72 0.63
# 6: pearson normal  0  1 15 15 -0.9 -0.4 0.7285717 0.74 0.75 0.77 0.72

# example extraction of plot
filled.contour(x = seq(-0.9,-0.6,0.1), 
               y = seq(-0.9,-0.6,0.1), 
               z = matrix(unlist(test3[(method=="pearson")&
                                         (dist=="normal")&
                                         (p1==0)&
                                         (p2==1)&
                                         (n1==15)&
                                         (n2==15)&
                                         (rho1<=-0.6)&
                                         (rho2<=-0.6),
                                       "gtv"]),
                          nrow=4,ncol=4, byrow = TRUE))

# run main simulation study
system.time(dt_normal[, tests:={
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


system.time(dt_gamma1[, tests:={
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


system.time(dt_gamma2[, tests:={
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


### scratch code for getting results processed from environment on work computer
# tmp.env <- new.env() # create a temporary environment
# load("r_power_work_partial_dt.RData", envir=tmp.env) # load workspace into temporary environment
# x <- get("dt", pos=tmp.env) # get the objects you need into your globalenv()
# #x <- tmp.env$x # equivalent to previous line
# rm(tmp.env) # remove the temporary environment to free up memory
# dt <- rbind(other[is.na(fz_nosim)==FALSE,],dt_test)


# convert to long
dt.long <- melt(dt, measure.vars = tests, variable.name = "test", value.name = "power")

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
          (dist=="normal"),list("fz_nosim" = round(mean(fz_nosim),2),
                               "fz" = round(mean(fz),2),
                               "zou" = round(mean(zou),2),
                               "gtv" = round(mean(gtv),2),
                               "slr" = round(mean(slr),2))
                               ,by=list(ratio)]

require(Publish)
dt.long[(method=="pearson")&(dist=="normal"),list(p50=round(quantile(power, .50, na.rm=TRUE),2),
                                                  p25=round(quantile(power, .025, na.rm=TRUE),2),
                                                  p75=round(quantile(power, .975, na.rm=TRUE),2)),by=test]
dt.long[(method=="pearson")&(dist=="gamma")&(p1=="1.5"),list(p50=round(quantile(power, .50, na.rm=TRUE),2),
                                                             p25=round(quantile(power, .25, na.rm=TRUE),2),
                                                             p75=round(quantile(power, .75, na.rm=TRUE),2)),by=test]
dt.long[(method=="pearson")&(dist=="gamma")&(p1=="1"),list(p50=round(quantile(power, .50, na.rm=TRUE),2),
                                                           p25=round(quantile(power, .025, na.rm=TRUE),2),
                                                           p75=round(quantile(power, .975, na.rm=TRUE),2)),by=test]


lparam1a <- 1
lparam1b <- 1
lparam2a <- 5
lparam2b <- 5
ln1 <- 60
ln2 <- 120
lalpha <- 0.05
lthreshold <- 0.8
lnsims <- 100
ltest <- "slr"
ldist <- "gamma"
lmethod <-"pearson"
lnames <- c("Mz twins","Dz twins")
corrxplot <- matrix(unlist(dt.long[(method==lmethod)&
                   (dist==ldist)&
                   (p1==lparam1a)&
                   (p2==lparam2a)&
                   (n1==ln1)&
                   (n2==ln2)&
                   (test==ltest),
                   power]),
                 nrow=19,ncol=19, byrow = TRUE)
test_lookup <- cbind(data.table("test" = c("fz_nosim","fz","gtv","slr","zou")),
                     data.table("name" = c("Fisher's z (no sim)","Fisher's Z","GTV","SLR","Zou's CI")))

filled.contour(x = corrs,y = corrs,z = corrxplot, nlevels = 10,
               xlim = c(-1,1), ylim = c(-1,1), zlim = c(0,1),
               plot.axes = {contour(x = corrs, y = corrs, z = corrxplot,
                                    levels = lthreshold, at = seq(-1, 1, 0.2), drawlabels = FALSE, axes = FALSE,
                                    add = TRUE, lwd = 3, col = "steelblue3");
                 abline(v = seq(-1, 1, 0.1), lwd = .5, col = "lightgray", lty = 2)
                 abline(h = seq(-1, 1, 0.1), lwd = .5, col = "lightgray", lty = 2)
                 axis(1, seq(-1,1,0.2))
                 axis(2, seq(-1,1,0.2))},
               plot.title = title(main = paste0(test_lookup[test==ltest,name]," test ~ ",
                                       ldist,"((",lparam1a,",",lparam2a,"),(",
                                                                          lparam1b,",",lparam2b,"))","\n",
                                                "Mz = ",ln1,
                                                ", Dz = ",ln2,
                                                ", alpha: ",lalpha, ", sims: ",lnsims),
                                  xlab = paste0("Correlation in ",lnames[1]),
                                  ylab = paste0("Correlation in ",lnames[2]), adj = 0),
               color.palette =  colorRampPalette(c("#f7fcf0","#525252")));
               arrows(0.63, 0.6, 0.845, 0.6, length = 0.14, lwd = 3, col = "steelblue3")    
               
# alternative plots 
library(ggplot2)
p <- ggplot(dt.long[(method=="pearson")&
                      (dist=="normal")&
                      (p1==0)&
                      (p2==1)&
                      (n1==15)&
                      (n2==15),])
p2 <- p + geom_point(aes(x = diff, 
                         y = power,
                         colour = test), size = 3)
p2
# could do: colour = interaction(ratio,n,sep="-",lex.order=TRUE)
p3 <- p + geom_smooth(aes(x = diff, y = power, colour = test))
p3

# plot of test by method
p <- ggplot(dt.long[(dist=="normal")&
                      (p1==0)&
                      (p2==1)&
                      (n1==15)&
                      (n2==15)&
                      (test=="fz"),])
# could do: colour = interaction(ratio,n,sep="-",lex.order=TRUE)
p3 <- p + geom_smooth(aes(x = diff, y = power, colour = method))
p3

# Plot power by N, holding ratio and correlations constant
p <- ggplot(dt.long[(method=="pearson")&
                      (dist=="normal")&
                      (p1==0)&
                      (p2==1)&
                      (ratio==1)&
                      (rho1==-0.9)&
                      (rho2==-0.6),])
p2 <- p + geom_point(aes(x = n, 
                         y = power,
                         colour = test), size = 3)
p2

# Plot power by ratio, holding n and correlations constant
p <- ggplot(dt.long[(method=="pearson")&
                      (dist=="normal")&
                      (p1==0)&
                      (p2==1)&
                      (rho1==-0.9)&
                      (rho2==-0.6),])
p <- p + geom_point(aes(x = ratio, 
                         y = power,
                         colour = test), size = 3)

p <- p + geom_tile(aes(fill=n))




dt.long[(method=="pearson")&
          (dist=="normal")&
          (p1==0)&
          (p2==1)&
          (n1==15)&
          (n2==15),mean(power),by=test]

dt.long[(method=="spearman")&
          (dist=="normal")&
          (p1==0)&
          (p2==1)&
          (n1==15)&
          (n2==15),mean(power),by=test]

dt[, lapply(.SD, sum, na.rm=TRUE), by=category ]

dt.long[(dist=="normal")&
          (p1==0)&
          (p2==1)&
          (ratio==1),
        mean(power),by=list(method,test,n)]


# plot power curve given parameters - Normal
um <- dt.long[(dist=="normal")&
                (method=="pearson")&
                (p1==0)&
                (p2==1)&
                (rho1==0.2)&
                (rho2==0.5)&
                (ratio==1)&
                (test=='zou'),
              power,by=list(n)]

threshold <- 0.8
# interpolate power curve
fit <- cbind(data.table("n" = min(um$n):max(um$n)),
             data.table("power" = splinefun(um$n, 
                                            um$power,
                                            method = "monoH.FC")(min(um$n):max(um$n))))

# find point of y which satisfied f(y) = desired power (e.g. 0.8)
cross <-fit[which.min(abs(threshold-fit$power)),] #find closest to Y intercept (0.8)
# plots results
plot(um,lwd=1)
lines(fit, lwd = 2)
abline(h=threshold)
abline(v=cross[[1]])

# plot power curve by N given parameters  - Normal
um <- dt.long[(dist=="normal")&
                (method=="pearson")&
                (p1==0)&
                (p2==1)&
                (rho1==0.2)&
                (rho2==0.5)&
                (ratio==0.5),
              power,by=list(n,test)]



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

cross <- data.table()
for(x in levels(um$test)) {
  cross <- rbind(cross, 
                 fit[(test %in% x)][which.min(abs(threshold-fit[(test %in% x),power]))])
}

cross$label <- c("FZ (no sim.)",
                 "FZ",
                 "GTV",
                 "SLR",
                 "Zou's CI")
cross[,"label":= paste0(label," (",n,")"),by=1:nrow(cross)]
cross <- cross[order(+rank(n))]

method <- "pearson"
dist <- "normal"
p <- c(0,1)
rho1 <- 0.2
rho2 <- 0.5
ratio <- 0.5
title <- paste0("Power to detect difference in ",method," correlations, by sample size\n",
                dist,"((",p[1],",",p[1],"),(",p[2],",",p[2],"))","\n",
                "rho: (",rho1,",",rho2,")",
                "; Mz to Dz ratio: ",ratio, "; sims: ",100)
ggplot(NULL, aes(x = n, y = power, colour = test, group = test))+ 
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
  scale_colour_discrete(name="Tests",
                        breaks=cross$test,
                        labels=cross$label)  +
  theme(panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  ggtitle(title)


# plot power curve by N given parameters - Gamma, subtle
um <- dt.long[(dist=="gamma")&
                (method=="pearson")&
                (p1==1.5)&
                (p2==0.09)&
                (rho1==0.2)&
                (rho2==0.5)&
                (ratio==1),
              power,by=list(n,test)]



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

cross <- data.table()
for(x in levels(um$test)) {
  cross <- rbind(cross, 
                 fit[(test %in% x)][which.min(abs(threshold-fit[(test %in% x),power]))])
}

cross$label <- c("FZ (no sim.)",
                 "FZ",
                 "GTV",
                 "SLR",
                 "Zou's CI")
cross[,"label":= paste0(label," (",n,")"),by=1:nrow(cross)]
cross <- cross[order(+rank(n))]

method <- "pearson"
dist <- "gamma"
p <- c(1.5,0.09)
rho1 <- 0.2
rho2 <- 0.5
ratio <- 1
title <- paste0("Power to detect difference in ",method," correlations, by sample size\n",
                dist,"((",p[1],",",p[1],"),(",p[2],",",p[2],"))","\n",
                "rho: (",rho1,",",rho2,")",
                "; Mz to Dz ratio: ",ratio, "; sims: ",100)
ggplot(NULL, aes(x = n, y = power, colour = test, group = test))+ 
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
  scale_colour_discrete(name="Tests",
                        breaks=cross$test,
                        labels=cross$label)  +
  theme(panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  ggitle(title)


# plot power curve by N given parameters - Gamma, extreme
um <- dt.long[(dist=="gamma")&
                (method=="pearson")&
                (p1==1)&
                (p2==5)&
                (rho1==0.2)&
                (rho2==0.5)&
                (ratio==0.5),
              power,by=list(n,test)]



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

cross <- data.table()
for(x in levels(um$test)) {
  cross <- rbind(cross, 
                 fit[(test %in% x)][which.min(abs(threshold-fit[(test %in% x),power]))])
}

cross$label <- c("FZ (no sim.)",
                    "FZ",
                    "GTV",
                    "SLR",
                    "Zou's CI")
cross[,"label":= paste0(label," (",n,")"),by=1:nrow(cross)]
cross <- cross[order(+rank(n))]


method <- "pearson"
dist <- "gamma"
p <- c(1,5)
rho1 <- 0.2
rho2 <- 0.5
ratio <- 0.5
title <- paste0("Power to detect difference in ",method," correlations, by sample size\n",
                dist,"((",p[1],",",p[1],"),(",p[2],",",p[2],"))","\n",
                "rho: (",rho1,",",rho2,")",
                "; Mz to Dz ratio: ",ratio, "; sims: ",100)

ggplot(NULL, aes(x = n, y = power, colour = test, group = test))+ 
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
  scale_colour_discrete(name="Tests",
                        breaks=cross$test,
                        labels=cross$label)  +
  theme(panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  ggtitle(title)



# # plot power curve by difference given parameters
# dist      <- "normal"
# method    <- "pearson"
# ratio     <- 1
# n         <- 240
# threshold <- 0.8
# 
# um <- dt.long[(dist==dist)&
#                 (method==method)&
#                 (p1==0)&
#                 (p2==1)&
#                 (ratio==ratio)&
#                 (n==n),
#               power,by=list(test,rho1,rho2,diff)]

um <- dt.long[(dist=="normal")&
                (method=="pearson")&
                (p1==0)&
                (p2==1)&
                (ratio==1)&
                (n==240),
              power,by=list(test,rho1,rho2,diff)]

um[,"z1":=atanh(rho1), by = 1:nrow(um)]
um[,"z2":=atanh(rho2), by = 1:nrow(um)]
um[,"diff1":=tanh(z1-z2), by = 1:nrow(um)]
um[,"diff2":=abs(tanh(z1-z2)), by = 1:nrow(um)]
# 
# # Plots to clarify why the use of tanh(atanh(difference)) is important
# # Power estimates are not monotonic here using just corr differences
# ggplot(um[(test %in% "fz")],aes(x = diff,y = power,group = rho1, colour = rho1))+geom_point()+theme_bw()
# 
# # Here we see that the differences are ordered when considered by tanh(atanh(diff))
# ggplot(um[(test %in% "fz")],aes(x = diff,y = diff1,group = rho1, colour = rho1))+geom_point()+theme_bw()
# 
# # Here we see that the differences are ordered when considered by and simplified by |tanh(atanh(diff))|
# ggplot(um[(test %in% "fz")],aes(x = diff,y = diff2,group = rho1, colour = rho1))+geom_point()+theme_bw()
# 
# # Here we see that power increases as diff increases, but not equally across series of (rho1,rho2) combinations
# ggplot(um[(test %in% "fz")],aes(x = diff,y = diff2,group = power, colour = power))+geom_point()+theme_bw()
# 
# 
# # Rho by Rho showing power; with squiggly line analogue of our filled conotur
# ggplot(um,aes(x = rho1,
#               y = rho2, 
#               group = power, 
#               colour = power))+
#   geom_point()+
#   geom_line()+
#   theme_bw() 
# 
# # Here, we see how the assymetrical tanh(atanh(diff)) works; fixed rho, fixed test
# ggplot(um[(test=="fz")&(rho1==0.5),],aes(x = diff1,
#               y = power, 
#               group = test, 
#               colour = test))+
#   geom_point()+
#   theme_bw() 
# 
# # power by diff1;  fixed rho, series by test
# ggplot(um[(rho1==0.5),],aes(x = diff1,
#                                          y = power, 
#                                          group = test, 
#                                          colour = test))+
#   geom_point()+
#   theme_bw() 
# 
# # power by diff2;  fixed rho, series by test
# ggplot(um[(rho1==0.5),],aes(x = diff2,
#                             y = power, 
#                             group = test, 
#                             colour = test))+
#   geom_point()+
#   theme_bw() 
# 
# # amazing abstract picture  -- how??!
# ggplot(um[(test %in% "fz")],aes(x = rho1,y = power,group = rho2, colour = rho2))+
#   geom_point()+
#   geom_line()+
#   theme_bw()
# 
# ggplot(um[(test %in% "fz")],aes(x = z1,y = power,group = z2, colour = z2))+
#   geom_point()+
#   geom_line()+
  theme_bw()

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
dist <- "normal"
p <- c(0,1)
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

ggplot(NULL, aes(x = diff1, y = power, colour = test, group = test))+ 
  scale_x_continuous(bquote("d* = |"~tanh(atanh(rho[1])-atanh(rho[2]))~"|"), 
                     breaks = seq(0,1,0.1),
                     limits = c(0,1)) +
  scale_y_continuous(bquote(Power~(1-beta)), 
                     breaks = seq(0,1,0.1),
                     limits = c(0,1)) +
  geom_point(data = um)  
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


# 
# ggplotly(p)
# 
# colors <- c('#4AC6B7', '#1972A4', '#965F8A', '#FF7070', '#C61951')
# colors2 <- c(rev(colors),colors[2:5])
# 
# p <- plot_ly(dt.long[(method=="pearson")&
#                   (dist=="normal")&
#                   (p1==0)&
#                   (p2==1)&
#                   (test=="zou"),], 
#              x = ~n1, y = ~n2, z = ~diff, color = ~power, size = ~n, colors = rev(colors),
#              marker = list(symbol = 'circle', sizemode = 'diameter'), sizes = c(1, 10),
#              text = ~paste('Power:', power, '<br>Dz correlation:', rho2, '<br>Mz correlation:', rho1,
#                            '<br>Twin pairs:', n)) %>%
#   layout(title = 'Power by correlation and group size',
#          scene = list(xaxis = list(title = 'Mz twins',
#                                    gridcolor = 'rgb(255, 255, 255)',
#                                    # range = c(2.003297660701705, 5.191505530708712),
#                                    zerolinewidth = 1,
#                                    ticklen = 5,
#                                    gridwidth = 2),
#                       yaxis = list(title = 'Dz twins',
#                                    gridcolor = 'rgb(255, 255, 255)',
#                                    # type = 'log',
#                                    # range = c(36.12621671352166, 91.72921793264332),
#                                    zerolinewidth = 1,
#                                    ticklen = 5,
#                                    gridwith = 2),
#                       zaxis = list(title = 'Difference in r',
#                                    gridcolor = 'rgb(255, 255, 255)',
#                                    zerolinewidth = 1,
#                                    ticklen = 5,
#                                    gridwith = 2)),
#          paper_bgcolor = 'rgb(243, 243, 243)',
#          plot_bgcolor = 'rgb(243, 243, 243)')
# ggplotly(p)
