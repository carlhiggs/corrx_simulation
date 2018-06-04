# R script for plots, ***following script corr_power_CH.r***
#  to simulate power for difference in correlations (Pearson, Spearman, and later ... ICC)
#  Pl
# Carl Higgs 2017

# may have to set working directory
setwd("C:/Users/Carl/OneDrive/Research/2 - BCA/Research project/bca_rp2/scripts")
require(compiler)
require(Rcpp)
require("simstudy")
sourceCpp('test.cpp')
require(data.table)
require(ggplot2)
library(parallel)

## The below function and entries all assume a data.table object containing simulation 
#  results has been produced in the preceding analysis stage and converted to long form
# This has a schema like the below (see other script for full details):
# > head(dt_s3.long)
#    method   dist p1 p2 n1  n2 rho1 rho2   ratio   n diff     test     power
# 1 pearson normal  0  1 15  15  0.2  0.5 1.00000  30  0.3 fz_nosim 0.1332765
# 2 pearson normal  0  1 15  30  0.2  0.5 0.50000  45  0.3 fz_nosim 0.1682679
# 3 pearson normal  0  1 15  60  0.2  0.5 0.25000  75  0.3 fz_nosim 0.1924844
# 4 pearson normal  0  1 15 120  0.2  0.5 0.12500 135  0.3 fz_nosim 0.2070783
# 5 pearson normal  0  1 15 240  0.2  0.5 0.06250 255  0.3 fz_nosim 0.2151483
# 6 pearson normal  0  1 15 480  0.2  0.5 0.03125 495  0.3 fz_nosim 0.2194003
# ...

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
    filled.contour(x = corrs,y = corrs,z = corrx_tmp, nlevels = 10,
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
    ## arrows for square plot
    arrows(0.6, 0.54, 0.77, 0.54, length = 0.14, lwd = 3, col = "steelblue3")
    ## arrows for wide rectangle plot
    # arrows(0.65, 0.54, 0.795, 0.54, length = 0.14, lwd = 3, col = "steelblue3")
    
    # finalise plot
    dev.off()
  }
  
  # produce contour comparison plot
  if(type=="contour2"){
    # make correlation-power matrix
    corrx_tmp <- data.long[(method==lmethod)&
                                           (dist==ldist)&
                                           (p1==lparam1a)&
                                           (p2==lparam2a)&
                                           (n1==ln1)&
                                           (n2==ln2),]
    # lookup table for test names
    test_lookup <- cbind(dt_s1.long[,mean(power),by=test],
                         data.table("label" = c("Fisher's z (no sim)","Fisher's Z","GTV","SLR","Zou's CI")))
    # order tests by smallest sample size estimate required to achieve power threshold
    test_lookup <- test_lookup[order(-rank(V1))]
    test_lookup[,"label":= paste0(label," (",round(V1,2),")"),by=1:nrow(cross)]
    
    
    # define plot title
    title <- paste0("Power to detect difference in ",method," correlations, by rho\n",
                    ldist,"((",lparam1a,",",lparam2a,"),(",lparam1b,",",lparam2b,"))","\n",
                    "Mz = ",ln1,", Dz = ",ln2,", alpha: ",alpha, ", sims: ",lnsims)
    # initialise output plot
    pdf(graph_out,width=gwidth,height=gheight)
    # make contour plot
    p <- ggplot(corrx_tmp,aes(x = rho1,y = rho2, z = power, colour = test, group = test))  +
      geom_contour(breaks = 0.8) + 
      scale_x_continuous(paste0("Correlation in ",lnames[2]),
                         breaks =seq(-.9, 0.9, 0.2),
                         limits = c(-.9,.9))  +
      scale_y_continuous(paste0("Correlation in ",lnames[2]),
                         breaks =seq(-.9, 0.9, 0.2),
                         limits = c(-.9,.9)) +
      scale_colour_discrete(name="Tests",
                            breaks=test_lookup$test,
                            labels=test_lookup$label)  +
      theme(panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black")) +
      ggtitle(title)  
    print(p)
    # finalise plot to pdf
    dev.off()
    # display plot on screen
    print(p)
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
    # finalise plot to pdf
    dev.off()
    # display plot on screen
    print(p)
  }
  if(type %in% c("diffpower","diffpowerabs")){
    ln = n1+n2
    # plot power curve by difference given parameters
    um <- data.long[(dist==ldist)&
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
    # finalise plot to pdf
    dev.off()
    # display plot on screen
    print(p)
  }
}

# Function to produce example distribution plots
dist_x_plot <- function(method = "pearson",dist = "normal",test="",n1 = 60,n2 = 120,
                        param1a ,param1b , param2a ,param2b ,
                        rho1 = 0.2, rho2 = 0.5, ratio = 0.5,
                        nsims = 100, names = c("MZ","DZ"),
                        graph_out = "distxplot.pdf",gwidth = 7,gheight=6.5) {
  a <- genCorGen(n1, nvars = 2, params1 = param1a, params2 = param2a, dist = dist, 
                 corMatrix = matrix(c(1, rho1, rho1, 1), ncol = 2), wide = TRUE)[,2:3]
  b <- genCorGen(n2, nvars = 2, params1 = param1b, params2 = param2b, dist = dist, 
                 corMatrix = matrix(c(1, rho2, rho2, 1), ncol = 2), wide = TRUE)[,2:3]
  r.a <- round(cor(a, method = method)[1,2],2)
  r.b <- round(cor(b, method = method)[1,2],2)
  a$Group <- names[1]
  b$Group <- names[2]
  c <- rbind(a,b)
  # initialise plot
  pdf(graph_out,width=gwidth,height=gheight)
  
  p <- ggplot(as.data.frame(c), aes_string('V1', 'V2')) +
    aes_string(colour = 'Group') +
    geom_point() + theme_bw(15) +
    theme(legend.position = c(1.04, 1.1),legend.text.align	 = 0)  +
    scale_color_manual(labels = c(bquote(paste("  60 ",.(names[1]),"  (r = ",.(r.a),")")), 
                                  bquote(paste("120 ",.(names[2]),"  (r = ",.(r.b),")"))),
                       values = c("MZ" = "#ef8a62","DZ" = "#67a9cf"))
  
  p <- ggExtra::ggMarginal(p,
                           type = 'density',
                           margins = 'both',
                           size = 5,
                           groupColour = TRUE,
                           groupFill = TRUE,
                           alpha = 0.4
  ) 
  
  
  print(p)
  # finalise plot to pdf
  dev.off()
  # display plot on screen
  print(p)
}

# SLR - Scenario 1pre - 100 sims 60:120
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

# SLR - Scenario 1  - norm - 1000 sims 60:120
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


# SLR - Scenario 2  - norm -  1000 sims 90:90
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

# SLR - Scenario 1  - gamma1 - 1000 sims 60:120
corrxplot(data.long  =  dt_s1.long,
          method     = "pearson",
          dist       =  "gamma",
          test       =  "slr",
          n1         =  60,
          n2         =  120,
          param1a    =  1.5,
          param1b    =  1.5,
          param2a    =  0.09,
          param2b    =  0.09,
          nsims      =  1000,
          names      =  c("Mz twins","Dz twins"),
          type       = "contour",
          graph_out  =  "../figs/corrx_contour_gamma1_sA1_n180_mzdz.5_slr_s1000.pdf",
          alpha      =  0.05,
          threshold  =  0.8) 


# SLR -  Scenario 2  - gamma1 -  1000 sims 90:90
corrxplot(data.long  =  dt_s2.long,
          method     = "pearson",
          dist       =  "gamma",
          test       =  "slr",
          n1         =  90,
          n2         =  90,
          param1a    =  1.5,
          param1b    =  1.5,
          param2a    =  .09,
          param2b    =  .09,
          nsims      =  1000,
          names      =  c("Mz twins","Dz twins"),
          type       = "contour",
          graph_out  =  "../figs/corrx_contour_gamma1_sA2_n180_mzdz1_slr_s1000.pdf",
          alpha      =  0.05,
          threshold  =  0.8) 

# SLR - Scenario 1  - gamma2 - 1000 sims 60:120
corrxplot(data.long  =  dt_s1.long,
          method     = "pearson",
          dist       =  "gamma",
          test       =  "slr",
          n1         =  60,
          n2         =  120,
          param1a    =  1,
          param1b    =  1,
          param2a    =  5,
          param2b    =  5,
          nsims      =  1000,
          names      =  c("Mz twins","Dz twins"),
          type       = "contour",
          graph_out  =  "../figs/corrx_contour_gamma2_sA1_n180_mzdz.5_slr_s1000.pdf",
          alpha      =  0.05,
          threshold  =  0.8) 


# SLR - Scenario 2  - gamma2 -  1000 sims 90:90 
corrxplot(data.long  =  dt_s2.long,
          method     = "pearson",
          dist       =  "gamma",
          test       =  "slr",
          n1         =  90,
          n2         =  90,
          param1a    =  1,
          param1b    =  1,
          param2a    =  5,
          param2b    =  5,
          nsims      =  1000,
          names      =  c("Mz twins","Dz twins"),
          type       = "contour",
          graph_out  =  "../figs/corrx_contour_gamma2_sA2_n180_mzdz1_slr_s1000.pdf",
          alpha      =  0.05,
          threshold  =  0.8) 

# gtv - Scenario 1pre - 100 sims 60:120
corrxplot(data.long  =  dt.long,
          method     = "pearson",
          dist       =  "normal",
          test       =  "gtv",
          n1         =  60,
          n2         =  120,
          param1a    =  0,
          param1b    =  0,
          param2a    =  1,
          param2b    =  1,
          nsims      =  100,
          names      =  c("Mz twins","Dz twins"),
          type       = "contour",
          graph_out  =  "../figs/corrx_contour_sA1pre_n180_mzdz.5_gtv_s100.pdf",
          alpha      =  0.05,
          threshold  =  0.8) 

# gtv - Scenario 1  - norm - 1000 sims 60:120
corrxplot(data.long  =  dt_s1.long,
          method     = "pearson",
          dist       =  "normal",
          test       =  "gtv",
          n1         =  60,
          n2         =  120,
          param1a    =  0,
          param1b    =  0,
          param2a    =  1,
          param2b    =  1,
          nsims      =  1000,
          names      =  c("Mz twins","Dz twins"),
          type       = "contour",
          graph_out  =  "../figs/corrx_contour_sA1_n180_mzdz.5_gtv_s1000.pdf",
          alpha      =  0.05,
          threshold  =  0.8) 


# gtv - Scenario 2  - norm -  1000 sims 90:90
corrxplot(data.long  =  dt_s2.long,
          method     = "pearson",
          dist       =  "normal",
          test       =  "gtv",
          n1         =  90,
          n2         =  90,
          param1a    =  0,
          param1b    =  0,
          param2a    =  1,
          param2b    =  1,
          nsims      =  1000,
          names      =  c("Mz twins","Dz twins"),
          type       = "contour",
          graph_out  =  "../figs/corrx_contour_sA2_n180_mzdz1_gtv_s1000.pdf",
          alpha      =  0.05,
          threshold  =  0.8) 

# gtv - Scenario 1  - gamma1 - 1000 sims 60:120
corrxplot(data.long  =  dt_s1.long,
          method     = "pearson",
          dist       =  "gamma",
          test       =  "gtv",
          n1         =  60,
          n2         =  120,
          param1a    =  1.5,
          param1b    =  1.5,
          param2a    =  0.09,
          param2b    =  0.09,
          nsims      =  1000,
          names      =  c("Mz twins","Dz twins"),
          type       = "contour",
          graph_out  =  "../figs/corrx_contour_gamma1_sA1_n180_mzdz.5_gtv_s1000.pdf",
          alpha      =  0.05,
          threshold  =  0.8) 


# gtv -  Scenario 2  - gamma1 -  1000 sims 90:90
corrxplot(data.long  =  dt_s2.long,
          method     = "pearson",
          dist       =  "gamma",
          test       =  "gtv",
          n1         =  90,
          n2         =  90,
          param1a    =  1.5,
          param1b    =  1.5,
          param2a    =  .09,
          param2b    =  .09,
          nsims      =  1000,
          names      =  c("Mz twins","Dz twins"),
          type       = "contour",
          graph_out  =  "../figs/corrx_contour_gamma1_sA2_n180_mzdz1_gtv_s1000.pdf",
          alpha      =  0.05,
          threshold  =  0.8) 

# gtv - Scenario 1  - gamma2 - 1000 sims 60:120
corrxplot(data.long  =  dt_s1.long,
          method     = "pearson",
          dist       =  "gamma",
          test       =  "gtv",
          n1         =  60,
          n2         =  120,
          param1a    =  1,
          param1b    =  1,
          param2a    =  5,
          param2b    =  5,
          nsims      =  1000,
          names      =  c("Mz twins","Dz twins"),
          type       = "contour",
          graph_out  =  "../figs/corrx_contour_gamma2_sA1_n180_mzdz.5_gtv_s1000.pdf",
          alpha      =  0.05,
          threshold  =  0.8) 


# gtv - Scenario 2  - gamma2 -  1000 sims 90:90 
corrxplot(data.long  =  dt_s2.long,
          method     = "pearson",
          dist       =  "gamma",
          test       =  "gtv",
          n1         =  90,
          n2         =  90,
          param1a    =  1,
          param1b    =  1,
          param2a    =  5,
          param2b    =  5,
          nsims      =  1000,
          names      =  c("Mz twins","Dz twins"),
          type       = "contour",
          graph_out  =  "../figs/corrx_contour_gamma2_sA2_n180_mzdz1_gtv_s1000.pdf",
          alpha      =  0.05,
          threshold  =  0.8) 



# Compare: Scenario 1pre - 100 sims 60:120
corrxplot(data.long  =  dt.long,
          method     = "pearson",
          dist       =  "normal",
          n1         =  60,
          n2         =  120,
          param1a    =  0,
          param1b    =  0,
          param2a    =  1,
          param2b    =  1,
          nsims      =  100,
          names      =  c("Mz twins","Dz twins"),
          type       = "contour2",
          graph_out  =  "../figs/corrx_contour_sA1pre_n180_mzdz.5_compare_s100.pdf",
          alpha      =  0.05,
          threshold  =  0.8) 

# Compare: Scenario 1  - norm - 1000 sims 60:120
corrxplot(data.long  =  dt_s1.long,
          method     = "pearson",
          dist       =  "normal",
          n1         =  60,
          n2         =  120,
          param1a    =  0,
          param1b    =  0,
          param2a    =  1,
          param2b    =  1,
          nsims      =  1000,
          names      =  c("Mz twins","Dz twins"),
          type       = "contour2",
          graph_out  =  "../figs/corrx_contour_sA1_n180_mzdz.5_compare_s1000.pdf",
          alpha      =  0.05,
          threshold  =  0.8) 


# Compare: Scenario 2  - norm -  1000 sims 90:90
corrxplot(data.long  =  dt_s2.long,
          method     = "pearson",
          dist       =  "normal",
          n1         =  90,
          n2         =  90,
          param1a    =  0,
          param1b    =  0,
          param2a    =  1,
          param2b    =  1,
          nsims      =  1000,
          names      =  c("Mz twins","Dz twins"),
          type       = "contour2",
          graph_out  =  "../figs/corrx_contour_sA2_n180_mzdz1_compare_s1000.pdf",
          alpha      =  0.05,
          threshold  =  0.8) 

# Compare: Scenario 1  - gamma1 - 1000 sims 60:120
corrxplot(data.long  =  dt_s1.long,
          method     = "pearson",
          dist       =  "gamma",
          n1         =  60,
          n2         =  120,
          param1a    =  1.5,
          param1b    =  1.5,
          param2a    =  0.09,
          param2b    =  0.09,
          nsims      =  1000,
          names      =  c("Mz twins","Dz twins"),
          type       = "contour2",
          graph_out  =  "../figs/corrx_contour_gamma1_sA1_n180_mzdz.5_compare_s1000.pdf",
          alpha      =  0.05,
          threshold  =  0.8) 


# Compare: Scenario 2  - gamma1 -  1000 sims 90:90
corrxplot(data.long  =  dt_s2.long,
          method     = "pearson",
          dist       =  "gamma",
          n1         =  90,
          n2         =  90,
          param1a    =  1.5,
          param1b    =  1.5,
          param2a    =  .09,
          param2b    =  .09,
          nsims      =  1000,
          names      =  c("Mz twins","Dz twins"),
          type       = "contour2",
          graph_out  =  "../figs/corrx_contour_gamma1_sA2_n180_mzdz1_compare_s1000.pdf",
          alpha      =  0.05,
          threshold  =  0.8) 

# Compare: Scenario 1  - gamma2 - 1000 sims 60:120
corrxplot(data.long  =  dt_s1.long,
          method     = "pearson",
          dist       =  "gamma",
          n1         =  60,
          n2         =  120,
          param1a    =  1,
          param1b    =  1,
          param2a    =  5,
          param2b    =  5,
          nsims      =  1000,
          names      =  c("Mz twins","Dz twins"),
          type       = "contour2",
          graph_out  =  "../figs/corrx_contour_gamma2_sA1_n180_mzdz.5_compare_s1000.pdf",
          alpha      =  0.05,
          threshold  =  0.8) 


# Compare: Scenario 2  - gamma2 -  1000 sims 90:90 
corrxplot(data.long  =  dt_s2.long,
          method     = "pearson",
          dist       =  "gamma",
          n1         =  90,
          n2         =  90,
          param1a    =  1,
          param1b    =  1,
          param2a    =  5,
          param2b    =  5,
          nsims      =  1000,
          names      =  c("Mz twins","Dz twins"),
          type       = "contour2",
          graph_out  =  "../figs/corrx_contour_gamma2_sA2_n180_mzdz1_compare_s1000.pdf",
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
          param1b    =  0,
          param2a    =  1,
          param2b    =  1,
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
          param1b    =  0,
          param2a    =  1,
          param2b    =  1,
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
          param1b    =  0,
          param2a    =  1,
          param2b    =  1,
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
          param1b    =  1.5,
          param2a    =  0.09,
          param2b    =  0.09,
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
          param1b    =  1.5,
          param2a    =  0.09,
          param2b    =  0.09,
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
          param1b    =  1.5,
          param2a    =  0.09,
          param2b    =  0.09,
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
          param1b    =  1,
          param2a    =  5,
          param2b    =  5,
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
          param1b    =  1,
          param2a    =  5,
          param2b    =  5,
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
          param1b    =  1,
          param2a    =  5,
          param2b    =  5,
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


# difference plot  - gamma1 - ratio 1:1
corrxplot(data.long  =  dt_s2.long,
          method     = "pearson",
          dist       = "gamma",
          param1a    =  1.5,
          param2a    =  0.09,
          ratio      =  1,
          n1         =  90,
          n2         =  90,
          nsims      =  1000,
          type       = "diffpowerabs",
          graph_out  =  "../figs/corrx_diffpower_mildskew_90_90_mzdz1_s1000.pdf",
          alpha      =  0.05,
          threshold  =  0.8) 


# difference plot  - gamma1 - ratio 0.5:1
corrxplot(data.long  =  dt_s1.long,
          method     = "pearson",
          dist       = "gamma",
          param1a    =  1.5,
          param2a    =  0.09,
          ratio      =  .5,
          n1         =  60,
          n2         =  120,
          nsims      =  1000,
          type       = "diffpowerabs",
          graph_out  =  "../figs/corrx_diffpower_mildskew_60_120_mzdz1_s1000.pdf",
          alpha      =  0.05,
          threshold  =  0.8) 


# difference plot  - gamma2 - ratio 1:1
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


# difference plot  - gamma2 - ratio 0.5:1
corrxplot(data.long  =  dt_s1.long,
          method     = "pearson",
          dist       = "gamma",
          param1a    =  1,
          param2a    =  5,
          ratio      =  0.5,
          n1         =  60,
          n2         =  120,
          nsims      =  1000,
          type       = "diffpowerabs",
          graph_out  =  "../figs/corrx_diffpower_extrskew_60_120_mzdz1_s1000.pdf",
          alpha      =  0.05,
          threshold  =  0.8) 


# Non-normal distribution simulation examples
# normal reference distribution
dist_x_plot(dist = "normal",
            n1 = 60,n2 = 120,
            param1a = 0,param1b = 0, param2a = 1,param2b = 1,
            rho1 = 0.2, rho2 = 0.5, 
            names = c("MZ","DZ"),
            graph_out = "../figs/distx_normal_60_120.pdf")


# mild skew distribution
dist_x_plot(dist = "gamma",
            n1 = 60,n2 = 120,
            param1a = 1.5,param1b = 1.5, param2a = .09,param2b = .09,
            rho1 = 0.2, rho2 = 0.5, 
            names = c("MZ","DZ"),
            graph_out = "../figs/distx_gamma_mildskew_60_120.pdf")

# extr skew distribution
dist_x_plot(dist = "gamma",
            n1 = 60,n2 = 120,
            param1a = 1,param1b = 1, param2a = 5,param2b = 5,
            rho1 = 0.2, rho2 = 0.5, 
            names = c("MZ","DZ"),
            graph_out = "../figs/distx_gamma_extrskew_60_120.pdf")


# gamma distribution example (mean/shape, and rate/dispersion)
# some positive skew (but distinctly non-normal)
gamma <- genCorGen(1000, nvars = 2, params1 = c(1.5,1.5), params2 = c(0.09,0.09),dist = "gamma",
                   corMatrix = matrix(c(1, 0.5, 0.5, 1), ncol = 2), wide = TRUE)[,2:3]
bivariate_distribution_gamma_m1.5_d0.09_n50 <- ggExtra::ggMarginal(data = as.data.frame(gamma), x = "V1", y = "V2")
plot(bivariate_distribution_gamma_m1.5_d0.09_n50)

# more gamma exploration, using our otherwise defaults
gamma_a <- genCorGen(20, nvars = 2, params1 = c(1,1), params2 = c(0.09,0.09),dist = "gamma",
                     corMatrix = matrix(c(1, 0.5, 0.5, 1), ncol = 2), wide = TRUE)[,2:3]
gamma_b <- genCorGen(90, nvars = 2, params1 = c(1,1), params2 = c(5,5),dist = "gamma",
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
