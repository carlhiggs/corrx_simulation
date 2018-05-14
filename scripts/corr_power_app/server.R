# install.packages("ggplot2")
# install.packages("scales")
require("shiny")
require("ggplot2")
require("scales")


# Mz and Dz Twin pairs correlations power plot
shinyServer(function(input, output) {
  mydata <- reactive({
    # Model Parameters:
      r1        <- input$r1
      r2        <- input$r2
      rdiff     <- r1 - r2
      maxn      <- input$maxn
      mzdz      <- input$mzdz
      # alpha     <- input$alpha
      sidedness <- 2
      # method    <- input$method
      method    <- "pearson"
      # Step 1: calculate Fisher's Z
      z1        <- atanh(r1)
      z2        <- atanh(r2)
      # Step    2: take difference
      zdiff     <- z1-z2
      # Step    3: calculate standard errmzor and test statistic
      tot_n         <- 15:maxn
      mzdz_vec <- seq(0,10,0.01)
      # n_p1         <- cbind(mz = tot_n*mzdz,dz = tot_n*(1-mzdz))
      n_p1         <- cbind(mz = tot_n/(mzdz+1)*mzdz,dz = tot_n/(mzdz+1))
      n_p2         <- cbind(mz = maxn/(mzdz_vec+1)*mzdz_vec,dz = maxn/(mzdz_vec+1))
      z_se_p1      <- sqrt(rowSums(1/(n_p1-3)))
      z_se_p2      <- sqrt(rowSums(1/(n_p2-3)))
      z_test_p1    <- zdiff/z_se_p1
      z_test_p2    <- zdiff/z_se_p2
      # optionally return p-value for observing diff at least this large under H0
      # z_p    <- sidedness*pnorm(-abs(z_test))
      z_ref    <- c(qnorm(1-0.1/2),qnorm(1-0.05/2))
      z_power1a <- 1-pnorm(z_ref[1] - abs(z_test_p1))
      z_power1b <- 1-pnorm(z_ref[2] - abs(z_test_p1))
      z_power1  <- rbind(cbind(z_power1a,z_ref[1]),cbind(z_power1b,z_ref[2]))
      z_power2a <- 1-pnorm(z_ref[1] - abs(z_test_p2))
      z_power2b <- 1-pnorm(z_ref[2] - abs(z_test_p2))
      z_power2  <- rbind(cbind(z_power2a,z_ref[1]),cbind(z_power2b,z_ref[2]))
      colnames(z_power1) <- colnames(z_power2) <- c("power","ref")
      
    # Collect and output results
    params1 <- paste0("MZ:DZ ratio: ",mzdz,"; rho_mz: ",r1,"; rho_dz: ",r2,"; delta: ",rdiff)
    params2 <- paste0("N: ",maxn,"; rho_mz: ",r1,"; rho_dz: ",r2,"; delta: ",rdiff)
    data1   <- cbind(n = tot_n,mzdz = mzdz,        mz = n_p1[,1],dz = n_p1[,2] , power1 = z_power1)
    data2   <- cbind(n = maxn, mzdz_vec = mzdz_vec,mz = n_p2[,1],dz = n_p2[,2] , power2 = z_power2)
    list(data1 =data1, 
         data2 = data2, 
         params1 = params1,
         params2 = params2 )
    })
  
  output$datatable <-
    renderTable({
      dt <- mydata()[["data"]]
      dt
    })
    

  output$graph1 <- renderPlot({
   p <- ggplot(as.data.frame(mydata()[["data1"]]),
        aes(x=n, y=power, group=as.character(round(ref,2))))
   p <- p +
     geom_line(aes(colour = as.character(round(ref,2))), size=1, alpha=.75) +
     ggtitle(paste0("Power estimate given parameters (",mydata()[["params1"]],")"))+
     scale_x_continuous(name="N")+
     scale_y_continuous(labels = comma, name="Power",limits = c(0,1), expand = c(0,0) ) + 
     labs(colour = "Normal ordinate")
   print(p)
  })
  
  output$graph2 <- renderPlot({
    p <- ggplot(as.data.frame(mydata()[["data2"]]),
                aes(x=mzdz_vec, y=power, group=as.character(round(ref,2))), log="x")
    p <- p +
      geom_line(aes(colour = as.character(round(ref,2))), size=1, alpha=.75) +
      ggtitle(paste0("Power estimate given parameters (",mydata()[["params2"]],")"))+
      scale_x_continuous(name="MZ:DZ ratio (log scale)", trans='log',breaks=c(seq(0,1,0.2),seq(2,10,2)))+
      scale_y_continuous(labels = comma, name="Power",limits = c(0,1), expand = c(0,0) ) + 
      labs(colour = "Normal ordinate")
    print(p)
  })
})