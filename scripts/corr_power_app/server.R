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
      maxn      <- input$maxn
      mzdz      <- input$mzdz
      alpha     <- 0.05
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
      n         <- cbind(mz = tot_n*mzdz,dz = tot_n*(1-mzdz))
      z_se      <- sqrt(rowSums(1/(n-3)))
      z_test    <- zdiff/z_se
      # optionally return p-value for observing diff at least this large under H0
      # z_p    <- sidedness*pnorm(-abs(z_test))
      z_ref   <- qnorm(1-alpha/sidedness)
      z_power <- 1-pnorm(z_ref - abs(z_test))
      
    # Table of results
    data <- cbind(n = tot_n,mzdz = mzdz,mz = n[,1],dz = n[,2] , power = z_power)
    list(data =data)
    })
  
  # output$datatable <- 
  #   renderTable({
  #     dt <- mydata()[["data"]]
  #     dt
  #   })
    

  output$graph1 <- renderPlot({
   p <- ggplot(as.data.frame(mydata()[["data"]]),
        aes(x=n, y=power))
   p <- p +
     geom_line(size=1, alpha=.75) +
     ggtitle("Power estimate")+
     scale_x_continuous(name="N")+
     scale_y_continuous(labels = comma, name="Power",limits = c(0,1), expand = c(0,0) )
   print(p)
  })
})