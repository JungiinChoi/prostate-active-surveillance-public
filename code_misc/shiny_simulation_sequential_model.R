library(shiny)
library(tidyverse)
ui <- fluidPage(
  titlePanel(""), # App title ----
  sidebarLayout( # Sidebar layout with input and output definitions ----
    sidebarPanel( # Sidebar panel for inputs ----
      numericInput(inputId = "xprob",
                  label = "X=1 probability",
                  value = 0.5),
      numericInput(inputId = "yprob1_x0",
                   label = "Y level 1 probability when x=0",
                   value = 0.3),
      numericInput(inputId = "yprob2_x0",
                   label = "Y level 2 probability when x=0",
                   value = 0.4),
      numericInput(inputId = "yprob3_x0",
                   label = "Y level 3 probability when x=0",
                   value = 0.2),
      numericInput(inputId = "yprob4_x0",
                   label = "Y level 4 probability when x=0",
                   value = 0.1),
      numericInput(inputId = "yprob1_x1",
                   label = "Y level 1 probability when x=1",
                   value = 0.1),
      numericInput(inputId = "yprob2_x1",
                   label = "Y level 2 probability when x=1",
                   value = 0.15),
      numericInput(inputId = "yprob3_x1",
                   label = "Y level 3 probability when x=1",
                   value = 0.25),
      numericInput(inputId = "yprob4_x1",
                   label = "Y level 4 probability when x=1",
                   value = 0.5)
    ),
    mainPanel(  # Main panel for displaying outputs ----
      plotOutput(outputId = "probPlot") # Output: Histogram ----
    )
  )
)
# ui <- fluidPage(
#   titlePanel(""), # App title ----
#   sidebarLayout( # Sidebar layout with input and output definitions ----
#                  sidebarPanel( # Sidebar panel for inputs ----
#                                numericInput(inputId = "xprob",
#                                             label = "x=1 probability",
#                                             value = 0.5),
#                                textInput(inputId = "yprob_x0",
#                                          placeholder = "Y distribution for x=0, ex, 0.3, 0.4, 0.2, 0.1",
#                                          label = "Y distribution for x=0"),
#                                textInput(inputId = "yprob_x1",
#                                          label = "Y distribution for x=1",
#                                          placeholder = "Y distribution for x=1, ex, 0.1, 0.15, 0.25, 0.5")),             
#                  mainPanel(  # Main panel for displaying outputs ----
#                              plotOutput(outputId = "probPlot") # Output: Histogram ----
#                  )
#   )
# )
server <- function(input, output) {
  # extract <- function(text) {
  #   text <- gsub(" ", "", text)
  #   split <- strsplit(text, ",", fixed = FALSE)[[1]]
  #   as.numeric(split)
  # }
  output$probPlot <- renderPlot({
    # yprob_x0 <- extract(input$yprob_x0)
    # yprob_x1 <- extract(input$yprob_x1)
    dat <- data.frame(id = 1:1000)
    set.seed(2022)
    dat$x <- rbinom(1000, 1, input$xprob)
    set.seed(2021)
    for(i in 1:length(dat$x)){
      dat$y[i] <- ifelse(dat$x[i] == 1, 
                         sample(c(1,2,3,4), 1000, T, c(input$yprob1_x1, input$yprob2_x1, 
                                                       input$yprob3_x1, input$yprob4_x1)),
                         sample(c(1,2,3,4), 1000, T, c(input$yprob1_x0, input$yprob2_x0, 
                                                       input$yprob3_x0, input$yprob4_x0)))
    }
    fit1<-glm(y_lev1 ~ x, family = binomial, 
              data = dat %>% 
                mutate(y_lev1 = factor(ifelse(y == 1, 1, 0), levels=c(0,1))))
    fit2<-glm(y_lev2 ~ x, family = binomial, 
              data = dat %>% 
                filter(y > 1) %>% 
                mutate(y_lev2 = factor(ifelse(y == 2, 1, 0), levels = c(0,1))))
    fit3<-glm(y_lev3 ~ x, family = binomial, 
              data = dat %>% 
                filter(y > 2) %>% 
                mutate(y_lev3 = factor(ifelse(y == 3, 1, 0), levels=c(0,1))))
    
    plot_x0 <- data.frame(id="x=0", x = 0)
    p_y1_given_x0 <- predict(fit1, newdata=plot_x0, type = "response")
    p_y2_given_x0ge1 <- predict(fit2, newdata=plot_x0, type = "response")
    p_y3_given_x0ge2 <- predict(fit3, newdata=plot_x0, type = "response")
    
    plot_x0$p1 <- p_y1_given_x0
    plot_x0$p2 <- (1-p_y1_given_x0) * p_y2_given_x0ge1
    plot_x0$p3 <- (1-p_y1_given_x0)*(1-p_y2_given_x0ge1) * p_y3_given_x0ge2
    plot_x0$p4 <- 1-plot_x0$p1-plot_x0$p2-plot_x0$p3
    
    plot_x1 <- data.frame(id="x=1", x = 1)
    p_y1_given_x1 <- predict(fit1, newdata=plot_x1, type = "response")
    p_y2_given_x1ge1 <- predict(fit2, newdata=plot_x1, type = "response")
    p_y3_given_x1ge2 <- predict(fit3, newdata=plot_x1, type = "response")
    
    plot_x1$p1 <- p_y1_given_x1
    plot_x1$p2 <- (1-p_y1_given_x1) * p_y2_given_x1ge1
    plot_x1$p3 <- (1-p_y1_given_x1) * (1-p_y2_given_x1ge1) * p_y3_given_x1ge2
    plot_x1$p4 <-1-plot_x1$p1-plot_x1$p2-plot_x1$p3
    
    
    p0 <- plot_x0 %>% 
      gather("level", "prob", p1:p4) 
    p1 <- plot_x1 %>% 
      gather("level", "prob", p1:p4) 
    
    dat_prob <- data.frame(id = rep("data", 4),
                           x = NA,
                            level = c("p1", "p2", "p3", "p4"),
                            prob = table(dat$y)/sum(table(dat$y)))
    dat_prob <- dat_prob[,c("id", "x","level", "prob.Freq")]
    colnames(dat_prob) <- c("id", "x","level", "prob")
    pplot <- rbind(dat_prob,p0, p1)
    ggplot(data = pplot)+
      geom_bar(aes(x = level, y = prob, fill = factor(id)),stat="identity",position = "dodge")+
      theme_classic()+ylim(c(0,1))+
      scale_fill_manual(name = "",
                          values = c("black","#00BFC4", "#F8766D"))
  })
}
shinyApp(ui = ui, server = server)
