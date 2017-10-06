#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(quantmod)
library(dygraphs)
library(Rsolnp)
library(lubridate)
library(highcharter)
library(foreach)
# library(parallel)
# cl <- makePSOCKcluster(2)

# Define server logic required to draw a histogram

shinyServer(function(input, output) {
    etf_returns <- eventReactive(input$go, {
        tickers <- unlist(strsplit(input$tickers,"\n"))
        sectors <- unlist(strsplit(input$names, "\n"))
        
        etf_yearly_returns <- function(ticker) {
            symbols <- getSymbols(ticker, src = "google", from =  floor_date(Sys.Date(), "years") - years(input$years))
            etf_prices <- do.call(merge, lapply(symbols, function(x) Cl(get(x))))
            
            etf_returns <- do.call(merge, lapply(etf_prices, 
                                                 function(x) periodReturn(x, 
                                                                          period = 'yearly', 
                                                                          type = 'log',
                                                                          subset = '1990::')))
            
            # Change the column names to the sector names from our dataframe above.
            
            colnames(etf_returns) <- sectors
            
            etf_returns*100
            
        }
        etf_yearly_returns(tickers)
    })
    
    solution <- eventReactive(input$go, {
        tickers <- unlist(strsplit(input$tickers,"\n"))
        sectors <- unlist(strsplit(input$names, "\n"))
        
        i <- length(tickers)
        m <- 11L
        
        z <- seq(-1, 1, length.out = m)
        
        # Pick a desired retrun mu0 and a risk sigma
        mu0 <- as.numeric(input$mu0)
        sigma <- as.numeric(input$sigma)
        
        # Average returns for eachs sector
        ret <- apply(etf_returns(), 2, mean)
        cov_ <- cov(etf_returns())
        
        # Initial values
        init_values <- rep(1/m, i*m)
        
        # Objective function to be minimized 
        obj <- function(Pi) {
            h <- ifelse(Pi == 0, 0, Pi*log(Pi))
            sum(h)
        }
        
        # Equality constraints
        eqfun <- function(Pi) {
            Pi_ <- matrix(Pi, nrow = i, ncol = m)
            z1 <- sum(Pi_ %*% z)
            z2 <- rowSums(Pi_)
            return(c(z1, z2))
        }
        
        # Values of the constraints
        eqB <- rep(1, i + 1)
        
        # Inequality Constraints
        ineqfun <- function(Pi) {
            Pi_ <- matrix(Pi, nrow = i, ncol = m)
            p <- as.vector(Pi_ %*% z)
            z1 <- p %*% ret
            z2 <- sqrt(p %*% cov_ %*% p)
            return(c(z1, z2))
        }
        
        ineqLB <- c(mu0, 0)
        ineqUB <- c(Inf, sigma)
        
        # Lower bound on the p (probabilities must be positive)
        
        LB <- rep(0, i*m)
        UB <- rep(Inf, i*m)
        
        # Solve the problem
        solution <- solnp(pars = init_values,
              fun = obj,
              eqfun = eqfun,
              eqB = eqB,
              ineqfun = ineqfun,
              ineqLB = ineqLB,
              ineqUB = ineqUB,
              LB = LB,
              control = list(trace = FALSE))
        
        p1 <- as.vector(matrix(solution$pars, nrow = i, ncol = m) %*% z)
        
        # Now solve the simple portfolio case
        eqfun <- function(p) {
            sum(p)
        }
        eqB <- 1
        ineqfun <- function(p) {
            z1 <- as.vector(p %*% ret)
            z2 <- as.vector(sqrt(p %*% cov_ %*% p))
            return(c(z1, z2))
        }
        LB <- rep(0, i)
        init_values <- rep(1/i, i)
        solution2 <- tryCatch(solnp(pars = init_values, fun = obj, eqfun = eqfun, eqB = eqB,
                           ineqfun = ineqfun, ineqLB = ineqLB, ineqUB = ineqUB, LB = LB,
                           control = list(trace = FALSE)),
                           error = function(e) {
                               warning(e)
                               return(rep(0,i))
                           })
        p2 <- solution2$pars
        cbind(p1, p2)
    })
    
    sectors <- eventReactive(input$go, {
        unlist(strsplit(input$names, "\n"))
    })
    
    tickers <- eventReactive(input$go, {
        unlist(strsplit(input$tickers,"\n"))
    })
    
    output$allocations <- renderHighchart({
        ret <- etf_returns()/100 + 1
        yearly_weighted_1 <- apply(ret,1, weighted.mean, solution()[,1])
        yearly_weighted_2 <- apply(ret,1, weighted.mean, solution()[,2])
        yearly_unweighted <- apply(ret,1, mean)
        mu_1 <- round((prod(yearly_weighted_1) - 1) * 100, 1)
        mu_2 <- round((prod(yearly_weighted_2) - 1) * 100, 1)
        mu_3 <- round((prod(yearly_unweighted) - 1) * 100, 1)
        
        ret_1 <- apply(etf_returns(),1, weighted.mean, solution()[,1])
        ret_2 <- apply(etf_returns(),1, weighted.mean, solution()[,2])
        ret_3 <- apply(etf_returns(),1, mean)
        
        # print(ret_2)
        
        s1 <- round((sd(ret_1)),2)
        s2 <- round((sd(ret_2)),2)
        s3 <- round((sd(ret_3)),2)
        
        txt <- paste0("Total Return over the period analyzed<br/>",
                      "Equally weighted portfolio Return: ", mu_3, "%,", 
                      " Portfolio with short sales Return: ", mu_1, "%,",
                      " Portfolio without short sales Return : ", mu_2, "%,",
                      "<br/>Equally weighted portfolio sd: ", s3, "%",
                      " Portfolio sd with short sales: ", s1, "%,",
                      " Portfolio sd without short sales: ", s2, "%"
                      )
        highchart() %>%
            hc_title(text = "Portfolio Allocation") %>%
            hc_subtitle(text = txt) %>%
            hc_add_series_labels_values(labels = sectors(), 
                                        values = round(solution()[,1]*100, 2), 
                                        type = "bar",
                                        name = "Allocation with shorts", 
                                        colorByPoint = FALSE,
                                        dataLabels = list(format =  "{point.name} {point.y}%")) %>%
            hc_add_series_labels_values(labels = sectors(), 
                                        values = round(solution()[,2]*100, 2), 
                                        type = "bar",
                                        name = "Allocation without shorts", 
                                        colorByPoint = FALSE,
                                        dataLabels = list(format =  "{point.name} {point.y}%")) %>%
            hc_xAxis(categories = sectors()) %>%
            hc_tooltip(pointFormat = "{point.y}%")
        
    })
    output$stock_returns <- renderDygraph({
        dygraph(cbind("Equally Weighted" = apply(etf_returns(), 1, mean),
                      "Info-Metrics with short sales" = apply(etf_returns(), 1, weighted.mean, solution()[,1]),
                      "Info-Metrics without short sales" = apply(etf_returns(), 1, weighted.mean, solution()[,2])),
                ylab = "Return per year", main = "Return per year") %>%
            dyOptions(drawPoints = TRUE, pointSize = 2)
    })
    
})