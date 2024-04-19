### shiny UI 

library(shiny)
library(shinythemes)
library(tidyverse)
library(geomtextpath)
library(DT)
library(ggtree)
library(cowplot)
library(scales)

ui <- shinyUI(navbarPage(
  #title="EVE Results Explorer",
  theme=shinytheme('flatly'),
  
  title = 'shinyEVE v1.3',
  
  tabsetPanel(tabPanel("Main",
    sidebarLayout(
      sidebarPanel(
        h3('Filter Results Data'),
        #select clade, gene, and temperature
        fluidRow(
          column(4, selectInput("iclade", "Select a Clade", choices = c("flexible", "hibernators"), selected='flexible')),
          column(4, uiOutput("gene_select")),
          column(4, radioButtons("itemp", "Select Temperature", choices = c("32C", "41C"), selected = '41C'))
        ),
        fluidRow(
          column(4, radioButtons('show_control', 'Show Control Data', choices = c('Yes', 'No'), selected = 'No')),
          column(4, radioButtons('idatalevel', 'Hypothesis Level', choices = c('exp.level', 'fold.change'), selected = 'exp.level')),
          column(4, radioButtons('show_thetas', 'Show Thetas', choices = c('Yes', 'No'), selected = 'No'))
        ),
        h3('Gene-Specific Stats'),
        plotOutput('sidebar_plot'),
        br(), br(),
        h3('Plot Results'),
        #refresh buttons
        fluidRow(
          column(4,actionButton("refresh_table", "Refresh Table")),
          column(4,actionButton("refresh_plot", "Plot Results")),
          column(4,div(tags$a(imageOutput('logo'),
                              href='https://github.com/dylan-unlv')))
        )
        
      ),
      mainPanel(
        plotOutput("plot"),
        #display table for gene exploring
        DT::dataTableOutput("res_table")
      )
    )
  ))
) )