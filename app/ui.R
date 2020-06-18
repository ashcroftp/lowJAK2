#' ui.R
#' Author: Peter Ashcroft, ETH Zurich

library(shiny)
library(shinyWidgets)

# Define UI for application that draws a histogram
shinyUI(
  fluidPage(
    withMathJax(),
    
    #' Slider style
    tags$head(
      tags$style(
        type = "text/css", 
        "label.control-label, .selectize-control.single { 
          display: table-cell; 
          text-align: left; 
          vertical-align: middle;
        }
        .form-group {
          display: table-row;
          width: 100%;
        }
        label.control-label {
          padding-right: 20px;
          width: 30%;
        }
        .selectize-control.single {
          padding-right: 20px;
          width: 100%;
        }"
      )
    ),
    
    # Application title
    titlePanel("MPN patients with low mutant JAK2 allele burden show late expansion restricted to erythroid and megakaryocytic lineages"),
    
    # Sidebar with a slider input for number of bins 
    fluidRow(
      column(10, offset = 1,
             addSpinner(plotOutput("plot"), spin = "circle", color = "#1b7cff")
      )
    ),
    
    hr(),
    h3("Controls"),
    fluidRow(
      column(12,
             actionButton("reset_to_default", "Reset to initial"),
             actionButton("reset_to_WT", "Reset to WT values"),
             actionButton("reset_to_P339", "Reset to P339"),
             actionButton("reset_to_P433", "Reset to P433"),
             actionButton("reset_to_P433_alt", "Reset to P433 - alt"),
             actionButton("reset_to_P328", "Reset to P328"),
             actionButton("reset_to_P328_alt", "Reset to P328 - alt"),
             actionButton("reset_to_late", "Reset to late expansion"),
             actionButton("reset_to_RS", "Reset to RS")
      )
    ),
    fluidRow(
      column(12,
             radioButtons("par.view", "Parameter view: ", choices = c("relative", "absolute"), inline = T)
      )
    ),
    
    hr(),
    h3("Mutant parameters"),
    
    fluidRow(
      column(3,
             h4("HSC mutant fraction"),
             sliderInput("mut.hsc", "\\(y_1^*\\)", min = 0, max = 1, value = default.mut$ic[[1]])
      )
    ),
    
    fluidRow(
      column(3,
             h4("Self-renewal probabilities, \\(a'\\)"),
             lapply(seq_len(n), function(i) {
               sliderInput(paste0("mut.sr", i), name.vec[[i]],
                           min = 0.0, max = 1, value = default.mut$sr[[i]])
             })
      ),
      column(3,
             h4("Division rates, \\(b'\\)"),
             lapply(seq_len(n), function(i) {
               sliderInput(paste0("mut.div", i), name.vec[[i]],
                           min = 0.1, max = 10, value = default.mut$div[[i]])
             })
             
      ),
      column(3,
             h4("Lineage bias, \\(c'\\)"),
             sliderInput("mut.c29", paste(name.vec[[2]], name.vec[[9]], sep = " \\(\\to\\) "), min = 0, max = 1, value = default.mut$C[[1]]),
             sliderInput("mut.c37", paste(name.vec[[3]], name.vec[[7]], sep = " \\(\\to\\) "), min = 0, max = 1, value = default.mut$C[[2]])
      ),
      column(3,
             h4("Death probabilities, \\(d'\\)"),
             lapply(seq_len(n), function(i) {
               sliderInput(paste0("mut.death", i), name.vec[[i]],
                           min = 0.0, max = 1, value = default.mut$death[[i]])
             })
      )
    ),
    
    hr(),
    h3("WT parameters"),
    fluidRow(
      column(3,
             h4("Self-renewal probabilities, \\(a\\)"),
             lapply(seq_len(n), function(i) {
               sliderInput(paste0("wt.sr", i), name.vec[[i]],
                           min = 0.0, max = 1, value = default.wt$sr[[i]])
             })
      ),
      column(3,
             h4("Division rates, \\(b\\)"),
             lapply(seq_len(n), function(i) {
               sliderInput(paste0("wt.div", i), name.vec[[i]],
                           min = 0.1, max = 10, value = default.wt$div[[i]])
             })
      ),
      column(3,
             h4("Lineage bias, \\(c\\)"),
             sliderInput("wt.c29", paste(name.vec[[2]], name.vec[[9]], sep = " \\(\\to\\) "), min = 0, max = 1, value = default.wt$C[[1]]),
             sliderInput("wt.c37", paste(name.vec[[3]], name.vec[[7]], sep = " \\(\\to\\) "), min = 0, max = 1, value = default.wt$C[[2]])
      ),
      column(3,
             h4("Death probabilities, \\(d\\)"),
             lapply(seq_len(n), function(i) {
               sliderInput(paste0("wt.death", i), name.vec[[i]],
                           min = 0.0, max = 1, value = default.wt$death[[i]])
             })
      )
    )
  )
)
