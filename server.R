#' server.R
#' Author: Peter Ashcroft, ETH Zurich

library(shiny)

source("R/basicFunctions.R")
source("R/simulateModel.R")
source("R/plotFunctions.R")

# Some values which are only used in server ----
y.pos <- c(1,1,0,-1,-1,-1,0,0,1,1)
colour.vec <- ifelse(seq_len(n) <= 3, "#FDB12B", ifelse(seq_len(n) <= 6, "#B97E42", ifelse(seq_len(n) <= 8, "#FA2D1C", "#919191")))
# adj <- createAdjacency(n, entries =
#                          list(
#                            c(1, 2, 1.0), #c(1, 4, 0.1),
#                            c(2, 3, 0.5), c(2, 9, 0.5),
#                            c(3, 4, 0.5), c(3, 7, 0.5),
#                            c(4, 5, 1.0), c(5, 6, 1.0),
#                            c(7, 8, 1.0),
#                            c(9, 10, 1.0)
#                          )
# )

# Function to update sliders upon reset ----
updateSliders <- function(values, session) {
  updateSliderInput(session, "mut.hsc", value = values$ic[[1]])
  for (i in seq_len(n)) {
    updateSliderInput(session, paste0("mut.div", i), value = values$div[[i]])
    updateSliderInput(session, paste0("mut.death", i), value = values$death[[i]])
    updateSliderInput(session, paste0("mut.sr", i), value = values$sr[[i]])
  }
  rm(i)
  updateSliderInput(session, "mut.c29", value = values$C[[1]])
  updateSliderInput(session, "mut.c37", value = values$C[[2]])

  #' Reset WT pars too
  for (i in seq_len(n)) {
    updateSliderInput(session, paste0("wt.div", i), value = default.wt$div[[i]])
    updateSliderInput(session, paste0("wt.death", i), value = default.wt$death[[i]])
    updateSliderInput(session, paste0("wt.sr", i), value = default.wt$sr[[i]])
  }
  rm(i)
  updateSliderInput(session, "wt.c29", value = default.wt$C[[1]])
  updateSliderInput(session, "wt.c37", value = default.wt$C[[2]])
}

# Define server logic ----
shinyServer(function(input, output, session) {
  #' Reset parameter values upon button press
  observeEvent(input$reset_to_default, updateSliders(default.mut, session))
  observeEvent(input$reset_to_WT, {
    params <- default.wt
    params$ic[[1]] <- 0.01
    updateSliders(params, session)
  })
  observeEvent(input$reset_to_late, updateSliders(late.mut, session))
  observeEvent(input$reset_to_RS, updateSliders(RS.mut, session))
  observeEvent(input$reset_to_P339, updateSliders(mut.P339, session))
  observeEvent(input$reset_to_P433, updateSliders(mut.P433, session))
  observeEvent(input$reset_to_P433_alt, updateSliders(mut.P433.alt, session))
  observeEvent(input$reset_to_P328, updateSliders(mut.P328, session))
  observeEvent(input$reset_to_P328_alt, updateSliders(mut.P328.alt, session))

  # observeEvent(input$reset_to_WT, {
  #     updateSliderInput(session, "mut.hsc", value = default.mut$ic)
  #     for (i in seq_len(n)) {
  #         updateSliderInput(session, paste0("mut.div", i), value = default.wt$div[[i]])
  #         updateSliderInput(session, paste0("mut.death", i), value = default.wt$death[[i]])
  #         updateSliderInput(session, paste0("mut.sr", i), value = default.wt$sr[[i]])
  #     }
  #     rm(i)
  #     updateSliderInput(session, "mut.c23", value = default.wt$C[[1]])
  #     updateSliderInput(session, "mut.c34", value = default.wt$C[[2]])
  # })
  #
  #
  # observeEvent(input$reset_to_4, {
  #     updateSliderInput(session, "mut.hsc", value = default.mut$ic)
  #     for (i in seq_len(n)) {
  #         updateSliderInput(session, paste0("mut.div", i), value = default.mut$div[[i]])
  #         updateSliderInput(session, paste0("mut.death", i), value = default.mut$death[[i]])
  #         updateSliderInput(session, paste0("mut.sr", i), value = default.mut$sr[[i]])
  #     }
  #     rm(i)
  #     updateSliderInput(session, "mut.c23", value = default.mut$C[[1]])
  #     updateSliderInput(session, "mut.c34", value = default.mut$C[[2]])
  # })
  #
  #
  # observeEvent(input$reset_to_late, {
  #     updateSliderInput(session, "mut.hsc", value = late.mut$ic)
  #     for (i in seq_len(n)) {
  #         updateSliderInput(session, paste0("mut.div", i), value = late.mut$div[[i]])
  #         updateSliderInput(session, paste0("mut.death", i), value = late.mut$death[[i]])
  #         updateSliderInput(session, paste0("mut.sr", i), value = late.mut$sr[[i]])
  #     }
  #     rm(i)
  #     updateSliderInput(session, "mut.c23", value = late.mut$C[[1]])
  #     updateSliderInput(session, "mut.c34", value = late.mut$C[[2]])
  # })
  #
  #
  # observeEvent(input$reset_to_RS, {
  #     updateSliderInput(session, "mut.hsc", value = RS.mut$ic)
  #     for (i in seq_len(n)) {
  #         updateSliderInput(session, paste0("mut.div", i), value = RS.mut$div[[i]])
  #         updateSliderInput(session, paste0("mut.death", i), value = RS.mut$death[[i]])
  #         updateSliderInput(session, paste0("mut.sr", i), value = RS.mut$sr[[i]])
  #     }
  #     rm(i)
  #     updateSliderInput(session, "mut.c23", value = RS.mut$C[[1]])
  #     updateSliderInput(session, "mut.c34", value = RS.mut$C[[2]])
  # })
  #
  # observeEvent(input$reset_to_S6A, {
  #     updateSliderInput(session, "mut.hsc", value = S6A.mut$ic)
  #     for (i in seq_len(n)) {
  #         updateSliderInput(session, paste0("mut.div", i), value = S6A.mut$div[[i]])
  #         updateSliderInput(session, paste0("mut.death", i), value = S6A.mut$death[[i]])
  #         updateSliderInput(session, paste0("mut.sr", i), value = S6A.mut$sr[[i]])
  #     }
  #     rm(i)
  #     updateSliderInput(session, "mut.c23", value = S6A.mut$C[[1]])
  #     updateSliderInput(session, "mut.c34", value = S6A.mut$C[[2]])
  # })


  output$plot <- renderPlot({

    #' Read WT parameters from input
    #' WT fraction
    wt.ic <- c(1 - input$mut.hsc, rep.int(0, n - 1))
    #' Division rates
    wt.div <- sapply(seq_len(n), function(i) input[[paste0("wt.div", i)]])
    #' Death probabilities
    wt.death <- sapply(seq_len(n), function(i) input[[paste0("wt.death", i)]])
    #' Self-renewal probabilities
    wt.sr <- sapply(seq_len(n), function(i) input[[paste0("wt.sr", i)]])
    #' Lineage bias
    wt.adj <- adj
    wt.adj[2,9] <- input$wt.c29
    wt.adj[2,3] <- 1 - wt.adj[2,9]
    wt.adj[3,7] <- input$wt.c37
    wt.adj[3,4] <- 1 - wt.adj[3,7]


    #' Read mutant parameters from input
    #' Mutant fraction
    mut.ic <- c(input$mut.hsc, rep.int(0, n - 1))
    #' Division rates
    mut.div <- sapply(seq_len(n), function(i) input[[paste0("mut.div", i)]])
    #' Death probabilities
    mut.death <- sapply(seq_len(n), function(i) input[[paste0("mut.death", i)]])
    #' Self-renewal probabilities
    mut.sr <- sapply(seq_len(n), function(i) input[[paste0("mut.sr", i)]])
    #' Lineage bias
    mut.adj <- wt.adj
    mut.adj[2,9] <- input$mut.c29
    mut.adj[2,3] <- 1 - mut.adj[2,9]
    mut.adj[3,7] <- input$mut.c37
    mut.adj[3,4] <- 1 - mut.adj[3,7]

    #' Steady state and clonal fraction
    steadyState <- numericSteadyState(
      wt.ic = wt.ic, mut.ic = mut.ic,
      wt.div = wt.div, mut.div = mut.div,
      wt.sr = wt.sr, mut.sr = mut.sr,
      wt.adj = wt.adj, mut.adj = mut.adj,
      wt.death = wt.death, mut.death = mut.death
    )
    clonal.fraction <- getClonalFraction(steadyState)

    if (input[["par.view"]] == "relative") {
      plot <- plotFractionParams(
        fraction.df = clonal.fraction,
        wt.div = wt.div, mut.div = mut.div,
        wt.sr = wt.sr, mut.sr = mut.sr,
        wt.adj = wt.adj, mut.adj = mut.adj,
        wt.death = wt.death, mut.death = mut.death,
        colour.vec = colour.vec,
        show.tree = T, y.pos = y.pos,
        compartment.names = name.vec, lab.offset = -0.35, y.lim = -1.75
      )
    }
    else{
      plot <- plotFractionParamsAbsolute(
        fraction.df = clonal.fraction,
        wt.div = wt.div, mut.div = mut.div,
        wt.sr = wt.sr, mut.sr = mut.sr,
        wt.adj = wt.adj, mut.adj = mut.adj,
        wt.death = wt.death, mut.death = mut.death,
        colour.vec = colour.vec,
        #alpha.vec = ifelse(seq_len(n) %in% c(2,5,8), 0.6, 1),
        show.tree = T, y.pos = y.pos, tree.height = 1.25,
        compartment.names = name.vec, lab.offset = -0.35, y.lim = -1.75
      )
    }
    return(plot)
  })

})
