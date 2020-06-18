#' hematopoiesis-figures.R
#' Author: Peter Ashcroft, ETH Zurich

#' Replicating specific experimental results.

#' Dependencies:
source("R/basicFunctions.R")
source("R/simulateModel.R")
source("R/plotFunctions.R")


# 0. Parameter and function definitions for this section ----
name.vec <- c(`1` = "HSC", `2` = "CMP", `3` = "MEP", `4` = "MkP", `5` = "Meg", `6` = "PLT",
              `7` = "EP", `8` = "RET", `9` = "CFU-G", `10` = "GRA")
n <- length(name.vec)
adj <- createAdjacency(n, entries =
                         list(
                           c(1, 2, 1.0), #c(1, 4, 0.1),
                           c(2, 3, 0.5), c(2, 9, 0.5),
                           c(3, 4, 0.5), c(3, 7, 0.5),
                           c(4, 5, 1.0), c(5, 6, 1.0),
                           c(7, 8, 1.0),
                           c(9, 10, 1.0)
                         )
)

#' Initial conditions
wt.ic <- c(0.99, rep.int(0, n - 1))
mut.ic <- c(0.01, rep.int(0, n - 1))
#' Wild-type cell parameters
wt.div <- c(1.0, 1.4, 1.8, 2.2, 2.6, 2.0, 2.2, 1.0, 2.8, 2.0)
wt.sr <- c(0.5, 0.3, 0.3, 0.3, 0.3, 0.0, 0.3, 0.0, 0.3, 0.0)
wt.adj <- adj
wt.death <- c(0, rep.int(0.2, n - 1))
wt.death[7] <- 0.6 # Higher due to lack of EPO

#' Save this information for reproducibility
#save(name.vec, n, adj, wt.ic, wt.div, wt.sr, wt.adj, wt.death, file = "zz-data/WT-params.RData")

#' Function to plot everything
plotEverything <- function(n, wt.ic, mut.ic, wt.div, mut.div, wt.sr, mut.sr, wt.adj, mut.adj, wt.death, mut.death, compartment.names) {
  #' Steady state and clonal fraction
  steadyState <- numericSteadyState(
    wt.ic = wt.ic, mut.ic = mut.ic,
    wt.div = wt.div, mut.div = mut.div,
    wt.sr = wt.sr, mut.sr = mut.sr,
    wt.adj = wt.adj, mut.adj = mut.adj,
    wt.death = wt.death, mut.death = mut.death
  )
  clonal.fraction <- getClonalFraction(steadyState)
  #' Plot definitions
  y.pos <- c(1,1,0,-1,-1,-1,0,0,1,1)
  #colour.vec <- ifelse(seq_len(n) <= 6, "#FDB12B", ifelse(seq_len(n) <= 8, "#FA2D1C", "#919191"))
  colour.vec <- ifelse(seq_len(n) <= 3, "#FDB12B", ifelse(seq_len(n) <= 6, "#B97E42", ifelse(seq_len(n) <= 8, "#FA2D1C", "#919191")))
  #' Plot the data
  plotFractionParams(
    fraction.df = clonal.fraction,
    wt.div = wt.div, mut.div = mut.div,
    wt.sr = wt.sr, mut.sr = mut.sr,
    wt.adj = wt.adj, mut.adj = mut.adj,
    wt.death = wt.death, mut.death = mut.death,
    colour.vec = colour.vec,
    #alpha.vec = ifelse(seq_len(n) %in% c(2,5,8), 0.6, 1),
    show.tree = T, y.pos = y.pos,
    compartment.names = name.vec, lab.offset = -0.35, y.lim = -1.75
  )
}


# 1. Platelet-biased pattern (P433) ----
## 1.1. Modify MEP lineage bias ----
#' Mutant cell parameters
mut.div <- wt.div
mut.div[1] <- 1.8 * wt.div[1]
mut.div[2] <- 1.8 * wt.div[2]
mut.div[3] <- 1.8 * wt.div[3]
mut.div[4] <- 1.8 * wt.div[4]
mut.div[5] <- 1.8 * wt.div[5]
mut.div[7] <- 2.0 * wt.div[7]
mut.div[9] <- 1.8 * wt.div[9]

mut.sr <- wt.sr
mut.sr[4] <- 0.49
mut.sr[5] <- 0.49
mut.sr[7] <- 0.46

mut.adj <- wt.adj
mut.adj[2,3] <- 0.7
mut.adj[2,9] <- 1 - mut.adj[2,3]
mut.adj[3,4] <- 0.95
mut.adj[3,7] <- 1 - mut.adj[3,4]

mut.death <- wt.death
mut.death[5] <- 0.1
mut.death[7] <- 0.1
#mut.death[8] <- 0

#' Save this information for reproducibility
#save(mut.ic, mut.div, mut.sr, mut.adj, mut.death, file = "zz-data/MUT-params-P433.RData")

plotEverything(
  n = n,
  wt.ic = wt.ic, mut.ic = mut.ic,
  wt.div = wt.div, mut.div = mut.div,
  wt.sr = wt.sr, mut.sr = mut.sr,
  wt.adj = wt.adj, mut.adj = mut.adj,
  wt.death = wt.death, mut.death = mut.death,
  compartment.names = name.vec
)
#savePlot(filename = "Figures/P433.pdf", width = 200, height.mult = 0.7)
plot.data <- createTable(
  wt.ic = wt.ic, mut.ic = mut.ic,
  wt.div = wt.div, mut.div = mut.div,
  wt.sr = wt.sr, mut.sr = mut.sr,
  wt.adj = wt.adj, mut.adj = mut.adj,
  wt.death = wt.death, mut.death = mut.death,
  compartment.names = name.vec
)
#save(plot.data, file = "zz-data/P433.RData")

## 1.2. No change to MEP lineage bias ----
#' Mutant cell parameters
mut.div <- wt.div
mut.div[1] <- 1.8 * wt.div[1]
mut.div[2] <- 1.8 * wt.div[2]
mut.div[3] <- 1.8 * wt.div[3]
mut.div[4] <- 1.8 * wt.div[4]
mut.div[5] <- 1.8 * wt.div[5]
mut.div[7] <- 2.0 * wt.div[7]
mut.div[9] <- 1.8 * wt.div[9]

mut.sr <- wt.sr
mut.sr[4] <- 0.49
mut.sr[5] <- 0.49

mut.adj <- wt.adj
mut.adj[2,3] <- 0.7
mut.adj[2,9] <- 1 - mut.adj[2,3]

mut.death <- wt.death
mut.death[4] <- 0.1
mut.death[5] <- 0.1
mut.death[7] <- 0.5

#' Save this information for reproducibility
#save(mut.ic, mut.div, mut.sr, mut.adj, mut.death, file = "zz-data/MUT-params-P433-alt.RData")

plotEverything(
  n = n,
  wt.ic = wt.ic, mut.ic = mut.ic,
  wt.div = wt.div, mut.div = mut.div,
  wt.sr = wt.sr, mut.sr = mut.sr,
  wt.adj = wt.adj, mut.adj = mut.adj,
  wt.death = wt.death, mut.death = mut.death,
  compartment.names = name.vec
)
#savePlot(filename = "Figures/P433-alt.pdf", width = 200, height.mult = 0.7)
plot.data <- createTable(
  wt.ic = wt.ic, mut.ic = mut.ic,
  wt.div = wt.div, mut.div = mut.div,
  wt.sr = wt.sr, mut.sr = mut.sr,
  wt.adj = wt.adj, mut.adj = mut.adj,
  wt.death = wt.death, mut.death = mut.death,
  compartment.names = name.vec
)
#save(plot.data, file = "zz-data/P433-alt.RData")

# 2. Red cell biased pattern (P328) ----
## 2.1. Modify MEP lineage bias ----
#' Mutant cell parameters
mut.div <- wt.div
mut.div[1] <- 1.8 * wt.div[1]
mut.div[2] <- 1.8 * wt.div[2]
mut.div[3] <- 1.8 * wt.div[3]
mut.div[4] <- 1.8 * wt.div[4]
mut.div[5] <- 1.8 * wt.div[5]
mut.div[7] <- 2.0 * wt.div[7]
mut.div[9] <- 1.8 * wt.div[9]

mut.sr <- wt.sr
mut.sr[5] <- 0.49
mut.sr[7] <- 0.46

mut.adj <- wt.adj
mut.adj[2,3] <- 0.7
mut.adj[2,9] <- 1 - mut.adj[2,3]
mut.adj[3,4] <- 0.4
mut.adj[3,7] <- 1 - mut.adj[3,4]

mut.death <- wt.death
mut.death[5] <- 0.1
mut.death[7] <- 0.1
#mut.death[8] <- 0

#' Save this information for reproducibility
#save(mut.ic, mut.div, mut.sr, mut.adj, mut.death, file = "zz-data/MUT-params-P328.RData")

plotEverything(
  n = n,
  wt.ic = wt.ic, mut.ic = mut.ic,
  wt.div = wt.div, mut.div = mut.div,
  wt.sr = wt.sr, mut.sr = mut.sr,
  wt.adj = wt.adj, mut.adj = mut.adj,
  wt.death = wt.death, mut.death = mut.death,
  compartment.names = name.vec
)
#savePlot(filename = "Figures/P328.pdf", width = 200, height.mult = 0.7)
plot.data <- createTable(
  wt.ic = wt.ic, mut.ic = mut.ic,
  wt.div = wt.div, mut.div = mut.div,
  wt.sr = wt.sr, mut.sr = mut.sr,
  wt.adj = wt.adj, mut.adj = mut.adj,
  wt.death = wt.death, mut.death = mut.death,
  compartment.names = name.vec
)
#save(plot.data, file = "zz-data/P328.RData")

## 2.2. No change to MEP lineage bias ----
#' Mutant cell parameters
mut.div <- wt.div
mut.div[1] <- 1.8 * wt.div[1]
mut.div[2] <- 1.8 * wt.div[2]
mut.div[3] <- 1.8 * wt.div[3]
mut.div[4] <- 1.8 * wt.div[4]
mut.div[5] <- 1.8 * wt.div[5]
mut.div[7] <- 2.0 * wt.div[7]
mut.div[9] <- 1.8 * wt.div[9]

mut.sr <- wt.sr
mut.sr[5] <- 0.49
mut.sr[7] <- 0.49

mut.adj <- wt.adj
mut.adj[2,3] <- 0.7
mut.adj[2,9] <- 1 - mut.adj[2,3]

mut.death <- wt.death
mut.death[5] <- 0.1
mut.death[7] <- 0.1

#' Save this information for reproducibility
#save(mut.ic, mut.div, mut.sr, mut.adj, mut.death, file = "zz-data/MUT-params-P328-alt.RData")

plotEverything(
  n = n,
  wt.ic = wt.ic, mut.ic = mut.ic,
  wt.div = wt.div, mut.div = mut.div,
  wt.sr = wt.sr, mut.sr = mut.sr,
  wt.adj = wt.adj, mut.adj = mut.adj,
  wt.death = wt.death, mut.death = mut.death,
  compartment.names = name.vec
)
#savePlot(filename = "Figures/P328-alt.pdf", width = 200, height.mult = 0.7)
plot.data <- createTable(
  wt.ic = wt.ic, mut.ic = mut.ic,
  wt.div = wt.div, mut.div = mut.div,
  wt.sr = wt.sr, mut.sr = mut.sr,
  wt.adj = wt.adj, mut.adj = mut.adj,
  wt.death = wt.death, mut.death = mut.death,
  compartment.names = name.vec
)
#save(plot.data, file = "zz-data/P328-alt.RData")

# 3. Platelet and red cell biased pattern (P339) ----

#' Mutant cell parameters
mut.div <- wt.div
mut.div[1] <- 1.8 * wt.div[1]
mut.div[2] <- 1.8 * wt.div[2]
mut.div[3] <- 1.8 * wt.div[3]
mut.div[4] <- 1.8 * wt.div[4]
mut.div[5] <- 1.8 * wt.div[5]
mut.div[7] <- 2.0 * wt.div[7]
mut.div[9] <- 1.8 * wt.div[9]

mut.sr <- wt.sr
mut.sr[5] <- 0.49
mut.sr[7] <- 0.46

mut.adj <- wt.adj
mut.adj[3,4] <- 1.4 * wt.adj[3,4]
mut.adj[3,7] <- 1 - mut.adj[3,4]

mut.death <- wt.death
mut.death[5] <- 0.1
mut.death[7] <- 0.1

#' Save this information for reproducibility
#save(mut.ic, mut.div, mut.sr, mut.adj, mut.death, file = "zz-data/MUT-params-P339.RData")

plotEverything(
  n = n,
  wt.ic = wt.ic, mut.ic = mut.ic,
  wt.div = wt.div, mut.div = mut.div,
  wt.sr = wt.sr, mut.sr = mut.sr,
  wt.adj = wt.adj, mut.adj = mut.adj,
  wt.death = wt.death, mut.death = mut.death,
  compartment.names = name.vec
)
#savePlot(filename = "Figures/P339.pdf", width = 200, height.mult = 0.7)
plot.data <- createTable(
  wt.ic = wt.ic, mut.ic = mut.ic,
  wt.div = wt.div, mut.div = mut.div,
  wt.sr = wt.sr, mut.sr = mut.sr,
  wt.adj = wt.adj, mut.adj = mut.adj,
  wt.death = wt.death, mut.death = mut.death,
  compartment.names = name.vec
)
#save(plot.data, file = "zz-data/P339.RData")

# 4. Save the plot data in a single excel file ----
plot.data.list <- sapply(c("P433", "P328", "P339"), function(name) {
  load(paste0("zz-Data/", name, ".RData"))
  return(plot.data)
}, simplify = F)
write.xlsx(plot.data.list, file = "zz-data/Patient-specific.xlsx")
