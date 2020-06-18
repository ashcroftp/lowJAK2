#' S9.R
#' Author: Peter Ashcroft, ETH Zurich

#' Six compartment model of general tissue in linear and branched structures

#' Dependencies:
source("R/basicFunctions.R")
source("R/simulateModel.R")
source("R/plotFunctions.R")

# 1. Linear population structure ----
#' How do changes in the mutant's birth rate, death probability,
#' or self-renewal probability affect the clonal burden.

#' WT parameters for this section
n <- 6
adj <- createAdjacency(n)
#' Initial conditions
wt.ic <- c(0.9, rep.int(0, n - 1))
mut.ic <- c(0.1, rep.int(0, n - 1))
#' Wild-type cell parameters
wt.div <- c(seq(from = 1.0, to = 3.0, length.out = n - 1), 1.0)
wt.sr <- c(0.5, rep.int(0.4, n - 2), 0.0)
wt.adj <- adj
wt.death <- c(0, rep.int(0.2, n - 1))

#' Function to plot the tree, parameters, and clonal fraction as a single plot
plotEverything <- function(n, wt.ic, mut.ic, wt.div, mut.div, wt.sr, mut.sr, wt.adj, mut.adj, wt.death, mut.death, show.var) {
  #' Steady state and clonal fraction
  steadyState <- numericSteadyState(
    wt.ic = wt.ic, mut.ic = mut.ic,
    wt.div = wt.div, mut.div = mut.div,
    wt.sr = wt.sr, mut.sr = mut.sr,
    wt.adj = wt.adj, mut.adj = mut.adj,
    wt.death = wt.death, mut.death = mut.death
  )
  clonal.fraction <- getClonalFraction(steadyState)
  #' Plot the data
  plotFractionParams(
    fraction.df = clonal.fraction,
    wt.div = wt.div, mut.div = mut.div,
    wt.sr = wt.sr, mut.sr = mut.sr,
    wt.adj = wt.adj, mut.adj = mut.adj,
    wt.death = wt.death, mut.death = mut.death,
    show.var = show.var,
    colour.vec = rep.int("#FA2D1C", n), alpha.vec = seq.int(0.2, 1, length.out = n),
    show.tree = T, y.pos = rep.int(1, n), tree.height = 0.2, y.lim = 0.5
  )
}


## 1.1. Birth rate ----
### 1.1.1. Faster dividing mutant cells in a non-HSC compartment ----
#' This results in a lower clonal fraction in the modified compartment.
#' All other fractions are preserved.

#' Mutant cell parameters
mut.div <- wt.div
mut.div[4] <- 2 * wt.div[4]
mut.sr <- wt.sr
mut.adj <- wt.adj
mut.death <- wt.death
#' Calculate steady state and plot
plotEverything(
  n = n,
  wt.ic = wt.ic, mut.ic = mut.ic,
  wt.div = wt.div, mut.div = mut.div,
  wt.sr = wt.sr, mut.sr = mut.sr,
  wt.adj = wt.adj, mut.adj = mut.adj,
  wt.death = wt.death, mut.death = mut.death,
  show.var = "div"
)
#savePlot(filename = "Figures/S9B-left.pdf", width = 60, height.mult = 1)
plot.data <- createTable(
  wt.ic = wt.ic, mut.ic = mut.ic,
  wt.div = wt.div, mut.div = mut.div,
  wt.sr = wt.sr, mut.sr = mut.sr,
  wt.adj = wt.adj, mut.adj = mut.adj,
  wt.death = wt.death, mut.death = mut.death
)
#save(plot.data, file = "zz-data/S9B-left.RData")

### 1.1.2. Faster dividing cells from HSC through to penultimate compartment ----
#' Mutant cell parameters
mut.div <- 2 * wt.div
mut.div[6] <- wt.div[6]
mut.sr <- wt.sr
mut.adj <- wt.adj
mut.death <- wt.death
#' Calculate steady state and plot
plotEverything(
  n = n,
  wt.ic = wt.ic, mut.ic = mut.ic,
  wt.div = wt.div, mut.div = mut.div,
  wt.sr = wt.sr, mut.sr = mut.sr,
  wt.adj = wt.adj, mut.adj = mut.adj,
  wt.death = wt.death, mut.death = mut.death,
  show.var = "div"
)
#savePlot(filename="Figures/S9B-right.pdf", width = 60, height.mult = 1)
plot.data <- createTable(
  wt.ic = wt.ic, mut.ic = mut.ic,
  wt.div = wt.div, mut.div = mut.div,
  wt.sr = wt.sr, mut.sr = mut.sr,
  wt.adj = wt.adj, mut.adj = mut.adj,
  wt.death = wt.death, mut.death = mut.death
)
#save(plot.data, file = "zz-data/S9B-right.RData")


## 1.2. Death rate ----
### 1.2.1. Less-death-prone cells early ----
#' All other fractions are preserved.
#' Mutant cell parameters
mut.div <- wt.div
mut.sr <- wt.sr
mut.adj <- wt.adj
mut.death <- wt.death
mut.death[2] <- 0.05
#' Calculate steady state and plot
plotEverything(
  n = n,
  wt.ic = wt.ic, mut.ic = mut.ic,
  wt.div = wt.div, mut.div = mut.div,
  wt.sr = wt.sr, mut.sr = mut.sr,
  wt.adj = wt.adj, mut.adj = mut.adj,
  wt.death = wt.death, mut.death = mut.death,
  show.var = "death"
)
#savePlot(filename="Figures/S9C-left.pdf", width = 60, height.mult = 1)
plot.data <- createTable(
  wt.ic = wt.ic, mut.ic = mut.ic,
  wt.div = wt.div, mut.div = mut.div,
  wt.sr = wt.sr, mut.sr = mut.sr,
  wt.adj = wt.adj, mut.adj = mut.adj,
  wt.death = wt.death, mut.death = mut.death
)
#save(plot.data, file = "zz-data/S9C-left.RData")

### 1.2.2. Less-death-prone cells late ----
#' All other fractions are preserved.
#' Mutant cell parameters
mut.div <- wt.div
mut.sr <- wt.sr
mut.adj <- wt.adj
mut.death <- wt.death
mut.death[n - 1] <- 0.05
#' Calculate steady state and plot
plotEverything(
  n = n,
  wt.ic = wt.ic, mut.ic = mut.ic,
  wt.div = wt.div, mut.div = mut.div,
  wt.sr = wt.sr, mut.sr = mut.sr,
  wt.adj = wt.adj, mut.adj = mut.adj,
  wt.death = wt.death, mut.death = mut.death,
  show.var = "death"
)
#savePlot(filename="Figures/S9C-right.pdf", width = 60, height.mult = 1)
plot.data <- createTable(
  wt.ic = wt.ic, mut.ic = mut.ic,
  wt.div = wt.div, mut.div = mut.div,
  wt.sr = wt.sr, mut.sr = mut.sr,
  wt.adj = wt.adj, mut.adj = mut.adj,
  wt.death = wt.death, mut.death = mut.death
)
#save(plot.data, file = "zz-data/S9C-right.RData")

## 1.3. Self-renewal ----
### 1.3.1. Increased self-renewal in a single  ----
#' Mutant cell parameters
mut.div <- wt.div
mut.sr <- wt.sr
mut.sr[3] <- 0.49
mut.adj <- wt.adj
mut.death <- wt.death
#' Calculate steady state and plot
plotEverything(
  n = n,
  wt.ic = wt.ic, mut.ic = mut.ic,
  wt.div = wt.div, mut.div = mut.div,
  wt.sr = wt.sr, mut.sr = mut.sr,
  wt.adj = wt.adj, mut.adj = mut.adj,
  wt.death = wt.death, mut.death = mut.death,
  show.var = "sr"
)
#savePlot(filename="Figures/S9D-left.pdf", width = 60, height.mult = 1)
plot.data <- createTable(
  wt.ic = wt.ic, mut.ic = mut.ic,
  wt.div = wt.div, mut.div = mut.div,
  wt.sr = wt.sr, mut.sr = mut.sr,
  wt.adj = wt.adj, mut.adj = mut.adj,
  wt.death = wt.death, mut.death = mut.death
)
#save(plot.data, file = "zz-data/S9D-left.RData")

### 1.3.2 Greater self-renewal in an amplifying compartment ----
#' Mutant cell parameters
mut.div <- wt.div
mut.sr <- wt.sr
mut.sr[n - 1] <- 0.49
mut.adj <- wt.adj
mut.death <- wt.death
#' Calculate steady state and plot
plotEverything(
  n = n,
  wt.ic = wt.ic, mut.ic = mut.ic,
  wt.div = wt.div, mut.div = mut.div,
  wt.sr = wt.sr, mut.sr = mut.sr,
  wt.adj = wt.adj, mut.adj = mut.adj,
  wt.death = wt.death, mut.death = mut.death,
  show.var = "sr"
)
#savePlot(filename="Figures/S9D-right.pdf", width = 60, height.mult = 1)
plot.data <- createTable(
  wt.ic = wt.ic, mut.ic = mut.ic,
  wt.div = wt.div, mut.div = mut.div,
  wt.sr = wt.sr, mut.sr = mut.sr,
  wt.adj = wt.adj, mut.adj = mut.adj,
  wt.death = wt.death, mut.death = mut.death
)
#save(plot.data, file = "zz-data/S9D-right.RData")

# 2. Branched population structure ----
#' How do changes in the mutant's differentiation bias affect the clonal burden?

#' Parameter defintions for this section
n <- 6
adj <- createAdjacency(
  n, entries = list(
    c(1,2,1),
    c(2,3,0.5), c(2,5,0.5),
    c(3,4,1),
    c(5,6,1)
  )
)
#' Initial conditions
wt.ic <- c(0.9, rep.int(0, n - 1))
mut.ic <- c(0.1, rep.int(0, n - 1))
#' Wild-type cell parameters
wt.div <- seq(from = 1.0, to = 3, length.out = n)
wt.div[c(4,6)] <- wt.div[1]
wt.sr <- c(0.5, rep.int(0.4, n - 2), 0.0)
wt.adj <- adj
wt.death <- c(0, rep.int(0.1, n - 1))
#' Mutant cell parameters
#' Here we just modify the branching fraction
mut.div <- wt.div
mut.sr <- wt.sr
mut.adj <- wt.adj
mut.adj[2,3] <- 1.5 * wt.adj[2,3]
mut.adj[2,5] <- 1 - mut.adj[2,3]
mut.death <- wt.death
#' Steady state and clonal fraction
steadyState <- numericSteadyState(
  wt.ic = wt.ic, mut.ic = mut.ic,
  wt.div = wt.div, mut.div = mut.div,
  wt.sr = wt.sr, mut.sr = mut.sr,
  wt.adj = wt.adj, mut.adj = mut.adj,
  wt.death = wt.death, mut.death = mut.death
)
clonal.fraction <- getClonalFraction(steadyState)
#' Plot the data
plotFractionParams(
  fraction.df = clonal.fraction,
  wt.div = wt.div, mut.div = mut.div,
  wt.sr = wt.sr, mut.sr = mut.sr,
  wt.adj = wt.adj, mut.adj = mut.adj,
  wt.death = wt.death, mut.death = mut.death,
  show.var = c("branch"),
  colour.vec = ifelse(seq_len(n) <= 4, "#FDB12B", "#FA2D1C"), alpha.vec = c(0.4, 0.55, 0.7, 1.0, 0.7, 1.0),
  show.tree = T, y.pos = c(1,1,0,0,1,1), tree.height = 0.25, y.lim = -0.5
)
#savePlot(filename="Figures/S9F.pdf", width = 60, height.mult = 1.3)
plot.data <- createTable(
  wt.ic = wt.ic, mut.ic = mut.ic,
  wt.div = wt.div, mut.div = mut.div,
  wt.sr = wt.sr, mut.sr = mut.sr,
  wt.adj = wt.adj, mut.adj = mut.adj,
  wt.death = wt.death, mut.death = mut.death
)
#save(plot.data, file = "zz-data/S9F.RData")

# 3. Save the plot data in a single excel file ----
plot.data.list <- sapply(c("S9B-left", "S9B-right", "S9C-left", "S9C-right", "S9D-left", "S9D-right", "S9F"), function(name) {
  load(paste0("zz-Data/", name, ".RData"))
  return(plot.data)
}, simplify = F)
write.xlsx(plot.data.list, file = "zz-data/S9.xlsx")

