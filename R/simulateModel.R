# simulateModel.R
# Author: Peter Ashcroft, ETH Zurich

## Dependencies ----
library(deSolve)
library(openxlsx)
#library(adaptivetau)

## Function list ----
#' - createAdjacency()
#' - odes()
#' - integrateODE()
#' - numericSteadyState()
#' - numericSteadyState2 () -- with amplification parameters
#' - getClonalFraction()
#' - createTable()

## Functions ----

#' Create an adjacency matrix of size n
#' No entries given creates a linear process
#' @param entries (list): each element is a vector c(i,j,\phi_{i,j})
createAdjacency <- function(n, entries = NULL) {
  # Create blank matrix
  adj <- matrix(0, nrow = n, ncol = n)
  
  # If no entries given, define a linear process
  if (is.null(entries)) {
    for (i in c(1:(n - 1))) adj[i, i + 1] <- 1.0
  }
  else {
    for (e in entries) adj[e[1], e[2]] <- e[3]
  }
  return(adj)
}


#' ODEs for WT and mutant compartments sizes
#' 
#' These will then be passed to the function \code{ode()}, part of the \code{deSolve} package.
#' 
#' @param t (numeric): Time argument
#' @param state (vector): State of population of dimension 2n. First half is WT, second half is Mut.
#' @param parm (list): A list of model parameters
#' Within the list parms, we have the following elements:
#' @param adj (matrix): The adjaceny matrix, whose elements are the branching
#' fractions from compartment \code{i} to compartment \code{j}.
#' The rows should sum to one.
#' @param divRate (vector): Rate at which cells divide in each compartment.
#' @param srProb (vector): Probability of self renewal.
#' Again dimensions are 2n, with first half WT rates, and second half Mut.
#' 
#' @return (list) List of infinitesimal changes to the size of each compartment
#' 
odes <- function(t, state, parms) {
  n <- length(state)/2
  wt.range <- c(1:n)
  mut.range <- c((n + 1):(2*n))
  
  #wt.influx <- 2 * ( (1 - parms$srProb[wt.range]) * parms$divRate[wt.range]) * state[wt.range]
  #mut.influx <- 2 * ( (1 - parms$srProb[mut.range]) * parms$divRate[mut.range]) * state[mut.range]
  
  wt.influx <- 2 * parms$divRate[wt.range] * (1 - parms$deathProb[wt.range]) * (1 - parms$srProb[wt.range]) * state[wt.range]
  mut.influx <- 2 * parms$divRate[mut.range] * (1 - parms$deathProb[mut.range]) * (1 - parms$srProb[mut.range]) * state[mut.range]
  
  # d.wt <- sapply(wt.range, function(i){
  #     sum( wt.influx * parms$wt.adj[,i] ) - (1 - 2 * parms$srProb[[i]]) * parms$divRate[[i]] * state[[i]]
  # })
  # 
  # d.mut <- sapply(mut.range, function(i){
  #     sum( mut.influx * parms$mut.adj[,(i-n)] ) - (1 - 2 * parms$srProb[[i]]) * parms$divRate[[i]] * state[[i]]
  # })
  d.wt <- sapply(wt.range, function(i){
    sum( wt.influx * parms$wt.adj[,i] ) - parms$divRate[[i]] * (1 - 2 * parms$srProb[[i]] * (1 - parms$deathProb[[i]])) * state[[i]]
  })
  
  d.mut <- sapply(mut.range, function(i){
    sum( mut.influx * parms$mut.adj[,(i - n)] ) - parms$divRate[[i]] * (1 - 2 * parms$srProb[[i]] * (1 - parms$deathProb[[i]])) * state[[i]]
  })
  return(list( c(d.wt, d.mut) ))
}


#' Integrate the combined ODEs for these parameters
integrateODE <- function(times, wt.ic, mut.ic, wt.div, mut.div, wt.sr, mut.sr, wt.adj, mut.adj, wt.death, mut.death) {
  # Build initial condition
  state <- c(wt.ic, mut.ic)
  # Build list of parameters
  param.list <- list(divRate = c(wt.div, mut.div), srProb = c(wt.sr, mut.sr), wt.adj = wt.adj, mut.adj = mut.adj, deathProb = c(wt.death, mut.death))
  # Solve
  sol <- ode(y = state, times = times, func = odes, parms = param.list)
  sol <- as.data.frame(sol)
  # Assign proper names to columns
  names(sol) <- c("time", createNames("WT", seq_along(wt.ic) ), createNames("Mut", seq_along(mut.ic) ) )
  return(sol)
}

#' Calculate the steady state size as the leading eigenvalue of each set of ODEs
numericSteadyState <- function(wt.ic, mut.ic, wt.div, mut.div, wt.sr, mut.sr, wt.adj, mut.adj, wt.death, mut.death) {
  # Number of compartments
  n <- length(wt.ic)
  
  ## WT steady state
  # Blank matrix
  mat <- matrix(0, nrow = n, ncol = n)
  # Populate matrix to describe the linear equations of the WT cells
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      mat[i,j] <- ifelse(i == j, (2 * wt.sr[[i]] * (1 - wt.death[[i]]) - 1) * wt.div[[i]], 0) + 2 * (1 - wt.sr[[j]]) * wt.div[[j]] * (1 - wt.death[[j]]) * wt.adj[j,i]
    }
  }
  # Find the eigenvalues and vectors
  eig <- eigen(mat)
  # Identify the zero which value
  i <- which(abs(eig$values) == min(abs(eig$values)))
  # Get corresponding eigenvector
  v <- eig$vectors[,i]
  # Scale by size of first compartment
  wt.ss <- v * (wt.ic[1] / v[1])
  # Add names
  names(wt.ss) <- createNames("WT", seq_len(n))
  
  
  ## Mutant steady state
  # Blank matrix
  mat <- matrix(0, nrow = n, ncol = n)
  # Populate matrix to describe the linear equations of the WT cells
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      mat[i,j] <- ifelse(i == j, (2 * mut.sr[[i]] * (1 - mut.death[[i]]) - 1) * mut.div[[i]], 0) + 2 * (1 - mut.sr[[j]]) * mut.div[[j]] * (1 - mut.death[[j]]) * mut.adj[j,i]
    }
  }
  # Find the eigenvalues and vectors
  eig <- eigen(mat)
  # Identify the zero which value
  i <- which(abs(eig$values) == min(abs(eig$values)))
  # Get corresponding eigenvector
  v <- eig$vectors[,i]
  # Scale by size of first compartment
  mut.ss <- v * (mut.ic[1] / v[1])
  # Add names
  names(mut.ss) <- createNames("Mut", seq_len(n))
  
  # Combine results
  ss <- c(wt.ss, mut.ss)
  return(ss)
}

#' Calculate the steady state size as the leading eigenvalue of each set of ODEs
numericSteadyState2 <- function(wt.ic, mut.ic, wt.div, mut.div, wt.sr, mut.sr, wt.adj, mut.adj, wt.amp, mut.amp) {
  # Number of compartments
  n <- length(wt.ic)
  
  ## WT steady state
  # Blank matrix
  mat <- matrix(0, nrow = n, ncol = n)
  # Populate matrix to describe the linear equations of the WT cells
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      mat[i,j] <- ifelse(i == j, (2 * wt.sr[[i]] - 1) * wt.div[[i]], 0) + 2 * wt.amp[[j]] * (1 - wt.sr[[j]]) * wt.div[[j]] * wt.adj[j,i]
    }
  }
  # Find the eigenvalues and vectors
  eig <- eigen(mat)
  # Identify the zero which value
  i <- which(abs(eig$values) == min(abs(eig$values)))
  # Get corresponding eigenvector
  v <- eig$vectors[,i]
  # Scale by size of first compartment
  wt.ss <- v * (wt.ic[1] / v[1])
  # Add names
  names(wt.ss) <- createNames("WT", seq_len(n))
  
  
  ## Mutant steady state
  # Blank matrix
  mat <- matrix(0, nrow = n, ncol = n)
  # Populate matrix to describe the linear equations of the WT cells
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      mat[i,j] <- ifelse(i == j, (2 * mut.sr[[i]] - 1) * mut.div[[i]], 0) + 2 * mut.amp[[j]]*(1 - mut.sr[[j]]) * mut.div[[j]] * mut.adj[j,i]
    }
  }
  # Find the eigenvalues and vectors
  eig <- eigen(mat)
  # Identify the zero which value
  i <- which(abs(eig$values) == min(abs(eig$values)))
  # Get corresponding eigenvector
  v <- eig$vectors[,i]
  # Scale by size of first compartment
  mut.ss <- v * (mut.ic[1] / v[1])
  # Add names
  names(mut.ss) <- createNames("Mut", seq_len(n))
  
  # Combine results
  ss <- c(wt.ss, mut.ss)
  return(ss)
}

#' Calculate the fraction of mutants per compartment
getClonalFraction <- function(state.vec) {
  id.vec <- unique(splitNames(names(state.vec))$id)
  DFapply(id.vec, function(id) {
    fraction <- state.vec[[paste("Mut", id, sep = ".")]] / (state.vec[[paste("WT", id, sep = ".")]] + state.vec[[paste("Mut", id, sep = ".")]])
    df <- data.frame(id = id, fraction = fraction)
    return(df)
  })
}

#' Create a table of variables and clonal fraction for each compartment
createTable <- function(wt.ic, mut.ic, wt.div, mut.div, wt.sr, mut.sr, wt.adj, mut.adj, wt.death, mut.death, compartment.names = NULL) {
  #' Build compartment names
  if (is.null(compartment.names)) compartment.names <- seq_along(wt.ic)
  
  #' Parameter differences
  getDiff <- function(wt, mut) {
    ifelse(wt == 0, 0, (mut - wt) / wt)
  }
  delta.sr <- getDiff(wt.sr, mut.sr)
  delta.div <- getDiff(wt.div, mut.div)
  delta.death <- getDiff(wt.death, mut.death)
  # delta.sr <- (mut.sr - wt.sr) / wt.sr
  # delta.div <- (mut.div - wt.div) / wt.div
  # delta.death <- (mut.death - wt.death) / wt.death
  
  #' Identify changes in branching rates
  wt.branch <- sapply(seq_len(nrow(wt.adj)), function(row.id) {
    #col.index <- which(wt.adj[row.id, ] > 0)[1]
    col.indices <- which(wt.adj[row.id, ] > 0)
    if (length(col.indices) == 2) {
      col.index <- col.indices[2]
      return(wt.adj[row.id, col.index])
    }
    else return(0)
  })
  mut.branch <- sapply(seq_len(nrow(mut.adj)), function(row.id) {
    #col.index <- which(wt.adj[row.id, ] > 0)[1]
    col.indices <- which(mut.adj[row.id, ] > 0)
    if (length(col.indices) == 2) {
      col.index <- col.indices[2]
      return(mut.adj[row.id, col.index])
    }
    else return(0)
  })
  delta.branch <- getDiff(wt.branch, mut.branch)
  #delta.branch <- (mut.branch - wt.branch) / wt.branch
  
  #' Compute the steady state and clonal fraction
  steadyState <- numericSteadyState(
    wt.ic = wt.ic, mut.ic = mut.ic,
    wt.div = wt.div, mut.div = mut.div,
    wt.sr = wt.sr, mut.sr = mut.sr,
    wt.adj = wt.adj, mut.adj = mut.adj,
    wt.death = wt.death, mut.death = mut.death
  )
  clonal.fraction <- getClonalFraction(steadyState)
  
  data.frame(
    compartment = compartment.names,
    a_WT = wt.sr, a_MUT = mut.sr, Delta.a = delta.sr,
    b_WT = wt.div, b_MUT = mut.div, Delta.b = delta.div,
    c_WT = wt.branch, c_MUT = mut.branch, Delta.c = delta.branch,
    d_WT = wt.death, d_MUT = mut.death, Delta.d = delta.death,
    x = steadyState[c(1:length(wt.ic))],
    y = steadyState[c(1:length(mut.ic)) + length(mut.ic)],
    z = clonal.fraction$fraction,
    row.names = NULL)
}
