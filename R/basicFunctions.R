#' basicFunctions.R
#' Author: Peter Ashcroft, ETH Zurich

#' Some basic functions that I will repeatedly use,
#' so I will define them here

## Dependencies ----
library(reshape2)

## Function list ----
#' - DFapply()
#' - createNames()
#' - splitNames()
#' - orderVariables()

## Functions ----
#' Sometimes we want to use lapply to bind dataframes,
#' so I define the function \code{DFapply} to do this for me
#' @param X (vector or list): What is to be looped over.
#' E.g. the list of dataframes we want to bind
#' @param FUN (function): First argument must be the value of \code{X}.
#' This should return a dataframe
#' @param ...: Further arguments to function
DFapply <- function(X, FUN, ...) {
  df.list <- lapply(X, FUN, ...)
  df <- do.call(rbind, df.list)
  return(df)
}

#' Create parameter/variable names
createNames <- function(name, compartment.names) {
  paste(name, compartment.names, sep = ".")
}

#' Opposite of createNames.
#' Here we split names into cell type and compartment id.
splitNames <- function(names) {
  DFapply(
    names,
    function(name) {
      split.name <- unlist(strsplit(name, split = ".", fixed = TRUE))
      df <- data.frame(cell.type = split.name[1], id = split.name[2])
      return(df)
    }
  )
}

#' The order in which we want variables to be plotted
orderVariables <- function(reverse = FALSE) {
  vars <- paste(c("WT", "Mut"), rep(c(1:12), each = 2), sep = ".")
  if (reverse) vars <- rev(vars)
  return(vars)
}
