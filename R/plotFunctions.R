#' plotFunctions.R
#' Author: Peter Ashcroft, ETH Zurich

#' Dependencies ----
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(scales)

#' Function list ----
#' - sciFormat()
#' - colourFunction()
#' - colourGradient()
#' - savePlot()
#' - plotODE()
#' - plotSS()
#' - plotFraction()
#' - plotTree()
#' - plotFractionParams()
#' - plotFractionParams2()
#' - plotFractionParamsAbsolute()

## Functions ----
#' Function for scientific notation on y-axis
#'
#' \code{sciFormat()}
#' @param l (string): Number to be converted to expoenential format
#' @return (string)
#' @examples sciFormat("1e-1")
sciFormat <- function(l) {
  l <- format(l, scientific = TRUE)
  ## Just plot 10^x
  l <- gsub("e","10^", gsub("^(.*)e", "e", l))
  # Return this as an expression
  parse(text = l)
}

#' Returns a vector of named colours as light/dark pairs for each compartment i
colourFunction <- function() {
  # Use the "paired" palette from ColorBrewer
  colours <- brewer.pal(12, "Paired")
  colours <- rep_len(colours, length.out = 24)
  # Add appropriate names
  names(colours) <- orderVariables()
  return(colours)
}

#'
colourGradient <- function(compartments, palette = "PuRd") {
  # Define all colours, select a nice distinguishable range, and then
  # sample the number we need
  colour.vec <- colorRampPalette(brewer.pal(9, palette)[c(3:9)])(length(compartments))
  # Add names
  names(colour.vec) <- compartments

  return(colour.vec)
}

#' My saving function to store figures with nice dimensions
#'
#' \code{savePlot}
#'
#' @param filename (string): name to be assigned to the stored figure.
#' @param plot (gg-object): the figure to be saved.
#' @param width (numeric): The width of the fiugure in millimeters.
#' @param height.mult (numeric): Ratio of height to width. Defaults to golden ratio.
#' @param ... : Options for ggsave()
#' @return NULL
#' @examples savePlot("./foo.pdf", plot, width=20, height.mult=2)
savePlot <- function(filename = "~/Desktop/foo.pdf", plot = last_plot(), width, height.mult = 2 / (1 + sqrt(5)), ...) {
  ggsave(filename = filename, plot = plot, width = width, height = width * height.mult, units = "mm", ...)
}

#' Plot the compartment sizes from ODE output
#'
plotODE <- function(data) {
  # Melt the data frame
  data.melt <- melt(data, id.vars = "time")
  # Arrange the values for pairwise comparison
  data.melt$variable <- factor(data.melt$variable, levels = orderVariables(reverse = TRUE))
  # Plot
  plot <- ggplot(data.melt, aes(x = time, y = value, colour = variable, fill = variable)) +
    geom_area(position = "stack") +
    scale_color_manual(values = colourFunction()) +
    scale_fill_manual(values = colourFunction()) +
    theme_bw()
  return(plot)
}

#' Plot the steady state compartment sizes
#'
#' @param SS (vector): A named numeric vector of the steady-state compartment sizes.
plotSS <- function(SS) {
  # Create a dataframe to feed to ggplot
  ss.data <- data.frame(variable = names(SS), value = SS, row.names = NULL)
  # Add column for grouping
  #ss.data$compartment <- splitNames(as.character(ss.data$variable))$id
  ss.data <- cbind(ss.data, splitNames(as.character(ss.data$variable)))
  # Arrange the values for pairwise comparison
  ss.data$variable <- factor(ss.data$variable, levels = orderVariables())
  # Plot
  plot <- ggplot(ss.data, aes(x = id, y = value, colour = variable, fill = variable, group = cell.type)) +
    geom_col(position = "dodge") +
    scale_color_manual(values = colourFunction()) +
    scale_fill_manual(values = colourFunction()) +
    theme_bw()
  return(plot)
}

#' Plot the clonal fraction per compartment
plotFraction <- function(fraction.df) {
  # Add mutant name colum to get the right colours
  fraction.df$mut <- paste("Mut", fraction.df$id, sep = ".")
  # Plot
  plot <- ggplot(fraction.df, aes(x = id, y = fraction, colour = mut, fill = mut)) +
    geom_col() +
    scale_color_manual(values = colourFunction()) +
    scale_fill_manual(values = colourFunction()) +
    theme_bw() +
    theme(legend.position = "None")
  return(plot)
}

#' Plot the population structure
plotTree <- function(adj, y.pos, show.diff = FALSE, colour.vec = NULL, alpha.vec = NULL, compartment.names = NULL, lab.offset = 0, y.lim = -1.5) {
  #' Default colour scheme
  if (is.null(colour.vec)) colour.vec <- colourGradient(seq_len(nrow(adj)))
  if (is.null(alpha.vec)) alpha.vec <- rep.int(1, nrow(adj))

  #' Construct the nodes dataframe
  nodes <- data.frame(
    id = factor(seq_len(nrow(adj))),
    x = factor(seq_len(nrow(adj))),
    y = y.pos,
    colour = colour.vec,
    alpha = alpha.vec,
    variable = factor("tree")
  )

  #' Construct the edges
  pos <- which(adj != 0, arr.ind = TRUE)
  edge.df <- data.frame(from = pos[, "row"], to = pos[, "col"], weight = adj[pos])
  edges <- data.frame(
    edge = factor(rep(seq_len(nrow(edge.df)), each = 2)),
    x = unlist(lapply(seq_len(nrow(edge.df)), function(i) c(nodes[nodes$id == edge.df[i, "from"], "x"], nodes[nodes$id == edge.df[i, "to"], "x"]))),
    y = unlist(lapply(seq_len(nrow(edge.df)), function(i) c(nodes[nodes$id == edge.df[i, "from"], "y"], nodes[nodes$id == edge.df[i, "to"], "y"]))),
    weight = unlist(lapply(seq_len(nrow(edge.df)), function(i) rep.int(edge.df[i, "weight"], 2)))
  )
  edges$weight <- sapply(edges$weight, function(w) ifelse(w == 1.0, 0.1, 3*w))

  #' Now build the tree
  tree <- ggplot(nodes, aes(x = x, y = y)) +
    # geom_line(data = edges, aes(group = edge, size = weight), show.legend = F) +
    # scale_size_continuous(range = c(0.1,1)) +
    geom_line(data = edges, aes(group = edge), show.legend = F) +
    geom_point(colour = "white", size = 6, show.legend = F) +
    geom_point(aes(colour = colour, alpha = alpha), size = 6, show.legend = F) +
    scale_colour_identity() +
    scale_alpha_identity() +
    facet_grid(variable ~ .) +
    #theme_void() +
    theme_cowplot() +
    theme(
      strip.background = element_rect(fill = "transparent"),
      strip.text.y = element_text(colour = "transparent"),
      axis.line = element_blank()
    ) +
    scale_x_discrete(name = "", breaks = NULL) +
    scale_y_continuous(name = "", breaks = NULL, limits = c(y.lim, 1.5))

  #' Add labels if necessary
  if (!is.null(compartment.names)) {
    dat <- nodes
    dat$name <- compartment.names[nodes$id]
    dat$y <- dat$y - 0.2
    tree <- tree +
      geom_text(data = dat, aes(label = name), color = "black", parse = F, angle = 0, hjust = "center", nudge_y = lab.offset)
  }
  else {
    tree <- tree + geom_text(aes(label = id), color = "white", size = 3)
  }

  #
  # if(show.diff) {
  #     diff <- (mut.adj[2,10]-wt.adj[2,10])/wt.adj[2,10]
  #     tree <- tree +
  #         geom_text(
  #             x=2.5,
  #             y=0.5,
  #             label=as.character(as.expression(substitute(Delta*a["2,10"] == d, list(d=round(diff, 2))))),
  #             parse=T,
  #             colour=ifelse(diff>0, "#0072B2", ifelse(diff<0, "#D55E00", "grey40"))
  #         )
  #
  #     diff <- (mut.adj[3,4]-wt.adj[3,4])/wt.adj[3,4]
  #     tree <- tree +
  #         geom_text(
  #             x=4.4,
  #             y=-0.5,
  #             label=as.character(as.expression(substitute(Delta*a["3,4"] == d, list(d=round(diff, 2))))),
  #             parse=T,
  #             colour=ifelse(diff>0, "#0072B2", ifelse(diff<0, "#D55E00", "grey40"))
  #         )
  # }
  return(tree)
}

#' Plot parameter changes, clonal fractions, and population structure on a single plot
plotFractionParams <- function(fraction.df, wt.div, mut.div, wt.sr, mut.sr, wt.adj, mut.adj, wt.death, mut.death, show.var = c("sr", "div", "branch", "death"), colour.vec = NULL, alpha.vec = NULL, show.tree = FALSE, y.pos = NULL, tree.height = 0.25, compartment.names = NULL, ...) {
  #' Default colour scheme
  if (is.null(colour.vec)) colour.vec <- colourGradient(seq_len(nrow(wt.adj)))
  if (is.null(names(colour.vec))) names(colour.vec) <- seq_len(nrow(wt.adj))
  if (is.null(alpha.vec)) alpha.vec <- rep.int(1, nrow(wt.adj))
  if (is.null(names(alpha.vec))) names(alpha.vec) <- seq_len(nrow(wt.adj))

  #' Load dataframe which we will extend
  plot.df <- fraction.df
  #' Find changes in division rates and self-renewal probabilities
  plot.df$div <- (mut.div - wt.div) / wt.div
  plot.df$sr <- (mut.sr - wt.sr) / wt.sr
  #' Identify changes in branching rates
  plot.df$branch <- sapply(seq_len(nrow(wt.adj)), function(row.id) {
    #' Next two lines use the original direction of \Delta c
    col.index <- which(wt.adj[row.id, ] > 0)[1]
    (mut.adj[row.id, col.index] - wt.adj[row.id, col.index]) / wt.adj[row.id, col.index]
    #' Next 6 lines flip the direction of \Delta c
    # col.indices <- which(wt.adj[row.id, ] > 0)
    # if (length(col.indices) == 2) {
    #   col.index <- col.indices[2]
    #   return((mut.adj[row.id, col.index] - wt.adj[row.id, col.index]) / wt.adj[row.id, col.index])
    # }
    # else return(0)
  })
  plot.df$death <- (mut.death - wt.death) / wt.death
  #' Melt the df for facet plotting
  melt.df <- melt(plot.df, id.vars = "id")
  melt.df$variable <- factor(melt.df$variable, levels = c("sr", "div", "branch", "death", "fraction"))
  #' Make sure we have no NaNs or NAs
  melt.df$value[!is.finite(melt.df$value)] <- 0

  #' Subset the data depending on what parameter changes we want to show
  melt.df <- subset(melt.df, subset = variable %in% c(show.var, "fraction"))
  #' Generate parseable names, and refactor the variables with these names as levels
  #var.names <- c(sr = "Delta*alpha", div = "Delta*beta", branch = "Delta*a")
  var.names <- c(sr = "Delta*a", div = "Delta*b", branch = "Delta*c", death = "Delta*d")
  melt.df$variable <- factor(melt.df$variable, levels = c(show.var, "fraction"))
  levels(melt.df$variable) <- c(var.names[show.var], "fraction")

  #' Colour as a function of id and value
  melt.df$colour <- ifelse(melt.df$value == 0, "grey40", colour.vec[melt.df$id])
  #' Alpha as a function of id and value
  melt.df$alpha <- ifelse(melt.df$value == 0, 1, alpha.vec[melt.df$id])
  #' Add a small perturbation to zero values so the bar can be seen
  melt.df$value[melt.df$value == 0] <- 0.02

  #' Plot the data
  plot <- ggplot(melt.df, aes(x = id, y = value)) +
    geom_col(aes(fill = colour, alpha = alpha)) +
    scale_fill_identity() +
    scale_alpha_identity() +
    facet_grid(variable ~ ., scales = "free_y", labeller = label_parsed) +
    geom_hline(data = data.frame(variable = var.names[show.var], y = 0), aes(yintercept = y), size = 0.5, linetype = "dotted") +
    scale_y_continuous(name = "", breaks = pretty_breaks(n = 2), expand = c(0, 0)) +
    geom_point(data = data.frame(variable = rep(var.names[show.var], each = 2), id = factor(1), value = rep.int(c(-1, 1), length(show.var))), alpha = 0) +
    geom_point(data = data.frame(variable = "fraction", id = factor(1), value = 0.2), alpha = 0) +
    #theme_classic() +
    theme_cowplot() +
    theme(panel.spacing = unit(1, "lines"))
  #' Do we need to add compartment names?
  if (is.null(compartment.names)) plot <- plot + scale_x_discrete(name = "compartment")
  else plot <- plot + scale_x_discrete(name = "compartment", labels = compartment.names)  #labels = parse(text = compartment.names))

  if (show.tree) {
    plot <- plot_grid(
      plotTree(adj = wt.adj, y.pos = y.pos, show.diff = F, colour.vec = colour.vec, alpha.vec = alpha.vec, compartment.names = compartment.names, ...) +
        theme(plot.margin = unit(c(0, 0, -0.2, 0), "npc")),
      plot +
        theme(plot.margin = unit(c(-0.01, 0.01, 0.01, 0.01), "npc")),
      ncol = 1, align = "v", axis = "l",
      rel_heights = c(tree.height, 1)
    )
  }

  return(plot)
}

#' Plot parameter changes, clonal fractions, and population structure on a single plot
plotFractionParams2 <- function(fraction.df, wt.div, mut.div, wt.sr, mut.sr, wt.adj, mut.adj, wt.amp, mut.amp, show.var = c("sr", "div", "amp"), colour.vec = NULL, alpha.vec = NULL, show.tree=FALSE, y.pos=NULL, tree.height = 0.25, compartment.names=NULL, ...) {
  #' Default colour scheme
  if (is.null(colour.vec)) colour.vec <- colourGradient(seq_len(nrow(wt.adj)))
  if (is.null(names(colour.vec))) names(colour.vec) <- seq_len(nrow(wt.adj))
  if (is.null(alpha.vec)) alpha.vec <- rep.int(1, nrow(wt.adj))
  if (is.null(names(alpha.vec))) names(alpha.vec) <- seq_len(nrow(wt.adj))

  #' Load dataframe which we will extend
  plot.df <- fraction.df
  #' Find changes in division rates and self-renewal probabilities
  plot.df$div <- (mut.div - wt.div) / wt.div
  plot.df$sr <- (mut.sr - wt.sr) / wt.sr
  plot.df$amp <- (mut.amp - wt.amp) / wt.amp
  #' Identify changes in branching rates
  plot.df$branch <- sapply(seq_len(nrow(wt.adj)), function(row.id) {
    col.index <- which(wt.adj[row.id, ] > 0)[1]
    (mut.adj[row.id, col.index] - wt.adj[row.id, col.index]) / wt.adj[row.id, col.index]
  })
  #' Melt the df for facet plotting
  melt.df <- melt(plot.df, id.vars = "id")
  melt.df$variable <- factor(melt.df$variable, levels = c("sr", "div", "branch", "amp", "fraction"))
  #' Make sure we have no NaNs or NAs
  melt.df$value[!is.finite(melt.df$value)] <- 0

  #' Subset the data depending on what parameter changes we want to show
  melt.df <- subset(melt.df, subset = variable %in% c(show.var, "fraction"))
  #' Generate parseable names, and refactor the variables with these names as levels
  #var.names <- c(sr = "Delta*alpha", div = "Delta*beta", branch = "Delta*a")
  var.names <- c(sr = "Delta*a", div = "Delta*b", branch = "Delta*c", amp = "Delta*h")
  melt.df$variable <- factor(melt.df$variable, levels = c(show.var, "fraction"))
  levels(melt.df$variable) <- c(var.names[show.var], "fraction")

  #' Colour as a function of id and value
  melt.df$colour <- ifelse(melt.df$value == 0, "grey40", colour.vec[melt.df$id])
  #' Alpha as a function of id and value
  melt.df$alpha <- ifelse(melt.df$value == 0, 1, alpha.vec[melt.df$id])
  #' Add a small perturbation to zero values so the bar can be seen
  melt.df$value[melt.df$value == 0] <- 0.02

  #' Plot the data
  plot <- ggplot(melt.df, aes(x = id, y = value)) +
    geom_col(aes(fill = colour, alpha = alpha)) +
    scale_fill_identity() +
    scale_alpha_identity() +
    facet_grid(variable ~ ., scales = "free_y", labeller = label_parsed) +
    geom_hline(data = data.frame(variable = var.names[show.var], y = 0), aes(yintercept = y), size = 0.5, linetype = "dotted") +
    scale_y_continuous(name = "", breaks = pretty_breaks(n = 2), expand = c(0, 0)) +
    geom_point(data = data.frame(variable = rep(var.names[show.var], each = 2), id = factor(1), value = rep.int(c(-1, 1), length(show.var))), alpha = 0) +
    geom_point(data = data.frame(variable = "fraction", id = factor(1), value = 0.2), alpha = 0) +
    theme(panel.spacing = unit(1, "lines"))
  #' Do we need to add compartment names?
  if (is.null(compartment.names)) plot <- plot + scale_x_discrete(name = "compartment")
  else plot <- plot + scale_x_discrete(name = "compartment", labels = compartment.names)  #labels = parse(text = compartment.names))

  if (show.tree) {
    plot <- plot_grid(
      plotTree(adj = wt.adj, y.pos = y.pos, show.diff = F, colour.vec = colour.vec, alpha.vec = alpha.vec, compartment.names = compartment.names, ...) +
        theme(plot.margin = unit(c(0, 0, -0.2, 0), "npc")),
      plot +
        theme(plot.margin = unit(c(-0.01, 0.01, 0.01, 0.01), "npc")),
      ncol = 1, align = "v", axis = "l",
      rel_heights = c(tree.height, 1)
    )
  }

  return(plot)
}

#' Plot parameter changes, clonal fractions, and population structure on a single plot, using absolute parameter values
plotFractionParamsAbsolute <- function(fraction.df, wt.div, mut.div, wt.sr, mut.sr, wt.adj, mut.adj, wt.death, mut.death, show.var = c("sr", "div", "branch", "death"), colour.vec = NULL, alpha.vec = NULL, show.tree = FALSE, y.pos = NULL, tree.height = 0.25, compartment.names = NULL, ...) {
  #' Default colour scheme
  if (is.null(colour.vec)) colour.vec <- colourGradient(seq_len(nrow(wt.adj)))
  if (is.null(names(colour.vec))) names(colour.vec) <- seq_len(nrow(wt.adj))
  if (is.null(alpha.vec)) alpha.vec <- rep.int(1, nrow(wt.adj))
  if (is.null(names(alpha.vec))) names(alpha.vec) <- seq_len(nrow(wt.adj))

  #' Dataframe of clonal fraction for plotting
  #' Colour as a function of id and value
  fraction.df$colour <- ifelse(fraction.df$fraction == 0, "grey40", colour.vec[fraction.df$id])
  #' Alpha as a function of id and value
  fraction.df$alpha <- ifelse(fraction.df$fraction == 0, 1, alpha.vec[fraction.df$id])
  #' Add a small perturbation to zero values so the bar can be seen
  fraction.df$fraction[fraction.df$fraction == 0] <- 0.02

  par.list = list(wt.div = wt.div, mut.div = mut.div,
                  wt.sr = wt.sr, mut.sr = mut.sr,
                  wt.death = wt.death, mut.death = mut.death)
  var.names <- c(sr = "a", div = "b", branch = "c", death = "d")
  #' Dataframe of parameters for plotting
  param.plots <- lapply(show.var, function(param) {
    if (param == "branch") {
      param.df <- rbind(
        data.frame(id = fraction.df$id, param = param, type = "wt", value = 0*fraction.df$fraction),
        data.frame(id = fraction.df$id, param = param, type = "mut", value = 0*fraction.df$fraction)
      )
      #' Identify changes in branching rates
      for (i in seq_len(nrow(wt.adj))) {
        j <- 1
        while (j < nrow(wt.adj)) {
          if (wt.adj[i,j] > 0 & wt.adj[i,j] < 1) {
            param.df$value[param.df$id == i & param.df$type == "wt"] <- wt.adj[i,j]
            param.df$value[param.df$id == i & param.df$type == "mut"] <- mut.adj[i,j]
            j <- nrow(wt.adj)
          }
          j <- j + 1
        }
      }
      rm(i,j)
      #' Colour and alpha value
      param.df$colour <- colour.vec[param.df$id]
      param.df$alpha <- alpha.vec[param.df$id]
    }
    else {
      param.df <- rbind(
        data.frame(id = fraction.df$id, param = param, type = "wt", value = par.list[[paste0("wt.", param)]]),
        data.frame(id = fraction.df$id, param = param, type = "mut", value = par.list[[paste0("mut.", param)]])
      )
      #' Colour as a function of id and value
      param.df$colour <- colour.vec[param.df$id] #ifelse(param.df$value == 0, "grey40", colour.vec[param.df$id])
      #' Alpha as a function of id and value
      param.df$alpha <- alpha.vec[param.df$id]#ifelse(param.df$value == 0, 1, alpha.vec[param.df$id])
      #' Add a small perturbation to zero values so the bar can be seen
      #param.df$value[param.df$value == 0] <- 0.02
    }
    #' Plot
    param.plot <- ggplot(param.df, aes(x = type, y = value)) +
      geom_line(aes(group = id, alpha = alpha), colour = "gray") +
      geom_point(aes(colour = colour, alpha = alpha), show.legend = F) +
      facet_grid(~ id) +
      scale_colour_identity() +
      scale_alpha_identity() +
      scale_y_continuous(name = parse(text = var.names[[param]]), limits = c(0, NA), breaks = pretty_breaks(n = 2)) +#, expand = c(0, 0)) +
      theme_cowplot() +
      theme(
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing = unit(0, "lines")
      )

    return(param.plot)
  })

  #' Plot the clonal fraction
  plot <- ggplot(fraction.df, aes(x = id, y = fraction)) +
    geom_col(aes(fill = colour, alpha = alpha)) +
    scale_fill_identity() +
    scale_alpha_identity() +
    scale_y_continuous(name = "fraction", limits = c(0, NA), breaks = pretty_breaks(n = 2)) +
    geom_point(data = data.frame(id = factor(1), fraction = 0.2), colour = "transparent", alpha = 0) +
    theme_cowplot()

  #' Do we need to add compartment names?
  if (is.null(compartment.names)) plot <- plot + scale_x_discrete(name = "compartment")
  else plot <- plot + scale_x_discrete(name = "compartment", labels = compartment.names)  #labels = parse(text = compartment.names))

  plot.list <- c(param.plots, list(plot))

  plot.heights <- c(rep.int(1, length(show.var)), 1.4)

  if (show.tree) {
    tree.plot <- plotTree(adj = wt.adj, y.pos = y.pos, show.diff = F, colour.vec = colour.vec, alpha.vec = alpha.vec, compartment.names = compartment.names, ...)
    tree.plot <- tree.plot + theme(plot.margin = unit(c(0, 0, -0.2, 0), "npc"))
    plot.list <- c(list(tree.plot), plot.list)
    plot.heights <- c(tree.height, plot.heights)
  }

  plot <- plot_grid(plotlist = plot.list, ncol = 1, align = "v", axis = "l", rel_heights = plot.heights)

  return(plot)
}
