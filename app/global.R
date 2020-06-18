#' global.R
#' Author: Peter Ashcroft, ETH Zurich

# Load default WT parameter values ----
WT <- new.env()
load("../zz-data/WT-params.RData", envir = WT)
#' Global defintions
name.vec <- WT$name.vec
n <- WT$n
adj <- WT$adj

default.wt <- list(
  ic = WT$wt.ic, #ic = 0.01, # This is the mutant frequency here
  div = WT$wt.div,
  sr = WT$wt.sr,
  C = c(WT$wt.adj[[2,9]], WT$wt.adj[[3,7]]),
  death = WT$wt.death
)

# Default mutant parameters (P339; platelet and red cell biased) ----
MUT <- new.env()
load("../zz-data/MUT-params-P339.RData", envir = MUT)

default.mut <- list(
  ic = MUT$mut.ic,
  div = MUT$mut.div,
  sr = MUT$mut.sr,
  C = c(MUT$mut.adj[[2,9]], MUT$mut.adj[[3,7]]),
  death = MUT$mut.death
)

# default.mut <- list(
#   ic = 0.01,
#   div = c(
#     1.8 * default.wt$div[1],
#     1.8 * default.wt$div[2],
#     1.8 * default.wt$div[3],
#     1.8 * default.wt$div[4],
#     1.8 * default.wt$div[5],
#     default.wt$div[6],
#     2.0 * default.wt$div[7],
#     default.wt$div[8],
#     1.8 * default.wt$div[9],
#     default.wt$div[10]
#   ),
#   sr = c(
#     default.wt$sr[1],
#     default.wt$sr[2],
#     default.wt$sr[3],
#     default.wt$sr[4],
#     0.49,
#     default.wt$sr[6],
#     0.46,
#     default.wt$sr[8],
#     default.wt$sr[9],
#     default.wt$sr[10]
#   ),
#   C = c(default.wt$C[1], 0.7),
#   death = c(
#     default.wt$death[1],
#     default.wt$death[2],
#     default.wt$death[3],
#     default.wt$death[4],
#     0.1,
#     default.wt$death[6],
#     0.1,
#     default.wt$death[8],
#     default.wt$death[9],
#     default.wt$death[10]
#   )
# )

# Mutant parameters with only late amplification ----
late.mut <- list(
  ic = c(0.01, rep.int(0, n - 1)),
  div = c(1.0, 1.4, 1.8, 4.3, 2.6, 2.0, 4.4, 1.0, 4.7, 2.0),
  sr = c(0.5, 0.3, 0.3, 0.49, 0.3, 0.0, 0.46, 0.0, 0.41, 0.0),
  C = c(0.5, 0.5),
  death = c(0.0, 0.2, 0.2, 0.0, 0.2, 0.2, 0.0, 0.2, 0.0, 0.2)
)

# Suggestion from Radek ----
RS.mut <- list(
  ic =  c(0.01, rep.int(0, n - 1)),
  div = c(2.0, 1.3, 1.5, 7.4, 1.0, 1.0, 10.0, 1.0, 5.0, 1.0),
  sr = c(0.5, 0.3, 0.3, 0.47, 0.0, 0.0, 0.46, 0.0, 0.3, 0.0),
  C = c(0.50, 0.75),
  death = c(0.0, 0.2, 0.2, 0.14, 0.01, 0.2, 0.04, 0.2, 0.2, 0.2)
)

# P339 (platelet and red cell biased) ----
MUT1 <- new.env()
load("../zz-data/MUT-params-P339.RData", envir = MUT1)

mut.P339 <- list(
  ic = MUT1$mut.ic,
  div = MUT1$mut.div,
  sr = MUT1$mut.sr,
  C = c(MUT1$mut.adj[[2,9]], MUT1$mut.adj[[3,7]]),
  death = MUT1$mut.death
)

# P433 (platelet biased) ----
MUT2 <- new.env()
load("../zz-data/MUT-params-P433.RData", envir = MUT2)

mut.P433 <- list(
  ic = MUT2$mut.ic,
  div = MUT2$mut.div,
  sr = MUT2$mut.sr,
  C = c(MUT2$mut.adj[[2,9]], MUT2$mut.adj[[3,7]]),
  death = MUT2$mut.death
)

# P433 (platelet biased) alternative ----
MUT2.alt <- new.env()
load("../zz-data/MUT-params-P433-alt.RData", envir = MUT2.alt)

mut.P433.alt <- list(
  ic = MUT2.alt$mut.ic,
  div = MUT2.alt$mut.div,
  sr = MUT2.alt$mut.sr,
  C = c(MUT2.alt$mut.adj[[2,9]], MUT2.alt$mut.adj[[3,7]]),
  death = MUT2.alt$mut.death
)


# P328 (red cell biased) ----
MUT3 <- new.env()
load("../zz-data/MUT-params-P328.RData", envir = MUT3)

mut.P328 <- list(
  ic = MUT3$mut.ic,
  div = MUT3$mut.div,
  sr = MUT3$mut.sr,
  C = c(MUT3$mut.adj[[2,9]], MUT3$mut.adj[[3,7]]),
  death = MUT3$mut.death
)

# P328 (red cell biased) alternative ----
MUT3.alt <- new.env()
load("../zz-data/MUT-params-P328-alt.RData", envir = MUT3.alt)

mut.P328.alt <- list(
  ic = MUT3.alt$mut.ic,
  div = MUT3.alt$mut.div,
  sr = MUT3.alt$mut.sr,
  C = c(MUT3.alt$mut.adj[[2,9]], MUT3.alt$mut.adj[[3,7]]),
  death = MUT3.alt$mut.death
)


