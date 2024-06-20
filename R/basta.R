#' @export
#'
stat.est <- function(x, order = 1) {
  invisible(tseries::garch(x, c(0, order), trace = FALSE)$coef)
}

#' BASTA Algorithm to Detect Multiple-Change-Points
#' 
#' If `x` is a piecewise stationary ARCH process, then `detection` returns the estimated
#' expectation of the function \eqn{g(.)} from the paper.
#' 
#' See the paper for additional details.
#'
#' @param x a numeric vector giving the time series to analyse
#' @param order order of the ARCH model
#' @param factor dampening factor `F` from the paper
#' @param c threshold constant `c` from the paper
#' @param epsilon epsilon from the paper
#' @param log logical, if true, then \eqn{U^{(4)}_t} is used; otherwise \eqn{U^{(3)}_t}
#' 
#' 
#' @return A list with the estimated expectation of the function \eqn{g(.)} and a vector of change-points.
#'
#' @export
#'
detection <- function(x, order = 1, factor = 8, c = 0.5, epsilon = 0, log = TRUE) {
  xx <- stat.resid(x, NULL, order, factor, epsilon)
  if (log) xx <- log(xx + epsilon)
  xx.buh <- best.unbal.haar(xx, max.inner.prod.p)
  n <- length(x)
  th <- c * n^(3 / 8)
  xx.buh.th <- bin.segm(xx.buh, th)
  xx.buh.th.rec <- reconstr(xx.buh.th)
  indicator <- sapply(
    unique(xx.buh.th.rec),
    function(x) min(which(xx.buh.th.rec == x)))
  if (length(indicator) < 2) stop("No change-points detected.")
  list(expected_g = xx.buh.th.rec, change_points = indicator[-1])
}

#' Create `xreg` For Structural Breaks Models
#' 
#' Create a matrix of dummy variables for use as external regressons.
#'
#' @param y dependent variable
#' @param demean logical, if `TRUE` then `y` will get demeaned prior to the BASTA algorithm
#' @param pulse logical, if `TRUE` then dummy variables equal one at detected change-point. If `FALSE` then dummy variables equal one at and after detected breakpoint.
#' @param ... additional arguments to pass to `detection`
#'
#' @return A matrix
#' @export
xreg_basta <- function(y, demean = FALSE, pulse = TRUE, ...) {
  y_ <- zoo::coredata(y)
  if (demean) y_ <- y_ - mean(y_)
  idx <- detection(y_, ...)$change_points

  res <- matrix(0L, nrow = length(y), ncol = length(idx))

  for (i in seq_along(idx)) {
    if (!pulse) {
      res[(idx[i]):(length(y)), i] <- 1L
    } else {
      res[idx[i], i] <- 1L
    }
  }

  res
}

###########################################################################
#
# The functions below first appeared in the R package "unbalhaar". However,
# these newer versions (below) must be used for the above functions to work
# correctly.
#
###########################################################################

best.unbal.haar <- function(x, criterion = max.inner.prod.p) {
  n <- length(x)
  if (n < 2) stop("Input vector too short")
  tree <- list(matrix(0, 5, 1))
  tree[[1]][1, 1] <- 1
  ipi <- inner.prod.iter(x)
  ind.max <- criterion(x)
  tree[[1]][3, 1] <- 1
  tree[[1]][4, 1] <- ind.max
  tree[[1]][5, 1] <- n
  tree[[1]][2, 1] <- ipi[ind.max]

  j <- 1

  while (sum(tree[[j]][5, ] - tree[[j]][3, ] - rep(1, dim(tree[[j]])[2]))) {
    tree <- c(tree, list(matrix(0, 5, 0)))
    no.parent.coeffs <- dim(tree[[j]])[2]
    no.child.coeffs <- 0
    for (i in 1:no.parent.coeffs) {
      if (tree[[j]][4, i] - tree[[j]][3, i] >= 1) {
        no.child.coeffs <- no.child.coeffs + 1
        tree[[j + 1]] <- matrix(c(tree[[j + 1]], matrix(0, 5, 1)), 5, no.child.coeffs)
        tree[[j + 1]][1, no.child.coeffs] <- 2 * tree[[j]][1, i] - 1
        ipi <- inner.prod.iter(x[(tree[[j]][3, i]):(tree[[j]][4, i])])
        ind.max <- criterion(x[(tree[[j]][3, i]):(tree[[j]][4, i])])
        tree[[j + 1]][2, no.child.coeffs] <- ipi[ind.max]
        tree[[j + 1]][3, no.child.coeffs] <- tree[[j]][3, i]
        tree[[j + 1]][5, no.child.coeffs] <- tree[[j]][4, i]
        tree[[j + 1]][4, no.child.coeffs] <- ind.max + tree[[j]][3, i] - 1
      }
      if (tree[[j]][5, i] - tree[[j]][4, i] >= 2) {
        no.child.coeffs <- no.child.coeffs + 1
        tree[[j + 1]] <- matrix(c(tree[[j + 1]], matrix(0, 5, 1)), 5, no.child.coeffs)
        tree[[j + 1]][1, no.child.coeffs] <- 2 * tree[[j]][1, i]
        ipi <- inner.prod.iter(x[(tree[[j]][4, i] + 1):(tree[[j]][5, i])])
        ind.max <- criterion(x[(tree[[j]][4, i] + 1):(tree[[j]][5, i])])
        tree[[j + 1]][2, no.child.coeffs] <- ipi[ind.max]
        tree[[j + 1]][3, no.child.coeffs] <- tree[[j]][4, i] + 1
        tree[[j + 1]][5, no.child.coeffs] <- tree[[j]][5, i]
        tree[[j + 1]][4, no.child.coeffs] <- ind.max + tree[[j]][4, i]
      }
    }
    j <- j + 1
  }

  smooth <- sum(x) / sqrt(n)

  z <- list(tree = tree, smooth = smooth)

  return(z)
}

#' @export
max.inner.prod.p <- function(x, p = 0.8) {
  ipi <- inner.prod.iter(x)
  n <- length(ipi)
  pip <- ipi[(1 + floor((1 - p) * n)):ceiling(p * n)]
  return(med(which(abs(pip) == max(abs(pip)))) + floor((1 - p) * n))
}

reconstr <- function(buh) {
  J <- length(buh$tree)
  n <- buh$tree[[1]][5, 1]
  rec <- rep(1 / sqrt(n) * buh$smooth, n)

  for (j in 1:J) {
    K <- dim(buh$tree[[j]])[2]
    for (k in 1:K) {
      rec[(buh$tree[[j]][3, k]):(buh$tree[[j]][5, k])] <-
        rec[(buh$tree[[j]][3, k]):(buh$tree[[j]][5, k])] + unbal.haar.vector(buh$tree[[j]][3:5, k]) * buh$tree[[j]][2, k]
    }
  }

  return(rec)
}

unbal.haar.vector <- function(a) {
  n <- a[3] - a[1] + 1
  m <- a[2] - a[1] + 1

  return(c(rep(sqrt(1 / m - 1 / n), m), rep(-1 / sqrt(n^2 / m - n), n - m)))
}
