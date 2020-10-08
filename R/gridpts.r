#' @importFrom tibble tibble
NULL
#' Grid points for group sequential design numerical integration
#'
#' Points and weights for Simpson's rule numerical integration from
#' p 349 - 350 of Jennison and Turnbull book.
#' This is not used for arbitrary integration, but for the canonical form of Jennison and Turnbull.
#' mu is computed elsewhere as drift parameter times sqrt of information.
#' Since this is a lower-level routine, no checking of input is done; calling routines should
#' ensure that input is correct.
#' Lower limit of integration can be \code{-Inf} and upper limit of integration can be \code{Inf}
#'
#' @details
#' Jennison and Turnbull (p 350) claim accuracy of \code{10E-6} with \code{r=16}.
#' The numerical integration grid spreads out at the tail to enable accurate tail probability calcuations.
#'
#'
#' @param r Integer, at least 2; default of 18 recommended by Jennison and Turnbull
#' @param mu Mean of normal distribution (scalar) under consideration
#' @param a lower limit of integration (scalar)
#' @param b upper limit of integration (scalar \code{> a})
#'
#' @return A \code{tibble} with grid points in \code{z} and numerical integration weights in \code{w}
#' @export
#'
#' @examples
#' # approximate variance of standard normal (i.e., 1)
#' gridpts() %>% summarise(var = sum(z^2 * w * dnorm(z)))
#'
#' # approximate probability above .95 quantile (i.e., .05)
#' gridpts(a = qnorm(.95), b = Inf) %>% summarise(p05 = sum(w * dnorm(z)))
gridpts <- function(r = 18, mu = 0, a = -Inf, b = Inf){
  # Define odd numbered grid points for real line
  x <- c(mu - 3 - 4 * log(r / (1:(r - 1))),
         mu - 3 + 3 * (0:(4 * r)) / 2 / r,
         mu + 3 + 4 * log(r / (r - 1):1)
         )
  # Trim points outside of [a, b] and include those points
  if (min(x) < a) x <- c(a, x[x > a])
  if (max(x) > b) x <- c(x[x < b], b)

  # Define even numbered grid points between the odd ones
  m <- length(x)
  y <- (x[2:m] + x[1:(m-1)]) / 2

  # Compute weights for odd numbered grid points
  i <- 2:(m-1)
  wodd <- c(x[2] - x[1],
            (x[i + 1] - x[i - 1]),
            x[m] - x[m - 1]) / 6

  weven <- 4 * (x[2:m] - x[1:(m-1)]) / 6

  # Now combine odd- and even-numbered grid points with their
  # corresponding weights
  z <- rep(0, 2*m - 1)
  z[2 * (1:m) - 1] <- x
  z[2 * (1:(m-1))] <- y
  w <- z
  w[2 * (1:m) - 1] <- wodd
  w[2 * (1:(m-1))] <- weven

  return(tibble::tibble(z=z, w=w))
}
