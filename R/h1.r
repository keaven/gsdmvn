#' @importFrom stats dnorm pnorm
#' @importFrom tibble tibble
NULL
#' Initialize numerical integration for group sequential design
#'
#' Compute grid points for first interim analysis in a group sequential design
#'
#' @param r Integer, at least 2; default of 18 recommended by Jennison and Turnbull
#' @param theta Drift parameter for first analysis
#' @param I Information at first analysis
#' @param a lower limit of integration (scalar)
#' @param b upper limit of integration (scalar \code{> a})
#'
#' @details Mean for standard normal distribution under consideration is \code{mu = theta * sqrt(I)}
#' @return A \code{tibble} with grid points in \code{z}, numerical integration weights in \code{w},
#' and a normal density with mean \code{mu = theta * sqrt{I}} and variance 1 times the weight in \code{w}.
#' @export
#'
#' @examples
#' # Replicate variance of 1, mean of 35
#' h1(theta = 5, I = 49) %>% summarise(mu = sum(z * h), var = sum((z - mu)^2 * h))
#'
#' # Replicate p-value of .0001 by numerical integration of tail
#' h1(a = qnorm(.9999)) %>% summarise(p = sum(h))
h1 <- function(r = 18, theta = 0, I = 1, a = -Inf, b = Inf){
  # fix for binding message
  z <- w <- h <- NULL
  # compute drift at analysis 1
  mu <- theta * sqrt(I);
  g <- gridpts(r, mu, a, b)
  # compute deviation from drift
  x <- g$z - mu
  # compute standard normal density, multiply by grid weight and return
  # values needed for numerical integration
  return(tibble::tibble(z = g$z, w = g$w, h = g$w * dnorm(x)))
}
