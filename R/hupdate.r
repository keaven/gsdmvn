#' @importFrom stats dnorm
#' @importFrom tibble tibble
NULL
#' Update numerical integration for group sequential design
#'
#' Update grid points for numerical integration from one analysis to the next
#'
#' @param r Integer, at least 2; default of 18 recommended by Jennison and Turnbull
#' @param theta Drift parameter for current analysis
#' @param I Information at current analysis
#' @param a lower limit of integration (scalar)
#' @param b upper limit of integration (scalar \code{> a})
#' @param thetam1  Drift parameter for previous analysis
#' @param Im1 Information at previous analysis
#' @param gm1 numerical integration grid from \code{h1()} or previous run of \code{hupdate()}
#'
#' @return A \code{tibble} with grid points in \code{z}, numerical integration weights in \code{w},
#' and a normal density with mean \code{mu = theta * sqrt{I}} and variance 1 times the weight in \code{w}.
#'
#' @examples
#' # 2nd analysis with no interim bound and drift 0 should have mean 0, variance 1
#' hupdate() %>% summarise(mu = sum(z * h), var = sum((z - mu)^2 * h))
#' @export
hupdate <- function(r = 18, theta = 0, I = 2, a = -Inf, b = Inf, thetam1 = 0, Im1 = 1, gm1 = h1()){
  # sqrt of change in information
  rtdelta <- sqrt(I - Im1)
  rtI <- sqrt(I)
  rtIm1 <- sqrt(Im1)
  g <- gridpts(r = r, mu = theta * rtI, a= a, b = b)
  # update integration
  mu <- theta * I - thetam1 * Im1
  h <- rep(0, length(g$z))
  for(i in seq_along(g$z)){
    x <- (g$z[i] * rtI - gm1$z * rtIm1 - mu ) / rtdelta
    h[i] <- sum(gm1$h * dnorm(x))
  }
  h <- h * g$w * rtI / rtdelta
  return(tibble::tibble(z = g$z, w = g$w, h = h))
}
