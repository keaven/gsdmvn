#' @importFrom dplyr summarize
#' @importFrom gsDesign gsDesign sfLDOF
#' @importFrom stats qnorm
NULL
#' Derive spending bound for group sequential boundary
#'
#' Computes one bound at a time based on spending under given distributional assumptions.
#' While user specifies \code{gs_spending_bound()} for use with other functions,
#' it is really not intended for use on its own.
#' Most important user specifications are made through a list provided to functions using \code{gs_spending_bound()}.
#' Function uses numerical integration and Newton-Raphson iteration to derive an individual bound for a group sequential
#' design that satisfies a targeted boundary crossing probability.
#' Algorithm has been modified from Chapter 19 of Jennison and Turnbull (2000).
#'
#' @param k analysis for which bound is to be computed
#' @param par a list with the following items: sf (class spending function),
#' timing (a vector containing values at which spending function is to be evaluated),
#' total_spend (total spend),
#' param (any parameters needed by the spending function)
#' @param hgm1 subdensity grid from h1 (k=2) or hupdate (k>2) for analysis k-1; if k=1, this is not used and may be NULL
#' @param theta natural parameter used for lower bound only spending;
#' represents average drift at each time of analysis at least up to analysis k;
#' upper bound spending is always set under null hypothesis (theta = 0)
#' @param info statistical information at all analyses, at least up to analysis k
#' @param efficacy TRUE (default) for efficacy bound, FALSE otherwise
#' @param test_bound a logical vector of the same length as \code{info} should indicate which analyses will have a bound
#' @param r  Integer, at least 2; default of 18 recommended by Jennison and Turnbull
#' @param tol Tolerance parameter for convergence (on Z-scale)
#'
#' @return returns a numeric bound (possibly infinite) or, upon failure, generates an error message.
#' @author Keaven Anderson \email{keaven\_anderson@@merck.}
#' @references Jennison C and Turnbull BW (2000), \emph{Group Sequential
#' Methods with Applications to Clinical Trials}. Boca Raton: Chapman and Hall.
#' @export
gs_spending_bound <- function(k = 1,
                              par = list(sf = gsDesign::sfLDOF,
                                              total_spend = 0.025,
                                              param = NULL,
                                              timing = NULL),
                              hgm1 = NULL,
                              theta = .1,
                              info = 1:3,
                              efficacy = TRUE,
                              test_bound = TRUE,
                              r = 18,
                              tol = 1e-6){
  if (is.null(par$timing)){timing <- info / (max(info))}else{timing <- par$timing}
  spend <- par$sf(alpha = par$total_spend, t = timing, param = par$param)$spend
  old_spend <- 0
  if (length(test_bound) == 1 && k > 1) test_bound <- rep(test_bound, k)
  for(i in 1:k){
    if (test_bound[i]){
      xx <- spend[i] - old_spend
      old_spend <- spend[i]
      spend[i] <- xx
    }else spend[i] <- 0
  }
  # Now just get spending for current bound
  spend <- spend[k]
  # lower bound
  if (!efficacy){
    if (length(theta) == 1) theta <- rep(theta, length(info))
      if (spend <= 0) return(-Inf) # If no spending, return -Inf for bound
    if (k == 1) return(qnorm(spend) + sqrt(info[1]) * theta[1]) # No need for iteration for first interim
    # Starting value
    a <- qnorm(spend) + sqrt(info[k]) * theta[k]
    # iterate to convergence; if not convergence after 10 iterations, return error
    for(iter in 0:10)
    {  # Get grid for rejection region
      hg <- hupdate(theta = theta[k], I =  info[k], a = -Inf, b = a, thetam1 = theta[k-1], Im1 = info[k-1], gm1 = hgm1, r = r)
      i <- nrow(hg)
      pik <- hg %>% summarise(sum(h)) %>% as.numeric() # pik is for lower bound crossing
      # Do first step in Newton-Raphson update to IA2 futility bound
      a_old <- a
      a <- a_old + (spend - pik) / hg$h[i] * hg$w[i] # analog to npebackground equation (8) for upper bound
      if (abs(a_old - a) < tol) return(a)
    }
    stop(paste("bound_update did not converge for lower bound calculation, analysis", k))
  }else{
    # upper bound
    if(spend <= 0) return(Inf)
    # Starting value
    b <- qnorm(spend, lower.tail = FALSE)
    if(k == 1) return(b) # No iteration needed for first bound
    for(iter in 0:10){
      # subdensity for final analysis in rejection region
      hg <- hupdate(theta = 0, I =  info[k], a = b, b = Inf, thetam1 = 0, Im1 = info[k-1], gm1 = hgm1)
      pik <- as.numeric(hg %>% summarise(sum(h))) # Probility of crossing final bound
      dpikdb <- hg$h[1] / hg$w[1] # Derivative of final bound crossing at b[k]
      b_old <- b
      b <- b - (spend - pik) / dpikdb # Newton-Raphson update
      if (abs(b - b_old) < tol) return(b)
    }
    stop(paste("bound_update did not converge for upper bound calculation, analysis", k))
  }
}
