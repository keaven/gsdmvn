#' @importFrom tibble tibble
#' @importFrom stats qnorm
NULL
#' Group sequential bound computation with non-constant effect
#'
#' \code{gs_power_npe()} derives group sequential bounds and boundary crossing probabilities for a design.
#' It allows a non-constant treatment effect over time, but also can be applied for the usual homogeneous effect size designs.
#' It requires treatment effect and statistical information at each analysis as well as a method of deriving bounds, such as spending.
#' The routine enables two things not available in the gsDesign package: 1) non-constant effect, 2) more flexibility in boundary selection.
#' For many applications, the non-proportional-hazards design function \code{gs_design_nph()} will be used; it calls this function.
#' Initial bound types supported are 1) spending bounds, 2) fixed bounds, and 3) Haybittle-Peto-like bounds.
#' The requirement is to have a boundary update method that can each bound without knowledge of future bounds.
#' As an example, bounds based on conditional power that require knowledge of all future bounds are not supported by this routine;
#' a more limited conditional power method will be demonstrated.
#' Boundary family designs Wang-Tsiatis designs including the original (non-spending-function-based) O'Brien-Fleming and Pocock designs
#' are not supported by \code{gs_power_npe()}.
#' @param theta natural parameter for group sequential design representing
#' expected incremental drift at all analyses; used for power calculation
#' @param theta1 natural parameter for alternate hypothesis, if needed for lower bound computation
#' @param info statistical information at all analyses for input \code{theta}
#' @param info0 statistical information under null hypothesis, if different than \code{info};
#' impacts null hypothesis bound calculation
#' @param info1 statistical information under hypothesis used for futility bound calculation if different from
#' \code{info}; impacts futility hypothesis bound calculation
#' @param binding indicator of whether futility bound is binding; default of FALSE is recommended
#' @param upper function to compute upper bound
#' @param lower function to compare lower bound
#' @param upar parameter to pass to upper
#' @param lpar parameter to pass to lower
#' @param test_upper indicator of which analyses should include an upper (efficacy) bound;
#' single value of TRUE (default)  indicates all analyses; otherwise,
#' a logical vector of the same length as \code{info} should indicate which analyses will have an efficacy bound
#' @param test_lower indicator of which analyses should include a lower bound;
#' single value of TRUE (default) indicates all analyses;
#' single value FALSE indicated no lower bound; otherwise,
#' a logical vector of the same length as \code{info} should indicate which analyses will have a lower bound
#' @param r  Integer, at least 2; default of 18 recommended by Jennison and Turnbull
#' @param tol Tolerance parameter for boundary convergence (on Z-scale)
#' @author Keaven Anderson \email{keaven\_anderson@@merck.}
#'
#' @export
#'
#' @examples
#'
#' library(gsDesign)
#' # Default (single analysis; Type I error controlled)
#' gs_power_npe(theta=0) %>% filter(Bound=="Upper")
#'
#' # Fixed bound
#' gs_power_npe(theta = c(.1, .2, .3), info = (1:3) * 40, info0 = (1:3) * 40,
#'               upper = gs_b, upar = gsDesign::gsDesign(k=3,sfu=gsDesign::sfLDOF)$upper$bound,
#'               lower = gs_b, lpar = c(-1, 0, 0))
#'
#' # Same fixed efficacy bounds, no futility bound (i.e., non-binding bound), null hypothesis
#' gs_power_npe(theta = rep(0,3), info = (1:3) * 40,
#' upar = gsDesign::gsDesign(k=3,sfu=gsDesign::sfLDOF)$upper$bound,
#' lpar = rep(-Inf, 3)) %>% filter(Bound=="Upper")
#'
#' # Fixed bound with futility only at analysis 1; efficacy only at analyses 2, 3
#' gs_power_npe(theta = c(.1, .2, .3), info = (1:3) * 40, upar = c(Inf, 3, 2), lpar = c(qnorm(.1), -Inf, -Inf))
#'
#' # Spending function bounds
#' # Lower spending based on non-zero effect
#' gs_power_npe(theta = c(.1, .2, .3), info = (1:3) * 40,
#'              upper = gs_spending_bound,
#'              upar = list(sf = gsDesign::sfLDOF, total_spend = 0.025, param = NULL, timing = NULL),
#'              lower = gs_spending_bound,
#'              lpar = list(sf = gsDesign::sfHSD, total_spend = 0.1, param = -1, timing = NULL))
#'
#' # Same bounds, but power under different theta
#' gs_power_npe(theta = c(.15, .25, .35), theta1 = c(.1, .2, .3), info = (1:3) * 40,
#'              upper = gs_spending_bound,
#'              upar = list(sf = gsDesign::sfLDOF, total_spend = 0.025, param = NULL, timing = NULL),
#'              lower = gs_spending_bound,
#'              lpar = list(sf = gsDesign::sfHSD, total_spend = 0.1, param = -1, timing = NULL))
#'
#' # Two-sided symmetric spend, O'Brien-Fleming spending
#' # Typically, 2-sided bounds are binding
#' xx <- gs_power_npe(theta = rep(0, 3), theta1 = rep(0, 3), info = (1:3) * 40,
#'                    upper = gs_spending_bound,
#'                    binding = TRUE,
#'                    upar = list(sf = gsDesign::sfLDOF, total_spend = 0.025, param = NULL, timing = NULL),
#'                    lower = gs_spending_bound,
#'                    lpar = list(sf = gsDesign::sfLDOF, total_spend = 0.025, param = NULL, timing = NULL))
#' xx
#'
#' # Re-use these bounds under alternate hypothesis
#' # Always use binding = TRUE for power calculations
#' upar <- (xx %>% filter(Bound=="Upper"))$Z
#' gs_power_npe(theta = c(.1, .2, .3), info = (1:3) * 40,
#'              binding = TRUE,
#'              upar = upar,
#'              lpar = -upar)
#'
gs_power_npe <- function(theta = .1, theta1 = NULL, info = 1, info1 = NULL, info0 = NULL,
                         binding = FALSE,
                         upper=gs_b, lower=gs_b, upar = qnorm(.975), lpar= -Inf,
                         test_upper = TRUE, test_lower = TRUE,
                         r = 18, tol = 1e-6){
  #######################################################################################
  # WRITE INPUT CHECK TESTS AND RETURN APPROPRIATE ERROR MESSAGES
  # theta should be a scalar or vector of real values; if vector, same length as info
  # info should be a scalar or vector of positive increasing values
  # info0 should be NULL or of the same form as info
  # test_upper and test_lower should be logical scalar or vector; if vector same length as info
  # END INPUT CHECKS
  #######################################################################################
  # SET UP PARAMETERS
  K <- length(info)
  if (is.null(info0)) info0 <- info
  if (is.null(info1)) info1 <- info
  if (length(info1) != length(info) || length(info0) != length(info)) stop("gs_power_npe: length of info, info0, info1 must be the same")
  if (length(theta) == 1 && K > 1) theta <- rep(theta, K)
  if (is.null(theta1)){theta1 <- theta}else if (length(theta1)==1) theta1 <- rep(theta1,K)
  if (length(test_upper) == 1 && K > 1) test_upper <- rep(test_upper, K)
  if (length(test_lower) == 1 && K > 1) test_lower <- rep(test_lower, K)
  a <- rep(-Inf, K)
  b <- rep(Inf, K)
  hgm1_0 <- NULL
  hgm1_1 <- NULL
  hgm1 <- NULL
  upperProb <- rep(NA, K)
  lowerProb <- rep(NA, K)
  ######################################################################################
  # COMPUTE BOUNDS
  for(k in 1:K){
    # Lower bound update
    a[k] <- lower(k = k, par = lpar, hgm1 = hgm1_1, theta = theta1, info = info1, r = r, tol = tol, test_bound = test_lower,
                  efficacy = FALSE)
    # Upper bound update
    b[k] <- upper(k = k, par = upar, hgm1 = hgm1_0, info = info0, r = r, tol = tol, test_bound = test_upper)
    if(k==1){
      upperProb[1] <- if(b[1] < Inf) {pnorm(b[1], mean = sqrt(info[1]) * theta[1], lower.tail = FALSE)}else{0}
      lowerProb[1] <- if(a[1] > -Inf){pnorm(a[1], mean = sqrt(info[1]) * theta[1])}else{0}
      hgm1_0 <- h1(r = r, theta = 0,         I = info0[1], a = if(binding){a[1]}else{-Inf}, b = b[1])
      hgm1_1 <- h1(r = r, theta = theta1[1], I = info1[1], a = a[1], b = b[1])
      hgm1   <- h1(r = r, theta = theta[1],  I = info[1],  a = a[1], b = b[1])
    }else{
      # Cross upper bound
      upperProb[k] <- if(b[k]< Inf){
        hupdate(r = r, theta = theta[k], I = info[k], a = b[k], b = Inf,
                thetam1 = theta[k - 1], Im1 = info[k - 1], gm1 = hgm1) %>%
          summarise(sum(h)) %>% as.numeric()
      }else{0}
      # Cross lower bound
      lowerProb[k] <- if(a[k] > -Inf){
        hupdate(r = r, theta = theta[k], I = info[k], a = -Inf, b = a[k],
                thetam1 = theta[k - 1], Im1 = info[k - 1], gm1 = hgm1) %>%
          summarise(sum(h)) %>% as.numeric()
      }else{0}
      if(k < K){
        hgm1_0 <- hupdate(r = r, theta = 0,         I = info0[k], a = if(binding){a[k]}else{-Inf}, b = b[k],
                          thetam1 = 0,           Im1 = info0[k-1], gm1 = hgm1_0)
        hgm1_1 <- hupdate(r = r, theta = theta1[k], I = info1[k], a = a[k], b = b[k],
                          thetam1 = theta1[k-1], Im1 = info1[k-1], gm1 = hgm1_1)
        hgm1   <- hupdate(r = r, theta = theta[k],  I = info[k],  a = a[k], b = b[k],
                          thetam1 = theta[k-1],  Im1 = info[k-1],  gm1 = hgm1)
      }
    }
  }
  return(tibble::tibble(
    Analysis = rep(1:K, 2),
    Bound = c(rep("Upper", K), rep("Lower", K)),
    Z= c(b, a),
    Probability = c(cumsum(upperProb),
                    cumsum(lowerProb)),
    theta = rep(theta, 2),
    theta1 = rep(theta1, 2),
    info = rep(info, 2),
    info0 = rep(info0, 2),
    info1 = rep(info1, 2))
  )
}

