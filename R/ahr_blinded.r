#' Blinded estimation of average hazard ratio
#'
#' Based on blinded data and assumed hazard ratios in different intervals, compute
#' a blinded estimate of average hazard ratio (AHR) and corresponding estimate of statistical information.
#' This function is intended for use in computing futility bounds based on spending assuming
#' the input hazard ratio (hr) values for intervals specified here.
#'
#' @param Srv input survival object (see \code{Surv}); note that only 0=censored, 1=event for \code{Surv}
#' @param intervals Vector containing positive values indicating interval lengths where the
#' exponential rates are assumed.
#' Note that a final infinite interval is added if any events occur after the final interval
#' specified.
#' @param hr vector of hazard ratios assumed for each interval
#' @param ratio ratio of experimental to control randomization.
#'
#' @return A \code{tibble} with one row containing
#' `AHR` blinded average hazard ratio based on assumed period-specific hazard ratios input in `failRates`
#' and observed events in the corresponding intervals
#' `Events` total observed number of events, `info` statistical information based on Schoenfeld approximation,
#' and info0 (information under related null hypothesis) for each value of `totalDuration` input;
#' if `simple=FALSE`, `Stratum` and `t` (beginning of each constant HR period) are also returned
#' and `HR` is returned instead of `AHR`
#'
#' @examples
#' library(simtrial)
#' library(survival)
#' ahr_blinded(Srv = Surv(time = simtrial::Ex2delayedEffect$month,
#'                        event = simtrial::Ex2delayedEffect$evntd),
#'             intervals = c(4,100),
#'             hr = c(1, .55),
#'             ratio = 1)
#'
#' @export
ahr_blinded <- function (Srv = Surv(time = simtrial::Ex1delayedEffect$month,
                                    event = simtrial::Ex1delayedEffect$evntd),
                         intervals = array(3, 3),
                         hr = c(1,.6),
                         ratio = 1)
{   msg <- "hr must be a vector of positive numbers"
    if (!is.vector(hr, mode="numeric")) stop(msg)
    if (min(hr) <= 0) stop(msg)

    events <- simtrial::pwexpfit(Srv, intervals)[,3]
    nhr <- length(hr)
    nx <- length(events)
    # Add to hr if length shorter than intervals
    if (length(hr) < length(events)) hr <- c(hr, rep(hr[nhr], nx - nhr))

    # Compute blinded AHR
    theta <- sum(log(hr[1:nx]) * events) / sum(events)

    # Compute adjustment for information
    Qe <- ratio / (1 + ratio)

    return(tibble::tibble(Events = sum(events), AHR=exp(theta), theta = theta, info0 = sum(events) * (1-Qe) * Qe))
}
