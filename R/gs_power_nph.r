#' @importFrom tibble tibble
#' @importFrom gsDesign gsDesign sfLDOF
#' @importFrom stats qnorm
#' @importFrom dplyr select left_join
NULL
#' Group sequential design power under non-proportional hazards
#'
#' @param enrollRates enrollment rates
#' @param failRates failure and dropout rates
#' @param ratio Experimental:Control randomization ratio (not yet implemented)
#' @param events Targeted events at each analysis
#' @param maxEvents Final planned events
#' @param analysisTimes Not yet implemented
#' @param upper Function to compute upper bound
#' @param upar Parameter passed to \code{upper()}
#' @param lower Function to compute lower bound
#' @param lpar Parameter passed to \code{lower()}
#' @param r Control for grid size; normally leave at default of \code{r=18}
#'
#' @return a \code{tibble} with columns Analysis, Bound, Z, Probability, theta, Time, avehr, Events
#' @details Need to be added
#' @export
#'
#' @examples
#' library(gsDesign)
#' library(gsDesign2)
#' gs_power_nph() %>% filter(abs(Z) < Inf)
gs_power_nph <- function(enrollRates=tibble::tibble(Stratum="All",
                                                  duration=c(2,2,10),
                                                  rate=c(3,6,9)),
                       failRates=tibble::tibble(Stratum="All",
                                                duration=c(3,100),
                                                failRate=log(2)/c(9,18),
                                                hr=c(.9,.6),
                                                dropoutRate=rep(.001,2)),
                       ratio=1,
                       events = c(30, 40, 50),
                       analysisTimes = NULL,
                       maxEvents = 45,       # max planned events
                       upper = gs_b,
                       # Default is Lan-DeMets approximation of
                       upar = gsDesign(k=length(events), test.type=1,
                                       n.I=events, maxn.IPlan = maxEvents,
                                       sfu=sfLDOF, sfupar = NULL)$upper$bound,
                       lower = gs_b,
                       lpar = c(qnorm(.1), rep(-Inf, 2)), # Futility only at IA1
                       r = 18
){
  avehr <- gs_info_ahr(enrollRates = enrollRates,
                   failRates = failRates,
                   ratio = ratio,
                   events = events,
                   analysisTimes = analysisTimes
  )
  x <- gs_prob(info = avehr$info,
              theta = -log(avehr$AHR),
              upper = upper,
              upar = upar,
              lower = lower,
              lpar = lpar,
              r = r) %>% select(-c("info","theta"))
  return(left_join(x, avehr, by = "Analysis"))
}
