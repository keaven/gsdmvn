#' gs_b: Default boundary generation
#'
#' \code{gs_b()} is the simplest version of a function to be used with the \code{upper} and \code{lower}
#' arguments in \code{gs_prob()},
#' \code{gs_power_nph} and \code{gs_design_nph()};
#' it simply returns the vector input in the input vector \code{Z} or, if \code{k} is specified \code{par[k]j} is returned.
#' Note that if bounds need to change with changing information at analyses, \code{gs_b()} should not be used.
#' For instance, for spending function bounds use
#' @param par For \code{gs_b()}, this is just Z-values for the boundaries; can include infinite values
#' @param info Information at analyses; not used in this function; present as it is a standard parameter for other boundary computation routines
#' @param k is NULL (default), return \code{par}, else return \code{par[k]}
#' @param ... further arguments passed to or from other methods
#' @return returns the vector input \code{par} if \code{k} is NULL, otherwise, \code{par[k]}
#' @export
#'
#' @examples
#' # Simple: enter a vector of length 3 for bound
#' gs_b(par = 4:2)
#'
#' # 2nd element of par
#' gs_b(4:2, k = 2)
#'
#' # Generate an efficacy bound using a spending function
#' # Use Lan-DeMets spending approximation of O'Brien-Fleming bound
#' # as 50%, 75% and 100% of final spending
#' # Information fraction
#' IF <- c(.5, .75, 1)
#' gs_b(par = gsDesign::gsDesign(alpha = .025, k= length(IF), test.type = 1, sfu = gsDesign::sfLDOF, timing = IF)$upper$bound)
gs_b <- function(par = NULL, info = NULL, k = NULL,...){
  if(is.null(k)){return(par)}else return(par[k])
}
