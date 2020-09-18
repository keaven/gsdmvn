## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  dev = "png", 
  dev.args= list(type = "cairo")
)
library(gsdmvn)

## -----------------------------------------------------------------------------
analysisTimes <- c(18, 24, 30, 36)
enrollRates <- tibble::tibble(Stratum = "All",
                              duration = c(2, 2, 2, 6),
                              rate = c(8, 12, 16, 24))
failRates <- tibble::tibble(Stratum = "All",
                            duration=c(3,100),
                            failRate=log(2)/c(8,14),
                            hr=c(.9,.6),
                            dropoutRate=.001
                           )

## -----------------------------------------------------------------------------
library(gsDesign2)
xx <- gsDesign2::AHR(enrollRates = enrollRates, failRates = failRates, totalDuration = analysisTimes)
Events <- ceiling(xx$Events)
yy <- gs_info_ahr(enrollRates = enrollRates, failRates = failRates, events = Events)

## -----------------------------------------------------------------------------
zz <- gs_power_npe(theta = yy$theta, info = yy$info, info0 = yy$info0, 
             upper = gs_spending_bound, lower = gs_spending_bound,
             upar = list(sf = gsDesign::sfLDOF, total_spend = 0.025, param = NULL, timing = NULL),
             lpar = list(sf = gsDesign::sfLDOF, total_spend = 0.025, param = NULL, timing = NULL)
)
zz

## -----------------------------------------------------------------------------
K <- 4
minx <- ((qnorm(.025) / sqrt(zz$info0[K]) + qnorm(.1) / sqrt(zz$info[K])) / zz$theta[K])^2
minx

## -----------------------------------------------------------------------------
gs_power_npe(theta = yy$theta[K], info = yy$info[K] * minx, info0 = yy$info0[K] * minx, upar = qnorm(.975), lpar = -Inf) %>% 
  filter(Bound == "Upper")

## -----------------------------------------------------------------------------
zz <- gs_power_npe(theta = yy$theta, info = yy$info * minx, info0 = yy$info0 * minx, 
             upper = gs_spending_bound, lower = gs_spending_bound,
             upar = list(sf = gsDesign::sfLDOF, total_spend = 0.025, param = NULL, timing = NULL),
             lpar = list(sf = gsDesign::sfLDOF, total_spend = 0.025, param = NULL, timing = NULL)
)
zz

## -----------------------------------------------------------------------------
zz <- gs_power_npe(theta = yy$theta, info = yy$info * minx * 1.2, info0 = yy$info0 * minx * 1.2, 
             upper = gs_spending_bound, lower = gs_spending_bound,
             upar = list(sf = gsDesign::sfLDOF, total_spend = 0.025, param = NULL, timing = NULL),
             lpar = list(sf = gsDesign::sfLDOF, total_spend = 0.025, param = NULL, timing = NULL)
)
zz

## -----------------------------------------------------------------------------
theta <- yy$theta
info <- yy$info
info0 <- yy$info0
upper = gs_spending_bound
lower = gs_spending_bound
upar = list(sf = gsDesign::sfLDOF, total_spend = 0.025, param = NULL, timing = NULL)
lpar = list(sf = gsDesign::sfLDOF, total_spend = 0.025, param = NULL, timing = NULL)
alpha = .025
beta = .1
binding = FALSE
test_upper = TRUE
test_lower = TRUE
r = 18
tol = 1e-06

zz <- gs_design_npe(theta = yy$theta, info = yy$info, info0 = yy$info0, 
             upper = gs_spending_bound, lower = gs_spending_bound,
             upar = list(sf = gsDesign::sfLDOF, total_spend = 0.025, param = NULL, timing = NULL),
             lpar = list(sf = gsDesign::sfLDOF, total_spend = 0.025, param = NULL, timing = NULL)
)
zz

