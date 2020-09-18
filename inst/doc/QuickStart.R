## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  dev = "png", 
  dev.args= list(type = "cairo")
)

## ----setup, warning = FALSE---------------------------------------------------
library(gsdmvn)
library(gsDesign)
library(gsDesign2)
library(simtrial)
library(knitr)

## -----------------------------------------------------------------------------
enrollRates <- tibble::tibble(Stratum = "All", duration = c(2, 2, 2, 6), rate = (1:4)/4)
enrollRates

## -----------------------------------------------------------------------------
medianSurv <- 12
failRates = tibble::tibble(Stratum = "All",
                           duration = c(4, Inf),
                           failRate = log(2) / medianSurv,
                           hr = c(1, .6),
                           dropoutRate = .001)
failRates

## -----------------------------------------------------------------------------
# Type I error
alpha <- .025
design <-
   gs_design_ahr(enrollRates = enrollRates,
                 failRates = failRates,
                 alpha = alpha,
                 beta = .1, # Type II error = 1 - power
                 analysisTimes = 36, # Planned trial duration
                 IF = 1, # single analysis at information-fraction of 1
                 upar = qnorm(1 - alpha), # Final analysis bound
                 lpar = -Inf # No futility bound
               )

## -----------------------------------------------------------------------------
names(design)

## -----------------------------------------------------------------------------
design$enrollRates

## ----message=FALSE,warning=FALSE----------------------------------------------
design$bounds %>% kable(digits=c(0,0,1,0,0,3,4,3,3,2,2))

## -----------------------------------------------------------------------------
gsDesign::nEvents(hr=design$bounds$AHR)

## -----------------------------------------------------------------------------
  gs_design_ahr(enrollRates = enrollRates,
              failRates = failRates,
              alpha = alpha,
              beta = .1, # Type II error = 1 - power
              analysisTimes = 30, # Planned trial duration
              IF = 1, # single analysis at information-fraction of 1
              upar = qnorm(1 - alpha), # Final analysis bound
              lpar = -Inf # No futility bound
            )$bounds %>% kable(digits=c(0,0,1,0,0,3,4,3,3,2,2))

## ---- echo=FALSE, fig.width=6.5-----------------------------------------------
t <- (0:50)/50
plot(t, 2 - 2 * pnorm(qnorm(1-.0125)/sqrt(t)), type="l", ylab = "f(t)", xlab = "t")

## -----------------------------------------------------------------------------
b <- gsDesign::gsDesign(k=3, timing = c(0.5, 0.75, 1), test.type=1, alpha = 0.025, sfu = gsDesign::sfLDOF)$upper$bound
b

## -----------------------------------------------------------------------------
design1s <- gs_design_ahr(enrollRates = enrollRates,
                          failRates = failRates,
                          analysisTimes = 36, # Trial duration
                          upper = gs_spending_bound,
                          upar = list(sf = gsDesign::sfLDOF, total_spend = 0.025),
                          lpar = rep(-Inf, 3), # No futility bound
                          IF = c(.5, .75, 1)
                        )

## -----------------------------------------------------------------------------
design1s$bounds %>% filter(Bound == "Upper") %>% kable(digits=c(0,0,1,0,0,3,4,3,3,2,2))

## -----------------------------------------------------------------------------
gs_power_ahr(enrollRates = design1s$enrollRates,
             failRates = design1s$failRates %>% mutate(hr = 1),
             upper = gs_spending_bound,
             upar = list(sf = gsDesign::sfLDOF, total_spend = 0.025),
             lpar = rep(-Inf, 3), # No futility bound
             events = design1s$bound$Events[1:3]
) %>% filter(Bound == "Upper") %>% kable(digits=c(0,0,1,0,0,3,4,3,3,2,2))

## -----------------------------------------------------------------------------
b2 <- gsDesign::gsDesign(test.type = 2, sfu = sfLDOF, alpha = 0.025, timing = c(.5, .75, 1))$upper$bound
b2

## -----------------------------------------------------------------------------
design2ss <- gs_design_ahr(enrollRates = enrollRates,
                           failRates = failRates,
                           analysisTimes = 36, # Trial duration
                           IF = c(.5, .75, 1), # Information fraction at analyses
                           upper = gs_spending_bound,
                           upar = list(sf = gsDesign::sfLDOF, total_spend = 0.025),
                           lower = gs_spending_bound,
                           lpar = list(sf = gsDesign::sfLDOF, total_spend = 0.025),
                           h1_spending = FALSE
                        )

## ---- message=FALSE-----------------------------------------------------------
design2ss$bounds %>% kable(digits=c(0,0,1,0,0,3,4,3,3,2,2))

## ----fig.width=6.5------------------------------------------------------------
ggplot(data = design2ss$bound, aes(x=Events, y=Z, group = Bound)) + geom_line(aes(linetype=Bound)) + geom_point() +
  ggtitle("2-sided symmetric bounds with O'Brien-Fleming-like spending")

## -----------------------------------------------------------------------------
b2sa <- c(Inf, b) # Same efficacy bound as before
a2sa <- c(rep(qnorm(.05), 2), rep(-Inf, 2)) # Single futility analysis bound

design2sa <- gs_design_nph(enrollRates = enrollRates,
                        failRates = failRates,
                        analysisTimes = 36, # Trial duration
                        upar = b2sa,
                        lpar = a2sa, # Asymmetric 2-sided bound
                        IF = c(.3, .5, .75, 1)
                        )

## -----------------------------------------------------------------------------
design2sa$enrollRates %>% summarise(N = ceiling(sum(rate * duration)))

## -----------------------------------------------------------------------------
design2sa$bounds %>% filter(abs(Z) < Inf) %>% kable(digits=c(0,0,1,0,0,3,4,3,3,2,2))

## -----------------------------------------------------------------------------
events <- (design2sa$bounds %>% filter(Bound == "Upper"))$Events
gs_power_nph(enrollRates = design1s$enrollRates,
           failRates = design2sa$failRates %>% mutate(hr = 1),
           upar = b2sa,
           lpar = a2sa,
           events = events,
           maxEvents = max(events)) %>% 
# filter eliminates bounds that are infinite
  filter(abs(Z) < Inf)  %>% kable(digits=c(0,0,1,0,0,3,4,3,3,2,2))

## -----------------------------------------------------------------------------
fr <- simfix2simPWSurv(failRates = failRates)
nsim <- 200 # Number of trial simulations
simresult <- NULL
N <- ceiling(design2sa$enrollRates %>% summarize(N = sum(rate/2 * duration)))*2
K <- max(design2sa$bounds$Analysis)
events <- ceiling(sort(unique(design2sa$bounds$Events)))

for(i in 1:nsim){
  sim <- simPWSurv(n = as.numeric(N),
                   enrollRates = design2sa$enrollRates,
                   failRates = fr$failRates,
                   dropoutRates = fr$dropoutRates
                  )
  for(k in 1:K){
    Z <- sim %>% cutDataAtCount(events[k]) %>% # cut simulation for analysis at targeted events
         tensurv(txval = "Experimental") %>% tenFH(rg = tibble(rho = 0, gamma = 0))
    simresult <- rbind(simresult,
                       tibble(sim = i, k = k, Z = -Z$Z)) # Change sign for Z
  }
}

## -----------------------------------------------------------------------------
bds <- tibble::tibble(k = sort(unique(design2sa$bounds$Analysis)),
            upper = (design2sa$bounds %>% filter(Bound == "Upper"))$Z,
            lower = (design2sa$bounds %>% filter(Bound == "Lower"))$Z
)
trialsum <- simresult %>% 
            full_join(bds, by = "k") %>%
            filter(Z < lower | Z >= upper | k == K) %>%
            group_by(sim) %>%
            slice(1) %>% 
            ungroup()
trialsum %>% summarize(nsim = n(),
                       "Early futility (%)" = 100 * mean(Z < lower),
                       "Power (%)" = 100 * mean(Z >= upper))

## -----------------------------------------------------------------------------
xx <- trialsum %>% mutate(Positive = (Z >= upper))
table(xx$k, xx$Positive)

## -----------------------------------------------------------------------------
design1sPH <- gs_design_nph(enrollRates = enrollRates,
                        failRates = failRates %>% mutate(hr = .7),
                        analysisTimes = 36, # Trial duration
                        upar = b,
                        lpar = rep(-Inf, 3), # No futility bound
                        IF = c(.5, .75, 1)
                        )
design1sPH$enrollRates %>% kable()

## -----------------------------------------------------------------------------
design1sPH$enrollRates %>% tail(1) %>% select(N)

## ----warning=FALSE------------------------------------------------------------
x <- gsSurv(k = 3,
            test.type = 1,
            alpha = .025,
            beta = .1,
            timing = c(.5, .75, 1),
            sfu = sfLDOF,
            lambda = log(2) / medianSurv,
            hr = .7, 
            eta = .001, # dropout rate
            gamma = c(.25, .5, .75, 1),
            R = c(2, 2, 2, 6),
            minfup = 24,
            T = 36)
gsBoundSummary(x) %>% kable(row.names = FALSE)

## -----------------------------------------------------------------------------
x$gamma %>% kable()

