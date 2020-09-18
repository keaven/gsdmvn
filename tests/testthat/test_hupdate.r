library(gsdmvn)
context("Updated grid and weights for numerical integration")

test_that("Testing hupdate() vs known results from gsProbability", {
  x <- gsDesign::gsProbability(k=2, n.I=1:2, theta = 2, a= c(-.5,0), b=c(3,2))
  expect_lt(abs(x$lower$prob[2] -
    hupdate(theta = 2, I = 2, b = 0, thetam1 = 2, Im1 = 1, gm1 = h1(theta = 2, I = 1, a = -.5, b = 3)) %>%
      summarize(plower2 = sum(h)) %>% as.numeric()),
    1e-6
  )
  expect_lt(abs(x$upper$prob[2] -
                 # Compare second upper crossing
                 hupdate(theta = 2, I = 2, a = 2, thetam1 = 2, Im1 = 1, gm1 = h1(theta = 2, I = 1, a = -.5, b = 3)) %>%
                 summarize(plower2 = sum(h)) %>% as.numeric()),
            1e-6
  )
})
