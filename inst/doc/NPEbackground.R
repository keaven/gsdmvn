## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  dev = "png", 
  dev.args= list(type = "cairo")
)

## ---- message=FALSE, warning=FALSE, echo=FALSE--------------------------------
library(tibble)
library(dplyr)
library(knitr)
mt <- tibble::tribble(~Statistic,    ~Example,   ~"Expected value", ~Variance,
  "$\\hat\\theta_k$", "$\\bar X_k$", "$\\theta(t_k)$",  "$\\mathcal{I}_k^{-1}$",
  "$Z_k=\\sqrt{\\mathcal{I}_k}\\hat\\theta_k$","$\\sqrt{n_k}\\bar X_k$", "$\\sqrt{\\mathcal{I}_k}\\theta(t_k)$",  "$1$",
  "$B_k=\\sqrt{t_k}Z_k$","$\\sum_{i=1}^{n_k}X_i/\\sqrt N$", "$t_k\\sqrt{\\mathcal{I}_K}\\theta(t_k)=\\mathcal{I}_k\\theta(t_k)/\\sqrt{\\mathcal{I}_K}$",  "$t_k$"
)
mt %>% kable(escape=FALSE)

