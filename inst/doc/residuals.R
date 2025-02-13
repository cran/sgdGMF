## ----setup, include = FALSE---------------------------------------------------
options(rmarkdown.html_vignette.check_title = FALSE)

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----sgdgmf-------------------------------------------------------------------
library(sgdGMF)

## ----libraries----------------------------------------------------------------
library(ggplot2)
library(ggpubr)
library(reshape2)

## ----data---------------------------------------------------------------------
# install.packages("mvabund")
# data(antTraits, package = "mvabund")

load(url("https://raw.githubusercontent.com/cran/mvabund/master/data/antTraits.RData"))

Y = as.matrix(antTraits$abund)
X = as.matrix(antTraits$env[,-3])
Z = matrix(1, nrow = ncol(Y), ncol = 1)

n = nrow(Y)
m = ncol(Y)

## ----family-------------------------------------------------------------------
family = poisson()

## ----rank---------------------------------------------------------------------
ncomp = sgdgmf.rank(Y = Y, X = X, Z = Z, family = family)$ncomp
cat("Selected rank: ", ncomp)

## ----fit----------------------------------------------------------------------
gmf = sgdgmf.fit(Y, X, Z, ncomp = ncomp, family = family, method = "airwls")

## ----resid--------------------------------------------------------------------
res = residuals(gmf, spectrum = TRUE, ncomp = 20)

## ----plot, fig.width = 7, fig.height = 5--------------------------------------
ggpubr::ggarrange(
  plot(gmf, type = "res-idx"),
  plot(gmf, type = "res-fit"),
  plot(gmf, type = "hist"),
  plot(gmf, type = "qq"),
  nrow = 2, ncol = 2, align = "hv")

## ----spectrum, fig.width = 7, fig.height = 3----------------------------------
ggpubr::ggarrange(
  screeplot(gmf, cumulative = FALSE, proportion = TRUE),
  screeplot(gmf, cumulative = TRUE, proportion = TRUE),
  nrow = 1, ncol = 2, align = "hv")

## ----resid2, fig.width = 7, fig.height = 3.5----------------------------------
plt.dev = image(gmf, type = "deviance", resid = TRUE, symmetric = TRUE)
plt.prs = image(gmf, type = "pearson", resid = TRUE, symmetric = TRUE)

ggpubr::ggarrange(
  plt.dev + labs(x = "Species", y = "Environments", title = "Deviance residuals"), 
  plt.prs + labs(x = "Species", y = "Environments", title = "Pearson residuals"),
  nrow = 1, ncol = 2, common.legend = FALSE, legend = "bottom", align = "hv")

