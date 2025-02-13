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
ncomp = sgdgmf.rank(Y = Y, X = X, Z = Z, family = family, method = "onatski")$ncomp
cat("Selected rank: ", ncomp)

## ----crossval-----------------------------------------------------------------
# Uncomment to run cross-validation
# crossval = sgdgmf.cv(Y = Y, X = X, Z = Z, family = family, ncomps = seq(1,5,1),
#                      method = "sgd", sampling = "block", control.cv = list(refit = FALSE))

## ----fit----------------------------------------------------------------------
gmf = sgdgmf.fit(Y, X, Z, ncomp = ncomp, family = family, method = "sgd", sampling = "block")

## ----cor----------------------------------------------------------------------
yhat_glm = fitted(gmf, type = "response", partial = TRUE) # VGLM model without matrix factorization
yhat_gmf = fitted(gmf, type = "response", partial = FALSE) # complete GMF model

cat(" VGLM: ", 100 * round(cor(c(Y), c(yhat_glm)), 4), "\n",
    "  GMF: ", 100 * round(cor(c(Y), c(yhat_gmf)), 4), "\n", sep = "")

## ----hist, fig.width = 7, fig.height = 5--------------------------------------
plt.res.fit.glm  = plot(gmf, type = "res-fit", partial = TRUE)
plt.res.hist.glm = plot(gmf, type = "hist", partial = TRUE)
plt.res.fit.gmf  = plot(gmf, type = "res-fit", partial = FALSE)
plt.res.hist.gmf = plot(gmf, type = "hist", partial = FALSE)

ggpubr::ggarrange(
  plt.res.fit.glm + ggtitle("Residuals vs Fitted values (VGLM)"), 
  plt.res.hist.glm + ggtitle("Histogram of the residuals (VGLM)"), 
  plt.res.fit.gmf + ggtitle("Residuals vs Fitted values (GMF)"), 
  plt.res.hist.gmf + ggtitle("Histogram of the residuals (GMF)"), 
  nrow = 2, ncol = 2, align = "hv")


## ----spectrum, fig.width = 7, fig.height = 3----------------------------------
plt.eig.glm = screeplot(gmf, partial = TRUE) + ggtitle("Residual screeplot (VGLM)")
plt.eig.gmf = screeplot(gmf, partial = FALSE) + ggtitle("Residual screeplot (GMF)")

ggpubr::ggarrange(plt.eig.glm, plt.eig.gmf, nrow = 1, ncol = 2, align = "hv")

## ----pred, fig.width = 7, fig.height = 3.5------------------------------------
plt.ant = image(gmf, limits = range(c(Y)), type = "data")
plt.fit = image(gmf, limits = range(c(Y)), type = "response")

ggpubr::ggarrange(
  plt.ant + labs(x = "Species", y = "Environments", title = "Observed abundance"), 
  plt.fit + labs(x = "Species", y = "Environments", title = "Predicted abundance"), 
  nrow = 1, ncol = 2, common.legend = TRUE, legend = "bottom", align = "hv")

## ----resid2, fig.width = 7, fig.height = 3.5----------------------------------
plt.dev = image(gmf, type = "deviance", resid = TRUE, symmetric = TRUE)
plt.prs = image(gmf, type = "pearson", resid = TRUE, symmetric = TRUE)

ggpubr::ggarrange(
  plt.dev + labs(x = "Species", y = "Environments", title = "Deviance residuals"), 
  plt.prs + labs(x = "Species", y = "Environments", title = "Pearson residuals"), 
  nrow = 1, ncol = 2, common.legend = FALSE, legend = "bottom", align = "hv")

## ----scores, fig.width = 7, fig.height = 4------------------------------------
biplot(gmf, titles = c("Environments", "Species"))

## -----------------------------------------------------------------------------
# Scores
head(coef(gmf, type = "scores"))

# Loadings
head(coef(gmf, type = "loadings"))

