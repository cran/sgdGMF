## ----setup, include = FALSE---------------------------------------------------
options(rmarkdown.html_vignette.check_title = FALSE)

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----sgdgmf, results = FALSE, message = FALSE, warning = FALSE----------------
library(sgdGMF)

## ----libraries, results = FALSE, message = FALSE, warning = FALSE-------------
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

## ----airwls, echo = FALSE, out.width = '80%'----------------------------------
# knitr::include_graphics("../man/alg-AIRWLS.png")

## ----newton, echo = FALSE, out.width = '80%'----------------------------------
# knitr::include_graphics("../man/alg-Newton.png")

## ----bsgd, echo = FALSE, out.width = '80%'------------------------------------
# knitr::include_graphics("../man/alg-BSGD.png")

## ----csgd, echo = FALSE, out.width = '80%'------------------------------------
# knitr::include_graphics("../man/alg-CSGD.png")

## ----fit----------------------------------------------------------------------
gmf.airwls = sgdgmf.fit(Y, X, Z, ncomp = ncomp, family = family, method = "airwls")
gmf.newton = sgdgmf.fit(Y, X, Z, ncomp = ncomp, family = family, method = "newton")
gmf.sgd    = sgdgmf.fit(Y, X, Z, ncomp = ncomp, family = family, method = "sgd", sampling = "block")

## ----cor----------------------------------------------------------------------
cat(" AIRWLS: ", 100 * round(cor(c(Y), c(gmf.airwls$mu)), 4), "\n",
    " Newton: ", 100 * round(cor(c(Y), c(gmf.newton$mu)), 4), "\n",
    " SGD:    ", 100 * round(cor(c(Y), c(gmf.sgd$mu)), 4), sep = "")

## ----resid--------------------------------------------------------------------
res.airwls = residuals(gmf.airwls, partial = FALSE, spectrum = TRUE, ncomp = 20)
res.newton = residuals(gmf.newton, partial = FALSE, spectrum = TRUE, ncomp = 20)
res.sgd = residuals(gmf.sgd, partial = FALSE, spectrum = TRUE, ncomp = 20)

## ----hist, echo = TRUE, include = FALSE, fig.width = 6, fig.height = 5--------
get.factor = function (levels, nrep) factor(rep(levels, each = nrep), levels = levels)

df = data.frame(
  residuals = c(res.airwls$residuals, res.newton$residuals, res.sgd$residuals),
  fitted = c(gmf.airwls$mu, gmf.newton$mu, gmf.sgd$mu),
  model = get.factor(c("AIRWLS", "Newton", "SGD"), n)
)

plt.res.vs.fit = ggplot(data = df, map = aes(x = fitted, y = residuals)) +
  geom_point(alpha = 0.5) + facet_grid(cols = vars(model)) + 
  geom_hline(yintercept = 0, col = 2, lty = 2)

plt.res.hist = ggplot(data = df, map = aes(x = residuals, y = after_stat(density))) +
  geom_histogram(bins = 30) + facet_grid(cols = vars(model)) + 
  geom_vline(xintercept = 0, col = 2, lty = 2)

ggpubr::ggarrange(plt.res.vs.fit, plt.res.hist, nrow = 2, align = "v")


## ----spectrum, echo = TRUE, include = FALSE, fig.width = 6, fig.height = 2.5----
data.frame(
  components = rep(1:20, 3),
  variance = c(res.airwls$lambdas, res.newton$lambdas, res.sgd$lambdas),
  model = get.factor(c("AIRWLS", "Newton", "SGD"), 20)) |>
  ggplot(map = aes(x = components, y = variance)) +
  geom_col(color = "white") + facet_grid(cols = vars(model))

## ----pred, echo = TRUE, include = FALSE, fig.width = 6, fig.height = 5--------
colnames(gmf.airwls$mu) = colnames(Y)
colnames(gmf.newton$mu) = colnames(Y)
colnames(gmf.sgd$mu) = colnames(Y)

df = rbind(
  reshape2::melt(Y, varnames = c("environments", "species"), value.name = "abundance"),
  reshape2::melt(gmf.airwls$mu, varnames = c("environments", "species"), value.name = "abundance"),
  reshape2::melt(gmf.newton$mu, varnames = c("environments", "species"), value.name = "abundance"),
  reshape2::melt(gmf.sgd$mu, varnames = c("environments", "species"), value.name = "abundance"))

df$model = get.factor(c("Data", "AIRWLS", "Newton", "SGD"), n*m)

ggplot(data = df, map = aes(x = species, y = environments, fill = abundance)) + 
  geom_raster() + facet_wrap(vars(model), nrow = 2, ncol = 2) + scale_fill_viridis_c() + 
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank()) +
  labs(x = "Species", y = "Environmens", fill = "Abundance", title = "Ant abundance data")

## ----scores, echo = TRUE, include = FALSE, fig.width = 8, fig.height = 5------
pca.airwls = RSpectra::svds(tcrossprod(gmf.airwls$U, gmf.airwls$V), 2)
pca.newton = RSpectra::svds(tcrossprod(gmf.newton$U, gmf.newton$V), 2)
pca.sgd = RSpectra::svds(tcrossprod(gmf.sgd$U, gmf.sgd$V), 2)

plt.envs = data.frame(
  idx = rep(1:n, times = 3),
  pc1 = c(scale(pca.airwls$u[,1]), scale(pca.newton$u[,1]), scale(pca.sgd$u[,1])),
  pc2 = c(scale(pca.airwls$u[,2]), scale(pca.newton$u[,2]), scale(pca.sgd$u[,2])),
  model = get.factor(c("AIRWLS", "Newton", "SGD"), n)) |>
  ggplot(map = aes(x = pc1, y = pc2, color = idx, label = idx)) +
  geom_hline(yintercept = 0, lty = 2, color = "grey40") + 
  geom_vline(xintercept = 0, lty = 2, color = "grey40") + 
  geom_point() + facet_grid(cols = vars(model)) + 
  scale_color_viridis_c(option = "viridis") +
  geom_text(color = 1, size = 2.5, nudge_x = -0.2, nudge_y = +0.2) +
  labs(x = "PC 1", y = "PC 2", color = "Environments")

plt.species = data.frame(
  idx = rep(1:m, times = 3),
  pc1 = c(scale(pca.airwls$v[,1]), scale(pca.newton$v[,1]), scale(pca.sgd$v[,1])),
  pc2 = c(scale(pca.airwls$v[,2]), scale(pca.newton$v[,2]), scale(pca.sgd$v[,2])),
  model = get.factor(c("AIRWLS", "Newton", "SGD"), m)) |>
  ggplot(map = aes(x = pc1, y = pc2, color = idx, label = idx)) +
  geom_hline(yintercept = 0, lty = 2, color = "grey40") + 
  geom_vline(xintercept = 0, lty = 2, color = "grey40") + 
  geom_point() + facet_grid(cols = vars(model)) + 
  scale_color_viridis_c(option = "inferno") +
  geom_text(color = 1, size = 2.5, nudge_x = -0.2, nudge_y = +0.2) +
  labs(x = "PC 1", y = "PC 2", color = "Species")

ggpubr::ggarrange(plt.envs, plt.species, nrow = 2, ncol = 1, align = "hv")

