% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sgdGMF-class.R
\name{biplot.sgdgmf}
\alias{biplot.sgdgmf}
\title{Biplot of a GMF model}
\usage{
\method{biplot}{sgdgmf}(
  x,
  ...,
  choices = 1:2,
  arrange = TRUE,
  byrow = FALSE,
  normalize = FALSE,
  labels = NULL,
  palette = NULL,
  titles = c(NULL, NULL)
)
}
\arguments{
\item{x}{an object of class \code{sgdgmf}}

\item{...}{further arguments passed to or from other methods}

\item{choices}{a length 2 vector specifying the components to plot}

\item{arrange}{if \code{TRUE}, return a single plot with two panels}

\item{byrow}{if \code{TRUE}, the panels are arranged row-wise (if \code{arrange=TRUE})}

\item{normalize}{if \code{TRUE}, orthogonalizes the scores using SVD}

\item{labels}{a vector of labels which should be plotted}

\item{palette}{the color-palette which should be used}

\item{titles}{a 2-dimensional string vector containing the plot titles}
}
\value{
If \code{arrange=TRUE}, a single ggplot object with the selected biplots,
otherwise, a list of two ggplot objects showing the row and column latent variables.
}
\description{
Plot the observations on a two-dimensional projection determined by the
estimated score matrix
}
\examples{
\donttest{# Load the sgdGMF package
library(sgdGMF)

# Generate data from a Poisson model
data = sim.gmf.data(n = 100, m = 20, ncomp = 5, family = poisson())

# Fit a GMF model
gmf = sgdgmf.fit(data$Y, ncomp = 3, family = poisson())

# Get the biplot of a GMF model
biplot(gmf)
}
}
