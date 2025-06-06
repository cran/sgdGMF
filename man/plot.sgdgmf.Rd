% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sgdGMF-class.R
\name{plot.sgdgmf}
\alias{plot.sgdgmf}
\title{Plot diagnostics for a GMF model}
\usage{
\method{plot}{sgdgmf}(
  x,
  ...,
  type = c("res-idx", "res-fit", "std-fit", "hist", "qq", "ecdf"),
  resid = c("deviance", "pearson", "working", "response", "link"),
  subsample = FALSE,
  sample.size = 500,
  partial = FALSE,
  normalize = FALSE,
  fillna = FALSE
)
}
\arguments{
\item{x}{an object of class \code{sgdgmf}}

\item{...}{further arguments passed to or from other methods}

\item{type}{the type of plot which should be returned}

\item{resid}{the type of residuals which should be used}

\item{subsample}{if \code{TRUE}, computes the residuals over o small fraction of the data}

\item{sample.size}{the dimension of the sub-sample which should be used}

\item{partial}{if \code{TRUE}, computes the partial residuals}

\item{normalize}{if \code{TRUE}, standardizes the residuals column-by-column}

\item{fillna}{if \code{TRUE}, fills the \code{NA} values with \code{0}}
}
\value{
A ggplot object showing the selected diagnostic plot.
}
\description{
Plots (one of) six diagnostics to graphically analyze the marginal and conditional
distribution of the residuals of a GMF model. Currently, the following plots are
available: residuals against observation indices, residuals agains fitted values,
absolute square-root residuals against fitted values, histogram of the residuals,
residual QQ-plot, residual ECDF-plot.
}
\examples{
\donttest{# Load the sgdGMF package
library(sgdGMF)

# Generate data from a Poisson model
data = sim.gmf.data(n = 100, m = 20, ncomp = 5, family = poisson())

# Fit a GMF model
gmf = sgdgmf.fit(data$Y, ncomp = 3, family = poisson())

# Plot the residual-based GMF diagnostics
plot(gmf, type = "res-fit") # Residuals vs fitted values
plot(gmf, type = "std-fit") # Abs-sqrt-transformed residuals vs fitted values
plot(gmf, type = "qq") # Residual QQ-plot
plot(gmf, type = "hist") # Residual histogram
}
}
