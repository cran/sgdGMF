% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/initGMF-class.R
\name{print.initgmf}
\alias{print.initgmf}
\title{Print the fundamental characteristics of an initialized GMF}
\usage{
\method{print}{initgmf}(x, ...)
}
\arguments{
\item{x}{an object of class \code{initgmf}}

\item{...}{further arguments passed to or from other methods}
}
\value{
No return value, called only for printing.
}
\description{
Print some summary information of an initialized GMF model.
}
\examples{
# Load the sgdGMF package
library(sgdGMF)

# Generate data from a Poisson model
data = sim.gmf.data(n = 100, m = 20, ncomp = 5, family = poisson())

# Fit a GMF model with 3 latent factors
init = sgdgmf.init(data$Y, ncomp = 3, family = poisson())

# Print the GMF object
print(init)

}
