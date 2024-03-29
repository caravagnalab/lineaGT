% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getter_params.R
\name{get_mean}
\alias{get_mean}
\title{Extract the estimated mean parameters.}
\usage{
get_mean(x)
}
\arguments{
\item{x}{a mvnmm object.}
}
\value{
dataframe of mean parameters.
}
\description{
Returns a dataframe \code{KxT} with the estimated mean paramaters \code{mu_kt} per
clone \code{k} and dimension \code{t}.
}
\examples{
if (FALSE) get_mean(x)

}
