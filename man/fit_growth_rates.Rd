% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/growth_model.R
\name{fit_growth_rates}
\alias{fit_growth_rates}
\title{Infer growth rates for each clone and subclone.}
\usage{
fit_growth_rates(
  x,
  steps = 500,
  highlight = c(),
  timepoints_to_int = c(),
  growth_model = "exp.log",
  force = T,
  tree_score = 1,
  py_pkg = NULL,
  mutations = F
)
}
\arguments{
\item{x}{a mvnmm object.}

\item{steps}{maximum number of steps for inference.}

\item{highlight}{set of clusters to run the inference for.
If not specified, it will be run on all the clusters.}

\item{timepoints_to_int}{if the provided timepoints are not integers nor a timepoints-to-int list is
stored in \code{x}, a list mapping their values to integers is required.}

\item{growth_model}{string specifying the type of growth model to use, between \code{exp} and \code{log}
corresponding to exponential and logistic models, respectively.}

\item{force}{if the model has already been fitted, setting \code{force} to \code{FALSE} will keep the
computed rates. Setting \code{force} to \code{TRUE} will fit the model again for the specified clusters.}
}
\value{
A \code{mvnmm} object with the additional tibble \code{growth.rates} containing the estimated population genetics parameters.
}
\description{
Infer growth rates for each clone and subclone.
}
