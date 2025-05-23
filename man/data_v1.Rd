% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_v1.R
\docType{data}
\name{data_v1}
\alias{data_v1}
\title{An example of joined national forest inventory and site data that is used
within the MLFS}
\format{
A data frame with 11984 rows and 15 variables:
\describe{
  \item{plotID}{a unique identifier for plot}
  \item{treeID}{a unique identifier for tree}
  \item{year}{year in which plot was visited}
  \item{speciesGroup}{identifier for species group}
  \item{code}{status of a tree: 0 (normal), 1(harvested), 2(dead), 3 (ingrowth)}
  \item{species}{species name}
  \item{height}{tree height in meters}
  \item{crownHeight}{crown height in meters}
  \item{protected}{logical, 1 if protected, otherwise 0}
  \item{slope}{slope on a plot}
  \item{elevation}{plot elevation}
  \item{northness}{plot northness, 1 is north, 0 is south}
  \item{siteIndex}{a proxy for site index, higher value represents more productive sites}
  \item{BA}{basal area of individual trees in m2}
  \item{weight}{upscale weight to calculate hectare values}
}
}
\usage{
data_v1
}
\description{
This is simulated data that reassemble the national forest inventory and
simulated data. We use it to show how to run examples for some specific functions
}
\keyword{datasets}
