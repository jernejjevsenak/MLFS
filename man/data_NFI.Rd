% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_NFI.R
\docType{data}
\name{data_NFI}
\alias{data_NFI}
\title{An example of national forest inventory data}
\format{
A data frame with 11984 rows and 10 variables:
\describe{
  \item{plotID}{a unique identifier for plot}
  \item{treeID}{a unique identifier for tree}
  \item{year}{year in which plot was visited}
  \item{speciesGroup}{identifier for species group}
  \item{code}{status of a tree: 0 (normal), 1(harvested), 2(dead), 3 (ingrowth)}
  \item{DBH}{diameter at breast height in cm}
  \item{species}{species name}
  \item{height}{tree height in meters}
  \item{crownHeight}{crown height in meters}
  \item{protected}{logical, 1 if protected, otherwise 0}
}
}
\usage{
data_NFI
}
\description{
This is simulated data that reassemble the national forest inventory
}
\keyword{datasets}
