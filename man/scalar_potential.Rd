% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HodgePotentialHumanFlow-package.R
\name{scalar_potential}
\alias{scalar_potential}
\title{Scalar potential of origin-destionation matrix by Hodge-Kodaira decompostion}
\usage{
scalar_potential(od_table)
}
\arguments{
\item{od_table}{A data.frame with three columns: origin, dest, and trips}
}
\value{
a data.frame with two columns: geozone and negative potential
}
\description{
This function calculates scalar potential of human flow for a given origin-destionation matrix.
}
