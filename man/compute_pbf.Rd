% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Compute_PBF.R
\name{compute_pbf}
\alias{compute_pbf}
\title{compute_pbf}
\usage{
compute_pbf(vec_vals, barcode, k = 2)
}
\arguments{
\item{vec_vals}{values of x for which the pbf function must be computed}

\item{barcode}{matrix with the barcode information generated by compute_homology or compute_homology_multi functions.}

\item{k}{Parameter k for exponentiation.}
}
\value{
Returns a list with the PBF for each tested dimension.
}
\description{
compute_pbf
}
\examples{
\dontrun{
compute_pbf(vec_vals,barcode,2)
}
}
