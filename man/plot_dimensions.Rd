% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Plotting_Functions.R
\name{plot_dimensions}
\alias{plot_dimensions}
\title{plot_dimensions}
\usage{
plot_dimensions(barcode_df, resolution = 0.1)
}
\arguments{
\item{barcode_df}{Barcode to be plotted}

\item{resolution}{Resolution, interval width, for the construction of the Persistent Betty Function.}
}
\value{
Returns a plot of the barcodes and the Persistent Betti Functions.
}
\description{
Function that plots the barcode and the PBFs for a specific barcode.
}
\examples{
\dontrun{
plot_dimensions(barcode_df,0.1)
}
}
