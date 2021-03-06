% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Similarity_Functions.R
\name{sample_structures_b}
\alias{sample_structures_b}
\title{sample_structures_b}
\usage{
sample_structures_b(x, rounds = 2)
}
\arguments{
\item{x}{output form sample_structures_a}

\item{rounds}{Number of elements to be selected from each group.}
}
\value{

}
\description{
Sample n structures from each group. Uses sample_structures_a output to select n similarities from each group.
}
\examples{
\dontrun{
sample_structures_b(x,rounds = 2)
}
}
