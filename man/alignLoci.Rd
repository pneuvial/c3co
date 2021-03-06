% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/alignLoci.R
\name{alignLoci}
\alias{alignLoci}
\title{Aligns Multiple Copy-Number Data Objects}
\usage{
alignLoci(dat)
}
\arguments{
\item{dat}{A list of PSCN data.frame:s with required columns \code{chr} and
\code{pos}.}
}
\value{
A list of the same length as \code{dat} where all data.frame:s have
the exact same set of fields \code{chr} and \code{pos} (and in the same
order).
}
\description{
Expands multiple copy-number data.frame:s to have the exact same set of
(chr, pos) loci.  If new rows are injected, the other non-(chr, pos) fields
are populated with missing values.
}
