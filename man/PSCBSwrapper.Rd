% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PSCBSwrapper.R
\name{PSCBSwrapper}
\alias{PSCBSwrapper}
\title{Transforms and aligns PSCBS data followed by a joint segmentation}
\usage{
PSCBSwrapper(PSCBSdata, stat, align = FALSE)
}
\arguments{
\item{PSCBSdata}{A list that contains PSCBS data}

\item{stat}{\code{"TCN"} or \code{"C1C2"} parameter to segment the data.
If \code{stat == "TCN"}, the segmentation will be done on TCN only.}

\item{align}{If \code{TRUE}, the PSCBS list elements are expanded before
segmentation to have the exact same set of (chr, pos) loci, which is the
union of all elements. See \code{\link[=alignLoci]{alignLoci()}} for details.}
}
\value{
A list which contains the breakpoints by chromosome and also the
binning of TCN, C1, and C2.
}
\description{
Transforms and aligns PSCBS data followed by a joint segmentation
}
