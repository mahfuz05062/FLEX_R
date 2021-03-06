% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PostProcessing.R
\name{getSubsetOfCoAnnRemoveIDs}
\alias{getSubsetOfCoAnnRemoveIDs}
\title{Remove a complex (or complexes) to re-evaluate the data. This removes all the positive pairs associated with them complex (thereby effectively removing the complex).}
\usage{
getSubsetOfCoAnnRemoveIDs(data.standard, ids, replace = FALSE)
}
\arguments{
\item{data.standard}{input data: either co-annotation standard or output of CalculatePredictionAndTrueOnLibraryProfiles as a data.frame (need the 'ID' column for both)}

\item{ids}{the complex ids to remove}

\item{replace}{Way1 (replace = false): Remove the positive pairs associated to ids.
This reduces the size of the data a bit.
Way2 (replace = true): Convert the positive pairs to negatives. 
This will maintain the size (number of pairs) of the data.}
}
\description{
Remove a complex (or complexes) to re-evaluate the data. This removes all the positive pairs associated with them complex (thereby effectively removing the complex).
}
