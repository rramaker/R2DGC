% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FindProblemIons.R
\name{FindProblemIons}
\alias{FindProblemIons}
\title{This function scans over a range of potential ions and looks for ions that are not present in any peak or are common enough to decrease alignment quality.
Can use as input to the consensus align function to avoid including these ions during alignment to speed up processing time and improve alignments}
\usage{
FindProblemIons(inputFile, possibleIons = c(70:600), numCores = 1,
  absentIonThreshold = 0.01, commonIonThreshold = 2, plotData = T)
}
\arguments{
\item{inputFile}{The file path of a representative chromatof file to use in searching for ions to filter}

\item{possibleIons}{A numeric vector of possible ions to search. Make sure each ion listed is present in the input file. Defaults to 70 through 600.}

\item{numCores}{The number of cores to use for parallel processing. Defaults to 1}

\item{absentIonThreshold}{Numeric indicating the fraction of total ion intensity an ion has to greater than in at least one peak to not be filtered as an absent ion. Defaults to 0.01.}

\item{commonIonThreshold}{Numeric indicating the number of standard deviations below the mean an ion has to decrease the global metabolite similarity score to be filtered as a common ion. Defaults to 2.}

\item{plotData}{Boolean. If true, relative ion impact scores will be plotted.}
}
\value{
Two column data frame identifying filtered ions and the reason they were filtered (absent or common). If plotData is TRUE, plots common ion scores. Y-axis is the z-scored number of pairwise metabolite comparisons with a similarity score greater than 50. X-axis is the ion with the filtered ions labeled in red.
}
\description{
This function scans over a range of potential ions and looks for ions that are not present in any peak or are common enough to decrease alignment quality.
Can use as input to the consensus align function to avoid including these ions during alignment to speed up processing time and improve alignments
}
\examples{
FindProblemIons(inputFile=system.file("extdata", "SampleA.txt", package="R2DGC"),
    possibleIons = 70:78,plotData=FALSE)

}
