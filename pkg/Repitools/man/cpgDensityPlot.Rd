\name{cpgDensityPlot}
\alias{cpgDensityPlot}
\alias{cpgDensityPlot,GRangesList-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot the distribution of sequencing reads CpG densities.}
\description{Function to generate a plot of the distribution of sequencing reads CpG densities.}
\usage{
  \S4method{cpgDensityPlot}{GRangesList}(x, cols=rainbow(length(x)), xlim=c(0,20), lty = 1, lwd = 1, main="CpG Density Plot", verbose=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A \code{GRangesList} object of reads to plot CpG density of}
  \item{cols}{The line colour for each element of \code{x}}
  \item{xlim}{\code{xlim} parameter passed to \code{plot}.}
  \item{lty}{The line type for each element of \code{x}}
  \item{lwd}{The line width for each element of \code{x}}
  \item{main}{\code{main} parameter passed to \code{plot}}
  \item{verbose}{Print details of processing.}
  \item{...}{Arguments passed into \code{cpgDensityCalc}. \code{seq.len}
             and \code{organism} are required.}
}
\details{
  See \code{cpgDensityCalc} for details of options for calculating the CpG density.
}
\value{
  A plot is created. The data processed by \code{cpgDensityCalc} is invisibly returned. 
}
\author{Aaron Statham}
\examples{
if(require(BSgenome.Hsapiens.UCSC.hg18))
{
  readsFile <- system.file("data","samplesList.RData", package="Repitools")
  load(readsFile)  # GRangesList of reads 'samplesList'

  cpgDensityPlot(samplesList, seq.len=300, organism=Hsapiens, lwd=4, verbose=TRUE)
}

}
