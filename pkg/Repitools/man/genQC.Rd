\name{genQC}
\alias{genQC}
\alias{genQC,character-method}

\title{Plot Quality Checking Information for Sequencing Data}
\description{A series of quality control plots for sequencing data are made.}
\usage{
  \S4method{genQC}{character}(paths, expt)
}
\arguments{
  \item{paths}{A vector of character strings, each containing an absolute path
               to an Rdata file of a \code{\linkS4class{SequenceQC}} object.}
  \item{expt}{The names of the experiments which the lanes are about.}
}
\value{
  The function is called for its output. The output is multiple pages, so the
  pdf device should be called before this function is.
}
\details{
  The paths can be named, in which case this gives the names of the lanes to use
  in the plotting. Otherwise the lanes will be given the names \code{"Lane 1"},
  \code{"Lane 2"}, ..., \code{"Lane n"}.
}
\author{Dario Strbenac}
\references{
  FastQC: \url{http://www.bioinformatics.bbsrc.ac.uk/projects/fastqc/}
}
\examples{
	# Not given.
}


