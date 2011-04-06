\name{findClusters}
\alias{findClusters}
\title{Find Clusters Epigenetically Modified Genes}
\description{Given a table of gene positions that has a score column, genes will
  first be sorted into positional order and consecutive windows of high or low
  scores will be reported.
}
\usage{
  findClusters(stats, score.col, w.size = 5, n.med, n.consec, cut.samps, maxFDR = 0.05,
               trend = c("down", "up"), n.perm = 100, getFDRs = FALSE, verbose = TRUE)
}
\arguments{
  \item{stats}{A \code{data.frame} with (at least) column \code{chr}, and a column
               of scores. Genes must be sorted in positional order.}
  \item{score.col}{A number that gives the column in \code{stats} which contains
                   the scores.}
  \item{w.size}{The number of consecutive genes to consider windows over.}
  \item{n.med}{Minimum number of genes in a window, that have median score
               centred around them above a cutoff.}
  \item{n.consec}{Minimum cluster size.}
  \item{cut.samps}{A vector of score cutoffs to calculate the FDR at.}
  \item{maxFDR}{The highest FDR level still deemed to be significant.}
  \item{trend}{Whether the clusters must have all positive scores (enrichment),
               or all negative scores (depletion).}
  \item{n.perm}{How many random tables to generate to use in the FDR calculations.}
  \item{getFDRs}{If TRUE, will also return the table of FDRs at a variety of score
                 cutoffs, from which the score cutoff for calling clusters was chosen.}
  \item{verbose}{Whether to print progress of computations.}
}
\details{
  First, the median over a window of size \code{w.size} is calculated in a rolling
  window and then associated with the middle gene of the window. Windows are again
  run over the genes, and the gene at the centre of the window is significant if
  there are also at least \code{n.med} genes with representative medians
  above the score cutoff, in the window that surrounds it. These marker genes
  are extended outwards, for as long as the score has the same sign. The
  order of the \code{stats} rows is randomised, and this process in done for
  every randomisation.

  The procedure for calling clusters is done at a range of score cutoffs.
  The first score cutoff to give an FDR below \code{maxFDR} is chosen as the
  cutoff to use, and clusters are then called based on this cutoff.
}
\value{
  If \code{getFDRs} is FALSE, then only the \code{stats} table, with an
  additional column, \code{cluster}. If \code{getFDRs} is TRUE, then a list with
  elements :
    \item{table}{The table \code{stats} with the additional column \code{cluster}.}
    \item{FDR}{The table of score cutoffs tried, and their FDRs.}
}
\author{Dario Strbenac, Aaron Statham}
\references{Saul Bert, in preparation}
\examples{
	# Not run. Time consuming.
	#chrs <- paste("chr", c(1:22, "X", "Y"), sep = "")
	#starts <- sample(1:10000000, 25000, replace = TRUE)
	#ends <- starts + 10000
	#strands <- rep('+', 25000)
	#genes <- data.frame(chr = sample(chrs, 25000, replace = TRUE),
        #             start = starts, end = ends, strand = '+', t = rnorm(100, 0, 2))
	#findClusters(genes, 2, 5, 2, 3, seq(1, 10, 0.5), trend = "up", n.perm = 10)
}
