\name{annotationBlocksCounts}
\alias{annotationBlocksCounts}
\alias{annotationBlocksCounts,ANY,data.frame-method}
\alias{annotationBlocksCounts,ANY,GRanges-method}
\title{Counts the number of sequencing reads within supplied genomic blocks.}
\description{
  Counts reads inside blocks.
}
\section{Usage}{
  \describe{
    \code{annotationBlocksCounts(x, anno, ...)}
  }
}
\section{Arguments}{
  \describe{
    \item{x:}{A character vector of BAM paths, a \code{GenomeDataList},
             \code{GRangesList}, or \code{GRanges} object.}
    \item{anno:}{A set of genomic features to make windows around a reference point of
               theirs. Either a \code{data.frame} with (at least) colums \code{chr},
               \code{start}, and \code{end}, or a \code{GRanges} object.}
    \item{seq.len:}{If sequencing reads need to be extended, the fragment size to be
                 used. Default: NULL (no extension).}
    \item{verbose:}{Whether to print progress. Default: TRUE.}
    \item{...:}{Parameters described above, that are not used in the top-level
               error-checking stage, but are passed further into a private
               function that uses them in its processing.}
  }
}
\value{
  A \code{matrix} of counts is returned, one column per sample and one row
  per row of genomic features supplied.
}
\author{Aaron Statham}
\seealso{\code{\link{annotationCounts}}, \code{\link{genomeBlocks}}}
\examples{#See the manual}
