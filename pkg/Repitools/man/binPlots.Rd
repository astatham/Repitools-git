\name{binPlots}
\alias{binPlots}
\alias{binPlots,GRangesList-method}
\alias{binPlots,AffymetrixCelSet-method}
\alias{binPlots,matrix-method}
\title{Create line plots of averaged signal across a promoter}
\description{
Using a specified ordering of genes, they are split into multiple bins.  In each bin, the signal across is summarized and displayed visually.  The plots are created in an array, one row for each column of the data matrix and one column for each ordering specified.
}
\usage{
  \S4method{binPlots}{GRangesList}(x, anno=NULL, design=NULL, up=7500, down=2500, by=100, bw=300, lib.size="lane", seq.len=NULL, verbose=FALSE, Acutoff=NULL, ...)
  \S4method{binPlots}{AffymetrixCelSet}(x, probe.map=NULL, anno=NULL, up=7500, down=2500, by=100, bw=300, log2.adj=TRUE, verbose=FALSE, ...)   
  \S4method{binPlots}{matrix}(x, lookup.table=NULL, ordering=NULL, ord.label="", plot.type = c("line", "heatmap", "terrain", "boxplot"), nbins = 10, cols = NULL, lwd = 3, lty = 1, same.scale = TRUE, symm.scale = FALSE, verbose = FALSE, remove.zeros=TRUE, use.mean=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{\code{matrix}, or \code{GRangesList} object, or \code{AffymetrixCelSet} object of data to use}
  \item{probe.map}{A pre-made mapping of probes to a anno}
  \item{anno}{\code{data.frame} giving the gene coordinates or regions of interest. The column names of this are \code{chr}, \code{name}, \code{start}, \code{end}, \code{strand}.}
  \item{lookup.table}{\code{matrix}, e.g. from \code{makeWindowLookupTable}}
  \item{design}{Optional for GDL data - design matrix specifying the contrast to compute (i.e. the samples to use and what differences to take)}
  \item{ordering}{\code{data.frame} of either numeric or factor variables, with the same number of rows as \code{lookup.table}}
  \item{ord.label}{Character string that describes what type of data the ordering is. e.g. "log2 expression". Used to label relevant plot axis.}
  \item{plot.type}{type of plot}
  \item{nbins}{number of bins}
  \item{cols}{vector of colours}
  \item{lwd}{line width (either scalar or vector)}
  \item{lty}{line type (either scalar or vector)}
  \item{same.scale}{should the scale on all plots be the same}
  \item{symm.scale}{should the scale on plots be symmetrical ie centred around 0}
  \item{verbose}{Print details of processing.}
  \item{remove.zeros}{Remove zero values.}
  \item{use.mean}{Whether to use mean instead of median within bins (useful for seq data).}
  \item{up}{number of bases upstream to plot}
  \item{down}{number of bases upstream to plot}
  \item{by}{how many bases between data points}
  \item{bw}{how many bases up and downstream to be counted}
  \item{log2.adj}{whether to take log2 of array intensities.}
  \item{lib.size}{Either the string \code{"ref"} to normalise to a reference or \code{"lane"} to multiply each library to be the same size.}
  \item{seq.len}{If sequencing reads need to be extended, the fragment size to be used}
  \item{Acutoff}{If \code{lib.size} is \code{"ref"}, this argument must be provided. Otherwise, it must not. Cutoff on the "A" values to take, before calculating trimmed mean.}
  \item{...}{further arguments to send to \code{bxp}}
}
\details{
If \code{plot.type="line"}, a line is plotted for each bin across the promoter.
  
If \code{plot.type="heatmap"}, a series of bins are plotted as a heatmap.  This can be useful to display a larger number of bins.

If \code{plot.type="terrain"}, a series of bins are plotted as a 3d-terrain map.  This can be useful to display a larger number of bins.

The option \code{plot.type="boxplot"} is really only suitable for 1- or 2-group plots.  In this case, instead of plotting the summarized data, it plots the whole distribution of signal across the whole promoter.
}
\value{
  Either a single- or multiple-panel figure.
}
\author{Mark Robinson}
\examples{
  annoFile <- system.file("data","chr21genes.csv", package="Repitools")
  annoDF <- read.csv(annoFile)

  readsFile <- system.file("data","samplesList.RData", package="Repitools")
  load(readsFile)  # GRangesList of reads 'samplesList'
  
  exprFile <- system.file("data","expr.RData", package="Repitools")
  load(exprFile)  # matrix of differential expression 'expr'

  # design matrix
  des <- matrix( c(0,1,0,-1), ncol=1, dimnames=list(names(samplesList),"PC-Norm") )

  binPlots(samplesList, annoDF, design=des, verbose=TRUE, seq.len=300, nbins=4,
           ordering=expr, ord.label="expression",plot.type="line")

  binPlots(samplesList, annoDF, design=des, verbose=TRUE, seq.len=300, nbins=8,
           ordering=expr, ord.label="expression",plot.type="heatmap")
}
