pkgname <- "Repitools"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('Repitools')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("BAM2GRanges")
### * BAM2GRanges

flush(stderr()); flush(stdout())

### Name: BAM2GenomicRanges
### Title: Read in a (list of) BAM file(s) into a GRanges(List) object.
### Aliases: BAM2GRanges BAM2GRangesList BAM2GRanges,character-method
###   BAM2GRangesList,character-method

### ** Examples

  tiny.BAM <- system.file("extdata", "ex1.bam", package = "Rsamtools")
  if(length(tiny.BAM) > 0)
    print(BAM2GRanges(tiny.BAM))



cleanEx()
nameEx("GDL2GRL")
### * GDL2GRL

flush(stderr()); flush(stdout())

### Name: GDL2GRL
### Title: Utility function to covert a GenomeDataList object into
###   GRangesList objects.
### Aliases: GDL2GRL GDL2GRL,GenomeDataList-method

### ** Examples

    require(BSgenome)
    gdl <- GenomeDataList(list(
                         GenomeData(list(
                                         chr1 = list(`-` = c(100, 200), `+` = c(800, 1000)),
                                         chr2 = list(`-` = c(450, 550), `+` = c(1500, 7500))
                                         )
                                    ),
                         GenomeData(list(
                                         chr1 = list(`-` = c(300, 700), `+` = c(850, 900)),
                                         chr2 = list(`-` = c(125, 250), `+` = c(500, 750))
                                         )
                                    )
                                )
                            )
    GDL2GRL(gdl)



cleanEx()
nameEx("annoDF2GR")
### * annoDF2GR

flush(stderr()); flush(stdout())

### Name: annoDF2GR
### Title: Convert a 'data.frame' to a 'GRanges'.
### Aliases: annoDF2GR annoDF2GR,data.frame-method

### ** Examples

  df <- data.frame(chr = c("chr1", "chr3", "chr7", "chr22"),
                   start = seq(1000, 4000, 1000),
                   end = seq(1500, 4500, 1000),
                   t = c(3.11, 0.93, 2.28, -0.18),
                   gc = c("High", "High", "Low", "High"))

  annoDF2GR(df)



cleanEx()
nameEx("annoGR2DF")
### * annoGR2DF

flush(stderr()); flush(stdout())

### Name: annoGR2DF
### Title: Convert an annotated 'GRanges' to a 'data.frame'.
### Aliases: annoGR2DF annoGR2DF,GRanges-method

### ** Examples

  chrs <- c("chr1", "chr3", "chr7", "chr22")
  starts <- seq(1000, 4000, 1000)
  ends <- seq(1500, 4500, 1000)
  t <- c(3.11, 0.93, 2.28, -0.18)
  gc <- c("High", "High", "Low", "High")
  gr <- GRanges(chrs, IRanges(starts, ends), strand = '*', t, gc)

  annoGR2DF(gr)



cleanEx()
nameEx("annotationBlocksCounts")
### * annotationBlocksCounts

flush(stderr()); flush(stdout())

### Name: annotationBlocksCounts
### Title: Counts the number of sequencing reads within supplied genomic
###   blocks.
### Aliases: annotationBlocksCounts
###   annotationBlocksCounts,ANY,data.frame-method
###   annotationBlocksCounts,ANY,GRanges-method

### ** Examples
#See the manual


cleanEx()
nameEx("annotationBlocksLookup")
### * annotationBlocksLookup

flush(stderr()); flush(stdout())

### Name: annotationBlocksLookup
### Title: Forms a mapping between probe locations and chromosomal blocks
###   (regions).
### Aliases: annotationBlocksLookup
###   annotationBlocksLookup,data.frame,data.frame-method
###   annotationBlocksLookup,data.frame,GRanges-method

### ** Examples

# create example set of probes and gene start sites
probeTab <- data.frame(position=seq(1000,3000,by=200), chr="chrX", strand="+")
genes <- data.frame(chr="chrX", start=c(2100,2200), end=c(2500, 2400), strand=c("+","-"))
rownames(genes) <- paste("gene",1:2,sep="")

# Call annotationLookup() and look at output
annotationBlocksLookup(probeTab, genes)



cleanEx()
nameEx("annotationCounts")
### * annotationCounts

flush(stderr()); flush(stdout())

### Name: annotationCounts
### Title: Counts the number of sequencing reads surrounding supplied
###   annotations
### Aliases: annotationCounts annotationCounts,ANY,data.frame-method
###   annotationCounts,ANY,GRanges-method

### ** Examples

    #See the manual



cleanEx()
nameEx("annotationLookup")
### * annotationLookup

flush(stderr()); flush(stdout())

### Name: annotationLookup
### Title: Forms a mapping between probes on a tiling array and windows
###   surrounding the TSSs of genes.
### Aliases: annotationLookup annotationLookup,data.frame,data.frame-method
###   annotationLookup,data.frame,GRanges-method

### ** Examples


# create example set of probes and gene start sites
probes <- data.frame(position=seq(1000, 3000, by = 200), chr = "chrX", strand = '-')
genes <- data.frame(chr = "chrX", start=c(2100, 1000), end = c(3000, 2200),
                    strand=c("+","-"))
rownames(genes) <- paste("gene", 1:2, sep = '')

# Call annotationLookup() and look at output
annotationLookup(probes, genes, 500, 500)



cleanEx()
nameEx("binPlots")
### * binPlots

flush(stderr()); flush(stdout())

### Name: binPlots
### Title: Create line plots of averaged signal across a promoter
### Aliases: binPlots binPlots,GenomeDataList-method
###   binPlots,GRangesList-method binPlots,AffymetrixCelSet-method
###   binPlots,matrix-method

### ** Examples

  annoFile <- system.file("data","chr21genes.csv", package="Repitools")
  annoDF <- read.csv(annoFile)

  readsFile <- system.file("data","samplesList.RData", package="Repitools")
  load(readsFile)  # GRangesList of reads 'samplesList'
  
  exprFile <- system.file("data","expr.RData", package="Repitools")
  load(exprFile)  # matrix of differential expression 'expr'

  # design matrix
  des <- matrix( c(0,1,0,-1), ncol=1, dimnames=list(names(samplesList),"PC-Norm") )

  binPlots(samplesList, annoDF, design=des, verbose=TRUE, seqLen=300, nbins=4,
           ordering=expr, ordLabel="expression",plotType="line")

  binPlots(samplesList, annoDF, design=des, verbose=TRUE, seqLen=300, nbins=8,
           ordering=expr, ordLabel="expression",plotType="heatmap")



cleanEx()
nameEx("blocksStats")
### * blocksStats

flush(stderr()); flush(stdout())

### Name: blocksStats
### Title: Calculate statistics for regions in the genome
### Aliases: blocksStats blocksStats,ANY,data.frame-method
###   blocksStats,ANY,GRanges-method

### ** Examples

  intensities <- matrix(c(6.8, 6.5, 6.7, 6.7, 6.9,
                          8.8, 9.0, 9.1, 8.0, 8.9), ncol = 2)
  colnames(intensities) <- c("Normal", "Cancer")
  d.matrix <- matrix(c(-1, 1))
  colnames(d.matrix) <- "Cancer-Normal"
  probe.anno <- data.frame(chr = rep("chr1", 5),
                           position = c(4000, 5100, 6000, 7000, 8000), 
                           index = 1:5)
  anno <- GRanges("chr1", IRanges(7500, 10000), '+', name = "Gene 1")
  blocksStats(intensities, anno, 2500, 2500, probe.anno, log2.adj = FALSE, design = d.matrix)



cleanEx()
nameEx("checkProbes")
### * checkProbes

flush(stderr()); flush(stdout())

### Name: checkProbes
### Title: Check Probe Specificity for Some Regions
### Aliases: checkProbes checkProbes,data.frame,data.frame-method
###   checkProbes,GRanges,GRanges-method

### ** Examples

	p.table <- data.frame(name = c("probeA", "probeB", "probeC", "probeC", "probeC"),
			    strand = c('+', '-', '+', '-', '-'),
                               chr = c("chr1", "chr2", "chr1", "chr2", "chr2"),
                             start = c(20, 276, 101, 101, 151),
                               end = c(44, 300, 125, 125, 175))
	r.table <- data.frame(name = c("gene1", "gene2", "gene3"),
                               chr = c("chr1", "chr2", "chr2"),
                            strand = c('+', '-', '+'),
                             start = c(20, 500, 75),
                               end = c(200, 800, 400))
	pdf("tmp.pdf", height = 6, width = 14)
	checkProbes(r.table, p.table, lwd = 4, col = "blue")
	dev.off()



cleanEx()
nameEx("clusterPlots")
### * clusterPlots

flush(stderr()); flush(stdout())

### Name: clusterPlots
### Title: Visualisation of tables of feature coverages.
### Aliases: clusterPlots clusterPlots,ClusteredCoverageList-method
###   clusterPlots,ScoresList-method

### ** Examples

  data.path <- system.file("data", package = "Repitools")
  load(file.path(data.path, "samplesList.RData"))
  samplesList <- samplesList[1:2]
  load(file.path(data.path, "expr.RData"))
  anno <- read.csv(file.path(data.path, "chr21genes.csv"), stringsAsFactors = FALSE)

  fs <- featureScores(samplesList, anno, up = 2000, down = 1000, freq = 500, s.width = 500)
  clusterPlots(fs, function(x) sqrt(x), expr = as.numeric(expr), plot.type = "heatmap")



cleanEx()
nameEx("cpgDensityCalc")
### * cpgDensityCalc

flush(stderr()); flush(stdout())

### Name: cpgDensityCalc
### Title: Calculate CpG Density in a Window
### Aliases: cpgDensityCalc cpgDensityCalc,GenomeDataList-method
###   cpgDensityCalc,GRangesList-method cpgDensityCalc,GRanges-method
###   cpgDensityCalc,data.frame-method

### ** Examples

require(BSgenome.Hsapiens.UCSC.hg18)
TSSTable <- data.frame(chr=paste("chr",c(1,2),sep=""), position=c(100000,200000))
cpgDensityCalc(TSSTable, organism=Hsapiens, window=600)



cleanEx()
nameEx("cpgDensityPlot")
### * cpgDensityPlot

flush(stderr()); flush(stdout())

### Name: cpgDensityPlot
### Title: Plot the distribution of sequencing reads CpG densities.
### Aliases: cpgDensityPlot cpgDensityPlot,GenomeDataList-method
###   cpgDensityPlot,GRangesList-method

### ** Examples

library(BSgenome.Hsapiens.UCSC.hg18)

readsFile <- system.file("data","samplesList.RData", package="Repitools")
load(readsFile)  # GRangesList of reads 'samplesList'

cpgDensityPlot( samplesList, seqLen=300, organism=Hsapiens, lwd=4, verbose=TRUE)



cleanEx()
nameEx("doSeqStats")
### * doSeqStats

flush(stderr()); flush(stdout())

### Name: doSeqStats
### Title: Calculate Statistics for Sequencing Data
### Aliases: doSeqStats doSeqStats,GenomeDataList-method
###   doSeqStats,GRangesList-method

### ** Examples

	# See user guide.



cleanEx()
nameEx("enrichmentCalc")
### * enrichmentCalc

flush(stderr()); flush(stdout())

### Name: enrichmentCalc
### Title: Calculate sequencing enrichment
### Aliases: enrichmentCalc enrichmentCalc,GenomeDataList-method
###   enrichmentCalc,GRangesList-method enrichmentCalc,GRanges-method

### ** Examples

library(BSgenome.Hsapiens.UCSC.hg18)

readsFile <- system.file("data","samplesList.RData", package="Repitools")
load(readsFile)  # GRangesList of reads 'samplesList'

tc <- enrichmentCalc(samplesList, seqLen=300, organism=Hsapiens)



cleanEx()
nameEx("enrichmentPlot")
### * enrichmentPlot

flush(stderr()); flush(stdout())

### Name: enrichmentPlot
### Title: Plot the distribution of sequencing enrichment.
### Aliases: enrichmentPlot enrichmentPlot,GenomeDataList-method
###   enrichmentPlot,GRangesList-method

### ** Examples

library(BSgenome.Hsapiens.UCSC.hg18)

## Not run: 
##D readsFile <- system.file("data","samplesList.Rdata", package="Repitools")
##D load(readsFile)  # GRangesList of reads 'samplesList'
##D 
##D enrichmentPlot(samplesList, seq.len=300, organism=Hsapiens, total.lib.size=FALSE)
## End(Not run)



cleanEx()
nameEx("featureBlocks")
### * featureBlocks

flush(stderr()); flush(stdout())

### Name: featureBlocks
### Title: Make windows for distances around a reference point.
### Aliases: featureBlocks featureBlocks,data.frame-method
###   featureBlocks,GRanges-method

### ** Examples

  genes <- data.frame(chr = c("chr1", "chr3", "chr7", "chr22"),
                   start = seq(1000, 4000, 1000),
                   end = seq(1500, 4500, 1000),
                   strand = c('+', '-', '-', '+'))

  featureBlocks(genes, 500, 500)



cleanEx()
nameEx("featureScores")
### * featureScores

flush(stderr()); flush(stdout())

### Name: featureScores
### Title: Get scores at regular sample points around genomic features.
### Aliases: featureScores featureScores,ANY,data.frame-method
###   featureScores,ANY,GRanges-method

### ** Examples

  data.path <- system.file("data", package = "Repitools")
  load(file.path(data.path, "samplesList.RData"))
  samplesList <- samplesList[1:2]
  load(file.path(data.path, "expr.RData"))
  anno <- read.csv(file.path(data.path, "chr21genes.csv"), stringsAsFactors = FALSE)

  fs <- featureScores(samplesList, anno, up = 2000, down = 1000, freq = 500, s.width = 500)



cleanEx()
nameEx("findClusters")
### * findClusters

flush(stderr()); flush(stdout())

### Name: findClusters
### Title: Find Clusters Epigenetically Modified Genes
### Aliases: findClusters

### ** Examples

  chrs <- sample(paste("chr", c(1:5), sep = ""), 500, replace = TRUE)
  starts <- sample(1:10000000, 500, replace = TRUE)
  ends <- starts + 10000
  genes <- data.frame(chr = chrs, start = starts, end = ends, strand = '+')
  genes <- genes[order(genes$chr, genes$start), ]
  genes$t.stat = rnorm(500, 0, 2)
  genes$t.stat[21:30] = rnorm(10, 4, 1)
  findClusters(genes, 5, 5, 2, 3, seq(1, 10, 1), trend = "up", n.perm = 2)



cleanEx()
nameEx("gcContentCalc")
### * gcContentCalc

flush(stderr()); flush(stdout())

### Name: gcContentCalc
### Title: Calculate The gcContent of a Region
### Aliases: gcContentCalc gcContentCalc,GRanges-method
###   gcContentCalc,data.frame-method

### ** Examples

require(BSgenome.Hsapiens.UCSC.hg18)
TSSTable <- data.frame(chr = paste("chr", c(1,2), sep = ""), position = c(100000, 200000))
gcContentCalc(TSSTable, 200, organism=Hsapiens)



cleanEx()
nameEx("genQC")
### * genQC

flush(stderr()); flush(stdout())

### Name: genQC
### Title: Plot Quality Checking Information for Sequencing Data
### Aliases: genQC genQC,character-method genQC,SequenceQCSet-method

### ** Examples

  ## Not run: 
##D     qc.files <- list.files(qc.dir, "QC.*RData", full.names = TRUE)
##D     genQC(qc.files, "My Simple Experiment")
##D   
## End(Not run)



cleanEx()
nameEx("genomeBlocks")
### * genomeBlocks

flush(stderr()); flush(stdout())

### Name: genomeBlocks
### Title: Creates bins across a genome.
### Aliases: genomeBlocks genomeBlocks,numeric-method
###   genomeBlocks,BSgenome-method

### ** Examples

  chr.lengths <- c(800, 200, 200)
  names(chr.lengths) <- c("chr1", "chr2", "chr3")
  genomeBlocks(chr.lengths, width = 200)



cleanEx()
nameEx("getCN")
### * getCN

flush(stderr()); flush(stdout())

### Name: getCN
### Title: Calculate Copy Number and Map To Enriched Regions
### Aliases: getCN getCN,data.frame,data.frame-method
###   getCN,GRanges,GRanges-method

### ** Examples

  inputs <- data.frame(chr = c("chr1", "chr1", "chr1", "chr2", "chr2"),
                     start = c(1, 50001, 100001, 1, 10001),
                       end = c(50000, 100000, 150000, 10000, 20000))
  enriched <- data.frame(chr = c("chr1", "chr1", "chr1", "chr1", "chr1", "chr2", "chr2"),
                       start = c(1, 1001, 2001, 10000001, 10005001, 9001, 9501),
                         end = c(1000, 2000, 3000, 10001000, 10006000, 9500, 10000))
  counts <- matrix(c(25, 39, 3, 10, 22, 29, 38, 5, 19, 31), nrow = 5)
  colnames(counts) <- c("Control", "Treatment")
  getCN(inputs, enriched, counts, p.method = "perm")



cleanEx()
nameEx("getProbePositionsDf")
### * getProbePositionsDf

flush(stderr()); flush(stdout())

### Name: getProbePositionsDf
### Title: Translate Affymetrix probe information in a table.
### Aliases: getProbePositionsDf
###   getProbePositionsDf,AffymetrixCdfFile-method

### ** Examples

## not run
# probePositions <- getProbePositionsDf(cdfU)



cleanEx()
nameEx("loadPairFile")
### * loadPairFile

flush(stderr()); flush(stdout())

### Name: loadPairFile
### Title: A routine to read Nimblegen tiling array intensities
### Aliases: loadPairFile

### ** Examples

# Not run
#
## Read in the NDF file 
# ndfAll <- processNDF("080310_HG18_chr7RSFS_AS_ChIP.ndf")
#
## Subset the NDF to only probes against chromosomes
# ndf <- ndfAll[grep("^chr", ndfAll$chr),]
#
## Read in a pair file using the chromosome only NDF
# arrayIntensity <- loadPairFile("Pairs/Array1_532.pair", ndf)
#



cleanEx()
nameEx("loadSampleDirectory")
### * loadSampleDirectory

flush(stderr()); flush(stdout())

### Name: loadSampleDirectory
### Title: A routine to read Nimblegen tiling array intensities
### Aliases: loadSampleDirectory

### ** Examples

# Not run
#
## Read in the NDF file 
# ndfAll <- processNDF("080310_HG18_chr7RSFS_AS_ChIP.ndf")
#
## Subset the NDF to only probes against chromosomes
# ndf <- ndfAll[grep("^chr", ndfAll$chr),]
#
## Read in a directory of pair files, returning both the Cy3 and Cy5 fluorescence in separate columns
# arrayIntensities <- loadSampleDirectory("Arrays", ndf, what="Cy3andCy5")
#



cleanEx()
nameEx("makeWindowLookupTable")
### * makeWindowLookupTable

flush(stderr()); flush(stdout())

### Name: makeWindowLookupTable
### Title: Using the output of 'annotationLookup', create a tabular storage
###   of the indices
### Aliases: makeWindowLookupTable

### ** Examples


# create example set of probes and gene start sites
probeTab <- data.frame(position=seq(1000,3000,by=200), chr="chrX", strand = '-')
genes <- data.frame(chr="chrX", start=c(2100, 1000), end = c(3000, 2200), strand=c("+","-"))
rownames(genes) <- paste("gene",1:2,sep="")

# Call annotationLookup() and look at output
aL <- annotationLookup(probeTab, genes, 500, 500)
print(aL)

# Store the results of annotationLookup() in a convenient tabular format
lookupTab <- makeWindowLookupTable(aL$indexes, aL$offsets, starts=seq(-400,200,by=200), ends=seq(-200,400,by=200))
print(lookupTab)




cleanEx()
nameEx("mappabilityCalc")
### * mappabilityCalc

flush(stderr()); flush(stdout())

### Name: mappabilityCalc
### Title: Calculate The Mappability of a Region
### Aliases: mappabilityCalc mappabilityCalc,GRanges-method
###   mappabilityCalc,data.frame-method

### ** Examples

# require(BSgenome.Hsapiens36bp.UCSC.hg18mappability)
# TSSTable <- data.frame(chr = paste("chr", c(1,2), sep = ""), position = c(100000, 200000))
# mappabilityCalc(TSSTable, 200, organism=Hsapiens36bp)



cleanEx()
nameEx("mergeReplicates")
### * mergeReplicates

flush(stderr()); flush(stdout())

### Name: mergeReplicates
### Title: Merge GRanges that are of replicate experiments.
### Aliases: mergeReplicates mergeReplicates,GenomeDataList-method
###   mergeReplicates,GRangesList-method

### ** Examples

  grl <- GRangesList(GRanges("chr1", IRanges(5, 10)),
                     GRanges("chr18", IRanges(25, 50)),
                     GRanges("chr22", IRanges(1, 100)))
  antibody <- c("MeDIP", "MeDIP", "H3K4me3")
  mergeReplicates(grl, antibody)



cleanEx()
nameEx("multiHeatmap")
### * multiHeatmap

flush(stderr()); flush(stdout())

### Name: multiHeatmap
### Title: Superfigure plots
### Aliases: multiHeatmap

### ** Examples

library(gplots)

cL <- NULL
br <- seq(-3,3,length=101)
col <- colorpanel(low="blue",mid="grey",high="red",n=100)
cL[[1]] <- list(breaks=br,colors=col)
br <- seq(-2,2,length=101)
col <- colorpanel(low="green",mid="black",high="red",n=100)
cL[[2]] <- list(breaks=br,colors=col)
br <- seq(0,20,length=101)
col <- colorpanel(low="black",mid="grey",high="white",n=100)
cL[[3]] <- list(breaks=br,colors=col)

testD <- list(matrix(runif(400),nrow=20),matrix(rnorm(100),nrow=20),matrix(rpois(100,lambda=10),nrow=20))
colnames(testD[[1]]) <- letters[1:20]
rownames(testD[[1]]) <- paste("row",1:20,sep="")

multiHeatmap(testD,cL,xspace=1)



cleanEx()
nameEx("plotClusters")
### * plotClusters

flush(stderr()); flush(stdout())

### Name: plotClusters
### Title: Plot Scores of Cluster Regions
### Aliases: plotClusters plotClusters,data.frame-method
###   plotClusters,GRanges-method

### ** Examples

    g.summary <- GRanges("chr1",
                         IRanges(seq(1000, 10000, 1000), width = 100),
                         rep(c('+', '-'), 5),
                         `t-statistic` = rnorm(10, 8, 2),
                         cluster = c(0, 0, 0, 0, 0, 1, 1, 1, 1, 0),
                         name = paste("Gene", 1:10))
    plotClusters(g.summary, 1, 0, ylim = c(4, 12), lwd = 5)



cleanEx()
nameEx("processNDF")
### * processNDF

flush(stderr()); flush(stdout())

### Name: processNDF
### Title: Reads in a Nimblegen microarray design file (NDF)
### Aliases: processNDF

### ** Examples

# Not run
#
## Read in the NDF file 
# ndfAll <- processNDF("080310_HG18_chr7RSFS_AS_ChIP.ndf")
#
## Subset the NDF to only probes against chromosomes
# ndf <- ndfAll[grep("^chr", ndfAll$chr),]



cleanEx()
nameEx("regionStats")
### * regionStats

flush(stderr()); flush(stdout())

### Name: regionStats
### Title: Find Regions of significance in microarray data
### Aliases: regionStats regionStats,AffymetrixCelSet-method
###   show,RegionStats-method RegionStats-class class:RegionStats
###   regionStats,matrix-method

### ** Examples

## Not run: 
##D library(Repitools)
##D library(aroma.affymetrix)
##D 
##D # assumes appropriate files are at annotationData/chipTypes/Hs_PromPR_v02/
##D cdf <- AffymetrixCdfFile$byChipType("Hs_PromPR_v02",verbose=-20)
##D cdfU <- getUniqueCdf(cdf,verbose=-20)
##D 
##D # assumes appropriate files are at rawData/experiment/Hs_PromPR_v02/
##D cs <- AffymetrixCelSet$byName("experiment",cdf=cdf,verbose=-20)
##D mn <- MatNormalization(cs)
##D csMN <- process(mn,verbose=-50)
##D csMNU <- convertToUnique(csMN,verbose=-20)
##D 
##D #> getNames(cs)
##D # [1] "samp1"  "samp2"  "samp3"  "samp4"
##D 
##D design <- matrix( c(1,-1,rep(0,length(cs)-2)), ncol=1, dimnames=list(getNames(cs),"elut5_L-P") )
##D 
##D # just get indices of chr7 here
##D ind <- getCellIndices(cdfU, unit = indexOf(cdfU, "chr7F"), unlist = TRUE, useNames = FALSE)
##D 
##D regs <- regionStats(csMNU, design, ind = ind, probeWindow = 500, verbose = TRUE)
## End(Not run)



cleanEx()
nameEx("sequenceCalc")
### * sequenceCalc

flush(stderr()); flush(stdout())

### Name: sequenceCalc
### Title: Find occurences of a DNA pattern
### Aliases: sequenceCalc sequenceCalc,GRanges-method
###   sequenceCalc,data.frame-method

### ** Examples

require(BSgenome.Hsapiens.UCSC.hg18)
TSSTable <- data.frame(chr=paste("chr",c(1,2),sep=""), position=c(100000,200000))
sequenceCalc(TSSTable, 600, organism=Hsapiens, pattern=DNAString("CG"))



cleanEx()
nameEx("significancePlots")
### * significancePlots

flush(stderr()); flush(stdout())

### Name: significancePlots
### Title: Create line plots of averaged signal across a promoter compared
###   to random sampling
### Aliases: significancePlots significancePlots,GenomeDataList-method
###   significancePlots,GRangesList-method
###   significancePlots,AffymetrixCelSet-method
###   significancePlots,matrix-method

### ** Examples

  # See examples in manual.



cleanEx()
nameEx("writeWig")
### * writeWig

flush(stderr()); flush(stdout())

### Name: writeWig
### Title: Writes sequencing data out into a wiggle files
### Aliases: writeWig writeWig,AffymetrixCelSet-method
###   writeWig,GenomeDataList-method writeWig,GRangesList-method

### ** Examples

#See examples in the manual



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
