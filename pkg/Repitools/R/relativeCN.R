setGeneric("relativeCN", function(input.windows, ip.windows, ...){standardGeneric("relativeCN")})

setMethod("relativeCN", c("GRanges", "GRanges"),
    function(input.windows, ip.windows, input.counts = NULL, gc.params = NULL, ..., verbose = TRUE)
{
    if(is.null(input.counts))
        stop("Matrix of counts for input types not given.")

    if(length(input.windows) != nrow(input.counts))
        stop("Rows of counts differ to rows of input regions.\n")
    
    require(GenomicRanges)
    require(DNAcopy)

    if(!is.null(gc.params)) # Do mappability / GC bias adjustment on counts.
        input.scores <- gcMappabilityAdjust(input.windows, input.counts, pc.params, verbose)
    else
        input.scores <- input.counts

    Mvalues <- log2((input.scores[, 2] / sum(input.scores[, 2])) /
                    (input.scores[, 1] / sum(input.scores[, 1])))
    cn <- CNA(chrom = as.character(seqnames(input.windows)),
             maploc = as.numeric(start(input.windows)),
          data.type = "logratio",
           genomdat = Mvalues,
           sampleid = paste(colnames(input.scores)[2], "/",
                            colnames(input.scores)[1], "Fold Change"))
    totals <- colSums(input.scores)
    wts <- ((totals[2] - input.scores[, 2]) / (input.scores[, 2] * totals[2])
          + (totals[1] - input.scores[, 1]) / (input.scores[, 1] * totals[1]))^-1
    non0 <- which(wts > 0)
    cn <- cn[non0, ]
    wts <- wts[non0]
    if(verbose == TRUE) message("Segmenting and smoothing.")
    cn <- segment(smooth.CNA(cn), weights = wts, ..., verbose = 0)
    # Extend CNV region to the end of the interval, since all positions are starts.
    cn$out[, "loc.end"] <- cn$out[, "loc.end"] + width(input.windows)[1]
    CNV.windows <- GRanges(cn$out[, "chrom"],
                           IRanges(cn$out[, "loc.start"], cn$out[, "loc.end"]))
	
    if(verbose == TRUE) message("Mapping copy number to IP windows.")
    map <- findOverlaps(ip.windows, CNV.windows, select = "first")
							   
    relative.cn <- 2^cn$out[map, "seg.mean"]
    # If CN ratio not called, assume it is 1.
    relative.cn[is.na(relative.cn)] <- 1
    names(relative.cn) <- .getNames(ip.windows)

    relative.cn
})

setMethod("relativeCN", c("data.frame", "data.frame"),
    function(input.windows, ip.windows, input.counts = NULL, gc.params = NULL, ..., verbose = TRUE)
{
    relativeCN(annoDF2GR(input.windows), annoDF2GR(ip.windows), input.counts, gc.params,
               ..., verbose = verbose)
})
