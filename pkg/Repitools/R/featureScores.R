setGeneric("featureScores", signature = c("x", "anno"), function(x, anno, ...)
                                           {standardGeneric("featureScores")})
setGeneric(".featureScores", signature = c("x", "y"), function(x, y, ...)
                                           {standardGeneric(".featureScores")})

setClassUnion(".SequencingData", c("character", "GenomeDataList", "GRanges",
                                  "GRangesList"))

setMethod(".featureScores", c("GRanges", ".CoverageSamples"),
    function(x, y, anno, up, down, dist, freq, s.width, verbose)
{
    # Unpack variables in y.
    pos.labels <- y@pos.labels
    cvg.samps <- y@cvg.samps
    max.out <- y@max.out
    chr.ord <- y@chr.ord
    anno.chr <- y@anno.chr
    old.ord <- y@old.ord

    # Only use sequencing data on annotated chromosomes.
    x <- x[seqnames(x) %in% seqlevels(anno)]
    seqlevels(x) <- seqlevels(anno)

    # Qualitatively near identical to running mean smoothing.
    if(verbose == TRUE) message("Extending all reads to smoothing width.")
    seqlengths(x) <- rep(NA, length(seqlengths(x)))
    x <- resize(x, s.width)

    # Infer chromosome end positions from feature annotations.
    # This means the user doesn't have to have a BSgenome object of their samples'
    # genome installed on their computer.
    max.anno <- seqapply(split(ranges(anno.chr), seqnames(anno.chr)),
                         function(z) max(end(z)))
    max.reads <- suppressWarnings(
                 seqapply(split(ranges(x), seqnames(x)), function(z) max(end(z))))
    max.reads <- max.reads[names(max.anno)]
    chr.lens <- mseqapply(max, max.anno, max.reads)
    gc()

    # Remove reads on chromosomes not in annotation.
    # Shift all reads and sampling positions right by the out-of-bounds upper bound,
    # and then adjust chr lengths accordingly.
    x <- shift(x, max.out)
    seqlengths(x) <- as.numeric(chr.lens) + 2 * max.out
    gc()

    # Get coverage.
    if(verbose == TRUE) message("Calculating coverage.")
    # Scale coverages for total reads.
    cvg <- coverage(x) / length(x)
    samp.chr <- as(cvg.samps, "RangesList")
    gc()

    # Do sampling, per chromosome.
    if(verbose == TRUE) message("Sampling coverage.")
    cvg.mat <- lapply(names(samp.chr), function(z)
	           {
		       inds <- samp.chr[[z]]
                       chr.mat <- cvg[[z]][inds, drop = TRUE]
		       matrix(chr.mat, ncol = length(pos.labels), byrow = TRUE)		   	
		    })
    cvg.mat <- do.call(rbind, cvg.mat)
    	
    colnames(cvg.mat) <- pos.labels
    if("name" %in% names(elementMetadata(anno.chr))) 
	rownames(cvg.mat) <- elementMetadata(anno.chr)[, "name"]
    else
	rownames(cvg.mat) <- 1:nrow(cvg.mat)
    cvg.mat <- cvg.mat[old.ord, ]

    new("ScoresList", names = "Undefined", scores = list(cvg.mat), anno = anno,
         up = up, down = down, dist = dist, freq = freq, s.width = s.width)
})

setMethod(".featureScores", c("GRangesList", ".CoverageSamples"),
    function(x, y, anno, up, down, dist, freq, s.width, verbose)
{
    if(length(s.width) == 1)
        s.width <- rep(s.width, length(x))
    scores <- mapply(function(z, i)
	           {
                        if(verbose == TRUE && !is.null(names(x)))
                            message("Processing sample ", names(x)[i])
		   	.featureScores(z, y, anno, up, down, dist,
                                        freq, s.width[i], verbose)
		   }, x, IntegerList(as.list(1:length(x))), SIMPLIFY = FALSE)

    if(!is.null(names(x)))
	names <- names(x)
    else
	names <- unname(sapply(scores, names))
    new("ScoresList", names = names, anno = anno, scores = unname(sapply(scores, tables)),
                            up = up, down = down, dist = dist, freq = freq,
                            s.width = s.width)
})

setMethod(".featureScores", c("character", ".CoverageSamples"),
    function(x, y, anno, up, down, dist, freq, s.width, verbose)
{
    if(length(s.width) == 1)
        s.width <- rep(s.width, length(x))

    scores <- mapply(function(y, i)
	           {
                        if(verbose == TRUE && !is.null(names(x)))
                            message("Processing sample ", names(x)[i])
		   	readsGR <- BAM2GRanges(y)
		   	.featureScores(readsGR, y, anno, up, down, dist, freq,
                                        s.width[i], verbose)
		   }, x, 1:length(x), SIMPLIFY = FALSE)
    if(!is.null(names(x)))
	names <- x
    else
	names <- unname(sapply(scores, names))
    new("ScoresList", names = names, anno = anno, scores = unname(sapply(scores, tables)),
                            up = up, down = down, dist = dist, freq = freq,
                            s.width = s.width)
})

setMethod(".featureScores", c("GenomeDataList", ".CoverageSamples"),
    function(x, y, ...)
{
    .featureScores(GDL2GRL(x), y, ...)
})

setMethod(".featureScores", c(".SequencingData", "GRanges"),
    function(x, y, up, down, dist = c("percent", "base"), freq, s.width,
             verbose = TRUE)
{
    if(is.null(s.width))
        stop("Mandatory argument 's.width' not provided.")

    dist <- match.arg(dist)

    str <- strand(y)
    st <- start(y)
    en <- end(y)	
    wd <- width(y) 
    pos <- str == '+'

    if(verbose == TRUE) message("Calculating sampling positions.")	
    if(all(str == '*'))
        ref.points <- round((st + en) / 2)
    else
        ref.points <- ifelse(pos, st, en)

    if(dist == "percent")
    {
        starts = as.numeric(ifelse(pos, ref.points - wd * up/100,
                                        ref.points - wd * down/100))
        ends = as.numeric(ifelse(pos, ref.points + wd * down/100,
                                        ref.points + wd * up/100))
    } else {
	starts = as.numeric(ifelse(pos, ref.points - up,
                                        ref.points - down))
	ends = as.numeric(ifelse(pos, ref.points + down,
                                      ref.points + up))
    }

    cov.winds <- GRanges(seqnames = seqnames(y),
                         IRanges(start = starts, end = ends),
                         strand = str)

    posns <- seq(-up, down, freq)
    if(dist == "percent")
	pos.labels <- paste(posns, '%')
    else
	pos.labels <- posns
    n.pos <- length(posns)

    # Make ranges for each sample point.
    cvg.samps <- rep(cov.winds, each = n.pos)
    if(dist == "percent")
        gap.size <- width(cvg.samps) / (n.pos - 1)
    else
        gap.size <- rep(freq, length(cvg.samps))

    ranges(cvg.samps) <- IRanges(start = as.numeric(
                              ifelse(strand(cvg.samps) %in% c('+', '*'), 
                                     start(cvg.samps) + 0:(n.pos - 1) * gap.size,
                                     end(cvg.samps) - 0:(n.pos - 1) * gap.size)
                                                   ),
                                     width = 1)

    # Find upper bound of how far a sampling position could be
    # outside of a chromosome.
    max.out <- ifelse(dist == "percent", max(width(y)) * 
                                         max(abs(up), abs(down)) / 100,
                                         max(abs(up), abs(down)))
    cvg.samps <- shift(cvg.samps, max.out)

    # Find order to get back from RangesList order to original order.
    chr.ord <- order(as.character(seqnames(y)))
    anno.chr <- y[chr.ord]
    old.ord <- order(chr.ord)

    samp.info <- new(".CoverageSamples", pos.labels = pos.labels, cvg.samps = cvg.samps,
                     max.out = max.out, chr.ord = chr.ord, anno.chr = anno.chr,
                     old.ord = old.ord)

    .featureScores(x, samp.info, y, up, down, dist, freq, s.width, verbose)
})

setMethod(".featureScores", c("matrix", "GRanges"),
    function(x, y, up, down, p.anno, mapping, freq, log2.adj = TRUE, verbose = TRUE)
{
    if("index" %in% colnames(x)) p.inds <- x$index else p.inds <- 1:nrow(x)    
    ind.col <- colnames(p.anno) == "index"
    mapping <- annotationLookup(p.anno[, !ind.col], y, up, down, verbose)
    p.used <- unique(unlist(mapping$indexes, use.names = FALSE))
    p.anno <- p.anno[p.used, ]
    mapping <- annotationLookup(p.anno[, !ind.col], y, up, down, verbose)
    
    intens <- x[p.anno$index, ]
    if(log2.adj == TRUE) intens <- log2(intens)

    points.probes <- makeWindowLookupTable(mapping$indexes, mapping$offsets,
                     starts = seq(-up, down - freq, freq),
                     ends = seq(-up + freq, down, freq))

    points.intens <- lapply(1:ncol(x), function(z)
                     {
                         .scoreIntensity(points.probes, intens[, z], returnMatrix = TRUE)
                     })    

    new("ScoresList", names = colnames(x), anno = y, scores = points.intens,
                            up = up, down = down, dist = NULL, freq = freq,
                            s.width = NULL)
})

setMethod(".featureScores", c("AffymetrixCelSet", "GRanges"),
    function(x, y, p.anno = NULL, mapping = NULL, ...)
{
    if(is.null(mapping) && is.null(p.anno))
        p.anno <- getProbePositionsDf(getCdf(x), verbose = verbose)

    intens <- extractMatrix(x, cells = p.anno$index, verbose = verbose)
    p.anno$index <- 1:nrow(p.anno)

    .featureScores(intens, y, p.anno = p.anno, mapping = mapping, ...)
})

setMethod("featureScores", c("ANY", "GRanges"), function(x, anno,
           up, down, ...)
{
    invisible(.validate(anno, up, down))
    .featureScores(x, anno, up = up, down = down, ...)
})

setMethod("featureScores", c("ANY", "data.frame"),
    function(x, anno, ...)
{
    featureScores(x, .annoDF2GR(anno), ...)
})
