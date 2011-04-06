# Common short functions used by multiple function in Repitools.

setGeneric(".validate", signature = c("anno"), function(anno, up, down)
                                        {standardGeneric(".validate")})

setMethod(".validate", "GRanges", function(anno, up, down)
{
    if(is.null(anno))
        stop("Mandatory argument 'anno' not provided.")
    if(is.null(up))
        stop("Mandatory argument 'up' not provided.")
    if(is.null(down))
        stop("Mandatory argument 'down' not provided.")

    str <- strand(anno)

    if(-up > down)
	stop("'up' and 'down' window boundaries cross over each other.\n")
	
    if(any(str == '*') && any(str %in% c('+', '-')))
	stop("Annotation contains mixed feature types.")

    # For unstranded features.
    if(any(str %in% '*') && up != down)
	stop("Different upstream and downstream window distances don't make
              sense for unstranded features.\n")
    NULL
})

.getNames <- function(annoGR)
{
    if("name" %in% names(elementMetadata(annoGR)))
        return(elementMetadata(annoGR)[, "name"])
    else
        return(1:length(annoGR))
}

.makeBlocks <- function(anno, up, down)
{
    str <- as.character(strand(anno))
    f.st <- start(anno)
    f.end <- end(anno)
    if('*' %in% str) # All are unstranded.
    {
        starts <- round((f.st + f.end) / 2) - up
        ends <- round((f.st + f.end) / 2) + down

    } else { # All stranded features.
        starts <- ifelse(str == '+', f.st - up, f.end - down)
        ends <- ifelse(str == '+', f.st + down, f.end + up)
    }
    f.names <- .getNames(anno)

    GRanges(seqnames(anno), IRanges(starts, ends), name = f.names)
}

setGeneric(".annoDF2GR", signature = "anno", function(anno, ...)
                                {standardGeneric(".annoDF2GR")})

setMethod(".annoDF2GR", "data.frame", function(anno)
{
    col.missing <- setdiff(c("chr", "start", "end"), colnames(anno))
    if(length(col.missing) > 0)
	stop("Columns ", paste(col.missing, collapse = ", "),
            " of annotation are not present.")

    extra <- !colnames(anno) %in% c("chr", "start", "end", "strand")
    DF <- DataFrame(anno[, extra, drop = FALSE])    

    GRanges(anno$chr,
	    IRanges(anno$start, anno$end),
            if("strand" %in% colnames(anno)) anno$strand else '*',
            DF)
})

setGeneric(".annoGR2DF", signature = "anno", function(anno, ...)
                                {standardGeneric(".annoGR2DF")})

setMethod(".annoGR2DF", "GRanges", function(anno)
{
    annoDF <- as.data.frame(anno)
    colnames(annoDF)[1] <- "chr"
    if('*' %in% annoDF$strand)
        annoDF <- annoDF[, -match("strand", colnames(annoDF))]

    annoDF
})

