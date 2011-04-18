setGeneric("genomeBlocks", function(genome, chrs = names(genome), width,
            spacing = width) {standardGeneric("genomeBlocks")})

setMethod("genomeBlocks", "numeric", function(genome, chrs = names(genome), width,
           spacing = width)
{
    require(GenomicRanges)

    chr.windows <- lapply(chrs, function(x)
                          GRanges(seqnames = names(genome[x]),
                                 ranges = IRanges(start = seq.int(spacing / 2,
                                                                genome[x],
                                                                spacing)
                                                        - width / 2 + 1,
                                                 end = seq.int(spacing / 2,
                                                             genome[x],
                                                             spacing)
                                                      + width / 2)))
    suppressWarnings(do.call(c, chr.windows))
})

setMethod("genomeBlocks", "BSgenome",
    function(genome, chrs = seqnames(genome), width, spacing = width)
{
    require(BSgenome)
    
    chr.lengths <- seqlengths(genome)[chrs]
    genomeBlocks(chr.lengths, chrs = chrs, width = width, spacing = spacing)
})
