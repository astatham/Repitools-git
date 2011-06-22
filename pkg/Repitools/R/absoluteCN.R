setGeneric("absoluteCN", function(input.windows, ...){standardGeneric("absoluteCN")})

setMethod("absoluteCN", c("GRanges"),
    function(input.windows, input.counts = NULL, gc.params = NULL, verbose = TRUE)
{
    n.bins <- gc.params@n.bins

    # Find which counting windows have sufficient mappability.
    input.win.mappability <- mappabilityCalc(input.windows, gc.params@mappability) * 100
    mappable <- input.win.mappability > gc.params@min.mappability
    input.windows <- input.windows[mappable, ]
    input.counts <- input.counts[mappable, ]
    input.counts <- input.counts * 100 / input.win.mappability[mappable]
    abs.CN <- apply(input.counts, 2, function(x) x / median(x))

    # Get the GC content of windows.
    gc <- gcContentCalc(input.windows, gc.params@genome)

    # Break GC content into bins, and find mode of bins, also including adjacent bins.
    gc.range <- range(gc)
    bins <- seq(gc.range[1], gc.range[2], length.out = n.bins) # Midpoint of each bin.
    mode.gcs <- apply(abs.CN, 2, function(x)
    {
        modes <- rep(NA, n.bins)
        for(index in 2:(gc.params@n.bins-1))
        {
            in.bins <- which(gc >= bins[index - 1] & gc < bins[index + 1])
            count.dens <- density(x[in.bins])
            modes[index] <- as.numeric(count.dens$x[which.max(count.dens$y)])
        }
        modes
    })

    # Find expected copy number for each window, based on its GC content.
    model.CN <- apply(mode.gcs, 2, function(x)
    {
        counts.model <- lm(x ~ poly(bins, gc.params@poly.degree))
        predict(counts.model, data.frame(bins = gc))
    })

    # Adjust the real counts by dividing by expected counts.
    t(t(abs.CN / model.CN) * gc.params@ploidy)
})

setMethod("absoluteCN", c("data.frame"),
    function(input.windows, ...)
{
    absoluteCN(annoDF2GR(input.windows), ...)
})
