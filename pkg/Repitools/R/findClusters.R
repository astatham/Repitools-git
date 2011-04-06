.makeClusters <- function(scores.chr, w.size, n.med, n.consec, cut, count = FALSE,
                          verbose)
{
    mergeOverlaps <- function(query, subject)
    {
        candidates <- IRanges(slice(coverage(c(query, subject)), lower = 1))
        candidates[countOverlaps(candidates, query) > 0]
    }

    clusterScores <- function(x)
    {
        med.scores <- rollmedian(x, k = w.size, na.pad = TRUE)
        med.cut <- med.scores > cut
        med.cut[is.na(med.cut)] <- FALSE
        med.cut.n <- rollmean(as.integer(med.cut), k = w.size, na.pad = TRUE) > (n.med / w.size)
        med.cut.n[is.na(med.cut.n)] <- FALSE
        score.cut <- x > 0

        # Window extension.
        mediansCut <- med.cut & med.cut.n
                
        cl.candidates <- mergeOverlaps(IRanges(mediansCut), IRanges(score.cut))
                
	# Minimum consecutive genes above score cutoff.
        cl.final <- cl.candidates[width(cl.candidates) >= n.consec]
		
        if(count == TRUE)
	    length(cl.final)
	else
	    as.numeric(coverage(cl.final, width = length(x)))
    }

    if(count == TRUE)
    	apply(scores.chr, 2, clusterScores)
    else
    	clusterScores(scores.chr[, 1])
}

findClusters <- function(stats, score.col, w.size = 5, n.med, n.consec, cut.samps,
                         maxFDR = 0.05, trend = c("down", "up"), n.perm = 100,
                         getFDRs = FALSE, verbose = TRUE)
{
    require(IRanges)
    require(zoo)
    trend <- match.arg(trend)
    
    # Do check in case users pass in some rows with NA scores.
    stats <- stats[!is.na(stats[, score.col]), ]

    # Simplifies processing in .makeClusters.
    if(trend == "down")
        scores <- -stats[, score.col]
    else
        scores <- stats[, score.col]
    cut.abs <- abs(cut.samps)

    perms <- 1:n.perm
    perm.scores <- sapply(perms, function(x) scores[sample(nrow(stats))])
    # Column 1 : Real Scores. Columns following : Permuted scores.
    chr.scores <- split(data.frame(scores, perm.scores), stats$chr)

    # Find the number of clusters at each cutoff, in the real and permuted data.
    clusts <- lapply(cut.abs, function(x){
                     if(verbose == TRUE) message("Counting clusters at cutoff ", x)
                     n.clusts.chrs <- lapply(chr.scores, .makeClusters, w.size,
                                             n.med, n.consec, x, TRUE, verbose)
                     n.clusts <- colSums(do.call(rbind, n.clusts.chrs))
                     
                     # Element 1 : Number of clusters in real data.
                     # Element 2 : Average number of clusters in permuted scores.
                     return(c(n.clusts[1], round(mean(n.clusts[2:length(n.clusts)]))))
                     })

    # Find the FDR of each cutoff tried.
    allFDR <- sapply(clusts, function(x)
                    {
                        FDR <- x[2] / x[1]
                        if((is.nan(FDR))) # 0 / 0 case
                            FDR = 0
                        FDR
                    })

    FDRtable <- data.frame(cutoff = cut.samps, FDR = allFDR)
    best.cut <- cut.abs[match(TRUE, allFDR < maxFDR)]
    
    if(verbose == TRUE)
        message("Using the cutoff ", ifelse(trend == "up", best.cut, -best.cut),
                " for a FDR of < ", maxFDR)

    in.clust <- unlist(lapply(chr.scores,
                              function(x) .makeClusters(x, w.size, n.med, n.consec,
                                                        best.cut)
                      ))

    # Join adjoining clusters, keeping consecutive indexing.
    cl.index = 1
    for(i in 2:length(in.clust))
    {
    	if(in.clust[i - 1] == 0 && in.clust[i] > 0)
	{
		in.clust[i] = cl.index
	} else if(in.clust[i - 1] > 0 && in.clust[i] > 0)
	{
		in.clust[i] = cl.index
	} else if((in.clust[i - 1] > 0 && in.clust[i] == 0))
	{
		cl.index = cl.index + 1
	}
    }

    stats$cluster <- in.clust

    if(getFDRs == TRUE)
    	list(table = stats, FDRs = FDRtable)
    else
    	stats
}
