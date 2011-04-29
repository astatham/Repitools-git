setGeneric("profilePlots", signature = "x", function(x, ...){standardGeneric("profilePlots")})

setMethod("profilePlots", "GRangesList", function(x, anno, design=NULL, up=7500, down=2500, by=100, bw=300, total.lib.size=TRUE, seq.len=NULL, verbose=FALSE, ...) {
        require(GenomicRanges)
	anno$position <- ifelse(anno$strand=="+", anno$start, anno$end)
	rownames(anno) <- anno$name
	blockPos <- seq.int(-up, down, by)
	if (verbose) message("made blockPos\n")
	annoBlocks <- data.frame(chr=rep(anno$chr, each=length(blockPos)),
		start=rep(anno$position-bw, each=length(blockPos)),
		end=rep(anno$position+bw, each=length(blockPos)),
		strand=rep(anno$strand, each=length(blockPos)))
	annoBlocks$start[annoBlocks$strand=="+"] <- annoBlocks$start[annoBlocks$strand=="+"] + blockPos
	annoBlocks$end[annoBlocks$strand=="+"] <- annoBlocks$end[annoBlocks$strand=="+"] + blockPos
	annoBlocks$start[annoBlocks$strand=="-"] <- annoBlocks$start[annoBlocks$strand=="-"] - blockPos
	annoBlocks$end[annoBlocks$strand=="-"] <- annoBlocks$end[annoBlocks$strand=="-"] - blockPos
	if (verbose) message("made annoBlocks\n")
	if (!is.null(design)) {
		stopifnot(all(design %in% c(-1,0,1)), nrow(design)==length(x))
		inUse <- !apply(design==0,1,all)
		design <- design[inUse,, drop=FALSE]
	} else inUse <- rep(TRUE, length(x))
	annoCounts <- annotationBlocksCounts(x[inUse], annoBlocks, seq.len, verbose)
	if (total.lib.size) {
		if (verbose) message("normalising to total library sizes\n")
		totalReads <- elementLengths(x[inUse])
		annoCounts <- t(t(annoCounts)/totalReads)*1000000
	}
	if (verbose) message("made annoCounts\n")
	if (!is.null(design)) {
		if (verbose) message("applying design matrix\n")
		design <- apply(design, 2, function(x) {
					x[x==1] <- 1/sum(x==1)
					x[x==-1] <- -1/sum(x==-1)
					return(x)
				})
		annoCounts <- annoCounts %*% design 
	}
	annoTable <- matrix(1:nrow(annoCounts), byrow=TRUE, ncol=length(blockPos), nrow=nrow(anno), dimnames=list(NULL, blockPos))
	if (verbose) message("made annoTable\n")
	profilePlots(annoCounts, annoTable, remove.zeros=FALSE, use.mean=TRUE, ...)
})

setMethod("profilePlots", "AffymetrixCelSet", function(x, probe.map=NULL, anno=NULL, up=7500, down=2500, by=100, bw=300, log2.adj=TRUE, verbose=FALSE, ...) {
        require(aroma.affymetrix)
			
	if (is.null(probe.map)) {
		if (is.null(anno)) stop("Either probe.map or anno must be supplied!")
		probePositions <- getProbePositionsDf( getCdf(x), verbose=verbose )
		anno$position <- ifelse(anno$strand=="+", anno$start, anno$end)
		rownames(anno) <- anno$name
		
		# run lookup twice.  first to get a list of smaller list of probes to use
		annot <- annotationLookup(probePositions, anno, up+bw, down+bw, verbose=verbose)
		pb <- unique(unlist(annot$indexes, use.names=FALSE))
		probePositions <- probePositions[pb,]
		annot <- annotationLookup(probePositions, anno, up+bw, down+bw, verbose=verbose)
		lookupT <- makeWindowLookupTable(annot$indexes, annot$offsets,
				starts = seq(-up-bw, down-bw, by), ends = seq(-up+bw, down+bw, by))
	} else {
		if (verbose) message("Using supplied probe.map\n")
		probePositions <- probe.map$probePositions
		lookupT <- probe.map$lookupT
	}
	
	dmM <- extractMatrix(x, cells = probePositions$index, verbose = verbose)
	if (log2.adj) dmM <- log2(dmM)
	
	profilePlots(dmM, lookupT, ...)
	invisible(list(lookupT=lookupT, probePositions=probePositions))
})


setMethod("profilePlots", "matrix", function(x, lookup.table, gene.list, titles=colnames(x), n.samples=1000, confidence=0.975, legend.plot="topleft", cols=rainbow(length(gene.list)), remove.zeros=TRUE, use.mean=FALSE, ...) {
	#Test gene.list for sanity
	for (i in 1:length(gene.list)) if (class(gene.list[[i]])=="logical") {
		if (length(gene.list[[i]])!=nrow(lookup.table)) 
		  stop("boolean gene.list element length must equal num of rows in lookup.table")
	} else if (class(gene.list[[i]])=="integer") {
		if(max(gene.list[[i]])>nrow(lookup.table)) 
		  stop("gene.list element value greater than num of rows in lookup.table") 
	} else stop("gene.list elements must a be boolean or integer vector")
	stopifnot(confidence>0.5, confidence<1)
	if (length(legend.plot)!=ncol(x)) if (length(legend.plot)!=1) stop("legend.plot must be either same length as columns in x or 1") else legend.plot <- rep(legend.plot, ncol(x))
	x.p <- as.numeric(colnames(lookup.table))
	#grab Intensities for all genes
	gene.list.max <- which.max(sapply(gene.list, FUN=function(u) if(class(u)=="logical") sum(u) else length(u)))
	for (i in 1:ncol(x)) {
		sMat <- .scoreIntensity(lookup.table, x[,i], removeZeros=remove.zeros, returnMatrix=TRUE, useMean=use.mean)
		sMat.gene.list <- lapply(gene.list, FUN=function(u) sMat[u, ])
		if (use.mean) trace.gene.list <- sapply(sMat.gene.list, FUN=function(u) apply(u, 2, mean, na.rm=TRUE)) else
		trace.gene.list <- sapply(sMat.gene.list, FUN=function(u) apply(u, 2, median, na.rm=TRUE))

		#choose n.samples random genelists
		inds <- lapply(seq_len(n.samples), FUN=function(u) sample(nrow(sMat), nrow(sMat.gene.list[[gene.list.max]])))
		if (use.mean) meds <- sapply(inds, FUN=function(u) apply(sMat[u,], 2, mean, na.rm=TRUE)) else
		meds <- sapply(inds, FUN=function(u) apply(sMat[u,], 2, median, na.rm=TRUE))
		meds.conf <- apply(meds, 1, quantile, p=c(1-confidence, 0.5, confidence))
	
		#plot tiem
		matplot(x.p, cbind(t(meds.conf), trace.gene.list), type="n", lty=c(2,1,2,1), lwd=c(1,3,1,3), xlab="Position relative to TSS", ylab="Signal", main=titles[i], ...)
		polygon(x=c(x.p, rev(x.p)), y=c(meds.conf[1,], rev(meds.conf[3,])), col="lightblue")
		matplot(x.p, cbind(t(meds.conf), trace.gene.list), type="l", lty=c(2,1,2,rep(1, length(gene.list))), lwd=c(1,3,1,rep(3, length(gene.list))), add=TRUE, col=c("blue","blue","blue",cols))
		if (!is.na(legend.plot[i])) legend(legend.plot[i], legend=names(gene.list), col=cols, lwd=3)
	}
})
