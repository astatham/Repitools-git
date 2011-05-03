setOldClass("AffymetrixCelSet")

setGeneric("binPlots", signature = "x", function(x, ...){standardGeneric("binPlots")})

setMethod("binPlots", "GRangesList", function(x, anno=NULL, design=NULL, up=7500, down=2500, by=100, bw=300, lib.size="lane", seq.len=NULL, verbose=FALSE, Acutoff=NULL, ...) {
        require(GenomicRanges)

        if(is.null(anno))
            stop("Annotation not given.")
	anno$position <- ifelse(anno$strand=="+", anno$start, anno$end)
	rownames(anno) <- anno$name
	if(lib.size == "ref" && is.null(Acutoff))
		stop("Must give value of Acutoff if using \"ref\" normalisation.\n")
	blockPos <- seq.int(-up, down, by)
	if (verbose) cat("made blockPos\n")
	annoBlocks <- data.frame(chr=rep(anno$chr, each=length(blockPos)),
                                 start=rep(anno$position-bw, each=length(blockPos)),
                                 end=rep(anno$position+bw, each=length(blockPos)),
                                 strand=rep(anno$strand, each=length(blockPos)))
	annoBlocks$start[annoBlocks$strand=="+"] <- annoBlocks$start[annoBlocks$strand=="+"] + blockPos
	annoBlocks$end[annoBlocks$strand=="+"] <- annoBlocks$end[annoBlocks$strand=="+"] + blockPos
	annoBlocks$start[annoBlocks$strand=="-"] <- annoBlocks$start[annoBlocks$strand=="-"] - blockPos
	annoBlocks$end[annoBlocks$strand=="-"] <- annoBlocks$end[annoBlocks$strand=="-"] - blockPos
	if (verbose) cat("made annoBlocks\n")
	if (!is.null(design)) {
		stopifnot(all(design %in% c(-1,0,1)), nrow(design)==length(x))
		inUse <- !apply(design==0,1,all)
		design <- design[inUse, , drop = FALSE]
	} else inUse <- rep(TRUE, length(x))
	annoCounts <- annotationBlocksCounts(x[inUse], annoBlocks, seq.len, verbose)

	totalReads <- elementLengths(x[inUse])
	if (lib.size == "ref") {
		if (verbose) cat("normalising to reference sample\n")
		annoCounts <- t(t(annoCounts) * calcNormFactors(annoCounts, Acutoff = Acutoff) * totalReads)
		
	} else { # lib.size = "lane"
		if (verbose) cat("normalising to total library sizes\n")
		annoCounts <- t(t(annoCounts)/totalReads)*1000000	
	}
	if (verbose) cat("made annoCounts\n")
	if (!is.null(design)) {
		if (verbose) cat("applying design matrix\n")
		design <- apply(design, 2, function(x) {
					x[x==1] <- 1/sum(x==1)
					x[x==-1] <- -1/sum(x==-1)
					return(x)
				})
		annoCounts <- annoCounts %*% design
		colnames(annoCounts) <- colnames(design)
	}
	annoTable <- matrix(1:nrow(annoCounts), byrow=TRUE, ncol=length(blockPos), nrow=nrow(anno), dimnames=list(NULL, blockPos))
	if (verbose) cat("made annoTable\n")
	binPlots(annoCounts, annoTable, remove.zeros=FALSE, use.mean=TRUE, ...)
})

setMethod("binPlots", "AffymetrixCelSet", function(x, probe.map=NULL, anno=NULL, up=7500, down=2500, by=100, bw=300, log2.adj=TRUE, verbose=FALSE, ...) {			
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
		if (verbose) cat("Using supplied probe.map\n")
		probePositions <- probe.map$probePositions
		lookupT <- probe.map$lookupT
	}
	
	dmM <- extractMatrix(x, cells = probePositions$index, verbose = verbose)
	if (log2.adj) dmM <- log2(dmM)

	binPlots(dmM, lookupT, ...)
	invisible(list(lookupT=lookupT, probePositions=probePositions))
})


setMethod("binPlots", "matrix", function(x, lookup.table=NULL, ordering=NULL, ord.label="", plot.type=c("line","heatmap","terrain","boxplot"), nbins=10, cols=NULL, lwd=3, lty=1, same.scale=TRUE, symm.scale=FALSE, verbose=FALSE, remove.zeros=TRUE, use.mean=FALSE, ...) {
    if(is.null(lookup.table))
        stop("Lookup table not given.")
    if(is.null(ordering))
        stop("Ordering not given.")    

  def.par <- par(no.readonly = TRUE) # save default, for resetting...
  plot.type <- match.arg(plot.type)
  if(!ncol(ordering) == ncol(x)) {
    if (!ncol(ordering) == 1)
      stop("ordering must have either 1 column or the same number of columns as dataMatrix.")
      orderingIndex <- rep(1, ncol(x))
  } else orderingIndex <- 1:ncol(x)

  if( is.null(cols) ) {
    require(gplots)
	if(plot.type=="line") {
	  cols <- colorpanel(nbins,"blue","green","red")
	} else {
	  cols <- colorpanel(64,"blue","white","red")
	}
  }
  
  label <- vector("list", ncol(ordering))
  .makeBins <- function(u) {
    if(class(u)=="numeric") {
      br <- quantile(u,p=(0:nbins)/nbins)
      list(breakpoints=br, intervals=cut(u,breaks=br))
	} else if(class(u)=="factor") {
	  nbins <- length(levels(u))
	  list(breakpoints=u, intervals=u)
	}
  }
  
  for(i in 1:ncol(ordering))
  {
	if(class(ordering[,i])=="numeric")
		label[[i]] <- " Order:"
	else if (class(ordering[,i])=="factor")
		label[[i]] <- " Factor:"
  }
  
  breaks <- apply(ordering, 2, .makeBins)
  if( plot.type %in% c("line","heatmap","terrain")) {
    intensScores <- array(NA,dim=c(ncol(x), ncol(lookup.table), nbins),
	                      dimnames=list(colnames(x),colnames(lookup.table),NULL))
  } else {
    intensScores <- vector("list",ncol(x))
	for(i in 1:length(intensScores))
		intensScores[[i]] <- vector("list", length(levels(breaks[[orderingIndex[i]]][["intervals"]])))
  }
  
  xval <- as.numeric(colnames(lookup.table))

  for(i in 1:ncol(x)) {
	if (verbose) cat(colnames(ordering)[orderingIndex[i]],": ",sep="")
	cutLevels <- levels( breaks[[orderingIndex[i]]][["intervals"]] )

	
	for(j in 1:length(cutLevels)){
		level <- cutLevels[j]
	    lookup.tableSubset <- lookup.table[breaks[[orderingIndex[i]]][["intervals"]]==level, ]
	  if( plot.type %in% c("line","heatmap","terrain")) {
	    intensScores[i,,j] <- .scoreIntensity(lookup.tableSubset, intensities=x[,i], minProbes=2, removeZeros=remove.zeros, useMean=use.mean)	
      } else {
		d <- .scoreIntensity(lookup.tableSubset, intensities=x[,i], minProbes=2, returnMatrix=TRUE, removeZeros=remove.zeros, useMean=use.mean)
		intensScores[[i]][[j]] <- boxplot(as.data.frame(d), plot=FALSE)
	  }
	}
  }

    if ( plot.type %in% c("line","heatmap","terrain")) if (same.scale) {
      rng <- range(intensScores, na.rm=TRUE)
      if (symm.scale) rng <- c(-max(abs(rng)),max(abs(rng)))
    }

  for(i in 1:ncol(x)) {
  cutLevels <- levels( breaks[[orderingIndex[i]]][["intervals"]] )

	  
  if(plot.type=="boxplot") {
	iS <- intensScores[[i]]
	n <- length(iS)
	df <- diff(xval)[1]
	for(j in 1:length(iS)) {
	  xvals <- as.numeric(iS[[j]]$names)
	  bxp(iS[[j]], at=xval+(j-1)*df/n,pars=list(boxwex=.7*df/n,medcol=cols[j],boxcol=cols[j],whiskcol=cols[j],outcol=cols[j]),
		  add=(j>1),show.names=(j==1),xlim=c(xvals[1]-df/2,xvals[length(xvals)]+df/2), ...)
	}
  } else {
	dm <- intensScores[i,,]
        if (!same.scale) {
          rng <- range(dm, na.rm=TRUE)
          if (symm.scale) rng <- c(-max(abs(rng)),max(abs(rng)))
        }

	titName <- paste("Signal:", colnames(x)[i], label[orderingIndex[i]], colnames(ordering)[orderingIndex[i]], sep="")
	if(plot.type=="line")
	{
		  layout(rbind(c(1, 2)), widths=c(3,2))
		  par(oma = c(0, 0, 2, 0))
		  par(mai=c(1.02,0.90,0.82,0))
		  matplot(xval,dm,type="l",col=cols,lty=lty,lwd=lwd,xlab="Position relative to TSS",ylab="Signal",ylim=rng)
		  par(mai=c(1.02,0.05,0.82,0))
		  plot.new()
		  legend(x="top", title = ord.label, col=cols, lty = 1, legend=cutLevels)
		  if (verbose) print(cols)
		  intervals <- breaks[[orderingIndex[i]]][["intervals"]]
		  if (verbose) print(intervals)
		  mtext(titName, line = 0.5, outer = TRUE)
	} else if(plot.type=="heatmap") {
		  layout(rbind(c(1,2,3)), widths=c(1,3,1))
		  par(mai=c(1.02,0.50,0.82,0.05))
		  par(oma = c(0, 0, 0, 0))
		  image(y=seq(1/nbins/2, 1-(1/nbins/2), 1/nbins),z=rbind(1:nbins), col=cols,axes=F, xlab="Signal Intensity", ylab = NA)
		  axis(2, at=(0:nbins)/nbins, labels=format(seq(rng[1], rng[2], length.out=nbins+1), digits=1))
		  par(mai=c(1.02,0.05,0.82,0.05))
		  image(xval,1:nbins,dm,xlab="Position relative to TSS", yaxt="n", ylab="Bin",col=cols,zlim=rng)
		  par(mai=c(1.02,0.05,0.82,0.50))
		  breakpoints <- breaks[[orderingIndex[i]]][["breakpoints"]]
		  plot(x=breakpoints,y=0:nbins, type="l", yaxt="n", lwd=3, xlab=ord.label, yaxs="i")
		  par(oma = c(0, 0, 2, 0))
		  mtext(titName, line = 0, outer = TRUE)
		} else if(plot.type=="terrain") {
		  layout(1)
		  par(oma = c(0, 0, 2, 0))
  		  dm.avg <- (dm[-1, -1] + dm[-1, -(ncol(dm) - 1)] +
             		dm[-(nrow(dm) -1), -1] + dm[-(nrow(dm) -1), -(ncol(dm) - 1)]) / 4

  		  this.cols = cols[cut(dm.avg, breaks = seq(rng[1], rng[2], length.out=length(cols)), include.lowest = T)] 
		  persp(xval, 1:nbins, dm, xlab="Position relative to TSS", yaxt="n", ylab=ord.label, col=this.cols, zlim=rng, theta=-25, phi=20, d=1.5, border=NA, ticktype="detailed", zlab="Signal")
		  mtext(titName, line = 0, outer = TRUE)
		}

	  }
	  par(def.par)#- reset to default
    }
})
