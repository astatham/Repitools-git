setGeneric("clusterPlots", function(c.list, ...){standardGeneric("clusterPlots")})

setMethod("clusterPlots", "ClusteredCoverageList",
	function(c.list, plot.ord = 1:length(c.list), plot.type = c("line", "heatmap"),
                 symm.scale = NULL, cols = NULL, t.name = NULL,
                 verbose = TRUE, ...)
{
    c.list <- c.list[plot.ord]
    plot.type <- match.arg(plot.type)
    if(is.null(symm.scale))
    {
        if(plot.type == "line")
            symm.scale <- rep(FALSE, length(unique(clusters(c.list))))
        else
            symm.scale <- rep(FALSE, length(c.list))
    } else if(length(symm.scale) == 1 && length(c.list) > 1) {
        if(plot.type == "line")
            symm.scale <- rep(symm.scale, length(unique(clusters(c.list))))
        else
            symm.scale <- rep(symm.scale, length(c.list))
    }

    n.marks <- length(c.list)

    if(is.null(cols) == TRUE)
    {
	require(gplots)
	if(plot.type == "line")
	{
	    cols <- colorpanel(n.marks, "blue", "green", "red")
	} else {
	    cols <- colorpanel(100, "blue", "white", "red")
	}
    }

    cvgs <- tables(c.list)

    # Get median expression for each cluster. Find ascending order.
    expr <- c.list@expr
    expr.name <- c.list@expr.name
    cl.id <- c.list@cluster.id
    n.clusters <- length(unique(cl.id))
    if(!is.null(expr))
    {
        cl.expr <- tapply(expr, factor(cl.id), median, na.rm = TRUE)
        cl.ord <- order(cl.expr)
    } else {
        cl.expr <- numeric()
        cl.ord <- 1:n.clusters
    }

    # Get x-axis pos and labels.
    pos.labels <- colnames(cvgs[[1]])
    pos <- as.integer(gsub('%', '', pos.labels)) # Get raw position if labels have
                                                 # percentage signs.
    
    if(verbose) message("Generating plot.")
    if(plot.type == "line")
    {
	# Group each cluster from all epigenetic marks.

	profiles <- list() 
	for(i in 1:n.clusters)
	    profiles[[i]] <- sapply(cvgs, function(x)
                                    colMeans(x[cl.id == cl.ord[i], , drop = FALSE]))

	# Plot the lineplots.
	invisible(mapply(function(x, y, z)
	{
            y.max = max(unlist(x)) * 1.1
            if(z) y.min = -y.max else y.min = 0
	    layout(rbind(1:2), widths=c(3, 1))
	    par(mai=c(1, 1, 0.8, 0.1))
	    matplot(pos, x, ylim = c(y.min, y.max), type = 'l', lty = 1, col = cols,
                    xlab = "Relative Position", ylab = "Smoothed Coverage", yaxt = 'n',
                    xaxs = 'i', yaxs = 'i',
                    main = paste("Within Cluster Coverage", if(!is.na(y)) paste(" (Median Expression :",
                    round(y, 2)), ")", sep = ''))
	    axis(2, at = c(y.min, (y.min + y.max) / 2, y.max), labels = c("None", "Medium", "High"))

	    plot.new()
	    par(mai=c(1,0.1,0.2,0.1))
	    legend("topleft", legend = names(c.list), title = "Mark", col = cols, lwd = 2)
	}, profiles, cl.expr[cl.ord], symm.scale))

    } else { # Plot a heatmap

        if(is.null(c.list@.old.ranges))
            ranges <- lapply(cvgs, range)
        else
            ranges <- c.list@.old.ranges

        plot.ranges <- mapply(function(x, y, z)
                             {
                                if(z)
                                    c(-max(abs(y)), max(abs(y)))
                                else
                                    range(x)
                             }, cvgs, ranges, symm.scale, SIMPLIFY = FALSE)
	
        par(oma = c(1, 1, 3, 1))
	sort.data <- c.list@sort.data
	sort.name <- c.list@sort.name
	# Get order of all features next.
	if(length(sort.data) == 0)
	    ord <- order(factor(cl.id, levels = cl.ord))
	else
	    ord <- order(factor(cl.id, levels = cl.ord), sort.data)
	
	# Re-arrange the ChIP and expression data and vector of sort data.
	cvgs <- lapply(cvgs, function(x) x[ord, ])
	if(!is.null(expr)) expr <- expr[ord]
	if(length(sort.data) > 0) sort.data <- sort.data[ord]

	# Plot heatmap.
        extras <- sum(!is.null(expr), !is.null(sort.data))
	switch(as.character(extras),
            `0` = layout(rbind(1:(n.marks + 1)), widths=c(1, rep(3, n.marks))),
	    `1` = layout(rbind(1:(n.marks + 2)), widths=c(1, rep(3, n.marks), 2)),
	    `2` = layout(rbind(1:(n.marks + 3)), widths=c(1, rep(3, n.marks), 2, 1)))
	par(mai=c(1.02,0.50,0.82,0.05))
  	
	n.bins = length(cols)
	image(y=seq(1/n.bins/2, 1-(1/n.bins/2), 1/n.bins), z=rbind(1:n.bins),
              col = cols, axes = FALSE, xlab = "Read Coverage", ylab = NA)
	axis(2, at=c(0, 0.5, 1), labels=c("None", "Medium", "High"))

	par(mai=c(1.02,0.05,0.82,0.05))
	mapply(function(x, y, z)
	{
	    image(pos, 1:nrow(x), t(x), zlim = z, xlab="Relative Position",
                  xaxt = "n", yaxt = "n", ylab = "Feature", col = cols, main = y)
	    axis(1, pos, labels = pos.labels)		

	    # Add lines delimiting the cluster boundaries.
	    bounds <- cumsum(table(cl.id)[cl.ord])[-n.clusters]
	    abline(h = bounds, lwd = 2)
	}, cvgs, names(c.list), plot.ranges)

	par(mai=c(1.02,0.05,0.82,0.50))
        if(!is.null(expr))
	    plot(expr, y = 1:length(expr), yaxs = 'i', xlab = expr.name, ylab = NA,
                 yaxt = "n", ...)
	if(!is.null(sort.data)) plot(sort.data, y = 1:length(sort.data), yaxs = 'i',
                     xlab = sort.name, ylab = NA, yaxt = "n", ...)	

        if(!is.null(t.name))
	    mtext(t.name, line = 0, font = 2, cex = 1.5, outer = TRUE)
    }
})

setMethod("clusterPlots", "ScoresList", function(c.list, scale = function(x) x,
          cap.q = 0.95, cap.type = c("sep", "all"), n.clusters = 5,
          plot.ord = 1:length(c.list), expr = NULL, expr.name = NULL, sort.data = NULL, sort.name = NULL,
          plot.type = c("line", "heatmap"), cols = NULL, t.name = NULL, verbose = TRUE,
          ...)
{
    c.list <- c.list[plot.ord]
    plot.type <- match.arg(plot.type)
    cap.type <- match.arg(cap.type)
    cvgs <- tables(c.list)

    if(!is.null(expr) && nrow(cvgs[[1]]) != length(expr))
	stop("Number of features does not match length of expression vector.\n")
    if(!is.null(sort.data) && is.null(sort.name))
    	stop("'sort.data' provided, 'sort.name' is not.")
    if(!is.null(sort.data) && !is.null(expr) && length(expr) != length(sort.data))
	stop("'sort.data' length not the same as number of features in coverage
               matrices or expression data.\n")

    # Precision sometimes means 0 is represented as very small negative numbers.
    cvgs <- lapply(cvgs, function(x) {x[x < 0] = 0; x})

    cvgs <- lapply(cvgs, scale)
    if(cap.type == "all") max.cvg <- quantile(do.call(cbind, cvgs), cap.q)

    # Find the maximum score allowable, then make any scores bigger be the
    # maximum score.
    cvgs <- lapply(cvgs, function(x)
    {
	if(cap.type == "sep") max.cvg <- quantile(x, cap.q)
	x[x > max.cvg] = max.cvg
	return(x)
    })

    # Do the k-means clustering for all marks together.
    set.seed(100)
    if(verbose) message("Doing k-means clustering.")
    all <- do.call(cbind, cvgs)
    cl.id <- kmeans(all, n.clusters, iter.max = 100)$cluster

    # Make sure order of IDs correspond to increasing expression, if present.
    if(!is.null(expr))
    {
        cl.expr <- tapply(expr, factor(cl.id), median, na.rm = TRUE)
        cl.ord <- order(cl.expr)
        cl.id <- sapply(cl.id, function(x) match(x, cl.ord))
    }

    ccl <- ClusteredCoverageList(c.list, scores = cvgs, cluster.id = cl.id, expr = expr,
                                 expr.name = expr.name, sort.data = sort.data,
                                 sort.name = sort.name)
    clusterPlots(ccl, plot.type = plot.type, cols = cols,
                 t.name = t.name, verbose = verbose, ...)
    invisible(ccl)
})
