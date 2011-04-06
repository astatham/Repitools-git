setClass("ScoresList", representation(
                                        names = "character",
					scores = "list", # list of matrices.
					anno = "GRanges",
					up = "numeric",
					down = "numeric",
					dist = "ANY", # character or NULL.
					freq = "numeric",
					s.width = "ANY" # character or NULL.
					))

setMethod("names", "ScoresList", function(x) x@names)
setGeneric("tables", function(x) {standardGeneric("tables")})
setMethod("tables", "ScoresList", function(x) x@scores)
setMethod("length", "ScoresList", function(x) length(x@scores))

setMethod("show", "ScoresList",
    function(object)
    {
        if(!is.null(object@dist))
	    distLabel <- ifelse(object@dist == "percent", '%', "bases")
        else
            distLabel <- "bases"
	cat("An object of class 'ScoresList'.\n")
	cat("Tables: ", paste(object@names, collapse = ", "), ".\n", sep = '')
	cat("Features:\n")
	print(object@anno)
	cat("Region:",  paste(object@up, distLabel, "up to", object@down,
                         distLabel, "down.\n"))
        if(!is.null(object@s.width))
        {
	    cat("Smoothing:", paste(object@s.width, collapse = ", "), "bases.\n")
	    cat("Sampling : ", object@freq, ' ', distLabel, ".\n\n",  sep = '')
        } else {
            cat("Window Width : ", object@freq, ' ', distLabel, ".\n\n",  sep = '')
        }        
    })

setMethod("[", "ScoresList",
    function(x, i)
    {
	new("ScoresList", names = x@names[i], anno = x@anno, scores = x@scores[i],
	                    up = x@up, down = x@down, dist = x@dist,
			    freq = x@freq, s.width = x@s.width[i])
	})

setReplaceMethod("names", "ScoresList",
    function(x, value)
    {
	if(length(value) != length(x@names))
	    stop("New mark name(s) are a different length to previous mark
                  name(s).\n")
	x@names <- value
	x
    }
)

setClass("ClusteredCoverageList", representation(
                                    cluster.id = "numeric",
				    expr = "numeric",
				    sort.data = "ANY",
				    sort.name = "ANY"),
                                 prototype(
				    sort.data = NULL,
				    sort.name = NULL),
                       contains = "ScoresList")

setMethod("show", "ClusteredCoverageList",
    function(object)
    {
	distLabel <- ifelse(object@dist == "percent", '%', "bases")
	cat("An object of class 'ClusteredCoverageList'.\n")
	cat("Tables: ", paste(object@names, collapse = ", "), ".\n", sep = '')
	cat("Region: ",  paste(object@up, distLabel, "up to", object@down,
	    distLabel, "down.\n"))
	cat("Features:\n")
	print(object@anno)
	cat("Smoothing:", paste(object@s.width, collapse = ", "), "bases.\n")
	cat("Sampling: ", object@freq, ' ', distLabel, ".\n",  sep = '')
	cat("Feature Expressions:", paste(paste(head(object@expr),
	    collapse = ", "), ", ...\n", sep = ''))
	cat("Feature Clusters:", paste(paste(head(object@cluster.id),
	    collapse = ", "), ", ...\n", sep = ''))
	if(!is.null(object@sort.data))
	    cat("Within Cluster Sorting: By ", object@sort.name, ". ",
	    paste(paste(head(object@sort.data), collapse = ", "), ", ...\n", sep = ''),
	    sep = '')		
    })

# Constructor
setGeneric("ClusteredCoverageList", function(x, ...)
           {standardGeneric("ClusteredCoverageList")})
setMethod("ClusteredCoverageList", "ScoresList",
    function(x, scores = tables(x), expr, cluster.id, sort.data = NULL,
             sort.name = NULL)
{
	new("ClusteredCoverageList", names = x@names, scores = scores, anno = x@anno,
	    up = x@up, down = x@down, dist = x@dist,
	    freq = x@freq, s.width = x@s.width, cluster.id = cluster.id,
	    expr = expr, sort.name = sort.name, sort.data = sort.data)
})

setMethod("[", "ClusteredCoverageList",
    function(x, i)
    {
	new("ClusteredCoverageList", names = x@names[i], scores = x@scores[i],
	    anno = x@anno, up = x@up, down = x@down, dist = x@dist,
	    freq = x@freq, s.width = x@s.width[i], cluster.id = x@cluster.id,
	    expr = x@expr, sort.data = x@sort.data, sort.name = x@sort.name)
    })

# A collection of variables that describe where the sampling will happen.
# An S4 class, so that they are created once, then dispatched on when
# required for multiple samples.

setClass(".CoverageSamples",
         representation(
                        pos.labels = "ANY", # character or numeric.
                        cvg.samps = "GRanges",
                        max.out = "numeric",
                        chr.ord = "numeric",
                        anno.chr = "GRanges",
                        old.ord = "numeric"
			))

# container for output of regionStats()    
setClass("RegionStats",representation("list"))

setMethod("show", "RegionStats",function(object) {
  cat("Object of class 'RegionStats'.\n")
  cat("Results for: ", paste(names(object$regions),collapse=" "), "\n")
  cat("Names:", paste(names(object),collapse=" "), "\n")
})

# container for output of ChromaBlocks()
setClass("ChromaResults",
    representation(
        blocks="GRanges", 
        regions="RangesList",
        FDRTable="matrix",
        cutoff="numeric"
    )
)

setMethod("show", "ChromaResults", function(object) {
  cat("Object of class 'ChromaResults'.\n")
  cat(sum(sapply(object@regions, length)), "regions found with using a cutoff of", object@cutoff, "\n")
})

#ChromaResults Generics
setGeneric("blocks", function(x) standardGeneric("blocks"))
setGeneric("regions", function(x) standardGeneric("regions"))
setGeneric("FDRTable", function(x) standardGeneric("FDRTable"))
setGeneric("cutoff", function(x) standardGeneric("cutoff"))

#ChromaResults Accessors
setMethod("blocks", "ChromaResults", function(x) x@blocks)
setMethod("regions", "ChromaResults", function(x) x@regions)
setMethod("FDRTable", "ChromaResults", function(x) x@FDRTable)
setMethod("cutoff", "ChromaResults", function(x) x@cutoff)

