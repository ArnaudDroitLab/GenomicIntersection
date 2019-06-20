
# Function to peproces the indices parameters in the intersect_ and union_
# functions.
preprocess_indices <- function(x, indices) {
    if(is(indices, "character")) {
        indices = names(x) %in% indices
    } else if(is(indices, "numeric")) {
        logical_indices = rep(FALSE, length(x))
        logical_indices[indices] = TRUE
        indices = logical_indices
    }
    stopifnot(is(indices, "logical"))
    
    indices
}

#' Calculate the intersection of a set of elements within a 
#' \linkS4class{GenomicOverlap} object.
#'
#' @param x A \linkS4class{GenomicOverlap} object.
#' @param indices The names of the elements from \code{x} whose intersection 
#'                should be found.
#' @param exclusive If \code{TRUE}, the returned intersection will exclude those 
#'                  ranges where the elements from \code{indices} are not the
#'                  only ones present.
#' @return A vector of the numeric indices of the element of \code{regions(x)} 
#'         which fit the given criteria.
#' @export
intersect_indices <- function(x, indices=names(x), exclusive=FALSE) {
    stopifnot(is(x, "GenomicOverlap"))
   
    indices = preprocess_indices(indices)
    stopifnot(is(exclusive, "logical"))
    
    has.factor = apply(intersect_matrix(x)[, indices, drop=FALSE] >= 1, 1, all)

    if(!exclusive) {
        no.others = rep(TRUE, nrow(intersect_matrix(x)))
    } else {
        no.others  = apply(intersect_matrix(x)[, !indices, drop=FALSE] == 0, 1, all)
    }

    return(has.factor & no.others)
}

#' Calculate the intersection of a set of elements within a 
#' \linkS4class{GenomicOverlap} object.
#'
#' @param x A \linkS4class{GenomicOverlap} object.
#' @param indices The names of the elements from \code{x} whose intersection 
#'                should be found.
#' @param exclusive If \code{TRUE}, the returned intersection will exclude those 
#'                  ranges where the elements from \code{indices} are not the
#'                  only ones present.
#' @return A \code{GRanges} objects containing the ranges from \code{regions(x)} 
#'         which fit the given criteria.
#' @export
intersect_regions <- function(x, indices=names(x), exclusive=TRUE) {
    res_indices = overlaps(x, indices, exclusive)
    regions(x)[res_indices]
}

#' Calculate the union of a set of elements within a 
#' \linkS4class{GenomicOverlap} object.
#'
#' @param x A \linkS4class{GenomicOverlap} object.
#' @param indices The names of the elements from \code{x} whose union 
#'                should be found.
#' @return A vector of the numeric indices of the element of \code{regions(x)} 
#'         which fit the given criteria.
#' @export
union_indices <- function(x, indices=names(x)) {
    stopifnot(is(x, "GenomicOverlap"))
    
    indices = preprocess_indices(indices)
    stopifnot(is(exclusive, "logical"))
    
    has.factor = apply(intersect_matrix(x)[, indices, drop=FALSE] >= 1, 1, any)
    
    return(has.factor)
}

#' Calculate the union of a set of elements within a 
#' \linkS4class{GenomicOverlap} object.
#'
#' @param x A \linkS4class{GenomicOverlap} object.
#' @param indices The names of the elements from \code{x} whose union 
#'                should be found.
#' @return A \code{GRanges} objects containing the ranges from \code{regions(x)} 
#'         which fit the given criteria.
#' @export
union_regions <- function(x, indices=names(x)) {
    res_indices = union_indices(x, indices, exclusive)
    regions(x)[res_indices]
}

#' Determines which regions form a "consensus" from all input regions.
#'
#' @param x A \linkS4class{GenomicOverlap} object.
#' @param consensus_threshold The fraction of input regions which must have
#'                            a given region for it to be selected.
#' @return A vector of the numeric indices of the element of \code{regions(x)} 
#'         which fit the given criteria.
#' @export
consensus_indices <- function(x, consensus_threshold) {
    stopifnot(is(x, "GenomicOverlap"))
    stopifnot(is(consensus_threshold, "numeric"))
    
    return((rowSums(x@matrix) / length(x)) > consensus_threshold)
}

#' Determines which regions form a "consensus" from all input regions.
#'
#' @param x A \linkS4class{GenomicOverlap} object.
#' @param consensus_threshold The fraction of input regions which must have
#'                            a given region for it to be selected.
#' @return A \code{GRanges} objects containing the ranges from \code{regions(x)} 
#'         which fit the given criteria.
#' @export
consensus_regions <- function(x, consensus_threshold) {
    stopifnot(is(x, "GenomicOverlap"))
    res_indices = consensus_indices(x, consensus_threshold)
    
    return(x@regions[res_indices])
}


# Utility function which aggregates one meta-data clumn from grl and
# returns a data-frame with aggregated values for each element of all_regions.
import_column <- function(grl, all_regions, col_name, aggregate_func) {
    # Apply import logic to each GRangesList item
    gr_col = lapply(grl, function(x) {
        
        # Handle the empty-list case by specifying a default all-NA result.
        col.values.vec = rep(NA, length(all_regions))
        if(length(x) > 0) {
            # Map final regions to initial regions
            indices <- GenomicRanges::findOverlaps(all_regions, x)
            
            # Get column values in a data-frame along with their target
            # index so we can identify which ones have a many-to-one
            # mapping.
            col.values.df <- data.frame(from = S4Vectors::from(indices), 
                                        value = mcols(x)[[col_name]][S4Vectors::to(indices)])
    
            # Apply the summarizing function to all groups of values
            # with the same target region.
            col.values.df <- S4Vectors::aggregate(value~from,
                                                  data = col.values.df, 
                                                  FUN = aggregate_func)
            
            # Reorder everything in a vector so we can reassign it to the
            # GRangesList mcols object.
            col.values.vec[col.values.df$from] = col.values.df$value
        }
        
        return(col.values.vec)
    })
    
    col_df = do.call(cbind, gr_col)
    colnames(col_df) <- paste0(col_name, ".", names(grl))
    
    return(col_df)
}

#' Create a \linkS4class{GenomicOverlap} object.
#'
#' Given a \linkS4class{GRangesList} object, determines which regions overlap 
#' with each others. The input regions are "flattened", and all overlapping
#' regions (within and across the elements of the input \code{grl} parameter) 
#' are mapped to a new region whose \code{start} and \code{end} are the minimum 
#' \code{start} and maximum \code{end} of the initial overlapping regions.
#'
#' A matrix indicating which input regions fall within these new mapped regions 
#' is then produced.
#'
#' @param grl The \linkS4class{GRangesList} object whose elements need to be 
#'            overlapped with each others.
#' @param import_spec A list of columns which should be imported into
#'                    the resulting object. Each element should be named after
#'                    a column in \code{mcols(grl)}, and should contain a 
#'                    function to be used to aggregate multiple values of that 
#'                    column.
#' @return An object of class \code{GenomicOverlap}.
#' @importFrom GenomicRanges reduce
#' @importFrom GenomicRanges mcols
#' @importMethodsFrom GenomicRanges countOverlaps findOverlaps
#' @importMethodsFrom S4Vectors from to aggregate
#' @export
GenomicOverlap <- function(grl, import_spec=list()) {
    # Validate input parameters.
    stopifnot(is(grl, "GRangesList"))
    stopifnot(is(keep.signal, "logical"))
    stopifnot(all(names(import_spec) %in% names(mcols(grl))))
    
    # Flatten the GRangesList so we can get a list of all possible regions.
    all_regions = GenomicRanges::reduce(unlist(grl))

    # Build a matrix to hold the results.
    overlap.matrix <- matrix(0, nrow=length(all_regions), ncol=length(grl))

    # Loop over all ranges, intersecting them with the flattened list of all possible regions.
    for(i in 1:length(grl)) {
        overlap.matrix[,i] <- GenomicRanges::countOverlaps(all_regions, grl[[i]], type="any")
    }
    colnames(overlap.matrix) <-  names(grl)

    # Import all columns in import list.
    for(col_name in names(import_spec)) {
        col_df = import_column(grl, all_regions, col_name, import_spec[[col_name]])
        
        mcols(all_regions) <- cbind(mcols(all_regions), col_df)
    }

    new("GenomicOverlap",
        regions=all_regions, 
        matrix=overlap.matrix)
}

setClass("GenomicOverlap",
         slots=list(regions="GRanges", 
                    matrix="matrix"))

setMethod("names",
          c(x="GenomicOverlap"),
          function(x) {
            colnames(x@matrix)
          })

setMethod("names<-",
          c(x="GenomicOverlap", value="character"),
          function(x, value) {
            colnames(x@matrix) <- value
            x
          })
          
setMethod("length",
          c(x="GenomicOverlap"),
          function(x) {
            ncol(x@matrix)
          })   
          
#' Returns the combined regions from a \linkS4class{GenomicOverlap} object.
#'
#' @param x The \linkS4class{GenomicOverlap} object.
#' @return A \code{GRanges} object representing the combined regions.
#' @export
regions <- function(x) {
    stopifnot(is(x, "GenomicOverlap"))
    x@regions
}

#' Returns the intersection matrix from a \linkS4class{GenomicOverlap} object.
#'
#' @param x The \linkS4class{GenomicOverlap} object.
#' @return A matrix with has as many columns as the number of items in the 
#'         initial \linkS4class{GRangesList}, and as many rows as the number of 
#'         combined regions. Its values indicate how many regions from the 
#'         initial \linkS4class{GRangesList} element map to each individual
#'         combined range. 
#' @export
intersect_matrix <- function(x) {
    stopifnot(is(x, "GenomicOverlap"))
    x@matrix
}


#' Calculate the pairwise overlaps of all factors within an \code{intersect.object}.
#'
#' @param intersect.object An intersect object returned by \code{\link{build_intersect}}.
#' @param filename A name for prefix for the names of the output files where the
#'   annotations will be written. Pass \code{NULL} to skip saving to a file.
#' @return A matrix containing the pairwise overlaps of all factors in
#' \code{intersect.object}. The row's factor is used as a denominator. Therefore, the
#    matrix is not symmetric.
#' @export
pairwise_overlap <- function(intersect.object, filename=NULL) {
    overlap.percentage <- matrix(0, nrow=intersect.object$Length, ncol=intersect.object$Length,
                                 dimnames=list(intersect.object$Names, intersect.object$Names))

    # Compare factors two by two.
    for(i in 1:intersect.object$Length) {
        for(j in 1:intersect.object$Length) {
            i.vector = intersect.object$Matrix[,i] >= 1
            j.vector = intersect.object$Matrix[,j] >= 1
            overlap.percentage[i,j] = sum(i.vector & j.vector) / sum(i.vector)
        }
    }

    if(!is.null(filename)) {
        write.table(overlap.percentage, file=filename, sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE)
    }

    return(overlap.percentage)
}

#' Projects the ranges in query into the ranges in target.
#'
#' Returns the ranges in target overlapping the ranges in query, adjusting
#' their boundaries so that only the overlapping parts of the target ranges 
#' are returned.
#'
#' @param query The ranges to be projected.
#' @param target The ranges to be projected against.
#'
#' @return The projection of the query ranges on the target ranges.
#' @export
project_ranges <- function(query, target) {
    hits = findOverlaps(target, query)

    ranges.df = data.frame(seqname=seqnames(target)[queryHits(hits)],
                            start=pmax(start(query)[subjectHits(hits)], start(target)[queryHits(hits)]),
                            end=pmin(end(query)[subjectHits(hits)], end(target)[queryHits(hits)]),
                            strand=strand(target)[queryHits(hits)])
    
    ranges.df = cbind(ranges.df, mcols(target)[queryHits(hits),], mcols(query)[subjectHits(hits),])
    colnames(ranges.df) = c(colnames(ranges.df)[1:4], colnames(mcols(target)), colnames(mcols(query)))
    
    return(GRanges(ranges.df))
}


#' Collapses a list of genomic ranges into a single set of unique, 
#' non-overlapping ranges.
#'
#' Ranges are prioritized in the input list order. So, if the first element
#' of the list (A) covers the range 1-10, and the second element (B) covers 
#' the range 5-15, then the resulting ranges will have a 1-10 range named A,
#' and a 11-15 range named 'B'.
#'
#' @param gr.list The ranges to be collapsed.
#'
#' @return The collapsed regions.
#' @export
collapse_regions <- function(gr.list) {
  # Resulting regions.
  collapsed.regions = list()
  
  # Keep track of the ranges that have already been assigned.
  combined.regions = GenomicRanges::GRanges()
  for(region.group in names(gr.list)) {
    # The ranges assigned to this element are all the specified ranges,
    # minus any range that has already been assigned.
    collapsed.regions[[region.group]] = GenomicRanges::setdiff(gr.list[[region.group]], combined.regions)
    collapsed.regions[[region.group]]$name = region.group
    
    # Add the newly assigned ranges to the set of assigned ranges.
    combined.regions = GenomicRanges::union(combined.regions, collapsed.regions[[region.group]])
  }

  # Return a single set of ranges.
  return(unlist(GenomicRanges::GRangesList(collapsed.regions)))
}

