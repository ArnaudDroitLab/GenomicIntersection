#' GenomicOverlap objects represent the overlaps of multiple GRanges.
#'
#' Some details about this class and my plans for it in the body.
#'
#' @slot regions A GRanges object representing the combined regions.
#' @slot matrix A matrix representing which elements overlap with the combined 
#'              regions.
#'
#' @name GenomicOverlap-class
#' @rdname GenomicOverlap-class
#' @export
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
    res_indices = union_indices(x, indices)
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

#' Calculate the pairwise overlaps of all factors within a 
#' \linkS4class{GenomicOverlap} object.
#'
#' @param x A \linkS4class{GenomicOverlap} object.
#' @return A matrix containing the pairwise overlaps of all elements in
#' \code{x}. The row's element is used as the denominator. Therefore, the
#' matrix is not symmetric.
#'
#' @export
pairwise_overlap <- function(x) {
    stopifnot(is(x, "GenomicOverlap"))
    overlap <- matrix(0, nrow=length(x), ncol=length(x),
                      dimnames=list(names(x), names(x)))

    # Compare factors two by two.
    for(i in seq_along(x)) {
        for(j in seq_along(x)) {
            i_vector = intersect_matrix(x)[,i] >= 1
            j_vector = intersect_matrix(x)[,j] >= 1
            overlap[i,j] = sum(i_vector & j_vector) / sum(i_vector)
        }
    }

    return(overlap)
}

