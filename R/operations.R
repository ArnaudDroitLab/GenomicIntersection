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
    hits = GenomicRanges::findOverlaps(target, query)

    new_start = pmax(start(query)[subjectHits(hits)],
                     start(target)[queryHits(hits)])
    new_end = pmin(end(query)[subjectHits(hits)], 
                   end(target)[queryHits(hits)])
    ranges_df = data.frame(seqname=seqnames(target)[queryHits(hits)],
                            start=new_start,
                            end=new_end,
                            strand=strand(target)[queryHits(hits)])
    
    ranges_df = cbind(ranges_df, 
                      mcols(target)[queryHits(hits),], 
                      mcols(query)[subjectHits(hits),])
    colnames(ranges_df) = c(colnames(ranges_df)[1:4], 
                            colnames(mcols(target)), 
                            colnames(mcols(query)))
    
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
#' @param grl The ranges to be collapsed.
#'
#' @return The collapsed regions.
#' @export
collapse_regions <- function(grl) {
  # Resulting regions.
  collapsed.regions = list()
  
  # Keep track of the ranges that have already been assigned.
  combined.regions = GenomicRanges::GRanges()
  for(region.group in names(grl)) {
    # The ranges assigned to this element are all the specified ranges,
    # minus any range that has already been assigned.
    collapsed.regions[[region.group]] = GenomicRanges::setdiff(grl[[region.group]], combined.regions)
    collapsed.regions[[region.group]]$name = region.group
    
    # Add the newly assigned ranges to the set of assigned ranges.
    combined.regions = GenomicRanges::union(combined.regions, collapsed.regions[[region.group]])
  }

  # Return a single set of ranges.
  return(unlist(GenomicRanges::GRangesList(collapsed.regions)))
}

