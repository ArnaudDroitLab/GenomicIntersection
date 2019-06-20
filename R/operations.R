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

#' Calculates enrichment ratios for query regions against a genome wide
#' partition of the genome.
#'
#' @param query_regions The regions whose enrichment ratios must be calculated.
#' @param genome_partition The genome partition indicating which part of the genome
#'    fall within which category. Each range should have a 'name' attribute
#'    indicating its category.
#' @param genome_order An optional ordering of the region types.
#' @return A data-frame containing the enrichment values.
#' @import GenomicRanges
#' @export
region_enrichment <- function(query_regions, genome_partition, genome_order=NULL) {
  # Project the query ranges into the genome ranges, so we can
  # know their repartition with basepair precision.
  in.query = project_ranges(query_regions, genome_partition)
  
  # Calculate total base-pair coverages for both the projected query 
  # and the target regions.
  all.region.types = sort(unique(genome_partition$name))
  coverages = matrix(0.0, ncol=2, nrow=length(all.region.types), 
                     dimnames=list(all.region.types, c("Query", "Genome")))
  for(region.type in all.region.types) {
      coverages[region.type, "Genome"] = sum(as.numeric(BiocGenerics::width(GenomicRanges::reduce(BiocGenerics::subset(genome_partition, name == region.type)))))
      coverages[region.type, "Query"]  = sum(as.numeric(BiocGenerics::width(GenomicRanges::reduce(BiocGenerics::subset(in.query, name == region.type)))))
  }

  # Transform the raw coverages into proportions.
  proportions = t(apply(coverages, 1, '/', apply(coverages, 2, sum)))
  
  # Build a data frame for output/plotting.
  enrichment.df = data.frame(QueryCoverage=coverages[,"Query"],
                             GenomeCoverage=coverages[,"Genome"],
                             QueryProportion=proportions[,"Query"],
                             GenomeProportion=proportions[,"Genome"],
                             Enrichment=log2(proportions[,"Query"] / proportions[,"Genome"]), 
                             RegionType=all.region.types)
  if(is.null(genome_order)) {
    genome_order = all.region.types
  }
  enrichment.df$RegionType = factor(enrichment.df$RegionType, levels=genome_order)
  
  return(enrichment.df)
}