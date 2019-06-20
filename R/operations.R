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
#' @importFrom S4Vectors subjectHits queryHits
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
    
    return(GRanges(ranges_df))
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
        coverages[region.type, "Genome"] = sum(as.numeric(BiocGenerics::width(GenomicRanges::reduce(genome_partition[genome_partition$name == region.type]))))
        coverages[region.type, "Query"]  = sum(as.numeric(BiocGenerics::width(GenomicRanges::reduce(in.query[in.query$name == region.type]))))
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

#' Performs region enrichment on a set of regions and returns summarized results.
#'
#' @param queries.regions A list of regions whose enrichment ratios must be calculated.
#' @param genome.regions The genome partition indicating which part of the genome
#'    fall within which category. Each range should have a 'name' attribute
#' @param factor.order An optional ordering of the region types for the produced plot.
#' @param file.prefix An optional file name prefix for tables and graphical representation.
#' @param plot.width The width of any resulting summary plot.
#' @param plot.height The height of any resulting summary plot.
#' @param individual.plots If true, produce individual plots as well as combined plots.
#' @return A list of summarized enrichment metrics.
#' @export
multiple_region_enrichment <- function(queries_regions, genome_partition, query.order=NULL,
                                       genome_order=NULL) {
    results=lapply(queries_regions, region_enrichment, 
                   genome_partition=genome_partition, genome_order=genome_order)
    
   
    # Summarize the results and return them.
    region_enrichment_summary(results, file.prefix, query.order=query.order,
                              genome.order=genome.order, plot.width=plot.width, plot.height=plot.height)
}

#' Performs a summary of region enrichment results.
#'
#' @param result.list a list of results returned by region_enrichment.
#' @param file.prefix An optional file name prefix for tables and graphical representation.
#' @param genome.regions The genome partition indicating which part of the genome
#'    fall within which category. Each range should have a 'name' attribute
#' @param factor.order An optional ordering of the region types for the produced plot.
#' @param plot.width The width of any resulting plot.
#' @param plot.height The height of any resulting plot.
#' @return A list of summarized enrichment metrics.
#' @importFrom reshape2 melt
#' @export
region_enrichment_summary <- function(result.list, file.prefix=NULL, query.order=NULL, genome.order=NULL) {
    # Put all of metrics into a single multidimensional array.
    metrics = c("QueryCoverage", "QueryProportion", "Enrichment")
    result.summary = array(dim=c(length(result.list), nrow(result.list[[1]]), 3), 
                           dimnames=list(names(result.list), rownames(result.list[[1]]), metrics))
    for(result in names(result.list)) {
        for(metric in metrics) {
            result.summary[result, ,metric] = result.list[[result]][,metric]
        }
    }

    # Add genomic/background information where appropriate. Enrichment is a ratio, so it cannot
    # have genomic/background data.
    results.data = list(Coverage=rbind(Genome=result.list[[1]]$GenomeCoverage, result.summary[,,"QueryCoverage"]),
                        Proportion=rbind(Genome=result.list[[1]]$GenomeProportion, result.summary[,,"QueryProportion"]),
                        Enrichment=result.summary[,,"Enrichment"])

    return(list(Data=results.data, Plots=results.plot))
}

#' @importFrom reshape2 melt
plot_multiple_region_enrichment <- function(enrichment_results, metric="Enrichment") {
    result.df = reshape2::melt(enrichment_results[[metric]], varnames=c("Query", "Category"))
    
    # Reorder queries
    if(is.null(query.order)) {
        query.order = rownames(results.data[[metric]])
    }
    result.df$Query = factor(result.df$Query, levels=query.order)

    # Reorder categories.
    if(is.null(genome.order)) {
        genome.order = colnames(results.data[[metric]])
    }
    result.df$Category = factor(result.df$Category, levels=rev(genome.order))
    
    result = ggplot(result.df, aes(x=Query, y=Category, fill=value)) +
        geom_tile() +
        theme_bw() +
        theme(axis.line = element_line(colour = "black"),
              axis.text = element_text(color="black"),
              axis.text.x = element_text(angle = 90, hjust = 1),
              axis.title = element_text(size=14),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank())
        

    # Change text labels depending on the type of data.
    if(metric=="Proportion") {
        result = result + geom_text(mapping=aes(label=sprintf("%.0f%%", value*100)))
    } else if(metric=="Enrichment") {
        result = result + geom_text(mapping=aes(label=sprintf("%.1f", value)))
    }
    
    # Change type of scale (two colors or three colors) depending on the type of data.
    if(metric=="Enrichment") {
        result = result + scale_fill_gradient2(low="dodgerblue", mid="white", high="red", name=metric)
    } else {
        result = result + scale_fill_gradient(low="white", high="red", name=metric)
    }
    
    
    ggsave(paste0(file.prefix, " all ", metric, ".pdf"), plot=result, width=plot.width, height=plot.height, limitsize=FALSE)
}