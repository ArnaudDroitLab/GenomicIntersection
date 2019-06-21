#' GenomicEnrichment objects represent the enrichment of a set of 
#' genomic regions over a genomic partition.
#'
#' @slot regions A GRangesList object representing the regions on which 
#'               enrichment should be computed.
#' @slot partition A GRanges object representing a partition of the genome
#'                 to be used as background for the enrichment analysis.
#' @slot enrichments A list of data-frames representing the enrichment results
#'                   for each elements of \code{regions}.
#' @name GenomicEnrichment-class
#' @rdname GenomicEnrichment-class
#' @export
setClass("GenomicEnrichment",
         slots=list(regions="GRangesList", 
                    partition="GRanges",
                    enrichments="list"))

#' Returns the names of the elements of a \linkS4class{GenomicEnrichment} object. 
#'
#' @param x The \linkS4class{GenomicEnrichment} object.
#' @return The names of the elements in \code{x}.
#' @export
setMethod("names",
          c(x="GenomicEnrichment"),
          function(x) {
            names(x@regions)
          })

#' Set names of the elements of a \linkS4class{GenomicEnrichment} object. 
#'
#' @param x The \linkS4class{GenomicEnrichment} object.
#' @param value The new names for the elements of the 
#'              \linkS4class{GenomicEnrichment} object.
#' @export
setMethod("names<-",
          c(x="GenomicEnrichment", value="character"),
          function(x, value) {
            names(x@regions) <- value
            names(x@enrichments) <- value
            x
          })

#' Returns the number of elements of a \linkS4class{GenomicEnrichment} object. 
#'
#' @param x The \linkS4class{GenomicEnrichment} object.
#' @return The number of elements in \code{x}.   
#' @export     
setMethod("length",
          c(x="GenomicEnrichment"),
          function(x) {
            length(x@regions)
          }) 


#' Calculates enrichment ratios for query regions against a genome wide
#' partition of the genome.
#'
#' @param query_regions The regions whose enrichment ratios must be calculated.
#' @param genome_partition The genome partition indicating which part of the 
#'    genome fall within which category. Each range should have a 'name' 
#'    attribute indicating its category.
#' @return A data-frame containing the enrichment values.
#' @import GenomicRanges
#' @importFrom BiocGenerics width
#' @export
single_region_enrichment <- function(query_regions, genome_partition) {
    stopifnot(is(query_regions, "GRanges"))
    stopifnot(is(genome_partition, "GRanges"))
    stopifnot("name" %in% names(mcols(genome_partition)))
    
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
    name_levels = all.region.types
    if(is(genome_partition$name, "factor")) {
        name_levels = levels(genome_partition$name)
    }
    enrichment.df$RegionType = factor(enrichment.df$RegionType, levels=name_levels)
    
    return(enrichment.df)
}

#' Performs region enrichment on a set of regions and returns summarized results.
#'
#' @param queries_regions A list of regions whose enrichment ratios must be 
#'                       calculated.
#' @param genome_partition The genome partition indicating which part of the 
#'    genome fall within which category. Each range should have a 'name' 
#'    attribute indicating its category.
#' @return An object of class \linkS4class{GenomicEnrichment}.
#' @export
GenomicEnrichment <- function(queries_regions, genome_partition) {
    results=lapply(queries_regions, single_region_enrichment, 
                   genome_partition=genome_partition)
                              
    methods::new("GenomicEnrichment",
                 regions = queries_regions,
                 partition = genome_partition,
                 enrichments = results)
}

# Utility function to extract columns from a GenomicEnrichment object and turn
# them into a matrix.
metric_matrix <- function(x, metric) {
    stopifnot(is(x, "GenomicEnrichment"))
    stopifnot(metric %in% c(c("QueryCoverage", "QueryProportion", "Enrichment")))
    
    results = do.call(cbind, lapply(x@enrichments, function(y) {y[,metric]}))
    
    rownames(results) = rownames(x@enrichments[[1]])
    
    results
}

# Utility function to extract columns from a GenomicEnrichment object and turn
# them into a data-frame including Genome-wide info, if relevant.
metric_df <- function(x, metric) {
    result_matrix = metric_matrix(x, metric)
    
    if(metric=="QueryCoverage") {
        results = cbind(Genome=x@enrichments[[1]]$GenomeCoverage,
                        result_matrix)
    } else if(metric=="QueryProportion") {
        results = cbind(Genome=x@enrichments[[1]]$GenomeProportion,
                        result_matrix)    
    } else {
        results = result_matrix
    }

    results
}

#' Returns a data-frame giving the coverages for all elements of a 
#' \linkS4class{GenomicEnrichment} object.
#'
#' @param x A \linkS4class{GenomicEnrichment} object.
#' @return A data-frame giving genome-wide and per-element coverages in 
#'         nucleotides.
#' @export
coverage_df <- function(x) {
    stopifnot(is(x, "GenomicEnrichment"))
    return(metric_df(x, "QueryCoverage"))
}

#' Plots the coverages for all elements of a 
#' \linkS4class{GenomicEnrichment} object.
#'
#' @param x A \linkS4class{GenomicEnrichment} object.
#' @return A plot of genome-wide and per-element coverages in 
#'         nucleotides.
#' @export
coverage_plot <- function(x) {
    plot_genomic_enrichment_metric(x, "Coverage")
}

#' Returns a data-frame giving the coverage proportions for all elements of a 
#' \linkS4class{GenomicEnrichment} object.
#'
#' @param x A \linkS4class{GenomicEnrichment} object.
#' @return A data-frame giving genome-wide and per-element coverage proportions 
#'         in nucleotides.
#' @export
proportion_df <- function(x) {
    stopifnot(is(x, "GenomicEnrichment"))
    return(metric_df(x, "QueryProportion"))
}

#' Plots the coverage proportions for all elements of a 
#' \linkS4class{GenomicEnrichment} object.
#'
#' @param x A \linkS4class{GenomicEnrichment} object.
#' @return A plot giving genome-wide and per-element coverage proportions 
#'         in nucleotides.
#' @export
proportion_plot <- function(x) {
    plot_genomic_enrichment_metric(x, "Proportion")
}

#' Returns a data-frame giving the enrichments for all elements of a 
#' \linkS4class{GenomicEnrichment} object.
#'
#' @param x A \linkS4class{GenomicEnrichment} object.
#' @return A data-frame giving per-element log2-enrichments of all partition 
#'         categories.
#' @export
enrichment_df <- function(x) {
    stopifnot(is(x, "GenomicEnrichment"))
    return(metric_df(x, "Enrichment"))
}

#' Plots the enrichments for all elements of a 
#' \linkS4class{GenomicEnrichment} object.
#'
#' @param x A \linkS4class{GenomicEnrichment} object.
#' @return A plot of per-element log2-enrichments of all partition 
#'         categories.
#' @export
enrichment_plot <- function(x) {
    plot_genomic_enrichment_metric(x, "Enrichment")
}

# Utility function to retrieve the levels/order of the partition categories.
partition_levels = function(x) {
    if(is(x@partition$name, "factor")) {
        return(levels(x@partition$name))
    } else {
        return(sort(unique(x@partition$name)))
    }
}

#' @import ggplot2
#' @importFrom reshape2 melt
plot_genomic_enrichment_metric <- function(x, metric="Enrichment") {
    if(metric=="Coverage") {
        data_df = coverage_df(x)
    } else if(metric=="Proportion") {
        data_df = proportion_df(x)
    } else if(metric=="Enrichment") {
        data_df = enrichment_df(x)
    } else {
        stop("Unrecognized metric ", metric)
    }
    
    plot_df = reshape2::melt(data_df, varnames=c("Category", "Query"))
    
    # Reorder queries and categories
    query_factors = names(x)
    if(metric != "Enrichment") {
        query_factors = c("Genome", query_factors)
    }
    plot_df$Query = factor(plot_df$Query, levels=query_factors)
    plot_df$Category = factor(plot_df$Category, levels=rev(partition_levels(x)))
    
    result = ggplot(plot_df, aes_string(x="Query", y="Category", fill="value")) +
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
        result = result + geom_text(mapping=aes_string(label='sprintf("%.0f%%", value*100)'))
    } else if(metric=="Enrichment") {
        result = result + geom_text(mapping=aes_string(label='sprintf("%.1f", value)'))
    }
    
    # Change type of scale (two colors or three colors) depending on the type of data.
    if(metric=="Enrichment") {
        result = result + scale_fill_gradient2(low="dodgerblue", mid="white", high="red", name=metric)
    } else {
        result = result + scale_fill_gradient(low="white", high="red", name=metric)
    }
    
    result
}