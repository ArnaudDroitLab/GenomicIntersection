#' GenomicEnrichment objects represent the enrichment of a set of 
#' genomic regions over a genomic partition.
#'
#' Some details about this class and my plans for it in the body.
#'
#' @slot regions A GRanges object representing the combined regions.
#' @slot partition A matrix representing which elements overlap with the combined 
#'              regions.
#' @slot enrichments
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
setMethod("length",
          c(x="GenomicEnrichment"),
          function(x) {
            length(x@regions)
          }) 


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
single_region_enrichment <- function(query_regions, genome_partition, genome_order=NULL) {
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
    results=lapply(queries_regions, single_region_enrichment, 
                   genome_partition=genome_partition, genome_order=genome_order)
    
   
    # Summarize the results and return them.
    region_enrichment_summary(results, query.order=query.order,
                              genome.order=genome.order)
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
region_enrichment_summary <- function(result.list, query.order=NULL, genome.order=NULL) {
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