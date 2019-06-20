#' Generates a venn diagram from an \code{intersect.object}.
#'
#' @param intersect.object An intersect object returned by \code{\link{build_intersect}}.
#' @param filename A filename for the resulting venn diagram. Pass \code{NULL}
#'   to skip saving to a file.
#' @return The grid object representing the venn.diagram.
#' @importFrom VennDiagram venn.diagram
#' @importFrom grid grid.draw
#' @export
plot_venn <- function(x, title=NULL) {
    stopifnot(is(x, "GenomicOverlap"))
    
    if(length(x) > 5) {
        stop("Cannot plot venn diagram of more than 5 groups!")
    }

    fill_vec = c("red", "yellow", "green", "blue", "orange")[seq_along(x)]
    grid_obj = VennDiagram::venn.diagram(intersect.object$List,
                                     fill=fill_vec,
                                     filename=NULL,
                                     print.mode="raw",
                                     main=title)
    grid.draw(grid_obj)
}


#' @import ggplot2
#' @export
plot_region_enrichment <- function(enrichment_df) {
    # Replace +/-Inf with NAs.
    enrichment_df$Enrichment[is.infinite(enrichment_df$Enrichment)] <- NA
    
    maxEnrich = max(abs(enrichment_df$Enrichment), na.rm=TRUE)
    ggplot(enrichment_df, aes(fill=Enrichment, y=RegionType, x="Network regions")) +
        geom_tile(color="black") + 
        geom_text(mapping=aes(label=sprintf("%.2f", enrichment_df$Enrichment))) +
        scale_fill_gradient2(low="dodgerblue", mid="white", high="red",
                             midpoint=0, limits=c(-maxEnrich, maxEnrich)) +
        labs(y="Region type", x=NULL) +
        theme(axis.line=element_blank(),
              axis.ticks=element_blank(),
              axis.title=element_blank())
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
multiple_region_enrichment <- function(queries.regions, genome.regions, query.order=NULL,
                                       genome.order=NULL, file.prefix=NULL, plot.width=7, plot.height=7,
                                       individual.plots=FALSE) {
    results=list()
    
    # Loop over all given query regions and perform enrichments.
    for(query in names(queries.regions)) { 
        # If we have an output prefix, figure out the name for the query-specific output.
        if(!is.null(file.prefix) && individual.plots) {
            file.out = paste0(file.prefix, " ", query, ".pdf")
        } else {
            file.out = NULL
        }

        results[[query]] = region_enrichment(queries.regions[[query]], genome.regions,
                                             genome_order=genome.order)
    }
    
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
region_enrichment_summary <- function(result.list, file.prefix=NULL, query.order=NULL, genome.order=NULL, plot.width=7, plot.height=7) {
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

    # If a file prefix was provided, write out the tables/plots.
    results.plot = list()
    if(!is.null(file.prefix)) {
        for(metric in names(results.data)) {
            write.table(results.data[[metric]], file=paste0(file.prefix, " ", metric, ".txt"), sep="\t", col.names=TRUE, row.names=TRUE)
            
            result.df = reshape2::melt(results.data[[metric]], varnames=c("Query", "Category"))
            
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
            
            results.plot[[metric]] = ggplot(result.df, aes(x=Query, y=Category, fill=value)) +
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
                results.plot[[metric]] = results.plot[[metric]] + geom_text(mapping=aes(label=sprintf("%.0f%%", value*100)))
            } else if(metric=="Enrichment") {
                results.plot[[metric]] = results.plot[[metric]] + geom_text(mapping=aes(label=sprintf("%.1f", value)))
            }
            
            # Change type fo scale (two colors or three colors) depending on the type of data.
            if(metric=="Enrichment") {
                results.plot[[metric]] = results.plot[[metric]] + scale_fill_gradient2(low="dodgerblue", mid="white", high="red", name=metric)
            } else {
                results.plot[[metric]] = results.plot[[metric]] + scale_fill_gradient(low="white", high="red", name=metric)
            }
            ggsave(paste0(file.prefix, " all ", metric, ".pdf"), plot=results.plot[[metric]], width=plot.width, height=plot.height, limitsize=FALSE)
        }
    }
    
    return(list(Data=results.data, Plots=results.plot))
}
