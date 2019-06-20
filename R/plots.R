#' Generates a venn diagram from a \linkS4class{GenomicOverlap} object.
#'
#' @param x A \linkS4class{GenomicOverlap} object.
#' @param title A title for the generated venn diagram.
#' @return The grid object representing the venn.diagram.
#' @importFrom VennDiagram venn.diagram
#' @importFrom grid grid.draw
#' @importFrom methods is
#' @export
plot_venn <- function(x, title=NULL) {
    stopifnot(is(x, "GenomicOverlap"))
    
    if(length(x) > 5) {
        stop("Cannot plot venn diagram of more than 5 groups!")
    }

    item_list = lapply(names(x), function(y) { union_indices(x, y) })

    fill_vec = c("red", "yellow", "green", "blue", "orange")[seq_along(x)]
    grid_obj = VennDiagram::venn.diagram(item_list,
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
    ggplot(enrichment_df, aes_string(fill="Enrichment", y="RegionType", x='"Network regions"')) +
        geom_tile(color="black") + 
        geom_text(mapping=aes(label=sprintf("%.2f", enrichment_df$Enrichment))) +
        scale_fill_gradient2(low="dodgerblue", mid="white", high="red",
                             midpoint=0, limits=c(-maxEnrich, maxEnrich)) +
        labs(y="Region type", x=NULL) +
        theme(axis.line=element_blank(),
              axis.ticks=element_blank(),
              axis.title=element_blank())
}
