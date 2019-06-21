#' Generates a venn diagram from a \linkS4class{GenomicOverlaps} object.
#'
#' @param x A \linkS4class{GenomicOverlaps} object.
#' @param title A title for the generated venn diagram.
#' @return The grid object representing the venn.diagram.
#' @importFrom VennDiagram venn.diagram
#' @importFrom grid grid.draw
#' @importFrom methods is
#' @export
plot_venn <- function(x, title=NULL) {
    stopifnot(is(x, "GenomicOverlaps"))
    
    if(length(x) > 5) {
        stop("Cannot plot venn diagram of more than 5 groups!")
    }

    item_list = lapply(names(x), function(y) { which(union_indices(x, y)) })
    names(item_list) = names(x)
    
    fill_vec = c("red", "yellow", "green", "blue", "orange")[seq_along(x)]
    grid_obj = VennDiagram::venn.diagram(item_list,
                                     fill=fill_vec,
                                     filename=NULL,
                                     print.mode="raw",
                                     main=title)
    grid.draw(grid_obj)
}
