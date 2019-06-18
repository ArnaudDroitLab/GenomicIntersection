#' Internal exclusive overlap.
#'
#' Returns the set of regions within \code{all.regions} where all
#' conditions/proteins at indices which.factors are present,
#' and no other factors are. Presence/absence at a given locus
#' is obtained from the \code{overlap.matrix} parameter.
#'
#' @param all.regions A \linkS4class{GRanges} object with the regions represented
#    by the rows of \code{overlap.matrix}.
#' @param overlap.matrix A matrix, with rows corresponding to regions
#'   in \code{all.rgions}, and columns corresponding to proteins/factors. A
#'   non-zero value in the matrix indicates the protein/factor of interest
#'   is present at this region.
#' @param which.factors Indices of columns where the factor should be present.
#' @return A \linkS4class{GRanges} object with the regions matching the given criteria.
exclusive.overlap.internal <- function(overlap.matrix, which.factors) {
    has.factor = rep(TRUE, nrow(overlap.matrix))
    if(sum(which.factors) != 0) {
        has.factor = apply(overlap.matrix[,  which.factors, drop=FALSE] >= 1, 1, all)
    }

    no.others = rep(TRUE, nrow(overlap.matrix))
    if(sum(!which.factors) != 0) {
        no.others  = apply(overlap.matrix[, !which.factors, drop=FALSE] == 0, 1, all)
    }

    return(has.factor & no.others)
}

#' Internal inclusive overlap.
#'
#' Returns the set of regions within \code{all.regions} where all
#' conditions/proteins at indices which.factors are present,
#' regardless of whether or not other factors are. Presence/absence
#' at a given locus is obtained from the overlap.matrix parameter.
#'
#' @param all.regions A \linkS4class{GRanges} object with the regions represented
#    by the rows of \code{overlap.matrix}.
#' @param overlap.matrix A matrix, with rows corresponding to regions
#'   in \code{all.rgions}, and columns corresponding to proteins/factors. A
#'   non-zero value in the matrix indicates the protein/factor of interest
#'   is present at this region.
#' @param which.factors Indices of columns where the factor should be present.
#' @return A \linkS4class{GRanges} object with the regions matching the given criteria.
inclusive.overlap.internal <- function(overlap.matrix, which.factors) {
    has.factor = apply(overlap.matrix[,  which.factors, drop=FALSE] >= 1, 1, all)
    return(has.factor)
}

#' Calculate an overlap of certain factors within an intersect.object.
#'
#' Given an \code{intersect.object} , finds all regions where all
#' factors described by either indices or names are present.
#' If both indices and names are \code{NULL}, the inner intersection
#' is calculated.
#'
#' @param intersect.object An \code{intersect.object} returned by \code{\link{build_intersect}}.
#'   by the rows of \code{overlap.matrix}.
#' @param indices A vector of indices into the factors of \code{intersect.object}.
#' @param names A vector of factor names in \code{intersect.object}.
#' @param exclusive If \code{TRUE}, a region will be returned if the factors in
#    indices or names are the ONLY the factors present at that region.
#' @return A \linkS4class{GRanges} object with the regions matching the given criteria.
#' @export
overlaps <- function(x, indices=names(x), exclusive=TRUE) {
    stopifnot(is(x, "GenomicIntersection"))
    
    if(is(indices, "character")) {
        indices = names(x) %in% indices
    } else if(is(indices, "numeric")) {
        logical_indices = rep(FALSE, length(x))
        logical_indices[indices] = TRUE
        indices = logical_indices
    }
    stopifnot(is(indices, "logical"))
    stopifnot(is(exclusive, "logical"))
    
    if(exclusive) {
        return(exclusive.overlap.internal(intersect_matrix(x), indices))
    } else {
        return(inclusive.overlap.internal(intersect_matrix(x), indices))
    }
}

overlaps_regions <- function(x, indices=names(x), exclusive=TRUE) {
    res_indices = overlaps(x, indices, exclusive)
    regions(x)[res_indices]
}

#' Given a \linkS4class{GRangesList} object, determine which items overlap each others.
#'
#' @param grl The \linkS4class{GRangesList} object whose elements need to be overlapped with
#' each others.
#' @param keep.signal Should the values of signal be kept?
#' @return A list with the following elements: \describe{
#' \item{Regions}{A \linkS4class{GRanges} object with all genomic ranges occupied by at least one item.
#'   All ranges are "flattened", so if two of the initial ranges overlapped each other
#'   imperfectly, they are combined into a single region spanning both of them.}
#' \item{Matrix}{A matrix, with \code{ncol=} the number of items in the initial \linkS4class{GRangesList} and \code{nrow=}
#'   the total number of loci, which is equal to the length of \code{Regions}. A value of 1 or more
#'   within the matrix indicates that the regions described by the column overlapped
#'   the region defined by the row.}
#' \item{List}{A list of \code{length(grl)} numeric vectors indicating which indices of \code{Regions} overlap
#'   with the given condition/protein. Useful to translate the regions into unique names
#'   for drawing venn diagrams.}
#' \item{Names}{The names of the initial grl items, corresponding to the column names of \code{Matrix} and the names
#'   of the element of \code{List}.}
#' \item{Length}{The number of items in the initial \linkS4class{GRangesList}, corresponding to the number of columns in \code{Matrix}
#'   and the number of elements in \code{List}}.}
#' @importFrom GenomicRanges reduce
#' @importFrom GenomicRanges mcols
#' @importMethodsFrom GenomicRanges countOverlaps findOverlaps
#' @export
GenomicIntersection <- function(grl, keep.signal = FALSE) {
    # Validate input parameters.
    stopifnot(is(grl, "GRangesList"))
    stopifnot(is(keep.signal, "logical"))
    
    # Flatten the GRangesList so we can get a list of all possible regions.
    all.regions = GenomicRanges::reduce(unlist(grl))

    # Build a matrix to hold the results.
    overlap.matrix <- matrix(0, nrow=length(all.regions), ncol=length(grl))

    if (keep.signal){
      signal.df <- data.frame(matrix(nrow = length(all.regions), ncol = length(grl)))
      colnames(signal.df) <- paste0("signal.", names(grl))
    }

    # Loop over all ranges, intersecting them with the flattened list of all possible regions.
    for(i in 1:length(grl)) {
        overlap.matrix[,i] <- GenomicRanges::countOverlaps(all.regions, grl[[i]], type="any")

        if (keep.signal){
          indices <- GenomicRanges::findOverlaps(all.regions, grl[[i]])
          signal.values.df <- data.frame(from = indices@from, signal = grl[[i]]@elementMetadata@listData$signalValue[indices@to])
          if (nrow(signal.values.df) != 0 & sum(!is.na(signal.values.df$signal)) != 0){
            signal.values.df <- aggregate(signal~from, data = signal.values.df, FUN = mean, na.rm = TRUE)
            signal.df[signal.values.df$from, i] <- signal.values.df$signal
          }
        }
    }
    colnames(overlap.matrix) <-  names(grl)

    if (keep.signal) mcols(all.regions) <- signal.df

    new("GenomicIntersection",
        regions=all.regions, 
        matrix=overlap.matrix)
}

setClass("GenomicIntersection",
         slots=list(regions="GRanges", 
                    matrix="matrix"))

setMethod("names",
          c(x="GenomicIntersection"),
          function(x) {
            colnames(x@matrix)
          })

setMethod("names<-",
          c(x="GenomicIntersection", value="character"),
          function(x, value) {
            colnames(x@matrix) <- value
            x
          })
          
setMethod("length",
          c(x="GenomicIntersection"),
          function(x) {
            ncol(x@matrix)
          })          

regions <- function(x) {
    stopifnot(is(x, "GenomicIntersection"))
    x@regions
}

intersect_matrix <- function(x) {
    stopifnot(is(x, "GenomicIntersection"))
    x@matrix
}

#' Generates a venn diagram from an \code{intersect.object}.
#'
#' @param intersect.object An intersect object returned by \code{\link{build_intersect}}.
#' @param filename A filename for the resulting venn diagram. Pass \code{NULL}
#'   to skip saving to a file.
#' @return The grid object representing the venn.diagram.
#' @importFrom VennDiagram venn.diagram
#' @export
intersect_venn_plot <- function(intersect.object, filename=NULL, title=NULL) {
    if(intersect.object$Length > 5) {
        stop("Cannot plot venn diagram of more than 5 groups!")
    }

    return(VennDiagram::venn.diagram(intersect.object$List,
                                     fill=c("red", "yellow", "green", "blue", "orange")[1:intersect.object$Length],
                                     filename=filename,
                                     print.mode="raw",
                                     main=title))
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

#' Calculates enrichment ratios for quer regions against a genome wide
#' partition of the genome.
#'
#' @param query.regions The regions whose enrichment ratios must be calculated.
#' @param genome.wide The genome partition indicating which part of the genome
#'    fall within which category. Each range should have a 'name' attribute
#'    indicating its category.
#' @param factor.order An optional ordering of the region types for the produced plot.
#' @param file.out An optional file name for a graphical representation of the enrichments.
#' @return A data-frame containing the enrichment values.
#' @export
#' @import GenomicRanges
#' @import ggplot2
region_enrichment <- function(query.regions, genome.wide, genome.order=NULL, file.out=NULL) {
  # Project the query ranges into the genome ranges, so we can
  # know their repartition with basepair precision.
  in.query = project_ranges(query.regions, genome.wide)
  
  # Calculate total base-pair coverages for both the projected query 
  # and the target regions.
  all.region.types = sort(unique(genome.wide$name))
  coverages = matrix(0.0, ncol=2, nrow=length(all.region.types), dimnames=list(all.region.types, c("Query", "Genome")))
  for(region.type in all.region.types) {
      coverages[region.type, "Genome"] = sum(as.numeric(width(reduce(BiocGenerics::subset(genome.wide, name == region.type)))))
      coverages[region.type, "Query"]  = sum(as.numeric(width(reduce(BiocGenerics::subset(in.query, name == region.type)))))
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
  if(is.null(genome.order)) {
    genome.order = all.region.types
  }
  enrichment.df$RegionType = factor(enrichment.df$RegionType, levels=rev(genome.order))
  
  # Plot the results.
  if(!is.null(file.out)) {
    # Replace +/-Inf with NAs.
    enrichment.df$Enrichment[is.infinite(enrichment.df$Enrichment)] <- NA
    
    maxEnrich = max(abs(enrichment.df$Enrichment), na.rm=TRUE)
    ggplot(enrichment.df, aes(fill=Enrichment, y=RegionType, x="Network regions")) +
        geom_tile(color="black") + 
        geom_text(mapping=aes(label=sprintf("%.2f", enrichment.df$Enrichment))) +
        scale_fill_gradient2(low="dodgerblue", mid="white", high="red", midpoint=0, limits=c(-maxEnrich, maxEnrich)) +
        labs(y="Region type", x=NULL) +
        theme(axis.line=element_blank(),
              axis.ticks=element_blank(),
              axis.title=element_blank())
    
    ggsave(file.out, width=7, height=7)
    
    write.table(enrichment.df, file=paste0(file.out, ".txt"), sep="\t", row.names=FALSE, col.names=TRUE)
  }
  
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
                                             genome.order=genome.order, file.out=file.out)
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

#' Obtains consensus regions from a set of regions.
#'
#' @param regions The regions to be consensus'ed.
#' @param keep.signal Whether or not the signal should be kept.
#'
#' @return The collapsed regions.
#' @export
region_consensus <- function(regions, keep.signal=TRUE, fake.signal=FALSE) {
    consensus = intersect_overlap(build_intersect(regions, keep.signal=keep.signal))
    if(keep.signal) {
        mcols(consensus) <- rowMeans(as.data.frame(mcols(consensus)), na.rm = TRUE)
        names(mcols(consensus)) <- "signalValue"
    }
    
    if(!keep.signal && fake.signal) {
        mcols(consensus)$signalValue = NA
    }
    
    return(consensus)
}