### {{{ --- Test setup ---

if(FALSE) {
    library( "RUnit" )
    library( "GenomicRanges" )
    library( "GenomicOperations" )
}

### }}}

test_list = GRangesList(List1=GRanges(c("1:1000-2000:+",
                                        "1:1900-3000:-",
                                        "1:4000-5000:+",
                                        "1:6000-7000:*")),
                        List2=GRanges(c("1:1000-1100:-",
                                        "1:1100-1200:+",
                                        "1:5900-6100:+")),
                        List3=GRanges(c("1:4500-4600:+",
                                        "1:8000-8900:+")))
                                        
test.GenomicIntersection <- function() {
    # Build the GenomicIntersection object.
    test = GenomicOverlaps(test_list)

    # Make sure its "dimensions" match the input GRangesList
    checkIdentical(length(test), length(test_list))
    checkIdentical(names(test), names(test_list))
    
    # Make sure the reduced regions are okay.
    reduced_regions = reduce(unlist(test_list))
    checkIdentical(regions(test), reduced_regions)
        
    # Test overlap calculations.
    overlap_1 = countOverlaps(reduced_regions, test_list$List1) > 0
    overlap_2 = countOverlaps(reduced_regions, test_list$List2) > 0
    overlap_3 = countOverlaps(reduced_regions, test_list$List3) > 0
    
    checkIdentical(intersect_indices(test, "List1", exclusive=FALSE), overlap_1)
    checkIdentical(intersect_indices(test, "List2", exclusive=FALSE), overlap_2)
    checkIdentical(intersect_indices(test, "List3", exclusive=FALSE), overlap_3)
    
    checkIdentical(intersect_indices(test, c("List1", "List2"), exclusive=FALSE), overlap_1 & overlap_2)
    checkIdentical(intersect_indices(test, c("List1", "List3"), exclusive=FALSE), overlap_1 & overlap_3)
    checkIdentical(intersect_indices(test, c("List2", "List3"), exclusive=FALSE), overlap_2 & overlap_3)
    
    checkIdentical(intersect_indices(test, c("List1", "List2", "List3"), exclusive=FALSE), overlap_1 & overlap_2 & overlap_3)    
    
    checkIdentical(intersect_indices(test, "List1", exclusive=TRUE), overlap_1 & !overlap_2 & !overlap_3)
    checkIdentical(intersect_indices(test, "List2", exclusive=TRUE), !overlap_1 & overlap_2 & !overlap_3)
    checkIdentical(intersect_indices(test, "List3", exclusive=TRUE), !overlap_1 & !overlap_2 & overlap_3)
    
    checkIdentical(intersect_indices(test, c("List1", "List2"), exclusive=TRUE), overlap_1 & overlap_2 & !overlap_3)
    checkIdentical(intersect_indices(test, c("List1", "List3"), exclusive=TRUE), overlap_1 & !overlap_2 & overlap_3)
    checkIdentical(intersect_indices(test, c("List2", "List3"), exclusive=TRUE), !overlap_1 & overlap_2 & overlap_3)
    
    checkIdentical(intersect_indices(test, c("List1", "List2", "List3"), exclusive=TRUE), overlap_1 & overlap_2 & overlap_3)        

    plot_venn(test)
}