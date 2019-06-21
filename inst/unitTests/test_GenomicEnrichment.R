if(FALSE) {
    library( "RUnit" )
    library( "GenomicRanges" )
    library( "GenomicOperations" )
}

test1 = rtracklayer::import(system.file("extData/enrichment_test_1.bed", package="GenomicOperations"))
test2 = rtracklayer::import(system.file("extData/enrichment_test_2.bed", package="GenomicOperations"))
partition = rtracklayer::import(system.file("extData/partition.bed", package="GenomicOperations"))

test.GenomicEnrichment <- function() {
    test_grl = GRangesList(Test1=test1, Test2=test2)
    test = GenomicEnrichment(test_grl, partition)

    # Test basic accessors.
    checkIdentical(names(test), names(test_grl))
    checkIdentical(length(test), length(test_grl))
    
    # Make sure the resulting data-frames have the right columns.
    checkIdentical(colnames(coverage_df(test)),
                   c("Genome", names(test_grl)))
    checkIdentical(colnames(proportion_df(test)),
                   c("Genome", names(test_grl)))                   
    checkIdentical(colnames(enrichment_df(test)),
                   names(test_grl))
                   
    # Make sure the numerical values match:
    # 1a. Genomic coverage
    tssA_coverage_g = sum(width(partition[partition$name=="1_TssA"]))               
    checkTrue(coverage_df(test)["1_TssA", "Genome"] == tssA_coverage_g)
    # 1b. Query coverage
    projected = project_ranges(test_grl[["Test1"]], partition)
    tssA_coverage_q = sum(width(projected[projected$name=="1_TssA"]))
    checkTrue(coverage_df(test)["1_TssA", "Test1"] == tssA_coverage_q)
    
    # 2a. Genomic proportion
    tssA_proportion_g = tssA_coverage_g / sum(width(partition))
    checkTrue((proportion_df(test)["1_TssA", "Genome"] - tssA_proportion_g) < 0.000001)
    # 2b. Query proportion
    tssA_proportion_q = tssA_coverage_q / sum(width(test_grl[["Test1"]]))
    checkTrue((proportion_df(test)["1_TssA", "Test1"] - tssA_proportion_q) < 0.000001)
    
    # 3. Enrichment
    tssA_enrichment = log2(tssA_proportion_q / tssA_proportion_g)
    checkTrue(enrichment_df(test)["1_TssA", "Test1"] - tssA_enrichment < 0.0000001)

    # Check that plotting works.
    print(coverage_plot(test))
    print(proportion_plot(test))
    print(enrichment_plot(test))
    
    
}
