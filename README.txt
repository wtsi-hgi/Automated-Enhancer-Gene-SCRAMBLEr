A tool to help the Genome Scramble project.

Gene metrics - First milestone is to compare different genes, providing differente relavent attributes
    So far the metrics used are:
    - The gene size
    - The expression of the gene in the cell line of interest
    - The mean expression of the gene in many cell lines
    - The standard deviation of the expression of the gene in many cell lines
    - The number of enhancers near the gene, within a defined search area
    - The density of enhancers withing the defined search window


    Metrics I would like to add are:
    - The number of clusters of enhancers within the defined search window
    - The density of clusters of enhancers within the defined search window
    - An outlier score for expression of the gene within the cell line of interest compared to all other cell lines
    - Nearby other genes
    - Nearby enhancers in 3D space

The config file can add the following settings:
    results_directory                       =   Path to directory to save figures and tables produced
    gene_annotation_reference               =   Path to gtf file with gene annotations of genome
    regulatory_elements_reference           =   Path to gff or bed file with annotations of regulatory elements
    cell_lines_expression_reference         =   Path to csv file with expression of each gene in each cell line
    cell_line_of_interest                   =   The cell line that will be compared against the others
    epigenetic_flags_of_interest            =   The flag in the bed file which will define the regulatory regions identified
    search_type                             =   "whole_gene" or "start_site" incidating whether to search downstream from the start or end of the gene
    search_within_gene                      =   Boolean value indicating whether to include the gene itself within the search
    upstream_search                         =   Number of bases to search upstream
    downstream_search                       =   Number of bases to search downstream
    relative_non_housekeeping_weight        =   Weight assigned to standard deviation of expression of gene in many cell types, when scoring genes, nominally 1
    relative_enhancer_count_weight          =   Weight assigned to number of enhancers found within search area, when scoring genes, nominally 1
    relative_enhancer_proportion_weight     =   Weight assigned to density of enhancer regions within search area, when scoring genes, nominally 1
    relative_cell_line_expression_weight    =   Weight assigned to expression of gene within cell line of interest, when scoring genes, nominally 1
    relative_gene_size_weight               =   Weight assigned to the size of the gene, when scoring genes, nominally 1
    kernel_size_type                        =   "absolute" or "relative", how to define the size of the kernel used for convolution of regulatory elements
    absolute_kernel_size                    =   Size of kernel used when convolving regulatory elements within search area, in bases
    relative_kernel_size                    =   Size of kernel used when convolving regulatory elements within search area, compared to the size of the gene which generated search area, recommended ~ 0.01 - 0.2
    kernel_shape                            =   Shape of kernel used when convoling regulatory elements within search area, "flat" or "guassian"
    sigmoidal_slope                         =   Further information soon
    sigmoidal_midpoint                      =   Further information soon

