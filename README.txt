A tool to help the Genome Scramble project.

The "config.json" file can be used to change the settings of the tool:

    results_directory

        The directory various results will be saved to, including plots and some tables of data.

    gene_prioritisation_report_directory

        The directory the gene prioritisation report is saved to, eventually this will likely be the same as the results directory.

    gene_annotation_reference

        The file used to fetch the gene annotation data, currently a gtf file from ensembl of all genes.

    regulatory_elements_reference

        The file used to fetch the epigenetic flags for regions along the genome, currently a bed file.

    general_expression_by_cell_line_reference_path

        The file used to fetch the expression of all genes in many cell line.

    specific_expression_by_cell_line_reference_path

        The file used to fetch the expression of all genes in the cell line we are specifically interested in.
    
    cell_line_of_interest

        The cell line we are interested in.

    chromosomes_of_interest
    
        An array of the chromosomes that we wish to look at with the tool.
    
    enhancer_epigenetic_flags_of_interest

        An array of the epigenetic flags which will be treated as enhancers (has only previously been tried with single value).

    quiescent_epigenetic_flags_of_interest

        An array of the epigenetic flags which will be treated as quiescent regions (has only previously been tried with single value).

    search_type

        Can be the string "whole_gene" or "start_site".
        When creating a search window, tells the tool whether upstream and downstream are in reference to the start and end of the gene, or just the start fo the gene.
        Suggested: "whole_gene", may cause errors with "start_site".

    search_within_gene

        Can be true or false.
        When creating a search window, tells the tool whether to include the gene itself within the search window.
        Suggested: true, may cause errors if false.

    upstream_search

        Number of bases upstream to create search window

    downstream_search

        Number of bases downstream to create search window

    relative_std_weight

        Weight given to standard deviation of expression of genes across cell lines in interest scoring.

    relative_anomalous_expression_weight

        Weight given to z score of gene expression in cell line of interest in comparison to other cell lines in interest scoring.

    relative_enhancer_count_weight

        Weight given to number of enhancers in search window in interest scoring.

    relative_enhancer_proportion_weight
    
        Weight given to proportion of bases which are flagged as part of enhancers within search window, in interest scoring.

    relative_cell_line_expression_weight

        Weight given to the expression of gene within cell line of interest, in interest scoring.

    relative_gene_size_weight

        Weight given to the size of the gene in interest scoring.
        Smaller genes score higher, larger genes asymptotically tend to 0.

    enhancer_kernel_shape

        Can be "guassian" or "flat".
        Shape of kernel used to convolve enhancer signal.

    enhancer_kernel_size_type

        Can be "relative" or "absolute".
        Tells the tool how to define the size of the convolution kernel.

    absolute_enhancer_kernel_size

        Size of the convolution kernel in bases.

    relative_enhancer_kernel_size

        Size of the convolution kernel as a proportion of the gene which created the search window being convolved.

    relative_enhancer_kernel_sigma

        Sigma of the convolution kernel as a proportion of the gene which created the search window being convolved.

    min_absolute_enhancer_cluster_width

        Depreciated, will remove.

    min_enhancer_cluster_prominence

        Depreciated, will remove.

    quiescent_kernel_shape

        Can be "guassian" or "flat".
        Shape of kernel used to convolve enhancer signal.
        Recommend keeping this the same as the enhancer version for now.

    quiescent_kernel_size_type

        Can be "relative" or "absolute".
        Tells the tool how to define the size of the convolution kernel.
        Recommend keeping this the same as the enhancer version for now.

    absolute_quiescent_kernel_size

        Size of the convolution kernel in bases.
        Recommend keeping this the same as the enhancer version for now.

    relative_quiescent_kernel_size

        Size of the convolution kernel as a proportion of the gene which created the search window being convolved.
        Recommend keeping this the same as the enhancer version for now.

    relative_quiescent_kernel_sigma

        Sigma of the convolution kernel as a proportion of the gene which created the search window being convolved.
        Recommend keeping this the same as the enhancer version for now.

    min_absolute_quiescent_cluster_width

        Depreciated, will remove.

    min_quiescent_cluster_prominence

        Depreciated, will remove.

    sigmoidal_slope

        Depreciated, will remove.

    sigmoidal_midpoint

        Depreciated, will remove.

    cell_line_specific_expression_threshold

        The threshold of cell line specific gene expression required by a gene to count as "interferring" and thus block the search window.
        Currently this is measured as a log2 of the gene expression.

    interferring_gene_overlaps

        true of false.
        Whether to allow any gene to be searched in which an interferring gene overlaps with the target gene.

    enhancer_convolution

        Once ranked by interest score, the tool will perform convolutions on this many genes.

    quiescent_convolution

        Whether to convolve quiescent regions as well as enhancers.
        Currently experimental.

    enhancer_convolution_weight

        Not yet in use.

    quiescent_convolution_weight

        Not yet in use.

    plateau_threshold

        Not yet in use.

