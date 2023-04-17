# AEG SCRAMBLE
## Introduction
The Automated Enhancer-Genome SCRAMBLEr is a tool to aid the automated identification of suitable regions to insert LoxPSym sites within a genome. The tool first finds genes that are suitable based on a few criteria and prioritises them. The locations of elements around selected genes are then used to identify suitable locations to insert LoxPSym sites.

The Genome SCRAMBLE project is run by Jonas Koeppel under the supervision of Leopold Parts, and will attempt to introduce controlled breaks between clusters of enhancers near genes and recombine the piece in a random manner to investigate phenotypic changes.

The milestones of this project are:
* Prioritise the genes of interest
* Identify regulatory elements around the gene
* Finding pegRNA insertion sites at suitable locations for LoxPSym sites
* Designing pegRNA sequences

## Configuration
Running the file requires a configuration JSON file. This file sets out a number of variables that can change how the program runs, as described below:
* results_directory - This is used to set the path of the directory you wish for figures and output BED and WIG files to be saved in.
* gene_prioritisation_report_directory - Similarly, this is used to set the path of the directory you wish the report of prioritised genes to be saved in.
* gene_annotation_reference - The path to the file the tool will use as its reference for gene annotations, this is in the format of a GTF file that is available from https://www.ensembl.org/Homo_sapiens/Info/Index
* regulatory_elements_reference - The path to the file the tool will use as its reference for regulatory element annotations, this is in the format of a BED file which splits the genome into regions and assigns each one from a set of epigenetic flags.
* general_expression_by_cell_line_reference_path - The path to the file the tool will use as its reference for expression of many genes in many cell lines, the format of this is a CSV file obtainable at https://sites.broadinstitute.org/ccle/datasets
* specific_expression_by_cell_line_reference_path - The path to the file the tool will use as its reference for expression of many genes in the cell line you are specifically interested in (this is the same information as in the general_expression_by_cell_line_reference_path, but more specfic to the cell line you are interested in).
* reference_genome - The path to the file the tool will use as its reference for the genome used to provide sequences of regions. The format is a FASTA file obtainable from https://www.ensembl.org/Homo_sapiens/Info/Index

* cell_line_of_interest - This is the cell line which the expression of which will be used to select genes, it should be a string in the same format as the name given in general_expression_by_cell_line_reference_path
* chromosomes_of_interest - An array of strings, each string is a chromosome for which data will be included, the strings here do not have "Chr" prepended.
* enhancer_epigenetic_flags_of_interest - An array of strings, each string is an epigenetic flag that will be associated with being an enhancer.
* quiescent_epigenetic_flags_of_interest - An array of strings, each string is an epigenetic flag that will be associated with being a quiescent region.

* search_type - When generating search windows for each gene, this tells the tool whether to count bases from just the start site of transcription, of from the start and end of the gene. The value should be "start_site" or "whole_gene respectively.
* upstream_search - When generating search windows for each gene, tells the tool how many bases upstream of the gene to search.
* downstream_search - When generating search windows for each gene, tells the tool how many bases downstream of the gene to search.

* std_hard_filter_max - A threshold to cut off genes with a standard deviation above this variable.
* std_hard_filter_min - A threshold to cut off genes with a standard deviation below this variable.
* anomalous_expression_hard_filter_max - A threshold to cut off genes with an anomalous z-score above this variable.
* anomalous_expression_hard_filter_min - A threshold to cut off genes with an anomalous z-score below this variable.
* enhancer_count_hard_filter_max - A threshold to cut off genes with a count of enhancers within the search window above this variable.
* enhancer_count_hard_filter_min - A threshold to cut off genes with a count of enhancers within the search window below this variable.
* enhancer_proportion_hard_filter_max - A threshold to cut off genes with a proportion of enhancers within the search window above this variable.
* enhancer_proportion_hard_filter_min - A threshold to cut off genes with a proportion of enhancers within the search window below this variable.
* cell_line_expression_hard_filter_max - A threshold to cut off genes with an expression within the cell line of interest above this variable.
* cell_line_expression_hard_filter_min - A threshold to cut off genes with an expression within the cell line of interest below this variable.
* gene_size_hard_filter_max - A threshold to cut off genes larger than this variable.
* gene_size_hard_filter_min - A threshold to cut off genes smaller than this variable.

* relative_std_weight - Weight given to standard deviation of expression of genes across cell lines, when interest scoring.
* relative_anomalous_expression_weight - Weight given to z score of gene expression in cell line of interest in comparison to other cell lines, when interest scoring.
* relative_enhancer_count_weight - Weight given to number of enhancers found within search window, when interest scoring.
* relative_enhancer_proportion_weight - Weight given to proportion of bases which are flagged as part of enhancers within search window, when interest scoring.
* relative_cell_line_expression_weight - Weight given to the expression of genes within the cell line of interest, when interest scoring.
* relative_gene_size_weight - Weight given to the size of the gene, when interest scoring. Smaller genes score higher, larger genes asymptotically tend to 0.

* enhancer_kernel_shape - Can be "guassian" or "flat"; shape of kernel used to convolve enhancer signal.
* enhancer_kernel_size_type - Can be "relative" or "absolute"; tells the tool how to define the size of the convolution kernel. Absolute not yet implemented.
* absolute_enhancer_kernel_size - Size of the convolution kernel in bases, used when enhancer_kernel_size_type is "absolute". Not yet implemented.
* absolute_enhancer_kernel_sigma - Sigma of the convolution kernel. Not yet implemented.
* relative_enhancer_kernel_size - Size of the convolution kernel as a proportion of the gene which created the search window being convolved, used when enhancer_kernel_size_type is "relative".
* relative_enhancer_kernel_sigma - Sigma of the convolution kernel as a proportion of the gene which created the search window being convolved.

* quiescent_kernel_shape - Can be "guassian" or "flat"; shape of kernel used to convolve quiescent signal.
* quiescent_kernel_size_type - Can be "relative" or "absolute"; tells the tool how to define the size of the convolution kernel. Absolute not yet implemented.
* absolute_quiescent_kernel_size - Size of the convolution kernel in bases, used when quiescent_kernel_size_type is "absolute". Not yet implemented.
* absolute_quiescent_kernel_sigma - Sigma of the convolution kernel. Not yet implemented.
* relative_quiescent_kernel_size - Size of the convolution kernel as a proportion of the gene which created the search window being convolved, used when quiescent_kernel_size_type is "relative".
* relative_quiescent_kernel_sigma - Sigma of the convolution kernel as a proportion of the gene which created the search window being convolved.

* cell_line_specific_expression_threshold - The threshold of cell line specific gene expression required by a gene to count as "interferring" and thus block a search window.
* interferring_gene_overlaps - When searching for elements that overlap gene search windows, this boolean value tells the tool whether to include those which overlap with the gene itself. Currently being revised.

* convolution_limit - Convolution applies to each gene in order of prioritisation according to its interest score; this will variable is how many genes convolution will be applied to as it is resource intensive.
* enhancer_convolution_weight - If enhancer and quiescent convolution are combined, this may weight the contibution of the enhancer convolution.
* quiescent_convolution_weight - If enhancer and quiescent convolution are combined, this may weight the contibution of the quiescent convolution.
* plateau_threshold - The threshold used when calling whether the signal of a convolution qualifies a base as being part of a plateau.
