import sys
import re
import json
import pandas as pd
import numpy as np

import find_metrics as fm

def read_config_file():
    
    #Generates global variables, reads in json file, assigns data to them based
    #on lines in config file.
    
    print("Reading configuration file...")

    global RESULTS_DIRECTORY
    global GENE_PRIORITISATION_REPORT_DIRECTORY
    global GENE_ANNOTATION_REFERENCE_PATH
    global REGULATORY_ELEMENTS_ANNOTATION_REFERENCE_PATH
    global GENERAL_EXPRESSION_BY_CELL_LINE_REFERENCE_PATH
    global SPECIFIC_EXPRESSION_BY_CELL_LINE_REFERENCE_PATH
    global REFERENCE_GENOME

    global CELL_LINE_OF_INTEREST
    global CHROMOSOMES_OF_INTEREST
    global ENHANCER_EPIGENETIC_FLAGS_OF_INTEREST
    global QUIESCENT_EPIGENETIC_FLAGS_OF_INTEREST

    global SEARCH_TYPE
    global SEARCH_WITHIN_GENE
    global UPSTREAM_SEARCH
    global DOWNSTREAM_SEARCH

    global STD_WEIGHT
    global ANOMALOUS_EXPRESSION_WEIGHT
    global ENHANCER_COUNT_WEIGHT
    global ENHANCER_PROPORTION_WEIGHT
    global CELL_LINE_EXPRESSION_WEIGHT
    global GENE_SIZE_WEIGHT

    global ENHANCER_KERNEL_SHAPE
    global ENHANCER_KERNEL_SIZE_TYPE
    global ABSOLUTE_ENHANCER_KERNEL_SIZE
    global RELATIVE_ENHANCER_KERNEL_SIZE
    global RELATIVE_ENHANCER_KERNEL_SIGMA
    global MIN_ABSOLUTE_ENHANCER_CLUSTER_WIDTH
    global MIN_ENHANCER_CLUSTER_PROMINENCE
    
    global QUIESCENT_KERNEL_SHAPE
    global QUIESCENT_KERNEL_SIZE_TYPE
    global ABSOLUTE_QUIESCENT_KERNEL_SIZE
    global RELATIVE_QUIESCENT_KERNEL_SIZE
    global RELATIVE_QUIESCENT_KERNEL_SIGMA
    global MIN_ABSOLUTE_QUIESCENT_CLUSTER_WIDTH
    global MIN_QUIESCENT_CLUSTER_PROMINENCE

    global SIGMOIDAL_SLOPE
    global SIGMOIDAL_MIDPOINT
    global CELL_LINE_SPECIFIC_EXPRESSION_THRESHOLD
    global INTERFERRING_GENE_OVERLAPS
    
    global ENHANCER_CONVOLUTION
    global QUIESCENT_CONVOLUTION
    global ENHANCER_CONVOLUTION_WEIGHT
    global QUIESCENT_CONVOLUTION_WEIGHT
    global PLATEAU_THRESHOLD
   
    try:
        
        with open(sys.argv[1], "r") as config_file:

            settings = json.load(config_file)

        RESULTS_DIRECTORY = settings["results_directory"]
        GENE_PRIORITISATION_REPORT_DIRECTORY = settings["gene_prioritisation_report_directory"]
        GENE_ANNOTATION_REFERENCE_PATH = settings["gene_annotation_reference"]
        REGULATORY_ELEMENTS_ANNOTATION_REFERENCE_PATH = settings["regulatory_elements_reference"]
        GENERAL_EXPRESSION_BY_CELL_LINE_REFERENCE_PATH = settings["general_expression_by_cell_line_reference_path"]
        SPECIFIC_EXPRESSION_BY_CELL_LINE_REFERENCE_PATH = settings["specific_expression_by_cell_line_reference_path"]
        REFERENCE_GENOME = settings["reference_genome"]

        CELL_LINE_OF_INTEREST = settings["cell_line_of_interest"]
        CHROMOSOMES_OF_INTEREST = settings["chromosomes_of_interest"]
        ENHANCER_EPIGENETIC_FLAGS_OF_INTEREST = settings["enhancer_epigenetic_flags_of_interest"]
        QUIESCENT_EPIGENETIC_FLAGS_OF_INTEREST = settings["quiescent_epigenetic_flags_of_interest"]

        SEARCH_TYPE = settings["search_type"]
        SEARCH_WITHIN_GENE = settings["search_within_gene"]
        UPSTREAM_SEARCH = settings["upstream_search"]
        DOWNSTREAM_SEARCH = settings["downstream_search"]

        STD_WEIGHT = settings["relative_std_weight"]
        ANOMALOUS_EXPRESSION_WEIGHT = settings["relative_anomalous_expression_weight"]
        ENHANCER_COUNT_WEIGHT = settings["relative_enhancer_count_weight"]
        ENHANCER_PROPORTION_WEIGHT = settings["relative_enhancer_proportion_weight"]
        CELL_LINE_EXPRESSION_WEIGHT = settings["relative_cell_line_expression_weight"]
        GENE_SIZE_WEIGHT = settings["relative_gene_size_weight"]

        ENHANCER_KERNEL_SHAPE = settings["enhancer_kernel_shape"]
        ENHANCER_KERNEL_SIZE_TYPE = settings["enhancer_kernel_size_type"]
        ABSOLUTE_ENHANCER_KERNEL_SIZE = settings["absolute_enhancer_kernel_size"]
        RELATIVE_ENHANCER_KERNEL_SIZE = settings["relative_enhancer_kernel_size"]
        RELATIVE_ENHANCER_KERNEL_SIGMA = settings["relative_enhancer_kernel_sigma"]
        MIN_ABSOLUTE_ENHANCER_CLUSTER_WIDTH = settings["min_absolute_enhancer_cluster_width"]
        MIN_ENHANCER_CLUSTER_PROMINENCE = settings["min_enhancer_cluster_prominence"]

        QUIESCENT_KERNEL_SHAPE = settings["quiescent_kernel_shape"]
        QUIESCENT_KERNEL_SIZE_TYPE = settings["quiescent_kernel_size_type"]
        ABSOLUTE_QUIESCENT_KERNEL_SIZE = settings["absolute_quiescent_kernel_size"]
        RELATIVE_QUIESCENT_KERNEL_SIZE = settings["relative_quiescent_kernel_size"]
        RELATIVE_QUIESCENT_KERNEL_SIGMA = settings["relative_quiescent_kernel_sigma"]
        MIN_ABSOLUTE_QUIESCENT_CLUSTER_WIDTH = settings["min_absolute_quiescent_cluster_width"]
        MIN_QUIESCENT_CLUSTER_PROMINENCE = settings["min_quiescent_cluster_prominence"]

        SIGMOIDAL_SLOPE = settings["sigmoidal_slope"]
        SIGMOIDAL_MIDPOINT = settings["sigmoidal_midpoint"]
        CELL_LINE_SPECIFIC_EXPRESSION_THRESHOLD = settings["cell_line_specific_expression_threshold"]
        INTERFERRING_GENE_OVERLAPS = settings["interferring_gene_overlaps"]

        ENHANCER_CONVOLUTION = settings["enhancer_convolution"]
        QUIESCENT_CONVOLUTION = settings["quiescent_convolution"]
        ENHANCER_CONVOLUTION_WEIGHT = settings["enhancer_convolution_weight"]
        QUIESCENT_CONVOLUTION_WEIGHT = settings["quiescent_convolution_weight"]
        PLATEAU_THRESHOLD = settings["plateau_threshold"]
            
    except:
        
        print("ERROR: Config file could not be read.")
        
def read_gene_annotations():
    
    #Assigns the gene annotations gtf file to a pandas dataframe.
    
    print("Reading gene annotations file...")

    try:
        gene_annotations = pd.read_csv(GENE_ANNOTATION_REFERENCE_PATH, sep = "\t",
                                       names = ["Chromosome", "Source", "Type", "Start", "End",
                                                "Score", "Strand", "Phase", "Attributes"],
                                       skiprows = 5,
                                       dtype = {"Chromosome" : str, "Source" : str, "Type" : str,
                                                "Start" : int, "End" : int, "Score" : str, "Strand" : str,
                                                "Phase" : str, "Attributes" : str}
                                       )
    
        return clean_genes(gene_annotations)
        
    except:
        print("ERROR: Gene annotations file could not be read.")
        
def read_general_expression_data():
    
    #Assigns the expression-for-many-cell-types csv file to a pandas dataframe.
        
    print("Reading general expression file...")

    try:
        general_expression_data = pd.read_csv(GENERAL_EXPRESSION_BY_CELL_LINE_REFERENCE_PATH).transpose()
        
        return clean_general_expression_data(general_expression_data)  
    
    except:
        print("ERROR: General expression data could not be read.")
        
def read_specific_expression_data():
    
    #Assigns the expression-for-cell-line-of-interest csv file to a pandas dataframe.
        
    print("Reading specific expression file...")

    try:
        specific_expression_data = pd.read_csv(SPECIFIC_EXPRESSION_BY_CELL_LINE_REFERENCE_PATH,
                                               sep = "\t",
                                               names = ["Gene_name", "Specific_gene_expression"],
                                               skiprows = 1)
        
        return clean_specific_expression_data(specific_expression_data)    
    
    except:
        print("ERROR: Specific expression data could not be read.")
         
def read_regulatory_elements():
    
    #Assigns the regulatory elements gff or bed file to a pandas dataframe.
    
    print("Reading regulatory elements file...")

    try:
        regulatory_elements = pd.read_csv(REGULATORY_ELEMENTS_ANNOTATION_REFERENCE_PATH, sep = "\t",
                                              names = ["Chromosome", "Start", "End", "Flag"])
        
        return regulatory_elements    
    
    except:
        print("ERROR: Could not read regulatory elements file.")
        
def clean_genes(gene_annotations):
    
    #Cleans genetic annotations dataframe by removing chromosomes that are not
    #of interest, removing annotations that are not for genes, finds the gene
    #biotype and removes non-protein coding genes extracts the gene name as an
    #attribute, drop irrelevant columns, drops dulpicate genes
    
    print("Cleaning gene data...")
    
    gene_annotations = gene_annotations.loc[gene_annotations["Chromosome"].isin(CHROMOSOMES_OF_INTEREST)]
    gene_annotations = gene_annotations.drop(gene_annotations[gene_annotations["Type"] != "gene"].index)
    gene_annotations["Gene_biotype"] = gene_annotations["Attributes"].apply(lambda x : re.findall("gene_biotype \"(.*?)\"", x)[0] if re.search("gene_name \"(.*?)\"", x) != None else "None")
    gene_annotations = gene_annotations.drop(gene_annotations[gene_annotations["Gene_biotype"] != "protein_coding"].index)
    gene_annotations["Gene_name"] = gene_annotations["Attributes"].apply(lambda x : re.findall("gene_name \"(.*?)\"", x)[0] if re.search("gene_name \"(.*?)\"", x) != None else "None")
    gene_annotations = gene_annotations.drop(["Source", "Type", "Score", "Phase", "Attributes", "Gene_biotype"], axis = 1)
    gene_annotations = gene_annotations.drop_duplicates(keep = False, subset = ["Gene_name"])
    gene_annotations = gene_annotations.rename(columns = {"Start" : "Gene_start", "End" : "Gene_end"})
    
    return gene_annotations

def clean_general_expression_data(general_expression_data):
    
    #Cleans the general expression data by dropping irrelevant non-numeric
    #columns, converting the first row to headers, setting gene name as index,
    #renaming columns, finding the means and standard deviation of expression
    #for each gene and the z-score of the cell line of interest's expression for
    #each gene, dropping irrelevant columns, and dropping duplicate genes
    
    print("Cleaning general expression data...")

    general_expression_data = general_expression_data.drop(["depmapID", "primary_disease"], axis = 0)
    general_expression_data.columns = general_expression_data.iloc[0]
    general_expression_data = general_expression_data.iloc[1:]
    general_expression_data = general_expression_data.reset_index().rename(columns = {"index" : "Gene_name"})
    general_expression_data = general_expression_data.rename(columns = {CELL_LINE_OF_INTEREST : "General_gene_expression"})
    general_expression_data = fm.find_mean(general_expression_data)
    general_expression_data = fm.find_std(general_expression_data)
    general_expression_data = fm.find_anomalous_score_of_gene_expression(general_expression_data)
    general_expression_data = general_expression_data[["Gene_name", "General_gene_expression", "Mean", "Std", "Anomalous_score"]]
    general_expression_data = general_expression_data.drop_duplicates(keep = False, subset = ["Gene_name"])
    
    return general_expression_data

def clean_specific_expression_data(specific_expression_data):
    
    #Cleans the expression data specific to the cell line of interest by turning
    #minus infinite strings into a floating point representation of negative
    #infinity, duplicate genes are dropped
    
    print("Cleaning specific expression data...")
    
    #specific_expression_data = specific_expression_data.drop(specific_expression_data[specific_expression_data["Specific_gene_expression"] == "-Inf"].index)
    specific_expression_data["Specific_gene_expression"] = specific_expression_data["Specific_gene_expression"].apply(lambda expression : 0 if expression == "-Inf" else pow(2, expression))
    #specific_expression_data["Specific_gene_expression"] = specific_expression_data["Specific_gene_expression"].apply(lambda expression : np.NINF if expression == "-Inf" else expression)
    specific_expression_data = specific_expression_data.drop_duplicates(keep = False, subset = ["Gene_name"])
    
    return specific_expression_data

def clean_regulatory_elements(regulatory_elements):
    
    #Clean the regulatory elements dataframe by removing "chr" from chromosome
    #strings, including only regulatory elements from within flags of interest,
    #separately creating a quiescent dataframe in the same
    #way
    
    print("Cleaning regulatory elements data...")
    
    regulatory_elements["Chromosome"] = regulatory_elements["Chromosome"].apply(lambda x : x[3:])
        
    enhancers = regulatory_elements[regulatory_elements["Flag"].isin(ENHANCER_EPIGENETIC_FLAGS_OF_INTEREST)]
    enhancers = enhancers.drop(["Flag"], axis = 1)

    quiescent_regions = regulatory_elements[regulatory_elements["Flag"].isin(QUIESCENT_EPIGENETIC_FLAGS_OF_INTEREST)]
    quiescent_regions = quiescent_regions.drop(["Flag"], axis = 1)
        
    return enhancers, quiescent_regions