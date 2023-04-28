import pandas as pd
import pyranges as pr
from sklearn.preprocessing import StandardScaler
import hashlib
from scipy import stats
import numpy as np

import data_initialisation as di
import region_convolutions as rc
import data_visualisation as dv

interesting_features = [
    "Std", "Anomalous_score", "Enhancer_count",
    "Enhancer_proportion", "Specific_gene_expression",
    "Gene_size", "Symmetry_ratio"
    ]

def find_mean(expression_data):
    
    # Adds the mean of gene expression to the given expression data frame,
    # giving a mean for each gene
    
    print("Finding mean of expression for each gene...")
    
    expression_data["Mean"] = expression_data.loc[
        :, expression_data.columns != "Gene_name"].mean(axis = 1)
    
    return expression_data
    
def find_std(expression_data):
    
    #Adds the standard deviation to the given expression dataframe, giving the
    #standard deviation across expression of each gene in all provided cell
    #types
    
    print("Finding standard deviation of expression for each gene...")
    
    expression_data["Std"] = expression_data.loc[
        :, expression_data.columns != "Gene_name"].std(axis = 1)
    
    return expression_data

def find_anomalous_score_of_gene_expression(expression_data):
    
    # Adds the z-score of each gene based on its expression
    # in the cell line of interest compared to all others
    
    print("Finding anomalies...")
    
    expression_data["Anomalous_score"] = expression_data.\
        apply(lambda gene : (gene["General_gene_expression"] - 
                             gene["Mean"]) / gene["Std"], axis = 1)

    return expression_data

def find_gene_sizes(genes):
    
    #Adds the size of each gene based on its start and end point
    
    print("Finding gene sizes...")
    
    genes["Gene_size"] = genes["Gene_end"] - genes["Gene_start"]
    
    return genes

def find_interferring_genes(gene_data):
    
    #For each gene, finds the nearest gene upstream and downstream
    
    print("Finding interferring genes...")
    
    #"Start" and "End" are generated for use with PyRanges
    #module which requires temporary columns "Start" and "End",
    #they are removed at the end of function
    
    gene_data["Start"] = gene_data["Gene_start"]
    gene_data["End"] = gene_data["Gene_end"]
    
    interferring_genes_search = pr.PyRanges(
        gene_data.loc[gene_data["Specific_gene_expression"] > \
            di.CELL_LINE_SPECIFIC_EXPRESSION_THRESHOLD])
    
    gene_search = pr.PyRanges(gene_data)
    genes_nearest_upstream_pr = gene_search.nearest(
        interferring_genes_search, how = "upstream", 
        suffix = "_upstream_interferring_gene", 
        overlap = di.INTERFERRING_GENE_OVERLAPS)
    
    genes_nearest_downstream_pr = gene_search.nearest(
        interferring_genes_search, how = "downstream", 
        suffix = "_downstream_interferring_gene", 
        overlap = di.INTERFERRING_GENE_OVERLAPS)

    genes_nearest_upstream = genes_nearest_upstream_pr.df
    genes_nearest_downstream = genes_nearest_downstream_pr.df
    
    gene_data = pd.merge(
        gene_data, 
        genes_nearest_upstream.loc[
            :, 
            ["Gene_name",
             "Start_upstream_interferring_gene",
             "End_upstream_interferring_gene",
             "Gene_name_upstream_interferring_gene"
            ]
        ],
        on = "Gene_name",
        how = "inner"
    )
    
    gene_data = pd.merge(
        gene_data, 
        genes_nearest_downstream.loc[
            :, 
            ["Gene_name",
             "Start_downstream_interferring_gene",
             "End_downstream_interferring_gene",
             "Gene_name_downstream_interferring_gene"
            ]
        ],
        on = "Gene_name",
        how = "inner"
    )
    
    if not di.INTERFERRING_GENE_OVERLAPS:
        gene_data = gene_data.loc[
            gene_data["End_upstream_interferring_gene"] < \
                gene_data["Gene_start"]]
        gene_data = gene_data.loc[
            gene_data["Start_downstream_interferring_gene"] > \
                gene_data["Gene_end"]]
        
    gene_data.drop(["Start", "End"], axis = 1) 
    
    return gene_data

def find_search_windows(genes):

    #Defines a search window for each gene, based on the number of bases
    #upstream and downstream specified, the starting point of the window, and
    #whether the gene itself is included. Search window is foreshorterned if an
    #interferring gene is within the window, or if the window would include
    #negative bases.

    print("Finding sites to define search window...") 
    
    if (di.SEARCH_TYPE == "whole_gene"):
        downstream_search_start = "Gene_end"
        
    elif (di.SEARCH_TYPE == "start_site"):
        downstream_search_start = 'Gene_start'
        
    else: 
        print("ERROR : Invalid search type.")
        
    genes["Search_window_start"] = genes.apply(
        lambda gene : gene["Gene_start"] - di.UPSTREAM_SEARCH 
        if gene["Strand"] == "+" 
        else gene["Gene_start"] - di.DOWNSTREAM_SEARCH, axis = 1)
    
    genes["Search_window_end"] = genes.apply(
        lambda gene : gene[downstream_search_start] + di.DOWNSTREAM_SEARCH 
        if gene["Strand"] == "+" 
        else gene[downstream_search_start] + di.UPSTREAM_SEARCH, axis = 1)
    
    genes["Search_window_start"] = genes.apply(
        lambda gene : 0 
        if gene["Gene_start"] < 0 
        else gene["Search_window_start"], axis = 1)
    
    genes["Search_window_start"] = genes.apply(
        lambda gene : gene["End_upstream_interferring_gene"] 
        if gene["Search_window_start"] < \
            gene["End_upstream_interferring_gene"] 
        else gene["Search_window_start"], axis = 1)
    
    genes["Search_window_end"] = genes.apply(
        lambda gene : gene["Start_downstream_interferring_gene"] 
        if gene["Search_window_end"] > \
            gene["Start_downstream_interferring_gene"] 
        else gene["Search_window_end"], axis = 1)
        
    genes["Search_window_size"] = \
        (genes["Search_window_end"] - genes["Search_window_start"])
                
    return genes
    
def find_element_overlaps_within_search_window(elements, genes):
    
    #PyRanges is used to find specified element type overlaps within search
    #window given for each gene. "Start" and "End" are generated for PyRanges
    #and removed subsequently. Overlaps are stored in a new dataframe
    
    print("Searching for overlapping elements...")
    
    genes["Start"] = genes["Search_window_start"]
    genes["End"] = genes["Search_window_end"]
    
    gene_search = pr.PyRanges(genes)
    elements_search = pr.PyRanges(elements)
    overlaps = gene_search.intersect(elements_search, strandedness = False)
    overlaps = overlaps.df
    
    genes.drop(["Start", "End"], axis = 1)
    
    return overlaps
    
def count_overlaps_per_gene(genes, overlaps, element_type):
    
    #Number of specified element overlaps are counted for each gene.
    
    print("Counting overlaps...")

    #overlaps.drop(["Start", "End"], axis = 1)
    genes = pd.merge(
        genes, 
        overlaps.groupby("Gene_name").size().reset_index(
        name = (element_type + "_count")), 
            on = "Gene_name", 
            how = "inner"
        )
    
    return genes
    
def find_nearby_enhancer_densities(gene_data, overlaps):

    # Density of specifed element overlaps within the given search window is
    # calculated for each gene.
    
    print("Finding proximal enhancer densities...")
    
    overlaps["Enhancer_proportion"] = (
        overlaps.loc[:, "End"] - overlaps.loc[:, "Start"]) /\
            overlaps.loc[:, "Search_window_size"]
    
    overlaps = overlaps.loc[
        :, 
        ["Gene_name", "Enhancer_proportion"]].groupby(
        ["Gene_name"], 
        as_index = False)["Enhancer_proportion"].sum().reset_index()
    
    gene_data = pd.merge(gene_data, overlaps, on = "Gene_name")

    return gene_data
    
def find_symmetry_of_elements(gene_data, overlaps):
    
    # For each gene calculates the number of elements on one side compared
    # to the other and produces a score of symmetry
    
    print("Assessing symmetry of elements near genes...")
    
    overlaps["Overlaps_upstream"] = overlaps\
        .apply(lambda overlap : True if overlap["End"] <= 
               overlap["Gene_start"] else False, axis = 1)
        
    overlaps_upstream = overlaps.loc[:, ["Gene_name", "Overlaps_upstream"]]\
        .groupby(["Gene_name"], as_index = False)["Overlaps_upstream"]\
            .sum().reset_index()
            
    gene_data = pd.merge(gene_data, overlaps_upstream, on = "Gene_name")
    
    gene_data["Symmetry_ratio"] = \
        1 - (2 * (np.absolute((gene_data["Overlaps_upstream"] /
                          gene_data["Enhancer_count"]) - 0.5)))

    gene_data.drop(["Overlaps_upstream"], axis = 1)
    
    return gene_data
    
def calculate_interest_score(gene_data):
    
    # Various attributes of each gene are scaled and normallised, before being
    # weighted and combined into an interest score.
    # The export_gene_scores_report function is called
    
    print("Scoring genes...")
    
    scaler = StandardScaler()
    scaled_genes = gene_data.loc[:, (["Gene_name"] + interesting_features)]
    scaled_genes.loc[:, interesting_features] = \
        scaler.fit_transform(scaled_genes.loc[:, interesting_features])
    
    for feature in interesting_features:
        gene_data["Z-" + feature] = stats.zscore(gene_data[feature])
    
    dv.compare_metrics(
        scaled_genes, 
        "Comparison of Metrics within Z-space",
        "metrics_comparison"
    )
    
    scaled_genes = scaled_genes.assign(
        Interest_score = (
            scaled_genes["Std"] * di.STD_WEIGHT +
            scaled_genes["Anomalous_score"] * di.ANOMALOUS_EXPRESSION_WEIGHT +
            scaled_genes["Enhancer_count"] * di.ENHANCER_COUNT_WEIGHT +
            scaled_genes["Enhancer_proportion"] * \
                di.ENHANCER_PROPORTION_WEIGHT +
            scaled_genes["Specific_gene_expression"] * \
                di.CELL_LINE_EXPRESSION_WEIGHT +
            (di.GENE_SIZE_WEIGHT * pow((2), (-scaled_genes["Gene_size"] * \
                di.GENE_SIZE_WEIGHT * di.GENE_SIZE_WEIGHT))) +
            scaled_genes["Symmetry_ratio"] * di.SYMMETRY_WEIGHT
        )
    ).sort_values("Interest_score", ascending=False)
    
    scaled_genes = scaled_genes.rename(
        columns = {
            "Std" : "Scaled_std",
            "Anomalous_score" : "Scaled_anomalous_score",
            "Enhancer_count" : "Scaled_enhancer_count",
            "Enhancer_proportion" : "Scaled_enhancer_proportion",
            "Specific_gene_expression" : "Scaled_specific_gene_expression",
            "Gene_size" : "Scaled_gene_size",
            "Symmetry_ratio" : "Scaled_symmetry_ratio"
        }
    )
    
    gene_data = pd.merge(
        gene_data, 
        scaled_genes.loc[:, [
            "Gene_name", 
            "Interest_score", 
            "Scaled_std",
            "Scaled_anomalous_score", 
            "Scaled_enhancer_count",
            "Scaled_enhancer_proportion",
            "Scaled_specific_gene_expression", 
            "Scaled_gene_size",
            "Scaled_symmetry_ratio"
            ]
        ],
        on = "Gene_name"
    )
    gene_data = iterate_through_hard_filters(gene_data)
    gene_data = gene_data.sort_values(
        "Interest_score", ascending = False).reset_index()
    
    for feature in interesting_features:
        gene_data["Z-" + feature] = stats.zscore(gene_data[feature])
    
    return gene_data

def iterate_through_hard_filters(gene_data):

    # Calls apply_hard_filter for each feature's min and max filter
    
    print("Applying hard filters...")
    
    max_filters = [
        di.STD_MAX, 
        di.ANOMALOUS_EXPRESSION_MAX, 
        di.ENHANCER_COUNT_MAX, 
        di.ENHANCER_PROPORTION_MAX, 
        di.CELL_LINE_EXPRESSION_MAX, 
        di.GENE_SIZE_MAX,
        di.SYMMETRY_MAX
    ]
    
    min_filters = [
        di.STD_MIN, 
        di.ANOMALOUS_EXPRESSION_MIN, 
        di.ENHANCER_COUNT_MIN, 
        di.ENHANCER_PROPORTION_MIN, 
        di.CELL_LINE_EXPRESSION_MIN, 
        di.GENE_SIZE_MIN,
        di.SYMMETRY_MIN
    ]
    
    for feature in interesting_features:
        gene_data = apply_hard_filter(
            gene_data,
            max_filters[interesting_features.index(feature)], 
                feature, 
                "max"
        )
        
        gene_data = apply_hard_filter(
            gene_data,
            min_filters[interesting_features.index(feature)],
            feature, "min"
        )
    
    return gene_data

def apply_hard_filter(gene_data, filter, feature, minmax):

    #Drops data above max filter or below min filter

    if minmax == "max":
        if filter is not False: gene_data = gene_data.drop(
                gene_data[gene_data[feature] > filter].index
            )
        
    elif minmax == "min":
        if filter is not False: gene_data = gene_data.drop(
                gene_data[gene_data[feature] < filter].index
            )
    
    else: 
        print("ERROR : Could not identify minmax.")
    
    return gene_data