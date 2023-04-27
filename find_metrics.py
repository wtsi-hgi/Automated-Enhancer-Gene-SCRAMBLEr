import numpy as np
import pandas as pd
import pyranges as pr
from sklearn.preprocessing import StandardScaler
import hashlib

import data_initialisation as di
import data_visualisation as dv

interesting_features = ["Std", "Anomalous_score", "Enhancer_count", "Enhancer_proportion", "Specific_gene_expression", "Gene_size"]

def find_gene_sizes(genes):
    
    #Adds the size of each gene based on its start and end point
    
    print("Finding gene sizes...")
    
    genes["Gene_size"] = genes["Gene_end"] - genes["Gene_start"]
    
    return genes

def find_interferring_genes(gene_data):
    
    #For each gene, finds the nearest gene upstream and downstream
    
    print("Finding interferring genes...")
    
    #"Start" and "End" are generated for use with PyRanges which requires the
    # columns to be named "Start" and "End".
    
    gene_data["Start"] = gene_data["Gene_start"]
    gene_data["End"] = gene_data["Gene_end"]
    
    interferring_genes_search = pr.PyRanges(
        gene_data.loc[gene_data["Specific_gene_expression"] > di.CELL_LINE_SPECIFIC_EXPRESSION_THRESHOLD])
    
    gene_search = pr.PyRanges(gene_data)

    nearest_gene_upstream = gene_search.nearest(
        interferring_genes_search, how = "upstream", 
        suffix = "_upstream_interferring_gene", 
        overlap = di.INTERFERRING_GENE_OVERLAPS)
    nearest_gene_downstream = gene_search.nearest(
        interferring_genes_search, how = "downstream", 
        suffix = "_downstream_interferring_gene", 
        overlap = di.INTERFERRING_GENE_OVERLAPS)
    
    gene_data = merge_nearest_genes_to_dataframe(
        gene_data, nearest_gene_upstream, nearest_gene_downstream)
    
    if not di.INTERFERRING_GENE_OVERLAPS:
        gene_data = gene_data.loc[
            gene_data["End_upstream_interferring_gene"] < gene_data["Gene_start"]]
        gene_data = gene_data.loc[
            gene_data["Start_downstream_interferring_gene"] > gene_data["Gene_end"]]
        
    gene_data.drop(["Start", "End"], axis = 1) 
    
    return gene_data

def merge_nearest_genes_to_dataframe(dataframe, upstream_genes, downstream_genes):
    nearest_gene_upstream = upstream_genes.df
    nearest_gene_downstream = downstream_genes.df

    dataframe = pd.merge(dataframe, 
        nearest_gene_upstream.loc[:, ["Gene_name", 
                                    "Start_upstream_interferring_gene", 
                                    "End_upstream_interferring_gene", 
                                    "Gene_name_upstream_interferring_gene"]], 
                                    on = "Gene_name", 
                                    how = "inner")
    dataframe = pd.merge(dataframe, 
        nearest_gene_downstream.loc[:, ["Gene_name", 
                                        "Start_downstream_interferring_gene", 
                                        "End_downstream_interferring_gene", 
                                        "Gene_name_downstream_interferring_gene"]], 
                                        on = "Gene_name", 
                                        how = "inner")
    
    return dataframe

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
        raise Exception("ERROR : Invalid search type.")
        
    genes["Search_window_start"] = np.where(
        genes["Strand"] == "+",
        genes["Gene_start"] - di.UPSTREAM_SEARCH,
        genes["Gene_start"] - di.DOWNSTREAM_SEARCH
    )
    
    genes["Search_window_end"] = np.where(
        genes["Strand"] == "+",
        genes[downstream_search_start] + di.DOWNSTREAM_SEARCH,
        genes[downstream_search_start] + di.UPSTREAM_SEARCH
    )
    
    genes["Search_window_start"] = np.where(
        genes["Gene_start"] < 0,
        0,
        genes["Search_window_start"]
    )
    
    genes["Search_window_start"] = np.where(
        genes["Search_window_start"] < genes["End_upstream_interferring_gene"],
        genes["End_upstream_interferring_gene"],
        genes["Search_window_start"]
    )
    
    genes["Search_window_end"] = np.where(
        genes["Search_window_end"] > genes["Start_downstream_interferring_gene"],
        genes["Start_downstream_interferring_gene"],
        genes["Search_window_end"]
    )    
        
    genes["Search_window_size"] = (genes["Search_window_end"] - genes["Search_window_start"])
                
    return genes

def find_mean(dataframe):
    
    #Adds the mean of gene expression to the given expression data frame, giving
    #a mean for each gene
    
    print("Finding mean of expression for each gene...")
    
    dataframe["Mean"] = dataframe.loc[:, dataframe.columns != "Gene_name"].mean(axis = 1)
    
    return dataframe
    
def find_std(dataframe):
    
    #Adds the standard deviation to the given expression dataframe, giving the
    #standard deviation across expression of each gene in all provided cell
    #types
    
    print("Finding standard deviation of expression for each gene...")
    
    dataframe["Std"] = dataframe.loc[:, dataframe.columns != "Gene_name"].std(axis = 1)
    
    return dataframe

def find_anomalous_score_of_gene_expression(dataframe):
    
    #Adds the z-score of each gene based on its expression in the cell line of interest compared to all others
    
    print("Finding anomalies...")
    
    dataframe["Anomalous_score"] = dataframe.apply(lambda gene : (gene["General_gene_expression"] - gene["Mean"]) / gene["Std"], axis = 1)

    return dataframe
    
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

    overlaps.drop(["Start", "End"], axis = 1)
    genes = pd.merge(genes, 
                     overlaps.groupby("Gene_name").size().reset_index(
        name = (element_type + "_count")), 
                     on = "Gene_name", 
                     how = "inner")
    
    return genes
    
def find_nearby_enhancer_densities(gene_data, overlaps):

    #Density of specifed element overlaps within the given search window is
    #calculated for each gene.
    
    print("Finding proximal enhancer densities...")
    
    overlaps["Enhancer_proportion"] = (overlaps.loc[:, "End"] - overlaps.loc[:, "Start"]) / overlaps.loc[:, "Search_window_size"]
    overlaps = overlaps.loc[:, ["Gene_name", "Enhancer_proportion"]].groupby(["Gene_name"], as_index = False)["Enhancer_proportion"].sum().reset_index()
    gene_data = pd.merge(gene_data, overlaps, on = "Gene_name")

    return gene_data
    
def calculate_interest_score(gene_data):
    
    #Various attributes of each gene are scaled and normallised, before being
    #weighted and combined into an interest score. The export_gene_scores_report
    #function is called
    
    print("Scoring genes...")
    
    scaler = StandardScaler()
    scaled_genes = gene_data.loc[:, (["Gene_name"] + interesting_features)]
    scaler.fit(scaled_genes.loc[:, interesting_features])
    scaled_genes.loc[:, interesting_features] = scaler.transform(scaled_genes[interesting_features])
    
    dv.compare_metrics(scaled_genes, "Comparison of Metrics within Z-space", "metrics_comparison")
    
    scaled_genes = scaled_genes.assign(
        Interest_score = (
            scaled_genes["Std"] * di.STD_WEIGHT +
            scaled_genes["Anomalous_score"] * di.ANOMALOUS_EXPRESSION_WEIGHT +
            scaled_genes["Enhancer_count"] * di.ENHANCER_COUNT_WEIGHT +
            scaled_genes["Enhancer_proportion"] * di.ENHANCER_PROPORTION_WEIGHT +
            scaled_genes["Specific_gene_expression"] * di.CELL_LINE_EXPRESSION_WEIGHT +
            pow((2 * di.GENE_SIZE_WEIGHT), (scaled_genes["Gene_size"] * di.GENE_SIZE_WEIGHT * di.GENE_SIZE_WEIGHT))
        )
    ).sort_values("Interest_score", ascending=False)
    
    gene_data = pd.merge(gene_data, scaled_genes.loc[:, ["Gene_name", "Interest_score"]], on = "Gene_name")
    gene_data = iterate_through_hard_filters(gene_data)
    gene_data = gene_data.sort_values("Interest_score", ascending = False).reset_index()
    
    export_gene_scores_report(gene_data)
    
    return gene_data

def iterate_through_hard_filters(gene_data):
    
    #Calls apply_hard_filter for each feature's min and max filter
    
    print("Applying hard filters...")
    
    max_filters = [di.STD_MAX, 
                   di.ANOMALOUS_EXPRESSION_MAX, 
                   di.ENHANCER_COUNT_MAX, 
                   di.ENHANCER_PROPORTION_MAX, 
                   di.CELL_LINE_EXPRESSION_MAX, 
                   di.GENE_SIZE_MAX]
    
    min_filters = [di.STD_MIN, 
                   di.ANOMALOUS_EXPRESSION_MIN, 
                   di.ENHANCER_COUNT_MIN, 
                   di.ENHANCER_PROPORTION_MIN, 
                   di.CELL_LINE_EXPRESSION_MIN, 
                   di.GENE_SIZE_MIN]
    
    for feature in interesting_features:
        
        gene_data = apply_hard_filter(gene_data, 
                                      max_filters[interesting_features.index(feature)], 
                                      feature, "max")
        
        gene_data = apply_hard_filter(gene_data, 
                                      min_filters[interesting_features.index(feature)], 
                                      feature, "min")
    
    return gene_data

def apply_hard_filter(gene_data, filter, feature, minmax):

    #Drops data above max filter or below min filter

    if minmax == "max":
        
        if filter is not False: gene_data = gene_data.drop(gene_data[gene_data[feature] > filter].index)
        
    elif minmax == "min":

        if filter is not False: gene_data = gene_data.drop(gene_data[gene_data[feature] < filter].index)
    
    else: 
        print("ERROR : Could not identify minmax.")
    
    return gene_data

def export_gene_scores_report(gene_data):
    
    #Idealy this will not read from file but from passed argument
    
    #Md5 checksum of config file is generated. Gene prioritisation report file
    #is created and checksum is included in name to differentiate different
    #configs. Report saved in given location.
    
    print("Exporting gene prioritisation report...")
    
    checksum = generate_config_checksum()
    
    with open("config.json", "r") as config:
        
        report_name = "gene_prioritisation_report_" + checksum.hexdigest() + ".txt"
        report = open((di.GENE_PRIORITISATION_REPORT_DIRECTORY + report_name), "w")
        report.write(config.read() + "\n")
        report.close()
        report = open((di.GENE_PRIORITISATION_REPORT_DIRECTORY + report_name), "a")
        gene_data.loc[:, (["Gene_name"] + ["Interest_score"] + interesting_features)].to_csv(
            (di.GENE_PRIORITISATION_REPORT_DIRECTORY + report_name), sep = "\t", index = True, mode = "a")
        report.close()
        
def generate_config_checksum():

    checksum = hashlib.md5()
    
    with open("config.json", "rb") as config:
        for chunk in iter(lambda: config.read(4096), b""):
            checksum.update(chunk)
        
    return checksum