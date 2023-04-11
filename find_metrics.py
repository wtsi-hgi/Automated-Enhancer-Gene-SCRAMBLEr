import pandas as pd
import pyranges as pr
from sklearn.preprocessing import StandardScaler
from sklearn.compose import ColumnTransformer
import hashlib
import sys

import data_initialisation as di
import region_convolutions as rc
import data_visualisation as dv
    
def find_mean(expression_data):
    
    #Adds the mean of gene expression to the given expression data frame, giving
    #a mean for each gene
    
    print("Finding mean of expression for each gene...")
    
    expression_data["Mean"] = expression_data.loc[:, expression_data.columns != "Gene_name"].mean(axis = 1)
    
    return expression_data
    
def find_std(expression_data):
    
    #Adds the standard deviation to the given expression dataframe, giving the
    #standard deviation across expression of each gene in all provided cell
    #types
    
    print("Finding standard deviation of expression for each gene...")
    
    expression_data["Std"] = expression_data.loc[:, expression_data.columns != "Gene_name"].std(axis = 1)
    
    return expression_data

def find_anomalous_score_of_gene_expression(expression_data):
    
    #Adds the z-score of each gene based on its expression in the cell line of interest compared to all others
    
    print("Finding anomalies...")
    
    expression_data["Anomalous_score"] = expression_data.apply(lambda gene : (gene["General_gene_expression"] - gene["Mean"]) / gene["Std"], axis = 1)

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
    #module, but are removed at the end of function
    gene_data["Start"] = gene_data["Gene_start"]
    gene_data["End"] = gene_data["Gene_end"]
    
    interferring_genes_pr = pr.PyRanges(gene_data.loc[gene_data["Specific_gene_expression"] > di.CELL_LINE_SPECIFIC_EXPRESSION_THRESHOLD])
    genes_pr = pr.PyRanges(gene_data)
    genes_nearest_upstream_pr = genes_pr.nearest(interferring_genes_pr, how = "upstream", suffix = "_upstream_interferring_gene", overlap = di.INTERFERRING_GENE_OVERLAPS)
    genes_nearest_downstream_pr = genes_pr.nearest(interferring_genes_pr, how = "downstream", suffix = "_downstream_interferring_gene", overlap = di.INTERFERRING_GENE_OVERLAPS)

    genes_nearest_upstream = genes_nearest_upstream_pr.df
    genes_nearest_downstream = genes_nearest_downstream_pr.df
    
    gene_data = pd.merge(gene_data, genes_nearest_upstream.loc[:, ["Gene_name", "Start_upstream_interferring_gene", "End_upstream_interferring_gene", "Gene_name_upstream_interferring_gene"]], on = "Gene_name", how = "inner")
    gene_data = pd.merge(gene_data, genes_nearest_downstream.loc[:, ["Gene_name", "Start_downstream_interferring_gene", "End_downstream_interferring_gene", "Gene_name_downstream_interferring_gene"]], on = "Gene_name", how = "inner")
    
    if di.INTERFERRING_GENE_OVERLAPS == False:
        
        gene_data = gene_data.loc[gene_data["End_upstream_interferring_gene"] < gene_data["Gene_start"]]
        gene_data = gene_data.loc[gene_data["Start_downstream_interferring_gene"] > gene_data["Gene_end"]]
        
    gene_data.drop(["Start", "End"], axis = 1)
        
    return gene_data

def find_search_windows(genes):

    #Defines a search window for each gene, based on the number of bases
    #upstream and downstream specified, the starting point of the window, and
    #whether the gene itself is included. Search window is foreshorterned if an
    #interferring gene is within the window, or if the window would include
    #negative bases.

    print("Finding sites to define search window...") 
    
    if (di.SEARCH_TYPE == "whole_gene" and di.SEARCH_WITHIN_GENE == True):
        
        genes["Search_window_start"] = genes.apply(lambda gene : gene["Gene_start"] - di.UPSTREAM_SEARCH if gene["Strand"] == "+" else gene["Gene_start"] - di.DOWNSTREAM_SEARCH, axis = 1)
        genes["Search_window_end"] = genes.apply(lambda gene : gene["Gene_end"] + di.DOWNSTREAM_SEARCH if gene["Strand"] == "+" else gene["Gene_end"] + di.UPSTREAM_SEARCH, axis = 1)
        genes["Search_window_start"] = genes.apply(lambda gene : 0 if gene["Gene_start"] < 0 else gene["Search_window_start"], axis = 1)
        genes["Search_window_start"] = genes.apply(lambda gene : gene["End_upstream_interferring_gene"] if gene["Search_window_start"] < gene["End_upstream_interferring_gene"] else gene["Search_window_start"], axis = 1)
        genes["Search_window_end"] = genes.apply(lambda gene : gene["Start_downstream_interferring_gene"] if gene["Search_window_end"] > gene["Start_downstream_interferring_gene"] else gene["Search_window_end"], axis = 1)
            
        genes["Search_window_size"] = (genes["Search_window_end"] - genes["Search_window_start"])
        
    elif (di.SEARCH_TYPE == "whole_gene" and di.SEARCH_WITHIN_GENE == False):

        search_upstreams = genes
        search_downstreams = genes
        search_upstreams["Upstream_search_window_end"] = search_upstreams["Gene_start"]
        search_upstreams["Upstream_search_window_start"] = search_upstreams.apply(lambda upstream : upstream["Upstream_search_window_end"] - di.UPSTREAM_SEARCH if upstream["Strand"] == "+" else upstream["Upstream_search_window_end"] - di.DOWNSTREAM_SEARCH, axis = 1)
        search_upstreams["Upstream_search_window_start"] = search_upstreams.apply(lambda upstream : 0 if upstream["Upstream_search_window_start"] < 0 else upstream["Upstream_search_window_start"], axis = 1)
        search_upstreams["Upstream_search_window_start"] = search_upstreams.apply(lambda upstream : upstream["End_upstream_interferring_gene"] if upstream["Upstream_search_window_start"] < upstream["End_upstream_interferring_gene"] else upstream["Upstream_search_window_start"], axis = 1)
        
        search_upstreams["Upstream_search_window_size"] = search_upstreams["Upstream_search_window_end"] - search_upstreams["Upstream_search_window_start"]
        
        search_downstreams["Downstream_search_window_start"] = search_downstreams["Gene_end"]
        search_downstreams["Downstream_search_window_end"] = search_downstreams.apply(lambda downstream : downstream["Downstream_search_window_start"] + di.DOWNSTREAM_SEARCH if downstream.Strand == "+" else downstream["Downstream_search_window_start"] + di.UPSTREAM_SEARCH, axis = 1)
        search_downstreams["Downstream_search_window_end"] = search_downstreams.apply(lambda downstream : downstream["Start_downstream_interferring_gene"] if downstream["Downstream_search_window_end"] > downstream["Start_downstream_interferrer"] else downstream["Downstream_search_window_end"], axis = 1)
        
        search_downstreams["Downstream_search_window_size"] = search_downstreams["Downstream_search_window_end"] - search_downstreams["Downstream_search_window_start"]
        
        genes = pd.merge(search_upstreams, search_downstreams)
        genes["Search_window_size"] = genes["Upstream_search_window_size"] + genes["Downstream_search_window_size"]
        
    elif (di.SEARCH_TYPE == "start_site" and di.SEARCH_WITHIN_GENE == True):
        
        genes["Search_window_start"] = genes.apply(lambda gene : gene["Gene_start"] - di.UPSTREAM_SEARCH if gene["Strand"] == "+" else gene["Gene_start"] - di.DOWNSTREAM_SEARCH, axis = 1)
        genes["Search_window_end"] = genes.apply(lambda gene : gene["Gene_start"] + di.DOWNSTREAM_SEARCH if gene["Strand"] == "+" else gene["Gene_start"] + di.UPSTREAM_SEARCH, axis = 1)
        genes["Search_window_start"] = genes.apply(lambda gene : 0 if gene.Start < 0 else gene.Start, axis = 1)
        genes["Search_window_start"] = genes.apply(lambda gene : gene["End_upstream_interferring_gene"] if gene["Search_window_start"] < gene["End_upstream_interferring_gene"] else gene["Search_window_start"], axis = 1)
        genes["Search_window_end"] = genes.apply(lambda gene : gene["Start_downstream_interferring_gene"] if gene["Search_window_end"] > gene["Start_downstream_interferring_gene"] else gene["Search_window_end"], axis = 1)
        
        genes["Search_window_size"] = (genes["Search_window_end"] - genes["Search_window_start"])
        
    elif (di.SEARCH_TYPE == "start_site" and di.SEARCH_WITHIN_GENE == False):
        
        print("Cannot currently perform this search type - search within gene combination")
        #genes_search["End"] = genes_search['Start']
        #genes_search["Start"] = genes_search.apply(lambda gene : gene.Start - upstream_search if gene.Strand == "+" else gene.Start - downstream_search, axis = 1)
        #genes_search["End"] = genes_search.apply(lambda gene : gene.End + downstream_search if gene.Strand == "+" else gene.End + upstream_search, axis = 1)
        #genes_search["Start"] = genes_search.apply(lambda gene : 0 if gene.Start < 0 else gene.Start, axis = 1)
        
    return genes
    
def find_element_overlaps_within_search_window(elements, genes):
    
    #PyRanges is used to find specified element type overlaps within search
    #window given for each gene. "Start" and "End" are generated for PyRanges
    #and removed subsequently. Overlaps are stored in a new dataframe
    
    print("Searching for overlapping elements...")
    
    genes["Start"] = genes["Search_window_start"]
    genes["End"] = genes["Search_window_end"]
    
    search_pr = pr.PyRanges(genes)
    elements_pr = pr.PyRanges(elements)
    overlaps = search_pr.intersect(elements_pr, strandedness = False)
    overlaps = overlaps.df
    
    genes.drop(["Start", "End"], axis = 1)
    
    return overlaps
    
def count_overlaps_per_gene(genes, overlaps, element_type):
    
    #Number of specified element overlaps are counted for each gene.
    
    print("Counting overlaps...")

    overlaps.drop(["Start", "End"], axis = 1)
    genes = pd.merge(genes, overlaps.groupby("Gene_name").size().reset_index(name = (element_type + "_count")), on = "Gene_name", how = "inner")
    
    return genes
    
def find_nearby_enhancer_densities(gene_data, overlaps):

    #Density of specifed element overlaps within the given search window is
    #calculated for each gene.
    
    print("Finding proximal enhancer densities...")
    
    overlaps["Enhancer_proportion"] = (overlaps.loc[:, "End"] - overlaps.loc[:, "Start"]) / overlaps.loc[:, "Search_window_size"]
    overlaps = overlaps.loc[:, ["Gene_name", "Enhancer_proportion"]].groupby(["Gene_name"], as_index = False)["Enhancer_proportion"].sum().reset_index()
    gene_data = pd.merge(gene_data, overlaps, on = "Gene_name")

    return gene_data
    
def calculate_interest_score(genes):
    
    #Various attributes of each gene are scaled and normallised, before being
    #weighted and combined into an interest score. The export_gene_scores_report
    #function is called
    
    print("Scoring genes...")
    
    scaler = StandardScaler()
    scaled_genes = genes.loc[:, ["Gene_name", "Std", "Anomalous_score", "Specific_gene_expression", "Enhancer_count", "Enhancer_proportion", "Gene_size"]]
    scaler.fit(scaled_genes.loc[:, ["Std", "Anomalous_score", "Specific_gene_expression", "Enhancer_count", "Enhancer_proportion", "Gene_size"]])
    scaled_genes.loc[:, ["Std", "Anomalous_score", "Specific_gene_expression", "Enhancer_count", "Enhancer_proportion", "Gene_size"]] = scaler.transform(scaled_genes[["Std", "Anomalous_score", "Specific_gene_expression", "Enhancer_count", "Enhancer_proportion", "Gene_size"]])
    
    dv.compare_metrics(scaled_genes, "Comparison of Metrics within Z-space", "metrics_comparison_zspace_keep_inf_log2")
    
    scaled_genes["Interest_score"] = 0
    scaled_genes.loc[:, "Interest_score"] = scaled_genes.loc[:, "Interest_score"] + (scaled_genes.loc[:, "Std"] * di.STD_WEIGHT)
    scaled_genes.loc[:, "Interest_score"] = scaled_genes.loc[:, "Interest_score"] + (scaled_genes.loc[:, "Anomalous_score"] * di.ANOMALOUS_EXPRESSION_WEIGHT)
    scaled_genes.loc[:, "Interest_score"] = scaled_genes.loc[:, "Interest_score"] + (scaled_genes.loc[:, "Enhancer_count"] * di.ENHANCER_COUNT_WEIGHT)
    scaled_genes.loc[:, "Interest_score"] = scaled_genes.loc[:, "Interest_score"] + (scaled_genes.loc[:, "Enhancer_proportion"] * di.ENHANCER_PROPORTION_WEIGHT)
    scaled_genes.loc[:, "Interest_score"] = scaled_genes.loc[:, "Interest_score"] + (scaled_genes.loc[:, "Specific_gene_expression"] * di.CELL_LINE_EXPRESSION_WEIGHT)
    scaled_genes.loc[:, "Interest_score"] = scaled_genes.loc[:, "Interest_score"] + (di.GENE_SIZE_WEIGHT / (scaled_genes.loc[:, "Gene_size"] + di.GENE_SIZE_WEIGHT))  
    scaled_genes = scaled_genes.sort_values("Interest_score", ascending = False)
    
    genes = pd.merge(genes, scaled_genes.loc[:, ["Gene_name", "Interest_score"]], on = "Gene_name")
    genes = genes.sort_values("Interest_score", ascending = False).reset_index()
    
    
    genes.loc[:, ["Gene_name", "Std", "Anomalous_score", "Specific_gene_expression", "Enhancer_count", "Enhancer_proportion", "Gene_size", "Interest_score"]].to_csv(di.RESULTS_DIRECTORY + "gene_scores.tsv", sep = "\t", index = True)
    
    export_gene_scores_report()
    
    return genes
    
def export_gene_scores_report():
    
    #Idealy this will not read from file but from passed argument
    
    #Md5 checksum of config file is generated. Gene prioritisation report file
    #is created and checksum is included in name to differentiate different
    #configs. Report saved in given location.
    
    print("Exporting gene prioritisation report...")
    
    hash_md5 = hashlib.md5()
    
    with open(sys.argv[1], "rb") as config:
        
        for chunk in iter(lambda: config.read(4096), b""):
            
            hash_md5.update(chunk)
            
    with open(sys.argv[1], "r") as config:
        
        report_name = "gene_prioritisation_report_" + hash_md5.hexdigest() + ".txt"
        report = open((di.GENE_PRIORITISATION_REPORT_DIRECTORY + report_name), "w")
        report.write(config.read() + "\n")
        
        with open((di.RESULTS_DIRECTORY + "gene_scores.tsv"), "r") as scores:
            
            report.write(scores.read())
            
        report.close()