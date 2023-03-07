import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import re
import pyranges
import sys
from scipy.ndimage import gaussian_filter1d
import math
from scipy.signal import find_peaks
from distutils.util import strtobool
import hashlib

def main():
    read_command_line()
    set_weights()
    read_genes(gene_annotation_reference)
    clean_genes()
    read_expression(cell_lines_expression_reference)
    clean_expression()
    read_regulatory_elements(regulatory_elements_reference)
    clean_regulatory_elements()
    find_gene_sizes()
    find_search_sites()
    find_search_site_sizes()
    merge_annotation_expression()
    filter_hap1_expression()
    find_nearby_genes()
    find_nearby_enhancers()
    count_nearby_enhancers()
    find_nearby_enhancer_densities()
    data_exploration()
    gene_scoring()
    export_gene_scores_report()
    enhancer_convolution()

def read_command_line():
    print("Reading command line...")
    global input_arguments
    
    global results_directory
    global gene_annotation_reference
    global regulatory_elements_reference
    global cell_lines_expression_reference
    global gene_prioritisation_report_directory
    
    global cell_line_of_interest
    global chromosomes_of_interest
    global epigenetic_flags_of_interest
    
    global search_type
    global search_within_gene
    global upstream_search
    global downstream_search
    
    global relative_non_housekeeping_weight
    global relative_enhancer_count_weight
    global relative_enhancer_proportion_weight
    global relative_cell_line_expression_weight
    global relative_gene_size_weight
    
    global kernel_size_type
    global absolute_kernel_size
    global relative_kernel_size
    global kernel_shape
    global relative_kernel_sigma
    
    global min_absolute_enhancer_cluster_width
    global min_enhancer_cluster_prominence
    
    global sigmoidal_slope
    global sigmoidal_midpoint
    global hap1_threshold
    global interferring_gene_overlaps
    
    input_arguments = sys.argv[1]
    with open(input_arguments) as input_file:
        results_directory = re.findall("(?<= = ).*", input_file.readline())[0]
        gene_annotation_reference = re.findall("(?<= = ).*", input_file.readline())[0]
        regulatory_elements_reference = re.findall("(?<= = ).*", input_file.readline())[0]
        cell_lines_expression_reference = re.findall("(?<= = ).*", input_file.readline())[0]
        gene_prioritisation_report_directory = re.findall("(?<= = ).*", input_file.readline())[0]
    
        cell_line_of_interest = re.findall("(?<= = ).*", input_file.readline())[0]
        chromosomes_of_interest = (re.findall("(?<= = ).*", input_file.readline())[0]).split(sep = ",")
        epigenetic_flags_of_interest = re.findall("(?<= = ).*", input_file.readline())[0]
    
        search_type = re.findall("(?<= = ).*", input_file.readline())[0]
        search_within_gene = strtobool(re.findall("(?<= = ).*", input_file.readline())[0])
        upstream_search = int(re.findall("(?<= = ).*", input_file.readline())[0])
        downstream_search = int(re.findall("(?<= = ).*", input_file.readline())[0])
    
        relative_non_housekeeping_weight = float(re.findall("(?<= = ).*", input_file.readline())[0])
        relative_enhancer_count_weight = float(re.findall("(?<= = ).*", input_file.readline())[0])
        relative_enhancer_proportion_weight = float(re.findall("(?<= = ).*", input_file.readline())[0])
        relative_cell_line_expression_weight = float(re.findall("(?<= = ).*", input_file.readline())[0])
        relative_gene_size_weight = float(re.findall("(?<= = ).*", input_file.readline())[0])
        
        kernel_size_type = re.findall("(?<= = ).*", input_file.readline())[0]
        absolute_kernel_size = float(re.findall("(?<= = ).*", input_file.readline())[0])
        relative_kernel_size = float(re.findall("(?<= = ).*", input_file.readline())[0])
        kernel_shape = re.findall("(?<= = ).*", input_file.readline())[0]
        relative_kernel_sigma = float(re.findall("(?<= = ).*", input_file.readline())[0])
        
        min_absolute_enhancer_cluster_width = float(re.findall("(?<= = ).*", input_file.readline())[0])
        min_enhancer_cluster_prominence = float(re.findall("(?<= = ).*", input_file.readline())[0])
        
        sigmoidal_slope = float(re.findall("(?<= = ).*", input_file.readline())[0])
        sigmoidal_midpoint = float(re.findall("(?<= = ).*", input_file.readline())[0])
        hap1_threshold = float(re.findall("(?<= = ).*", input_file.readline())[0])
        interferring_gene_overlaps = strtobool(re.findall("(?<= = ).*", input_file.readline())[0])
    
def set_weights():
    print("Setting weights...")
    global non_housekeeping_weight
    global enhancer_count_weight
    global enhancer_proportion_weight
    global cell_line_expression_weight
    global gene_size_weight
    
    non_housekeeping_weight = (relative_non_housekeeping_weight / 4)
    enhancer_count_weight = (relative_enhancer_count_weight / 250)
    enhancer_proportion_weight = (relative_enhancer_proportion_weight / 0.2)
    cell_line_expression_weight = (relative_cell_line_expression_weight / 14)
    gene_size_weight = (relative_gene_size_weight / 1)

def read_genes(genes_file):
    print("Reading gene annotations file...")
    global genes
    genes = pd.read_csv(genes_file, sep = "\t", names = ["Chromosome", "Source", "Type", "Start", "End", "Score", "Strand", "Phase", "Attributes"], skiprows = 5, dtype = {"Chromosome" : str, "Source" : str, "Type" : str, "Start" : int, "End" : int, "Score" : str, "Strand" : str, "Phase" : str, "Attributes" : str})
    
def read_expression(expression_file):
    print("Reading gene expression file...")
    global expression
    expression = pd.read_csv(expression_file).transpose()

def read_regulatory_elements(regulatory_file):
    print("Reading regulatory elements file...")
    global regulatory_elements
    if (regulatory_elements_reference[-3:] == "gff"):
        regulatory_elements = pd.read_csv(regulatory_file, sep = "\t", names = ["Chromosome", "Source", "Type", "Start", "End", "Score", "Strand", "Phase", "Attributes"])
    elif (regulatory_elements_reference[-3:] == "bed"):
        regulatory_elements = pd.read_csv(regulatory_file, sep = "\t", names = ["Chromosome", "Start", "End", "Flag"])
    else:
        print("Could not read regulatory elements file: " + regulatory_elements_reference)

def clean_genes():
    print("Cleaning gene data...")
    global genes
    genes = genes[genes["Chromosome"].isin(chromosomes_of_interest)]
    genes = genes.drop(genes[genes["Type"] != "gene"].index)
    genes["Gene_biotype"] = genes["Attributes"].apply(lambda x : re.findall("gene_biotype \"(.*?)\"", x)[0] if re.search("gene_name \"(.*?)\"", x) != None else "None")
    genes = genes.drop(genes[genes["Gene_biotype"] != "protein_coding"].index)
    genes["Gene_name"] = genes["Attributes"].apply(lambda x : re.findall("gene_name \"(.*?)\"", x)[0] if re.search("gene_name \"(.*?)\"", x) != None else "None")
    genes = genes.drop(["Source", "Type", "Score", "Phase", "Attributes", "Gene_biotype"], axis = 1)
    genes = genes.drop_duplicates(keep = False, subset = ["Gene_name"])

def clean_expression():
    print("Cleaning expression data...")
    global expression
    expression = expression.drop(["depmapID", "primary_disease"], axis = 0)
    expression.columns = expression.iloc[0]
    expression = expression.iloc[1:]
    expression = expression.reset_index().rename(columns = {"index" : "Gene_name"})
    expression["Mean"] = expression.loc[:, expression.columns != "Gene_name"].mean(axis = 1)
    expression["Std"] = expression.loc[:, expression.columns != "Gene_name"].std(axis = 1)
    expression = expression[["Gene_name", cell_line_of_interest, "Mean", "Std"]]
    expression = expression.drop_duplicates(keep = False, subset = ["Gene_name"])

def clean_regulatory_elements():
    print("Cleaning regulatory elements data...")
    global regulatory_elements
    if (regulatory_elements_reference[-3:] == "gff"):
        regulatory_elements = regulatory_elements.drop(regulatory_elements[regulatory_elements["Type"] != "enhancer"].index)
        #regulatory_elements["Enhancer_ID"] = regulatory_elements["Attributes"].apply(lambda x : re.findall("ID=enhancer:(.*?);", x)[0])
        regulatory_elements = regulatory_elements.drop(["Source", "Type", "Score", "Strand", "Phase", "Attributes"], axis = 1)
    elif (regulatory_elements_reference[-3:] == "bed"):
        regulatory_elements = regulatory_elements.drop(regulatory_elements[regulatory_elements["Flag"] != epigenetic_flags_of_interest].index)
        regulatory_elements = regulatory_elements.drop(["Flag"], axis = 1)
        regulatory_elements["Chromosome"] = regulatory_elements["Chromosome"].apply(lambda x : x[3:])
    else:
        print("Could not clean regulatory data.")

def merge_annotation_expression():
    print("Merging annotations and expression data...")
    global genes, expression
    #print(genes.Gene_name.duplicated().sum())
    #print(genes.loc[genes.duplicated(subset = ["Gene_name"]), :])
    #print(expression.Gene_name.duplicated().sum())
    genes = pd.merge(genes, expression, on = "Gene_name", how = "inner")
    #print(genes.Gene_name.duplicated().sum())
    #print(genes.loc[genes.duplicated(subset = ["Gene_name"], keep = False), :])

def find_gene_sizes():
    print("Finding gene sizes...")
    global genes
    genes["Gene_size"] = genes["End"] - genes["Start"]

def find_search_sites():
    global genes, genes_search
    genes_search = genes.loc[:,["Gene_name", "Chromosome", "Start", "End", "Strand"]]
    if (search_type == "whole_gene" and search_within_gene == "True"):
        print("Finding search sites from each end of the each gene, search site will contain the gene...")
        genes_search["Start"] = genes_search.apply(lambda gene : gene.Start - upstream_search if gene.Strand == "+" else gene.Start - downstream_search, axis = 1)
        genes_search["End"] = genes_search.apply(lambda gene : gene.End + downstream_search if gene.Strand == "+" else gene.End + upstream_search, axis = 1)
        genes_search["Start"] = genes_search.apply(lambda gene : 0 if gene.Start < 0 else gene.Start, axis = 1)
    elif (search_type == "whole_gene" and search_within_gene == "False"):
        print("Finding search sites from each end of the each gene, search site won't contain the gene...")
        search_upstreams = genes_search.loc[:, ["Gene_name", "Chromosome", "Start", "End", "Strand"]]
        search_downstreams = genes_search.loc[:, ["Gene_name", "Chromosome", "Start", "End", "Strand"]]
        search_upstreams["End"] = search_upstreams["Start"]
        search_upstreams["Start"] = search_upstreams.apply(lambda upstream : upstream.Start - upstream_search if upstream.Strand == "+" else upstream.Start - downstream_search, axis = 1)
        search_downstreams["Start"] = search_downstreams["End"]
        search_downstreams["End"] = search_downstreams.apply(lambda downstream : downstream.End + downstream_search if downstream.Strand == "+" else downstream.End + upstream_search, axis = 1)
        genes_search = pd.concat([search_upstreams, search_downstreams])
        genes_search["Start"] = genes_search.apply(lambda gene : 0 if gene.Start < 0 else gene.Start, axis = 1)
    elif (search_type == "start_site" and search_within_gene == "True"):
        print("Finding search sites from the start of each gene, search site will contain the gene...")
        genes_search["End"] = genes_search['Start']
        genes_search["Start"] = genes_search.apply(lambda gene : gene.Start - upstream_search if gene.Strand == "+" else gene.Start - downstream_search, axis = 1)
        genes_search["End"] = genes_search.apply(lambda gene : gene.End + downstream_search if gene.Strand == "+" else gene.End + upstream_search, axis = 1)
        genes_search["Start"] = genes_search.apply(lambda gene : 0 if gene.Start < 0 else gene.Start, axis = 1)
    elif (search_type == "start_site" and search_within_gene == "False"):
        print("Cannot currently perform this search type - search within gene combination")
        #genes_search["End"] = genes_search['Start']
        #genes_search["Start"] = genes_search.apply(lambda gene : gene.Start - upstream_search if gene.Strand == "+" else gene.Start - downstream_search, axis = 1)
        #genes_search["End"] = genes_search.apply(lambda gene : gene.End + downstream_search if gene.Strand == "+" else gene.End + upstream_search, axis = 1)
        #genes_search["Start"] = genes_search.apply(lambda gene : 0 if gene.Start < 0 else gene.Start, axis = 1)

def find_search_site_sizes():
    print("Finding search site sizes...")
    global genes_search
    genes_search["Search_size"] = (genes_search["End"] - genes_search["Start"]) * 2 #STOP-GAP SOLUTION
    #print("debug1")
    #print(genes_search)
    #print("debug2")
    #print(genes_search.groupby("Gene_name"))
    #print("debug3")
    #print(genes_search.groupby("Gene_name")["Search_size"])
    #print("debug4")
    #print(genes_search.groupby("Gene_name")["Search_size"].transform("sum"))
    #print("debug5")
    #genes_search = genes_search.set_index("Gene_name")
    #print(genes_search.groupby("Gene_name")["Search_size"].transform("sum").sort_values("Gene_name"))
    #print(genes.merge(genes_search.groupby("Gene_name")["Search_size"].transform("sum").sort_values("Gene_name"), on = "Gene_name"))
    #print(genes_search.assign(Search_size = genes_search.groupby()["Search_size"].transform("sum").sort_values("Gene_name")))

def find_nearby_enhancers():
    print("Searching for proximal enhancers...")
    global regulatory_elements, overlaps, genes_search, genes
    search_pr = pyranges.PyRanges(genes_search)
    regulatory_elements_pr = pyranges.PyRanges(regulatory_elements)
    overlaps = search_pr.intersect(regulatory_elements_pr, strandedness = False)
    overlaps = overlaps.df
    
def count_nearby_enhancers():
    print("Counting proximal enhancers...")
    global overlaps, genes
    genes = pd.merge(genes, overlaps.groupby("Gene_name").size().reset_index(name = "Enhancer_count"), on = "Gene_name", how = "inner")
    
def find_nearby_enhancer_densities():
    global genes
    print("Finding proximal enhancer densities...")
    overlaps["Enhancer_proportion"] = (overlaps["End"] - overlaps["Start"]) / overlaps["Search_size"]
    print(overlaps)
    genes = pd.merge(genes, overlaps.groupby(["Gene_name"], as_index = False)["Enhancer_proportion"].sum().reset_index(), on = "Gene_name", how = "inner")

def data_exploration():
    print("Visualising data...")
    global relative_expression, overlaps, genes
    relative_expression = expression.sort_values(by = ["Std"], ascending = False)
    #relative_expression = relative_expression[relative_expression["mean"] < relative_expression["HAP1"]]
    relative_expression["Difference"] = relative_expression[cell_line_of_interest] - relative_expression["Mean"]

    figure, axis = plt.subplots(1, 3, figsize = (18.5, 10.5))
    axis[0].scatter(relative_expression["Std"], relative_expression[cell_line_of_interest], s = 0.3, c = "red")
    axis[0].scatter(relative_expression["Std"], relative_expression["Mean"], s = 0.3, c = "blue")
    axis[0].set(xlabel = "Standard deviation of expression within cell lines", ylabel = "Normalised expression")

    axis[1].scatter(relative_expression["Difference"], relative_expression[cell_line_of_interest], s = 0.3, c = "red")
    axis[1].scatter(relative_expression["Difference"], relative_expression["Mean"], s = 0.3, c = "blue")
    axis[1].set(xlabel = "Difference between expression within " + cell_line_of_interest + " and mean expression", ylabel = "Normalised expression")

    axis[2].scatter(relative_expression["Std"], relative_expression["Difference"], s = 0.3, c = "green")
    axis[2].set(xlabel = "Difference between expression within " + cell_line_of_interest + " and mean expression", ylabel = "Normalised expression")

    plt.savefig(results_directory + "Significant_" + cell_line_of_interest + "_expression.png")
    plt.close(figure)
    
    figure, axis = plt.subplots(2, 2, figsize = (18.5, 10.5))
    axis[0][0].scatter(genes["Enhancer_count"], genes["Std"], s = 0.3, c = "purple")
    axis[0][0].set(xlabel = "Number of enhancers within search area", ylabel = "Standard deviation of expression within cell lines")
    
    axis[1][0].scatter(genes["Enhancer_count"], genes[cell_line_of_interest], s = 0.3, c = "red")
    axis[1][0].scatter(genes["Enhancer_count"] + 0.5, genes["Mean"], s = 0.3, c = "blue")
    axis[1][0].set(xlabel = "Number of enhancers within search area", ylabel = "Normalised expression")
    
    axis[0][1].scatter(genes["Enhancer_proportion"], genes["Std"], s = 0.3, c = "purple")
    axis[0][1].set(xlabel = "Proportion of search area which is enhancer", ylabel = "Standard deviation of expression within cell lines")
    
    axis[1][1].scatter(genes["Enhancer_proportion"], genes[cell_line_of_interest], s = 0.3, c = "red")
    axis[1][1].scatter(genes["Enhancer_proportion"], genes["Mean"], s = 0.3, c = "blue")
    axis[1][1].set(xlabel = "Proportion of search area which is enhancer", ylabel = "Normalised expression")
    
    plt.savefig(results_directory + "Enhancers_near_genes.png")
    plt.close(figure)

def gene_scoring():
    print("Scoring genes...")
    global genes
    genes["Interest_score"] = (genes.loc[:, "Std"] * non_housekeeping_weight) + \
        (genes.loc[:, "Enhancer_count"] * enhancer_count_weight) + \
        (genes.loc[:, "Enhancer_proportion"] * enhancer_proportion_weight) + \
        (genes.loc[:, "HAP1"] * cell_line_expression_weight) + \
        (gene_size_weight / (genes.loc[:, "Gene_size"] + gene_size_weight))
    genes = genes.sort_values("Interest_score", ascending = False)
    genes[["Gene_name", "Std", "HAP1", "Enhancer_count", "Enhancer_proportion", "Gene_size", "Interest_score"]].to_csv(results_directory + "gene_scores.tsv", sep = "\t", index = False)
    
def export_gene_scores_report():
    print("Exporting gene prioritisation report...")
    hash_md5 = hashlib.md5()
    with open(input_arguments, "rb") as config:
        for chunk in iter(lambda: config.read(4096), b""):
            hash_md5.update(chunk)
            
    with open(input_arguments, "r") as config:
        report_name = "gene_prioritisation_report_" + hash_md5.hexdigest() + ".txt"
        report = open((gene_prioritisation_report_directory + report_name), "w")
        report.write(config.read() + "\n")
        with open((results_directory + "gene_scores.tsv"), "r") as scores:
            report.write(scores.read())
        report.close()
    
def filter_hap1_expression():
    print("Filtering genes by HAP1 expression...")
    global genes
    genes["HAP1_normalised"] =  genes.apply(lambda gene : 1 - (1 / (1 + math.exp((sigmoidal_slope * gene.HAP1) - (sigmoidal_midpoint * sigmoidal_slope)))), axis = 1)
    
def find_nearby_genes():
    print("Finding nearby genes...")
    global genes_search, genes
    search_pr = pyranges.PyRanges(genes_search)
    genes_pr = pyranges.PyRanges(genes.loc[genes["HAP1_normalised"] > hap1_threshold])
    genes_nearest_upstream_pr = genes_pr.nearest(genes_pr, how = "upstream", suffix = "_upstream_interferrer", overlap = interferring_gene_overlaps)
    genes_nearest_downstream_pr = genes_pr.nearest(genes_pr, how = "downstream", suffix = "_downstream_interferrer", overlap = interferring_gene_overlaps)
    genes_nearest_upstream = genes_nearest_upstream_pr.df
    genes_nearest_downstream = genes_nearest_downstream_pr.df
    genes_search = pd.merge(genes_search, genes_nearest_upstream.loc[:, ["Gene_name", "Start_upstream_interferrer", "End_upstream_interferrer", "Gene_name_upstream_interferrer"]], on = "Gene_name", how = "inner")
    genes_search = pd.merge(genes_search, genes_nearest_downstream.loc[:, ["Gene_name", "Start_downstream_interferrer", "End_downstream_interferrer", "Gene_name_downstream_interferrer"]], on = "Gene_name", how = "inner")
    print(genes_search[["Gene_name", "Start", "End"]])
    genes_search = genes_search.loc[genes_search["End_upstream_interferrer"] < genes_search["Start_downstream_interferrer"]]
    print(genes_search[["Gene_name", "Start", "End"]])
    genes_search["Start"] = genes_search["End_upstream_interferrer"].astype(int)
    genes_search["End"] = genes_search["Start_downstream_interferrer"].astype(int)
    print(genes_search[["Gene_name", "Start", "End"]])

def enhancer_convolution():
    print("Convolving enhancers...")
    global genes_search, overlaps, gene_overlaps
    for index, gene in genes_search.iterrows():
        print("Convolving " + gene["Gene_name"])
        window = get_kernel(int((relative_kernel_size * (gene.End - gene.Start))), int(relative_kernel_sigma * (gene.End - gene.Start)))
        gene_basewise = np.zeros((gene.End - gene.Start), dtype = int)
        basewise = np.arange(gene.Start, gene.End)
        gene_specific_enhancer_overlaps = overlaps.loc[overlaps["Gene_name"] == gene["Gene_name"]]
        for index, overlap in gene_specific_enhancer_overlaps.iterrows():
            overlap_basewise = np.where(np.logical_and(overlap.Start <= basewise, basewise <= overlap.End), 1, 0)
            gene_basewise = np.where(overlap_basewise == 1, 1, gene_basewise)
        step_x = np.arange(gene.Start, gene.End)
        gene_convolution = np.convolve(window, gene_basewise)
        conv_x = np.arange((gene.Start - (len(window) // 2)), (gene.Start - (len(window) // 2) + len(gene_convolution)))
        
        peaks, _ = find_peaks(x = gene_convolution, width = min_absolute_enhancer_cluster_width, prominence = (min_enhancer_cluster_prominence, None))
        
        figure, axis = plt.subplots(2, 1, figsize = (18.5, 10.5))
        axis[0].plot(step_x, gene_basewise, c = "green")
        axis[0].set(xlabel = "Coordinate on chromosome " + gene["Chromosome"])
        axis[0].plot(conv_x, gene_convolution, c = "orange")
        axis[0].plot((peaks + (gene.Start - ((len(conv_x) - len(step_x)) / 2))), gene_convolution[peaks], "x")
        #axis[0].plot(conv_x, peaks, "x")
        #axis[0].set(xlabel = "Coordinate on chromosome " + gene["Chromosome"])
        
        plt.savefig(results_directory + gene["Gene_name"] + "_enhancer_convolution")
        plt.close()

def get_kernel(size, sigma):
    if kernel_shape == "flat":
        kernel = np.ones(size)
        return kernel
    elif kernel_shape == "guassian":
        kernel = np.zeros(size)
        np.put(kernel, (size // 2), 10)
        kernel = gaussian_filter1d(kernel, sigma)
        return kernel
    
if __name__ == "__main__":
    main()