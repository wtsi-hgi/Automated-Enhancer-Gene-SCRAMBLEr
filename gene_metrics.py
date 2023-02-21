import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import re
import pyranges
import sys
import sklearn as sk
from sklearn.tree import DecisionTreeClassifier
from scipy.ndimage import gaussian_filter1d
import math



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
    proximal_enhancer_search()
    count_proximal_enhancers()
    find_proximal_enhancer_densities()
    #find_nearby_genes()
    #data_exploration()
    gene_scoring()
    convolution()

def read_command_line():
    print("Reading command line...")
    global results_directory
    global gene_annotation_reference
    global regulatory_elements_reference
    global cell_lines_expression_reference
    
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
    global relative_kernel_sigma
    
    global sigmoidal_slope
    global sigmoidal_midpoint
    
    if len(sys.argv) == 2:
        input_arguments = sys.argv[1]
        with open(input_arguments) as input_file:
            results_directory = re.findall("(?<= = ).*", input_file.readline())[0]
            gene_annotation_reference = re.findall("(?<= = ).*", input_file.readline())[0]
            regulatory_elements_reference = re.findall("(?<= = ).*", input_file.readline())[0]
            cell_lines_expression_reference = re.findall("(?<= = ).*", input_file.readline())[0]
        
            cell_line_of_interest = re.findall("(?<= = ).*", input_file.readline())[0]
            chromosomes_of_interest = re.findall("(?<= = ).*", input_file.readline())[0]
            epigenetic_flags_of_interest = re.findall("(?<= = ).*", input_file.readline())[0]
        
            search_type = re.findall("(?<= = ).*", input_file.readline())[0]
            search_within_gene = re.findall("(?<= = ).*", input_file.readline())[0]
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
            relative_kernel_sigma = float(re.findall("(?<= = ).*", input_file.readline())[0])
            
            sigmoidal_slope = float(re.findall("(?<= = ).*", input_file.readline())[0])
            sigmoidal_midpoint = float(re.findall("(?<= = ).*", input_file.readline())[0])
        
    else:
        results_directory = sys.argv[1]
        gene_annotation_reference = sys.argv[2]
        regulatory_elements_reference = sys.argv[3]
        cell_lines_expression_reference = sys.argv[4]
        
        cell_line_of_interest = sys.argv[5]
        epigenetic_flags_of_interest = sys.argv[6]
        
        search_type = sys.argv[7]
        search_within_gene = sys.argv[8]
        upstream_search = int(sys.argv[9])
        downstream_search = int(sys.argv[10])
        
        relative_non_housekeeping_weight = float(sys.argv[11])
        relative_enhancer_count_weight = float(sys.argv[12])
        relative_enhancer_proportion_weight = float(sys.argv[13])
        relative_cell_line_expression_weight = float(sys.argv[14])
        relative_gene_size_weight = float(sys.argv[15])
        relative_kernel_size = float(sys.argv[16])
    
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
        regulatory_elements = regulatory_elements.drop(regulatory_elements[regulatory_elements["Flag"] != epigenetic_flag_of_interest].index)
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

#def reduce_search_sites_due_to_interferring_genes(gene, genes, genes_search):
#    interferring_genes = genes.loc[genes["Chromosome"] ==  gene["Chromosome"]]
#    search_start = genes_
#    search_end = genes_
#    interferring_genes = interferring_genes.loc[genes_search["Startintereferring_genes["Start"]]
#    genes.loc[genes["End"] <=  genes_search["Gene_name"]]

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

def proximal_enhancer_search():
    print("Searching for proximal enhancers...")
    global regulatory_elements, overlaps, genes_search, genes
    search_pr = pyranges.PyRanges(genes_search)
    regulatory_elements_pr = pyranges.PyRanges(regulatory_elements)
    overlaps = search_pr.intersect(regulatory_elements_pr, strandedness = False)
    overlaps = overlaps.df
    
def count_proximal_enhancers():
    print("Counting proximal enhancers...")
    global overlaps, genes
    genes = pd.merge(genes, overlaps.groupby("Gene_name").size().reset_index(name = "Enhancer_count"), on = "Gene_name", how = "inner")
    
def find_proximal_enhancer_densities():
    global genes
    print("Finding proximal enhancer densities...")
    overlaps["Enhancer_proportion"] = (overlaps["End"] - overlaps["Start"]) / overlaps["Search_size"]
    genes = pd.merge(genes, overlaps.groupby(["Gene_name"], as_index = False)["Enhancer_proportion"].sum().reset_index(), on = "Gene_name", how = "inner")

def data_exploration():
    print("Visualising data...")
    global relative_expression
    global overlaps
    global genes
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
    genes[["Gene_name", "Std", "HAP1", "Enhancer_count", "Enhancer_proportion", "Gene_size", "Interest_score"]].to_csv(results_directory + "gene_scoring.tsv", sep = "\t")
    
def transform_hap1_expression():
    global genes
    genes["HAP1_score"] = 1 - (1 / (1 + math.exp((sigmoidal_slope * genes["HAP1"]) - (sigmoidal_midpoint * sigmoidal_slope))))
    
def find_nearby_genes():
    print("Finding nearby genes...")
    global genes_search, genes
    search_pr = pyranges.PyRanges(genes_search)
    genes_pr = pyranges.PyRanges(genes)
    gene_overlaps = search_pr.intersect(genes_pr, strandedness = False)
    gene_overlaps = gene_overlaps.df
    print(gene_overlaps)
    print(genes)

def convolution():
    print("Convolving...")
    global genes_search, overlaps
    for index, gene in genes_search.iterrows():
        window = get_guassian_kernel(int((relative_kernel_size * (gene.End - gene.Start))), int(relative_kernel_sigma * (gene.End - gene.Start)))
        gene_basewise = np.zeros((gene.End - gene.Start), dtype = int)
        basewise = np.arange(gene.Start, gene.End)
        gene_overlaps = overlaps.loc[overlaps["Gene_name"] == gene["Gene_name"]]
        for index, overlap in gene_overlaps.iterrows():
            overlap_basewise = np.where(np.logical_and(overlap.Start <= basewise, basewise <= overlap.End), 1, 0)
            gene_basewise = np.where(overlap_basewise == 1, 1, gene_basewise)
        step_x = np.arange(gene.Start, gene.End)
        print("Starting convolution")
        gene_convolution = np.convolve(window, gene_basewise)
        print("Convolution finished")
        conv_x = np.arange((gene.Start - (len(window) // 2)), (gene.Start - (len(window) // 2) + len(gene_convolution)))
        
        figure, axis = plt.subplots(2, 1, figsize = (18.5, 10.5))
        axis[0].plot(step_x, gene_basewise, c = "green")
        axis[0].set(xlabel = "Coordinate on chromosome " + gene["Chromosome"])
        axis[1].plot(conv_x, gene_convolution, c = "orange")
        axis[1].set(xlabel = "Coordinate on chromosome " + gene["Chromosome"])
        
        plt.savefig(results_directory + gene["Gene_name"] + "_enhancer_convolution")
        plt.close()

def get_guassian_kernel(size, sigma):
    kernel = np.zeros(size)
    print(size)
    print(np.count_nonzero(kernel))
    np.put(kernel, (size // 2), 1)
    print(np.count_nonzero(kernel))
    kernel = gaussian_filter1d(kernel, sigma)
    print(np.count_nonzero(kernel))
    return kernel
    
    #genes_search["Basewise_0s"] = genes_search.apply(lambda gene_search_region : np.zeros(gene_search_region.Search_size, dtype = int), axis = 1)
    #overlaps["Basewise_1s"] = overlaps.apply(lambda overlap : np.ones((overlap.End - overlap.Start), dtype = int), axis = 1)
    #print(genes_search)
    #print(overlaps)
    #overlaps = overlaps.merge(genes_search, how = "left")
    #print(overlaps)
    
#def convolution():
    #print("Convolving...")
    #global overlaps, genes
    #print(genes)
    #print(overlaps)
    #overlaps = overlaps.drop(["Strand"], axis = 1)
    #overlaps["Base_labels"] = np.nan
    #overlaps = overlaps.sort_values(["Chromosome", "Start"]).head(100)
    #print(overlaps)
    #genes = genes.sort_values(["Chromosome", "Start"]).head(10)
    #print(genes)
    #for gene in genes["Gene_name"]:
    #    convolve_gene(gene, genes, overlaps)
    
    #enhancer_convolution = np.convolve(window, enhancer_ranges)
    
#def convolve_gene(gene, genes, genes_search):
    #base_labels = np.zeros((2, gene["Search_size"]))
    #np.put(base_labels, [0, ])
    #print(overlaps.loc[overlaps["Gene_name"] == gene])
    #print(overlaps.loc[overlaps["Gene_name"] == gene, "Base_labels"])
    #overlap_coordinate_ranges = np.empty
    #for index, overlap in overlaps.loc[overlaps["Gene_name"] == gene].iterrows():
    #    overlap_coordinate_ranges = np.append(overlap_coordinate_ranges, np.arange(overlap.Start, (overlap.End + 1), 1))
    #overlap_coordinate_ranges[1:]
    #print("Debug1")
    #print(genes_search.loc[genes_search["Gene_name"] == gene])
    #print(genes_search.loc[genes_search["Gene_name"] == gene]["Start"])
    #overlaps.loc[overlaps["Gene_name"] == gene, "Base_labels"] = overlaps.loc[overlaps["Gene_name"] == gene].apply(lambda overlap : np.arange(overlap.Start, overlap.End, 1), axis = 1)
    #gene_coordinate_range = np.arange(genes_search.loc[genes_search["Gene_name"] == gene].Start, genes_search.loc[genes_search["Gene_name"] == gene].End, 1)
    #print(overlaps.loc[overlaps["Gene_name"] == gene]["Base_labels"])
    #print(gene_coordinate_range)
    #print(overlap_coordinate_ranges)
    #print(np.isin(gene_coordinate_range, overlap_coordinate_ranges))
    #print(overlaps.loc[overlaps["Gene_name"] == gene])
    
    
    #for overlap in overlaps.loc[overlaps["Gene_name"] == gene]:
    #    overlap
    #gene["Base_labels"] = gene.apply(lambda overlap : [*range(gene.Start, gene.End, 1)], axis = 1)
    #gene = gene.drop(["Start", "End"], axis = 1)
    #enhancer_convolution = np.convolve(window, gene["Base_labels"].to_numpy())

if __name__ == "__main__":
    main()