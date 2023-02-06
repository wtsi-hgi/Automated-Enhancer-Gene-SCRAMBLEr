import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import ensembl_rest
import requests
import re
import pybedtools
import pyranges
import sys

server = "http://rest.ensembl.org"
overlap_region_ext = "/overlap/region/human/"
chromosomes_of_interest = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"]
non_housekeeping_weight = 1
absolute_enhancers_weight = 1
enhancer_density_weight = 1

def main():
    read_command_line()
    read_genes_pandas(gene_annotation_reference)
    #read_genes_bedtools()
    clean_genes_pandas()
    #clean_genes_pyranges()
    read_expression_pandas(cell_lines_expression_reference)
    clean_expression_pandas()
    read_regulatory_elements_pandas(regulatory_elements_reference)
    clean_regulatory_elements_pandas()
    merge_annotation_expression_pandas()
    #proximal_enhancer_search_api()
    proximal_enhancer_search_pandas()
    #print(genes)
    #print(regulatory_elements)
    data_exploration()

def read_command_line():
    global gene_annotation_reference
    global cell_lines_expression_reference
    global regulatory_elements_reference
    global upstream_search
    global downstream_search
    global cell_line_of_interest
    
    gene_annotation_reference = sys.argv[1]
    cell_lines_expression_reference = sys.argv[2]
    regulatory_elements_reference = sys.argv[3]
    upstream_search = int(sys.argv[4])
    downstream_search = int(sys.argv[5])
    cell_line_of_interest = sys.argv[6]

def read_genes_pandas(genes_file):
    global genes
    genes = pd.read_csv(genes_file, sep = "\t", names = ["Chromosome", "Source", "Type", "Start", "End", "Score", "Strand", "Phase", "Attributes"], skiprows = 5, dtype = {"Chromosome" : str, "Source" : str, "Type" : str, "Start" : int, "End" : int, "Score" : str, "Strand" : str, "Phase" : str, "Attributes" : str})
    
def read_genes_bedtools(genes_file):
    global genes
    genes = pybedtools.Bedtool(genes_file)
    
def read_expression_pandas(expression_file):
    global expression
    expression = pd.read_csv(expression_file).transpose()

def read_regulatory_elements_pandas(regulatory_file):
    global regulatory_elements
    regulatory_elements = pd.read_csv(regulatory_file, sep = "\t", names = ["Chromosome", "Source", "Type", "Start", "End", "Score", "Strand", "Phase", "Attributes"])

def clean_genes_pandas():
    global genes
    genes = genes[genes["Chromosome"].isin(chromosomes_of_interest)]
    genes = genes.drop(genes[genes["Type"] != "gene"].index)
    genes["Gene_biotype"] = genes["Attributes"].apply(lambda x : re.findall("gene_biotype \"(.*?)\"", x)[0] if re.search("gene_name \"(.*?)\"", x) != None else "None")
    genes = genes.drop(genes[genes["Gene_biotype"] != "protein_coding"].index)
    genes["Gene_name"] = genes["Attributes"].apply(lambda x : re.findall("gene_name \"(.*?)\"", x)[0] if re.search("gene_name \"(.*?)\"", x) != None else "None")
    genes = genes.drop(["Source", "Type", "Score", "Phase", "Attributes", "Gene_biotype"], axis = 1)
    
def clean_genes_pyranges():
    global genes
    global genes_pr
    genes = genes.dropna()
    genes["Start"] = genes["Start"].astype('int')
    genes["End"] = genes["End"].astype('int')
    genes_pr = pyranges.PyRanges(genes)

def clean_genes_bedtools():
    return None

def clean_expression_pandas():
    global expression
    expression = expression.drop(["depmapID", "primary_disease"], axis = 0)
    expression.columns = expression.iloc[0]
    expression = expression.iloc[1:]
    expression = expression.reset_index().rename(columns = {"index" : "Gene_name"})
    expression["Mean"] = expression.loc[:, expression.columns != "Gene_name"].mean(axis = 1)
    expression["Std"] = expression.loc[:, expression.columns != "Gene_name"].std(axis = 1)
    expression = expression[["Gene_name", cell_line_of_interest, "Mean", "Std"]]

def clean_regulatory_elements_pandas():
    global regulatory_elements
    regulatory_elements = regulatory_elements.drop(regulatory_elements[regulatory_elements["Type"] != "enhancer"].index)
    regulatory_elements["Enhancer_ID"] = regulatory_elements["Attributes"].apply(lambda x : re.findall("ID=enhancer:(.*?);", x)[0])
    regulatory_elements = regulatory_elements.drop(["Source", "Type", "Score", "Strand", "Phase", "Attributes"], axis = 1)

def merge_annotation_expression_pandas():
    global genes, expression
    genes = pd.merge(genes, expression, on = "Gene_name", how = "inner")

def proximal_enhancer_search_pandas():
    global genes
    global regulatory_elements
    global overlaps
    genes_search = genes
    genes_search["Start"] = genes_search["Start"].apply(lambda x : x - upstream_search)
    genes_search["End"] = genes_search["End"].apply(lambda x : x + downstream_search)
    genes_search["Search_size"] = genes_search["End"] - genes_search["Start"]
    gene_pr = pyranges.PyRanges(genes_search)
    regulatory_elements_pr = pyranges.PyRanges(regulatory_elements)
    overlaps = gene_pr.intersect(regulatory_elements_pr, strandedness = False)
    overlaps_count = overlaps.df.groupby("Gene_name").size().reset_index(name = "Enhancer_count")
    overlaps_size = overlaps.df
    overlaps_size["Enhancer_content"] = overlaps_size["End"] - overlaps_size["Start"]
    overlaps_size = overlaps_size.groupby(["Gene_name", "Search_size"], as_index = False)["Enhancer_content"].sum().reset_index()
    overlaps_size["Enhancer_proportion"] = (overlaps_size["Enhancer_content"] / overlaps_size["Search_size"])
    print("Overlapping regions")
    print(overlaps_size)
    overlaps_count = pd.merge(overlaps_count, overlaps_size, on = "Gene_name", how = "inner")
    genes = pd.merge(genes, overlaps_count, on = "Gene_name", how = "inner")

def proximal_enhancer_search_api():
    genes["search_start"] = genes["start"].apply(lambda x : x - 200000 if x > 200000 else 1)
    genes["search_end"] = genes["end"].apply(lambda x : x + 200000)
    genes["nearby_enhancer"] = 0
    for index, gene in genes.iterrows():
        for index, enhancer in regulatory_elements.iterrows():
            if enhancer["chromosome"] == gene["chromosome"] and (gene["search_start"] < enhancer["start"] < gene["search_end"] or gene["search_start"] < enhancer["end"] < gene["search_end"]):
                gene["nearby_enhancer"] += 1

def data_exploration():
    global relative_expression
    global overlaps
    global genes
    relative_expression = expression.sort_values(by = ["Std"], ascending = False)
    #relative_expression = relative_expression[relative_expression["mean"] < relative_expression["HAP1"]]
    relative_expression["Difference"] = relative_expression[cell_line_of_interest] - relative_expression["Mean"]
    #relative_expression["interest_metric"] = relative_expression["std"] * relative_expression["difference"]
    #relative_expression["interest_metric"] = relative_expression["std"]
    #relative_expression = relative_expression.sort_values(by = ["interest_metric"], ascending = False)

    figure, axis = plt.subplots(1, 3, figsize = (18.5, 10.5))
    axis[0].scatter(relative_expression["Std"], relative_expression[cell_line_of_interest], s = 0.3, c = "red")
    axis[0].scatter(relative_expression["Std"], relative_expression["Mean"], s = 0.3, c = "blue")
    axis[0].set(xlabel = "Standard deviation of expression within cell lines", ylabel = "Normalised expression")

    axis[1].scatter(relative_expression["Difference"], relative_expression[cell_line_of_interest], s = 0.3, c = "red")
    axis[1].scatter(relative_expression["Difference"], relative_expression["Mean"], s = 0.3, c = "blue")
    axis[1].set(xlabel = "Difference between expression within " + cell_line_of_interest + " and mean expression across cell lines", ylabel = "Normalised expression")

    axis[2].scatter(relative_expression["Std"], relative_expression["Difference"], s = 0.3, c = "green")
    axis[2].set(xlabel = "Difference between expression within " + cell_line_of_interest + " and mean expression across cell lines", ylabel = "Normalised expression")

    plt.savefig("Significant_" + cell_line_of_interest + "_expression.png")
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
    
    plt.savefig("Enhancers_near_genes.png")
    plt.close(figure)

def api_find_enhancers_near_genes():
    
    
    relative_expression = relative_expression.head(5)
    relative_expression["chromosome"] = relative_expression["gene_name"].apply(lambda x : ensembl_rest.symbol_lookup(species = "homo sapiens", symbol = x)["seq_region_name"])
    relative_expression["start_location"] = relative_expression["gene_name"].apply(lambda x : ensembl_rest.symbol_lookup(species = "homo sapiens", symbol = x)["start"])
    relative_expression["end_location"] = relative_expression["gene_name"].apply(lambda x : ensembl_rest.symbol_lookup(species = "homo sapiens", symbol = x)["end"])
    relative_expression["strand"] = relative_expression["gene_name"].apply(lambda x : ensembl_rest.symbol_lookup(species = "homo sapiens", symbol = x)["strand"])
    relative_expression["search_region"] =  relative_expression["chromosome"].astype(str) + ":" + (relative_expression["start_location"] - 200000).astype(str) + ".." + (relative_expression["end_location"] + 200000).astype(str) + ":" + relative_expression["strand"].astype(str)
    relative_expression["predicted_enhancers"] = relative_expression["search_region"].apply(lambda x : requests.get(server + overlap_region_ext + x + "?feature=regulatory").text.count("enhancer"))
    relative_expression["interest_metric"] = relative_expression["std"] * relative_expression["predicted_enhancers"]
    relative_expression = relative_expression.sort_values(by = ["interest_metric"], ascending = False)

    figure, axis = plt.subplots(1, 2, figsize = (18.5, 10.5))
    axis[0].scatter(relative_expression["std"], relative_expression["predicted_enhancers"], s = 0.5, c = "black")
    axis[0].set(xlabel = "Standard deviation of expression within cell lines", ylabel = "Number of enhancers within 200kb either side of gene")


    plt.savefig("Enhancers_near_genes.png")
    plt.close(figure)

    #print(CCLE_expression)
    #print(HAP1_expression)
    print(relative_expression)

def gene_scoring():
    global genes
    genes["Interest_score"] = genes["Std"] * non_housekeeping_weight + genes["Enhancer_count"] + absolute_enhancers_weight
    genes = genes.sort_values("Interest_score")

if __name__ == "__main__":
    main()