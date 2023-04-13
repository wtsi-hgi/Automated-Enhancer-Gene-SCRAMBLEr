import numpy as np
import pandas as pd
import pyranges as pr
from scipy.ndimage import gaussian_filter1d

import data_initialisation as di
import data_visualisation as dv
import sequence_seeking as ss

def define_step_function_of_element_overlaps_within_search_window(gene_data, overlaps, region_name):
    
    print("Generating step function for elements within search window...")
    
    gene_data["Search_window_start", "Search_window_end"] = gene_data["Search_window_start", "Search_window_end"].astype("int")

    gene_data[(region_name + "_step_function_x")] = [np.empty(0, dtype = float)] * len(gene_data)
    gene_data[(region_name + "_step_function_y")] = [np.empty(0, dtype = float)] * len(gene_data)

def convolution(genes_search, overlaps, region_name):
    
    #X and Y coordinates are generated for both step function and convolution
    #Empty columns must first be defined
    #For each gene, the get_kernel function is called, an array of zeroes is generated the size of the search window,
    #overlaps within the window are then iterated over, and sections of the array become ones at the relevant coordinates.
    #This step function is then convolved with the step function.
    #X coordinates for each are then retroactively found from lengths of each
    
    genes_search["Search_window_start"] = genes_search["Search_window_start"].astype("int")
    genes_search["Search_window_end"] = genes_search["Search_window_end"].astype("int")
    
    genes_search[(region_name + "_searched_coordinates")] = [np.empty(0, dtype = float)] * len(genes_search)
    genes_search[(region_name + "_step_function")] = [np.empty(0, dtype = float)] * len(genes_search)
    genes_search[(region_name + "_convolution")] = [np.empty(0, dtype = float)] * len(genes_search)
    genes_search[(region_name + "_convolved_coordinates")] = [np.empty(0, dtype = float)] * len(genes_search)

    genes_search = genes_search.sort_values("Interest_score", ascending = False).reset_index(drop = True)
    
    for index, gene in genes_search.head(di.ENHANCER_CONVOLUTION).iterrows():
        
        print("Convolving " + region_name.lower() + "s for " + gene["Gene_name"] + " (" + str(index + 1) + " of " + str(di.ENHANCER_CONVOLUTION) + ")...")
        
        kernel = get_kernel(di.ENHANCER_KERNEL_SHAPE, int((di.RELATIVE_ENHANCER_KERNEL_SIZE * (gene["Search_window_end"] - gene["Search_window_start"]))), int(di.RELATIVE_ENHANCER_KERNEL_SIGMA * (gene["Search_window_end"] - gene["Search_window_start"])))
        step_function_y = np.zeros((gene["Search_window_end"] - gene["Search_window_start"]), dtype = int)
        basewise = np.arange(gene["Search_window_start"], gene["Search_window_end"])
        gene_specific_overlaps = overlaps.loc[overlaps["Gene_name"] == gene["Gene_name"]]
        
        for overlap_index, overlap in gene_specific_overlaps.iterrows():
            
            overlap_basewise = np.where(np.logical_and(overlap["Start"] <= basewise, basewise <= overlap["End"]), 1, 0)
            step_function_y = np.where(overlap_basewise == 1, 1, step_function_y)
            
        step_function_x = ((np.arange(gene["Search_window_start"], gene["Search_window_end"])))
        convolution_y = np.convolve(kernel, step_function_y)
        convolution_x = (np.arange((gene["Search_window_start"] - (len(kernel) // 2)), (gene["Search_window_start"] - (len(kernel) // 2) + len(convolution_y))))
        
        genes_search.at[index, (region_name + "_searched_coordinates")]  = step_function_x
        genes_search.at[index, (region_name + "_step_function")] = step_function_y
        genes_search.at[index, (region_name + "_convolution")] = np.convolve(kernel, step_function_y)
        genes_search.at[index, (region_name + "_convolved_coordinates")] = convolution_x
        
        dv.plot_convolutions(gene, step_function_x, step_function_y, convolution_x, convolution_y)
        
    return genes_search

def get_kernel(kernel_shape, size, sigma):
    
    #Kernel is generated as numpy array depending on desired shape and size
    
    if kernel_shape == "flat":
        
        kernel = np.ones(size)
        
    
    elif kernel_shape == "guassian":
        
        kernel = np.zeros(size)
        np.put(kernel, (size // 2), 10)
        kernel = gaussian_filter1d(kernel, sigma)
        
    else:
        raise Exception("Kernel shape is neither Flat nor Guassian")
        
    return kernel

def combine_convolutions(enhancer_convolution, quiescent_convolution):
    
    #Eventually will be used to add convolutions together
    
    print("Merging convolutions...")

    enhancer_convolution = np.negative(enhancer_convolution)
    print(enhancer_convolution)
    combined_convolution = np.add((enhancer_convolution * di.ENHANCER_CONVOLUTION_WEIGHT), (quiescent_convolution * di.QUIESCENT_CONVOLUTION_WEIGHT))
    
    return combined_convolution
    
def export_convolutions(gene_data):
    
    #Coordinates of convolutions are exported to wig file, for each gene
    
    print("Exporting enhancer density convolutions to wig file...")
    
    for index, gene in gene_data.head(di.ENHANCER_CONVOLUTION).iterrows():
    
        with open((di.RESULTS_DIRECTORY + gene["Gene_name"] + "_convolutions.wiggle"), "w") as f:
            
            f.write("fixedStep chrom=chr" + gene["Chromosome"] + " start=" + str(gene["Enhancer_convolved_coordinates"][0]) +" step=1")
            f.write("\n")
    
        convolution_signal = pd.DataFrame({"Convolution_signal" : gene_data.loc[index, "Enhancer_convolution"]})
        convolution_signal.to_csv((di.RESULTS_DIRECTORY + gene["Gene_name"] + "_convolutions.wig"), sep = "\t", index = False, mode = "a", header = False)
    
def find_plateaus(gene_data):
    
    #find_plateaus takes convolved coordinates, and applies a threshold to
    #separate the search window into regions based on the y-value of each
    #convolved base.
    
    gene_data["Plateau_coordinates"] = ""
    gene_data["Plateau_starts"] = ""
    gene_data["Plateau_ends"] = ""
    
    gene_data = gene_data.sort_values("Interest_score", ascending = False)
    
    for index, gene in gene_data.head(di.ENHANCER_CONVOLUTION).iterrows():
        
        print("Finding plateaus for gene " + gene["Gene_name"] + " (" + str(index + 1) + " of " + str(di.ENHANCER_CONVOLUTION) + ")...")
        
        convolved_x = gene["Enhancer_convolved_coordinates"]
        convolved_y = gene["Enhancer_convolution"]
        convolved_y = np.append(convolved_y, 0)
        boolean_below_threshold = convolved_y < di.PLATEAU_THRESHOLD
        boolean_below_threshold = np.concatenate((boolean_below_threshold[:(gene["Gene_start"] - convolved_x[0])], np.full((gene["Gene_end"] -  gene["Gene_start"]), False), boolean_below_threshold[(gene["Gene_end"] - convolved_x[0]):]))
        boolean_below_threshold = (boolean_below_threshold[:-1] != boolean_below_threshold[1:])
        plateau_coordinates = convolved_x[boolean_below_threshold]
        plateau_coordinates = np.concatenate([[convolved_x[0]], plateau_coordinates, [convolved_x[-1]]])
        
        gene_data.at[index, "Plateau_coordinates"] = plateau_coordinates
        gene_data.at[index, "Plateau_starts"] = plateau_coordinates[::2]
        gene_data.at[index, "Plateau_ends"] = plateau_coordinates[1::2]
        
    return gene_data
    
def export_plateaus(gene_data):
    #export_plateaus saves plateaus associated with each gene as a bed file
    
    print("Exporting plateaus to bed file...")
    
    with open((di.RESULTS_DIRECTORY + "plateaus.bed"), "w") as f:
            f.write("Chromosome Start	End	Gene_name")
            f.write("\n")

    for index, gene in gene_data.head(di.ENHANCER_CONVOLUTION).iterrows():
        
        plateau_regions = pd.DataFrame({"Start" : gene_data.loc[index, "Plateau_starts"], "End" : gene_data.loc[index, "Plateau_ends"]})
        plateau_regions["Gene_name"] = gene["Gene_name"]
        plateau_regions["Chromosome"] = "chr" + gene["Chromosome"]
        plateau_regions["Strand"] = gene["Strand"]
        
        #plateau_regions = ss.find_fasta(plateau_regions)
        
        plateau_regions.to_csv((di.RESULTS_DIRECTORY + "plateaus.bed"), sep = "\t", index = False, columns = ["Chromosome", "Start", "End", "Gene_name"], mode = "a", header = False)