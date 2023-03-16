import numpy as np
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d

import data_initialisation as di

def enhancer_convolution(genes_search, overlaps):
    
    genes_search["Start"] = genes_search["Search_window_start"].astype("int")
    genes_search["End"] = genes_search["Search_window_end"].astype("int")
    
    genes_search["Enhancer_searched_coordinates"] = [np.empty(0, dtype = float)] * len(genes_search)
    genes_search["Enhancer_step_function"] = [np.empty(0, dtype = float)] * len(genes_search)
    genes_search["Enhancer_convolution"] = [np.empty(0, dtype = float)] * len(genes_search)
    genes_search["Enhancer_convolved_coordinates"] = [np.empty(0, dtype = float)] * len(genes_search)

    genes_search = genes_search.sort_values("Interest_score", ascending = False).reset_index(drop = True)
    
    for index, gene in genes_search.iterrows():
        
        print("Convolving enhancers for gene " + gene["Gene_name"] + " (" + str(index + 1) + " of " + str(di.ENHANCER_CONVOLUTION) + ")...")
        
        window = get_kernel(di.ENHANCER_KERNEL_SHAPE, int((di.RELATIVE_ENHANCER_KERNEL_SIZE * (gene.End - gene.Start))), int(di.RELATIVE_ENHANCER_KERNEL_SIGMA * (gene.End - gene.Start)))
        gene_basewise = np.zeros((gene.End - gene.Start), dtype = int)
        basewise = np.arange(gene.Start, gene.End)
        gene_specific_enhancer_overlaps = overlaps.loc[overlaps["Gene_name"] == gene["Gene_name"]]
        
        for overlap_index, overlap in gene_specific_enhancer_overlaps.iterrows():
            
            overlap_basewise = np.where(np.logical_and(overlap.Start <= basewise, basewise <= overlap.End), 1, 0)
            gene_basewise = np.where(overlap_basewise == 1, 1, gene_basewise)

        genes_search.loc[index, "Enhancer_searched_coordinates"] = ((np.arange(gene["Start"], gene["End"]))).tostring()
        genes_search.loc[index, "Enhancer_step_function"] = gene_basewise.tostring()
        convolution = np.convolve(window, gene_basewise)
        enhancer_convolution = convolution.tostring()
        genes_search.loc[index, "Enhancer_convolution"] = enhancer_convolution
        genes_search.loc[index, "Enhancer_convolved_coordinates"] = ((np.arange((gene["Start"] - (len(window) // 2)), (gene["Start"] - (len(window) // 2) + len(convolution))))).tostring()

        if index > di.ENHANCER_CONVOLUTION:
            break
        
    genes_search.drop(["Start", "End"], axis = 1)
    #print(genes_search)
        
    return genes_search

def quiescent_convolution(genes_search, overlaps):
    
    print("Convolving quiescent regions...")
    
    genes_search["Start"] = genes_search["Search_window_start"].astype("int")
    genes_search["End"] = genes_search["Search_window_end"].astype("int")
    
    genes_search["Quiescent_searched_coordinates"] = ""
    genes_search["Quiescent_convolution"] = ""
    genes_search["Quiescent_convolved_coordinates"] = ""

    genes_search = genes_search.sort_values("Interest_score", ascending = False)

    i = 1
    
    for index, gene in genes_search.iterrows():
        
        window = get_kernel(di.ENHANCER_KERNEL_SHAPE, int((di.RELATIVE_ENHANCER_KERNEL_SIZE * (gene.End - gene.Start))), int(di.RELATIVE_ENHANCER_KERNEL_SIGMA * (gene.End - gene.Start)))
        gene_basewise = np.zeros((gene.End - gene.Start), dtype = int)
        basewise = np.arange(gene.Start, gene.End)
        gene_specific_enhancer_overlaps = overlaps.loc[overlaps["Gene_name"] == gene["Gene_name"]]
        
        for overlap_index, overlap in gene_specific_enhancer_overlaps.iterrows():
            
            overlap_basewise = np.where(np.logical_and(overlap.Start <= basewise, basewise <= overlap.End), 1, 0)
            gene_basewise = np.where(overlap_basewise == 1, 1, gene_basewise)
            
        gene["Quiescent_searched_coordinates"] = np.arange(gene.Start, gene.End)
        gene["Quiescent_convolution"] = np.convolve(window, gene_basewise)
        gene["Quiescent_convolved_coordinates"] = np.arange((gene.Start - (len(window) // 2)), (gene.Start - (len(window) // 2) + len(gene["Quiescent_convolution"])))

        print("Convolved quiescent regions for gene " + str(index) + " of " + str(di.ENHANCER_CONVOLUTION) + "...")
        
        if i > di.ENHANCER_CONVOLUTION:
            break
        
    genes_search.drop(["Start", "End"], axis = 1)
        
    return genes_search
    
def combine_convolutions(enhancer_convolution, quiescent_convolution):
    
    print("Merging convolutions...")

    #enhancer_convolution = np.negative(enhancer_convolution)
    #print(enhancer_convolution)
    #recombination_convolution = np.add((enhancer_convolution * di.ENHANCER_CONVOLUTION_WEIGHT), (quiescent_convolution * di.QUIESCENT_CONVOLUTION_WEIGHT))
    
    return overall_convolution

def get_kernel(kernel_shape, size, sigma):
    
    if kernel_shape == "flat":
        
        kernel = np.ones(size)
        
        return kernel
    
    elif kernel_shape == "guassian":
        
        kernel = np.zeros(size)
        np.put(kernel, (size // 2), 10)
        kernel = gaussian_filter1d(kernel, sigma)
        
        return kernel
    
def find_plateaus(gene_data):
    
    gene_data["Plateau_starts"] = ""
    gene_data["Plateau_ends"] = ""
    
    gene_data = gene_data.sort_values("Interest_score", ascending = False)
    
    for index, gene in gene_data.iterrows():
        
        print("Finding plateaus for gene " + gene["Gene_name"] + " (" + str(index + 1) + " of " + str(di.ENHANCER_CONVOLUTION) + ")...")
        
        #print(gene)
        #print(gene["Enhancer_convolved_coordinates"])
        #convolved_x = np.fromstring(gene["Enhancer_convolved_coordinates"])
        #convolved_y = np.fromstring(gene["Enhancer_convolution"])
        convolved_x = np.frombuffer(gene["Enhancer_convolved_coordinates"], dtype = float)
        convolved_y = np.frombuffer(gene["Enhancer_convolution"], dtype = float)
        threshold = di.PLATEAU_THRESHOLD
        
        #convolved_y = np.negative((convolved_y))
        #threshold = threshold * -1
        
        pre_threshold_crossings = np.diff(convolved_y < threshold, append = False)
        pre_threshold_crossings = np.argwhere(pre_threshold_crossings)[:,0]
        pre_threshold_crossings = convolved_x[pre_threshold_crossings]
        
        post_threshold_crossings = np.diff(convolved_y < threshold, prepend = False)
        post_threshold_crossings = np.argwhere(post_threshold_crossings)
        post_threshold_crossings = convolved_x[post_threshold_crossings]
        
        plateau_starts = pre_threshold_crossings[::2]
        plateau_ends = post_threshold_crossings[1::2]
        
        gene_data.at[index, "Plateau_starts"] = plateau_starts
        gene_data.at[index, "Plateau_ends"] = plateau_ends
        
        if index > di.ENHANCER_CONVOLUTION:
            
            break
        
    return gene_data
    
    
    
    