import numpy as np
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d

import data_initialisation as di

def enhancer_convolution(genes_search, overlaps):
    
    print("Convolving enhancers...")
    
    genes_search["Start"] = genes_search["Search_window_start"].astype("int")
    genes_search["End"] = genes_search["Search_window_end"].astype("int")
    
    genes_search["Enhancer_searched_coordinates"] = ""
    genes_search["Enhancer_convolution"] = ""
    genes_search["Enhancer_convolved_coordinates"] = ""

    genes_search = genes_search.sort_values("Interest_score", ascending = False)

    i = 1
    
    for index, gene in genes_search.iterrows():
        
        window = get_kernel(di.ENHANCER_KERNEL_SHAPE, int((di.RELATIVE_ENHANCER_KERNEL_SIZE * (gene.End - gene.Start))), int(di.RELATIVE_ENHANCER_KERNEL_SIGMA * (gene.End - gene.Start)))
        gene_basewise = np.zeros((gene.End - gene.Start), dtype = int)
        basewise = np.arange(gene.Start, gene.End)
        gene_specific_enhancer_overlaps = overlaps.loc[overlaps["Gene_name"] == gene["Gene_name"]]
        
        for index, overlap in gene_specific_enhancer_overlaps.iterrows():
            
            overlap_basewise = np.where(np.logical_and(overlap.Start <= basewise, basewise <= overlap.End), 1, 0)
            gene_basewise = np.where(overlap_basewise == 1, 1, gene_basewise)
            
        gene["Enhancer_searched_coordinates"] = np.arange(gene.Start, gene.End)
        gene["Enhancer_convolution"] = np.convolve(window, gene_basewise)
        gene["Enhancer_convolved_coordinates"] = np.arange((gene.Start - (len(window) // 2)), (gene.Start - (len(window) // 2) + len(gene["Enhancer_convolution"])))

        print("Convolved enhancers for gene " + str(i) + " of " + str(di.ENHANCER_CONVOLUTION) + "...")
        i = i + 1
        
        if i > di.ENHANCER_CONVOLUTION:
            break
        
    genes_search.drop(["Start", "End"], axis = 1)
        
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
        
        for index, overlap in gene_specific_enhancer_overlaps.iterrows():
            
            overlap_basewise = np.where(np.logical_and(overlap.Start <= basewise, basewise <= overlap.End), 1, 0)
            gene_basewise = np.where(overlap_basewise == 1, 1, gene_basewise)
            
        gene["Quiescent_searched_coordinates"] = np.arange(gene.Start, gene.End)
        gene["Quiescent_convolution"] = np.convolve(window, gene_basewise)
        gene["Quiescent_convolved_coordinates"] = np.arange((gene.Start - (len(window) // 2)), (gene.Start - (len(window) // 2) + len(gene["Quiescent_convolution"])))

        print("Convolved quiescent regions for gene " + str(i) + " of " + str(di.ENHANCER_CONVOLUTION) + "...")
        i = i + 1
        
        if i > di.ENHANCER_CONVOLUTION:
            break
        
    genes_search.drop(["Start", "End"], axis = 1)
        
    return genes_search
    
def overlay_convolutions(enhancer_convolution, quiescent_convolution):
    
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
    
    print("Finding plateaus...")
    
    gene_data["plateau_starts"] = ""
    gene_data["plateau_ends"] = ""
    
    i = 1
    
    for index, gene in gene_data.iterrows():
        
        convolved_x = gene.loc[:, "Enhancer_convolved_coordinates"]
        convolved_y = gene.loc[:, "Enhancer_convolution"]
        threshold = di.PLATEAU_THRESHOLD
        
        convolved_y = np.negative(convolved_y)
        threshold = threshold * -1
        
        pre_threshold_crossings = np.diff(convolved_y < threshold, append = False)
        pre_threshold_crossings = np.argwhere(pre_threshold_crossings)[:,0]
        pre_threshold_crossings = convolved_x[pre_threshold_crossings]
        
        post_threshold_crossings = np.diff(convolved_y < threshold, prepend = False)
        post_threshold_crossings = np.argwhere(post_threshold_crossings)
        post_threshold_crossings = convolved_x[post_threshold_crossings]
        
        plateau_starts = pre_threshold_crossings[::2]
        plateau_ends = post_threshold_crossings[1::2]
        
        gene["plateau_starts"] = plateau_starts
        gene["plateau_ends"] = plateau_ends
        
        print("Found plateaus for gene " + str(i) + " of " + str(di.ENHANCER_CONVOLUTION) + "...")
        i = i + 1
        
        if i > di.ENHANCER_CONVOLUTION:
            break
        
    return gene_data
    
    
    
    