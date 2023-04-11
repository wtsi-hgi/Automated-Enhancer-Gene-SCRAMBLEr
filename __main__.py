import pandas as pd

import data_initialisation as di
import find_metrics as fm
import region_convolutions as rc
import data_visualisation as dv


def main():
    
    di.read_config_file()
    gene_annotations = di.read_gene_annotations()
    general_expression_data = di.read_general_expression_data()
    specific_expression_data = di.read_specific_expression_data()
    gene_data = pd.merge(gene_annotations, general_expression_data, on = "Gene_name", how = "inner")
    gene_data = pd.merge(gene_data, specific_expression_data, on = "Gene_name", how = "inner")
    
    del gene_annotations, general_expression_data, specific_expression_data
    
    gene_data = fm.find_gene_sizes(gene_data)
    gene_data = fm.find_interferring_genes(gene_data)
    gene_data = fm.find_search_windows(gene_data)
    
    regulatory_elements = di.read_regulatory_elements()
    enhancers, quiescent_regions = di.clean_regulatory_elements(regulatory_elements)
    del regulatory_elements
    
    enhancer_overlaps = fm.find_element_overlaps_within_search_window(enhancers, gene_data)
    if di.QUIESCENT_CONVOLUTION == True:
        quiescent_overlaps = fm.find_element_overlaps_within_search_window(quiescent_regions, gene_data)
    del enhancers, quiescent_regions
    
    gene_data = fm.count_overlaps_per_gene(gene_data, enhancer_overlaps, "Enhancer")
    gene_data = fm.find_nearby_enhancer_densities(gene_data, enhancer_overlaps)
    gene_data = fm.calculate_interest_score(gene_data)
    
    gene_data = rc.convolution(gene_data, enhancer_overlaps, "Enhancer")

    if di.QUIESCENT_CONVOLUTION == True:
        gene_data = rc.quiescent_convolution(gene_data, quiescent_overlaps)
        del quiescent_overlaps
    del enhancer_overlaps
    
    gene_data = rc.find_plateaus(gene_data)
    rc.export_convolutions(gene_data)
    rc.export_plateaus(gene_data)
    #dv.gene_report(gene_data)
    
    
    #enhancer_convolution, recombination_convolution = rc.overlay_convolutions(rc.enhancer_convolution(gene_data, enhancer_overlaps), rc.quiescent_convolution(gene_data, quiescent_overlaps)) 
    
if __name__ == "__main__":
    main()